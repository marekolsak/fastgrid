/*
    AutoGrid

    Copyright (C) 1989-2007, Garrett M. Morris, David S. Goodsell, Ruth Huey, Arthur J. Olson,
    All Rights Reserved.
    Copyright (C) 2008-2009, Marek Olsak (maraeo@gmail.com), All Rights Reserved.

    AutoGrid is a Trade Mark of The Scripps Research Institute.

    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation; either version 2
    of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
*/

#include "utils.h"
#include "programparameters.h"
#include "exceptions.h"
#include "pairwiseinteractionenergies.h"
#include "desolvexpfunc.h"
#include "bondvectors.h"
#include "inputdataloader.h"
#include "math.h"
#include <new>

// Useful macros

#define FOR_EACH_GRID_POINT(gridPos, outputIndex) \
    /* Z axis */ \
    for (int z = 0; z < input->numGridPoints[Z]; z++) \
    { \
        /* gridPos contains the current grid point. */ \
        double gridPosZ = (z - input->numGridPointsDiv2[Z]) * input->gridSpacing; \
        int outputIndexZBase = z * input->numGridPoints[X] * input->numGridPoints[Y]; \
\
        /* Y axis */ \
        for (int y = 0; y < input->numGridPoints[Y]; y++) \
        { \
            double gridPosY = (y - input->numGridPointsDiv2[Y]) * input->gridSpacing; \
            int outputIndexZYBase = outputIndexZBase + y * input->numGridPoints[X]; \
\
            /* X axis */ \
            for (int x = 0; x < input->numGridPoints[X]; x++) \
            { \
                double gridPos[XYZ] = {(x - input->numGridPointsDiv2[X]) * input->gridSpacing, gridPosY, gridPosZ}; \
                int outputIndex = outputIndexZYBase + x;

#define END_FOR() } } }

__forceinline double roundOutput(double a)
{
    if (fabs(a) < 0.0005)
        return 0;
    return round3dp(a);
}

void initCovalentMaps(const InputData *input, const GridMapList &gridmaps)
{
    // TODO: once covalent maps are supported, rewrite this function using NVIDIA CUDA

    for (int m = 0; m < gridmaps.getNumMaps(); m++)
        if (gridmaps[m].isCovalent)
#if defined(AG_OPENMP)
    #pragma AG_OPENMP_PARALLEL_FOR
#endif
            FOR_EACH_GRID_POINT(gridPos, outputIndex)
            {
                // Calculate the distance from the current grid point to the covalent attachment point
                double distance[XYZ];
                subtractVectors(distance, input->covalentPoint, gridPos);

                // Distance squared from current grid point to the covalent attachment point
                double rcovSq = lengthSquared(distance[X], distance[Y], distance[Z]);
                rcovSq = rcovSq * input->covHalfWidthSquaredInv;
                if (rcovSq < sq(APPROX_ZERO))
                    rcovSq = sq(APPROX_ZERO);
                double energy = input->covBarrier * (1 - exp(-0.69314718055994529 * rcovSq)); // -0.69314718055994529 = log(0.5)

                gridmaps[m].energies[outputIndex] = energy;
            }
            END_FOR();
}

void calculateElectrostaticMap(const InputData *input, const GridMapList &gridmaps, const ParameterLibrary &parameterLibrary)
{
    // TODO: rewrite this function using NVIDIA CUDA

#if defined(AG_OPENMP)
    #pragma AG_OPENMP_PARALLEL_FOR
#endif
    FOR_EACH_GRID_POINT(gridPos, outputIndex)
    {
        double &outputEnergy = gridmaps.getElectrostaticMap().energies[outputIndex];

        //  Do all Receptor (protein, DNA, etc.) atoms...
        for (int ia = 0; ia < input->numReceptorAtoms; ia++)
        {
            //  Get distance, r, from current grid point, c, to this receptor atom, input->receptorAtomCoord,
            double distance[XYZ];
            subtractVectors(distance, input->receptorAtomCoord[ia], gridPos);
            double rSq = lengthSquared(distance[X], distance[Y], distance[Z]);
            if (rSq < sq(APPROX_ZERO))
                rSq = sq(APPROX_ZERO);
            double invR = rsqrt(rSq);

            // apply the estat forcefield coefficient/weight here
            double tmp = input->charge[ia] * min(invR, 2.0) * parameterLibrary.coeff_estat;
            if (input->distDepDiel) // distDepDiel is a constant
            {
                // Distance-dependent dielectric...
                int indexR = min(lookup(1 / invR), MAX_DIST-1); // make sure lookup index is in the table
                outputEnergy += tmp * input->epsilon[indexR];
            }
            else
                // Constant dielectric...
                outputEnergy += tmp * input->invDielCal;
        }

        // Round to 3 decimal places
        outputEnergy = roundOutput(outputEnergy);
    }
    END_FOR();
}

void calculateFloatingGrid(const InputData *input, const GridMapList &gridmaps)
{
    // TODO: rewrite this function using NVIDIA CUDA

#if defined(AG_OPENMP)
    #pragma AG_OPENMP_PARALLEL_FOR
#endif
    FOR_EACH_GRID_POINT(gridPos, outputIndex)
    {
        double fgridInvRMin = 1.0 / BIG;

        //  Do all Receptor (protein, DNA, etc.) atoms...
        for (int ia = 0; ia < input->numReceptorAtoms; ia++)
        {
            //  Get distance, r, from current grid point, c, to this receptor atom, input->receptorAtomCoord,
            double distance[XYZ];
            subtractVectors(distance, input->receptorAtomCoord[ia], gridPos);
            double rSq = lengthSquared(distance[X], distance[Y], distance[Z]);
            if (rSq < sq(APPROX_ZERO))
                rSq = sq(APPROX_ZERO);

            // Calculate the so-called "Floating Grid"...
            fgridInvRMin = max(rsqrt(rSq), fgridInvRMin);
        }

        gridmaps.getFloatingGridMins()[outputIndex] = float(1 / fgridInvRMin);
    }
    END_FOR();
}

struct HBondInfo
{
    struct Info
    {
        double min, max;
        bool flag;
    } info[MAX_MAPS];

    __forceinline HBondInfo(int numAtomMaps)
    {
        Info def;
        def.min = 999999;
        def.max = -999999;
        def.flag = false;

        for (int m = 0; m < numAtomMaps; m++)
            info[m] = def;
    }

    __forceinline Info &operator [](int i)
    {
        return info[i];
    }

    __forceinline void insert(int mapIndex, double energy)
    {
        Info &i = info[mapIndex];
        i.min = min(i.min, energy);
        i.max = max(i.max, energy);
        i.flag = true;
    }
};

__forceinline int findClosestHBond(const InputData *input, const double *gridPos)
{
    int closestH = 0;
    double rminSq = sq(999999.0);
    for (int ia = 0; ia < input->numReceptorAtoms; ia++)
        if (input->hbond[ia] == DS || input->hbond[ia] == D1)
        {
            // DS or D1
            double distance[XYZ];
            subtractVectors(distance, input->receptorAtomCoord[ia], gridPos);
            double rSq = lengthSquared(distance[X], distance[Y], distance[Z]);
            if (rSq < rminSq)
            {
                rminSq = rSq;
                closestH = ia;
            }
        }           // Hydrogen test
    return closestH;
}

__forceinline void getHBondAngularFunction(const InputData *input, const BondVectors *bondVectors, int ia, int closestH, double (&distance)[XYZ],
                             double &racc, double &rdon, double &Hramp) // output
{
    racc = 1;
    rdon = 1;
    Hramp = 1;   // Hramp ramps in Hbond acceptor probes

    //  cosTheta = distance dot bondVectors->rvector == cos(angle) subtended.
    double cosTheta = -dotProduct(distance, bondVectors->rvector[ia]);

    // Calculate racc/rdon/Hramp
    if (input->hbond[ia] == D1)
    {
        // D1
        //  ia-th receptor atom = Hydrogen (4 = H)
        //  => receptor H-bond donor, OH or NH.
        //  calculate racc for H-bond ACCEPTOR PROBES at this grid pt.

        racc = 0;

        //  H->current-grid-pt vector < 90 degrees from
        //  N->H or O->H vector,
        if (cosTheta > 0)
        {
            //  racc = [cos(theta)]^2.0 for N-H
            //  racc = [cos(theta)]^4.0 for O-H
            racc = cosTheta;

            switch (bondVectors->rexp[ia]) // racc = pow(cosTheta, bondVectors->rexp[ia]);
            {
            case 4:
                racc = sq(racc);
            case 2:
                racc = sq(racc);
            }

            // NEW2 calculate dot product of bond vector with bond vector of best input->hbond
            if (ia != closestH)
            {
                double theta = angle(bondVectors->rvector[closestH], bondVectors->rvector[ia]);
                Hramp = 0.5 - 0.5 * cos(theta * (120.0 / 90.0));
            }   // ia test
            // END NEW2 calculate dot product of bond vector with bond vector of best input->hbond
        }
        // endif (input->atomType[ia] == hydrogen)
        // NEW Directional N acceptor
    }
    else if (input->hbond[ia] == A1)
    {
        // A1
        //  ia-th macromolecule atom = Nitrogen (4 = H)
        //  calculate rdon for H-bond Donor PROBES at this grid pt.

        //  H->current-grid-pt vector < 90 degrees from X->N vector
        rdon = 0;
        if (cosTheta > 0)
            rdon = sq(cosTheta); // for H->N
        // endif (input->atomType[ia] == nitrogen)
        // end NEW Directional N acceptor
    }
    else if (input->hbond[ia] == A2)
    {
        rdon = 0;

        if (bondVectors->disorder[ia])
        {
            // A2
            // cylindrically disordered hydroxyl

            racc = 0;
            double theta = acos_clamped(cosTheta);
            if (theta <= 1.24791 + (PI / 2))
            {
                // 1.24791 rad = 180 deg minus C-O-H bond angle, ** 108.5 deg
                rdon = sq(sq(cos(theta - 1.24791)));    // pow(.., 4)
                racc = rdon;
            }
        }
        else
        {
            // A2
            //  ia-th receptor atom = Oxygen
            //  => receptor H-bond acceptor, oxygen.

            // rdon expressions from Goodford
            if (cosTheta >= 0)
            {
                // ti is the angle in the lone pair plane, away from the
                // vector between the lone pairs,
                // calculated as (grid vector CROSS lone pair plane normal)
                // DOT C=O vector - 90 deg
                double cross[XYZ];
                crossProduct(cross, distance, bondVectors->rvector2[ia]);
                double rd2 = lengthSquared(cross[0], cross[1], cross[2]);
                if (rd2 < APPROX_ZERO)
                    rd2 = APPROX_ZERO;
                double inv_rd = rsqrt(rd2);

                double ti = fabs(acos_clamped(inv_rd * dotProduct(cross, bondVectors->rvector[ia])) - (PI / 2));
                // the 2.0*ti can be replaced by (ti + ti) in: rdon = (0.9 + 0.1*sin(2.0*ti))*cos(t0);
                rdon = 0.9 + 0.1 * sin(ti + ti);
                // 0.34202 = cos (100 deg)
            }
            else if (cosTheta >= -0.34202)
                rdon = 562.25 * cube(0.116978 - sq(cosTheta));

            // t0 is the angle out of the lone pair plane, calculated
            // as 90 deg - acos (vector to grid point DOT lone pair
            // plane normal)
            double t0 = (PI / 2) - angle(distance, bondVectors->rvector2[ia]);
            rdon *= cos(t0);

            // endif input->atomType == oxygen, not disordered
        }
    }
}

__forceinline void sumPairwiseInteractions(const InputData *input, const GridMapList &gridmaps, const PairwiseInteractionEnergies &energyLookup,
                             const DesolvExpFunc &desolvExpFunc, const BondVectors *bondVectors, HBondInfo &hbond,
                             int outputIndex, int m, int ia, int indexR, int hydrogen, double racc, double rdon, double Hramp)
{
    double &e = gridmaps[m].energies[outputIndex];
    double pwiEnergy = energyLookup(input->atomType[ia], indexR, m);

    if (gridmaps[m].isHBonder)
    {
        // PROBE forms H-bonds...

        // rsph ramps in angular dependence for distances with negative energy
        double rsph = clamp(pwiEnergy / 100, 0.0, 1.0);

        if ((gridmaps[m].hbond == AS || gridmaps[m].hbond == A2) && (input->hbond[ia] == DS || input->hbond[ia] == D1))
            // PROBE can be an H-BOND ACCEPTOR,
            e += (bondVectors->disorder[ia] ? energyLookup(hydrogen, max(0, indexR - 110), m) : pwiEnergy) *
                 Hramp * (racc + (1 - racc) * rsph);
        else if (gridmaps[m].hbond == A1 && (input->hbond[ia] == DS || input->hbond[ia] == D1))
            // A1 vs DS, D1
            hbond.insert(m, pwiEnergy * (racc + (1 - racc) * rsph));
        else if ((gridmaps[m].hbond == DS || gridmaps[m].hbond == D1) && input->hbond[ia] >= AS)
            // DS,D1 vs AS,A1,A2
            // PROBE is H-BOND DONOR,
            hbond.insert(m, pwiEnergy * (rdon + (1 - rdon) * rsph));
        else
            // hbonder PROBE-ia cannot form a H-bond...,
            e += pwiEnergy;
    }
    else
        // PROBE does not form H-bonds...,
        e += pwiEnergy;

    // add desolvation energy
    // forcefield desolv coefficient/weight in desolvExpFunc
    e += gridmaps[m].solparProbe * input->vol[ia] * desolvExpFunc(indexR) +
        (input->solpar[ia] + input->solparQ * fabs(input->charge[ia])) * gridmaps[m].volProbe * desolvExpFunc(indexR);
}

void calculateGridmaps(const InputData *input, const GridMapList &gridmaps, const ParameterLibrary &parameterLibrary,
                       const PairwiseInteractionEnergies &energyLookup, const DesolvExpFunc &desolvExpFunc, const BondVectors *bondVectors)
{
    int hydrogen = parameterLibrary.getAtomParameterRecIndex("HD");

    /*
        Ulozim si do mrizky pozice vsech receptor atoms.
        Sirka mrizky bude NBCUTOFF * 2, pak se staci divat jen do ctyr bunek. Indexovat se bude pomoci gridPos nebo {x,y,z}.
    */

#if defined(AG_OPENMP)
    #pragma AG_OPENMP_PARALLEL_FOR
#endif
    FOR_EACH_GRID_POINT(gridPos, outputIndex)
    {
        HBondInfo hbond(gridmaps.getNumAtomMaps());
        int closestH = findClosestHBond(input, gridPos);

        //  Do all Receptor (protein, DNA, etc.) atoms...
        for (int ia = 0; ia < input->numReceptorAtoms; ia++)
        {
            // If distance from grid point to atom ia is too large,
            // or if atom is a disordered hydrogen,
            //   add nothing to the grid-point's non-bond energy;
            //   just continue to next atom...

            //  distance[] = Unit vector from current grid pt to ia_th m/m atom.
            //  Get distance, r, from current grid point, c, to this receptor atom, input->receptorAtomCoord,
            double distance[XYZ];
            subtractVectors(distance, input->receptorAtomCoord[ia], gridPos);

            // rSq = |distance|^2
            double rSq = lengthSquared(distance[X], distance[Y], distance[Z]);

            if (rSq > sq(NBCUTOFF))
                continue;   // onto the next atom...
            if (input->atomType[ia] == hydrogen && bondVectors->disorder[ia])
                continue;   // onto the next atom...

            // Normalize the distance vector
            if (rSq < sq(APPROX_ZERO))
                rSq = sq(APPROX_ZERO);
            double invR = rsqrt(rSq);
            scalarProduct(distance, distance, invR);

            int indexR = min(lookup(1 / invR), MAX_DIST-1); // make sure lookup index is in the table

            double racc, rdon, Hramp;
            getHBondAngularFunction(input, bondVectors, ia, closestH, distance, racc, rdon, Hramp);

            // For each probe atom-type,
            // Sum pairwise interactions between each probe
            // at this grid point (gridPos[0:2])
            // and the current receptor atom, ia...
            for (int m = 0; m < gridmaps.getNumAtomMaps(); m++)
                if (!gridmaps[m].isCovalent)
                    sumPairwiseInteractions(input, gridmaps, energyLookup, desolvExpFunc, bondVectors, hbond, outputIndex, m, ia, indexR, hydrogen, racc, rdon, Hramp);

            gridmaps.getDesolvationMap().energies[outputIndex] += input->solparQ * input->vol[ia] * desolvExpFunc(indexR);
        } // ia loop, over all receptor atoms...

        for (int m = 0; m < gridmaps.getNumAtomMaps(); m++)
        {
            double &e = gridmaps[m].energies[outputIndex];
            if (hbond[m].flag)
                e += hbond[m].min + hbond[m].max;
            e = roundOutput(e);
        }
        double &e = gridmaps.getDesolvationMap().energies[outputIndex];
        e = roundOutput(e);
    }
    END_FOR();
}

// Function: Calculation of interaction energy grids for Autodock.
// Directional H_bonds from Goodford:
// Distance dependent dielectric after Mehler and Solmajer.
// Charge-based desolvation
// Copyright: (C) 2004, TSRI
//
// Authors: Garrett Matthew Morris, Ruth Huey, David S. Goodsell
//
// The Scripps Research Institute
// Department of Molecular Biology, MB5
// 10550 North Torrey Pines Road
// La Jolla, CA 92037-1000.
//
// e-mail: garrett@scripps.edu
// rhuey@scripps.edu
// goodsell@scripps.edu
//
// Helpful suggestions and advice:
// Arthur J. Olson
// Bruce Duncan, Yng Chen, Michael Pique, Victoria Roberts
// Lindy Lindstrom
//
// Inputs: Control file, receptor PDBQT file, parameter file
// Returns: Atomic affinity, desolvation and electrostatic grid maps.
void autogridMain(int argc, char **argv)
{
    // Get the time at the start of the run...
    tms tmsJobStart;
    Clock jobStart = times(&tmsJobStart);

    double versionNumber = 4.00;

    // Initialize the ProgramParameters object, which parses the command-line arguments
    ProgramParameters programParams(argc, argv);

    // Initialize the log file
    LogFile logFile(versionNumber, programParams.getProgramName(), programParams.getLogFilename());

    // Declaration of gridmaps, InputDataLoader::load takes care of their initialization
    GridMapList gridmaps(&logFile);

    // Initialization of free energy coefficients and atom parameters
    ParameterLibrary parameterLibrary(&logFile, programParams.getDebugLevel());

    // Reading in the grid parameter file
    InputDataLoader *inputDataLoader = new InputDataLoader(&logFile);
    inputDataLoader->load(programParams.getGridParameterFilename(), gridmaps, parameterLibrary);
    // TODO: shouldn't we put these out of the load function? :
    // - gridmaps initialization code
    // - initialization of atom parameters recIndex/mapIndex (in parameterLibrary)

    // Now we want to make the input data read-only
    const InputData *input = inputDataLoader;

    if (input->floatingGridFilename[0])
        gridmaps.enableFloatingGrid();

    // Inititializing arrays of output energies
    gridmaps.prepareGridmaps(input->numGridPointsPerMap);

    // TODO: add a smarter mechanism of checking for the available disk space, we need to know it as soon as possible.
    // the formerly implemented checks in the middle of calculations were done too late

    // Loading the parameter library from the file
    if (input->parameterLibraryFilename[0])
        parameterLibrary.load(input->parameterLibraryFilename);

    // Writing to AVS-readable gridmaps file (fld)
    saveAVSGridmapsFile(gridmaps, input, programParams, logFile);

#if defined(BOINCCOMPOUND)
    boinc_fraction_done(0.1);
#endif

    beginTimer(0);

    // Calculating the lookup table of the pairwise interaction energies
    PairwiseInteractionEnergies energyLookup;
    energyLookup.calculate(gridmaps, logFile, input->numReceptorTypes, input->receptorTypes, input->rSmooth);

    // Precalculating the exponential function for receptor and ligand desolvation
    DesolvExpFunc desolvExpFunc(parameterLibrary.coeff_desolv);

    // Calculating bond vectors for directional H-bonds
    BondVectors *bondVectors = new BondVectors(&logFile);
    bondVectors->calculate(input, parameterLibrary);

    endTimer(0);

    logFile.printFormatted("Beginning grid calculations.\n"
                           "\nCalculating %d grids over %d elements, around %d receptor atoms.\n\n"
                           "                    Percent   Estimated Time  Time/this plane\n"
                           "XY-plane  Z-coord   Done      Remaining       Real, User, System\n"
                           "            /Ang              /sec            /sec\n"
                           "________  ________  ________  ______________  __________________________\n\n",
                           gridmaps.getNumMapsInclFloatingGrid(), input->numGridPointsPerMap, input->numReceptorAtoms);

    // TODO: rewrite writing out progress in percents
    /* Former code:

        for (all z)
        {
            tms timesGridStart;
            Clock gridStartTime = times(&timesGridStart);

            for (all y) ...

            tms timerGridEnd;
            Clock gridEndTime = times(&timerGridEnd);
            logFile.printFormatted(" %6d   %8.3lf   %5.1lf%%   ", gridCoordZ, input->gridCornerMin[Z] + gridPosZ, ((z+1) * 100.0) / input->numGridPoints[Z]);
            logFile.printTimeInHMS((gridEndTime - gridStartTime) * (input->numGridPoints[Z] - z));
            logFile.print("  ");
            logFile.printExecutionTimes(gridStartTime, gridEndTime, &timesGridStart, &timerGridEnd);
        }
    */

    beginTimer(1);
    // Covalent Atom Types are not yet supported with the new AG4/AD4 atom typing mechanism...
    initCovalentMaps(input, gridmaps);
    endTimer(1);

    beginTimer(2);
    // Calculation of the electrostatic map
    calculateElectrostaticMap(input, gridmaps, parameterLibrary);
    endTimer(2);

    // Calculation of the atom maps and the desolvation map
    beginTimer(3);
    calculateGridmaps(input, gridmaps, parameterLibrary, energyLookup, desolvExpFunc, bondVectors);
    endTimer(3);

    // Calculate the so-called "floating grid"
    if (gridmaps.containsFloatingGrid())
    {
        beginTimer(4);
        calculateFloatingGrid(input, gridmaps);
        beginTimer(4);
    }

    delete bondVectors;

#if defined(BOINCCOMPOUND)
    boinc_fraction_done(0.9);
#endif

    // Save all gridmaps
    gridmaps.saveToFiles(input, programParams.getGridParameterFilename());

    delete inputDataLoader;

    // Writing out summary
    gridmaps.logSummary();

    // Get the time at the end of the run and print the difference
    tms tmsJobEnd;
    Clock jobEnd = times(&tmsJobEnd);
    logFile.printExecutionTimesInHMS(jobStart, jobEnd, &tmsJobStart, &tmsJobEnd);
}

int main(int argc, char **argv)
{
    try
    {
        // Initialize BOINC if needed
        boincInit();

        // AutoGrid's main function
        autogridMain(argc, argv);

        // This should not return if used
        boincDone();

        logTimers();
        return 0;
    }
    catch (ExitProgram &e)  // the ExitProgram exception is a replacement for C's exit function.
    {
        return e.getExitCode();
    }
    catch (std::bad_alloc &)
    {
        fprintf(stderr, "\n%s: FATAL ERROR: Not enough memory!\n", *argv);
        return 0xBADA110C; // BADALLOC
    }
}
