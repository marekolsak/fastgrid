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

//#define AG_OPENMP

// Useful macros

#define FOR_EACH_GRID_POINT(gridPos, outputIndex) \
    /* Z axis */ \
    for (int z = 0; z < input->numGridPoints[Z]; z++) \
    { \
        /* gridPos contains the current grid point. */ \
        double gridPosZ = (z - input->ne[Z]) * input->spacing; \
        int outputIndexZBase = z * input->numGridPoints[X] * input->numGridPoints[Y]; \
\
        /* Y axis */ \
        for (int y = 0; y < input->numGridPoints[Y]; y++) \
        { \
            double gridPosY = (y - input->ne[Y]) * input->spacing; \
            int outputIndexZYBase = outputIndexZBase +  y * input->numGridPoints[X]; \
\
            /* X axis */ \
            for (int x = 0; x < input->numGridPoints[X]; x++) \
            { \
                double gridPos[XYZ] = {(x - input->ne[X]) * input->spacing, gridPosY, gridPosZ}; \
                int outputIndex = outputIndexZYBase + x;

#define END_FOR() } } }

double roundOutput(double a)
{
    a = round3dp(a);
    if (fabs(a) < PRECISION)
        return 0;
    return a;
}

double getEnergyForCovalentMap(const InputData *input, const double *gridPos)
{
    // Calculate the distance from the current grid point to the covalent attachment point
    double distance[XYZ];
    for (int i = 0; i < XYZ; i++)
        distance[i] = input->covalentPoint[i] - gridPos[i];

    // Distance squared from current grid point to the covalent attachment point
    double rcovSq = hypotenuseSq(distance[X], distance[Y], distance[Z]);
    rcovSq = rcovSq / (input->covHalfWidth * input->covHalfWidth);
    if (rcovSq < APPROX_ZERO_SQ)
        rcovSq = APPROX_ZERO_SQ;
    return input->covBarrier * (1 - exp(-0.69314718055994529 * rcovSq)); // -0.69314718055994529 = log(0.5)
}

void calculateElectrostaticMap(const InputData *input, const GridMapList &gridmaps, const ParameterLibrary &parameterLibrary)
{
    // TODO: rewrite this function using NVIDIA CUDA

#if defined(AG_OPENMP)
    #pragma omp parallel for schedule(dynamic, 1)
#endif
    FOR_EACH_GRID_POINT(gridPos, outputIndex)
    {
        double &outputEnergy = gridmaps.getElectrostaticMap().energies[outputIndex];

        // Covalent Atom Types are not yet supported with the new AG4/AD4 atom typing mechanism...
        if (gridmaps.getElectrostaticMap().isCovalent)
            outputEnergy = getEnergyForCovalentMap(input, gridPos);

        //  Do all Receptor (protein, DNA, etc.) atoms...
        for (int ia = 0; ia < input->numReceptorAtoms; ia++)
        {
            //  Get distance, r, from current grid point, c, to this receptor atom, input->coord,
            double distance[XYZ];
            for (int i = 0; i < XYZ; i++)
                distance[i] = input->coord[ia][i] - gridPos[i];
            double rSq = hypotenuseSq(distance[X], distance[Y], distance[Z]);
            if (rSq < APPROX_ZERO_SQ)
                rSq = APPROX_ZERO_SQ;
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
#if defined(AG_OPENMP)
    #pragma omp parallel for schedule(dynamic, 1)
#endif
    FOR_EACH_GRID_POINT(gridPos, outputIndex)
    {
        double fgridInvRMin = 1.0 / BIG;

        //  Do all Receptor (protein, DNA, etc.) atoms...
        for (int ia = 0; ia < input->numReceptorAtoms; ia++)
        {
            //  Get distance, r, from current grid point, c, to this receptor atom, input->coord,
            double distance[XYZ];
            for (int i = 0; i < XYZ; i++)
                distance[i] = input->coord[ia][i] - gridPos[i];
            double rSq = hypotenuseSq(distance[X], distance[Y], distance[Z]);
            if (rSq < APPROX_ZERO_SQ)
                rSq = APPROX_ZERO_SQ;

            // Calculate the so-called "Floating Grid"...
            fgridInvRMin = max(rsqrt(rSq), fgridInvRMin);
        }

        gridmaps.getFloatingGridMins()[outputIndex] = float(1 / fgridInvRMin);
    }
    END_FOR();
}

struct HBondInfo
{
    double min[MAX_MAPS], max[MAX_MAPS];
    bool flag[MAX_MAPS];

    HBondInfo(int numAtomMaps)
    {
        for (int m = 0; m < numAtomMaps; m++)
        {
            min[m] = 999999;
            max[m] = -999999;
            flag[m] = false;
        }
    }
};

int findClosestHBond(const InputData *input, const double *gridPos)
{
    int closestH = 0;
    double rminSq = sq(999999.0);
    for (int ia = 0; ia < input->numReceptorAtoms; ia++)
        if ((input->hbond[ia] == DS) || (input->hbond[ia] == D1))
        {
            // DS or D1
            double distance[XYZ];
            for (int i = 0; i < XYZ; i++)
                distance[i] = input->coord[ia][i] - gridPos[i];
            double rSq = hypotenuseSq(distance[X], distance[Y], distance[Z]);
            if (rSq < rminSq)
            {
                rminSq = rSq;
                closestH = ia;
            }
        }           // Hydrogen test
    return closestH;
}

void calculateGridmaps(const InputData *input, const GridMapList &gridmaps, const ParameterLibrary &parameterLibrary,
                       const PairwiseInteractionEnergies &energyLookup, const DesolvExpFunc &desolvExpFunc, const BondVectors *bondVectors)
{
    int hydrogen = parameterLibrary.getAtomParameterRecIndex("HD");

#if defined(AG_OPENMP)
    #pragma omp parallel for schedule(dynamic, 1)
#endif
    FOR_EACH_GRID_POINT(gridPos, outputIndex)
    {
        // Covalent Atom Types are not yet supported with the new AG4/AD4 atom typing mechanism...
        for (int m = 0; m < gridmaps.getNumAtomMaps(); m++)
            if (gridmaps[m].isCovalent)
                gridmaps[m].energies[outputIndex] = getEnergyForCovalentMap(input, gridPos);
        if (gridmaps.getDesolvationMap().isCovalent)
            gridmaps.getDesolvationMap().energies[outputIndex] = getEnergyForCovalentMap(input, gridPos);

        HBondInfo hbond(gridmaps.getNumAtomMaps());
        int closestH = findClosestHBond(input, gridPos);

        //  Do all Receptor (protein, DNA, etc.) atoms...
        for (int ia = 0; ia < input->numReceptorAtoms; ia++)
        {
            //  Get distance, r, from current grid point, c, to this receptor atom, input->coord,
            double distance[XYZ];
            for (int i = 0; i < XYZ; i++)
                distance[i] = input->coord[ia][i] - gridPos[i];
            double rSq = hypotenuseSq(distance[X], distance[Y], distance[Z]);
            if (rSq < APPROX_ZERO_SQ)
                rSq = APPROX_ZERO_SQ;

            // If distance from grid point to atom ia is too large,
            // or if atom is a disordered hydrogen,
            //   add nothing to the grid-point's non-bond energy;
            //   just continue to next atom...

            if (rSq > sq(NBCUTOFF))
                continue;   // onto the next atom...
            if (input->atomType[ia] == hydrogen && bondVectors->disorder[ia])
                continue;   // onto the next atom...

            // Normalize the distance vector
            double invR = rsqrt(rSq);
            for (int i = 0; i < XYZ; i++)
                distance[i] *= invR;

            int indexR = min(lookup(1 / invR), MAX_DIST-1); // make sure lookup index is in the table

            double racc = 1;
            double rdon = 1;
            double Hramp = 1;   // Hramp ramps in Hbond acceptor probes

            // Calculate racc/rdon/Hramp
            if (input->hbond[ia] == D1)
            {
                // D1
                //  ia-th receptor atom = Hydrogen (4 = H)
                //  => receptor H-bond donor, OH or NH.
                //  calculate racc for H-bond ACCEPTOR PROBES at this grid pt.
                double cosTheta = 0;
                //  distance[] = Unit vector from current grid pt to ia_th m/m atom.
                //  cosTheta = d dot bondVectors->rvector == cos(angle) subtended.
                for (int i = 0; i < XYZ; i++)
                    cosTheta -= distance[i] * bondVectors->rvector[ia][i];

                if (cosTheta <= 0)
                    //  H->current-grid-pt vector >= 90 degrees from
                    //  N->H or O->H vector,
                    racc = 0;
                else
                {
                    //  racc = [cos(theta)]^2.0 for N-H
                    //  racc = [cos(theta)]^4.0 for O-H,
                    racc = cosTheta;

                    switch (bondVectors->rexp[ia])
                    {
                    case 4:
                        racc = sq(racc);
                    case 2:
                        racc = sq(racc);
                    }
                    // racc = pow(cosTheta, bondVectors->rexp[ia]);

                    // NEW2 calculate dot product of bond vector with bond vector of best input->hbond
                    if (ia == closestH)
                        Hramp = 1;
                    else
                    {
                        cosTheta = 0;
                        for (int i = 0; i < XYZ; i++)
                            cosTheta += bondVectors->rvector[closestH][i] * bondVectors->rvector[ia][i];
                        double theta = acos(clamp(cosTheta, -1.0, 1.0));
                        Hramp = 0.5 - 0.5 * cos(theta * (120.0 / 90.0));
                    }   // ia test
                    // END NEW2 calculate dot product of bond vector with bond vector of best input->hbond
                }
                // endif (input->atomType[ia] == hydrogen)
                // NEW Directional N acceptor
            }
            else if (input->hbond[ia] == A1)
            {           // A1
                //  ia-th macromolecule atom = Nitrogen (4 = H)
                //  calculate rdon for H-bond Donor PROBES at this grid pt.
                double cosTheta = 0;
                //  distance[] = Unit vector from current grid pt to ia_th m/m atom.
                //  cosTheta = d dot bondVectors->rvector == cos(angle) subtended.
                for (int i = 0; i < XYZ; i++)
                    cosTheta -= distance[i] * bondVectors->rvector[ia][i];

                if (cosTheta <= 0)
                    //  H->current-grid-pt vector >= 90 degrees from
                    //  X->N vector,
                    rdon = 0;
                else
                    //  racc = [cos(theta)]^2.0 for H->N
                    rdon = cosTheta * cosTheta;
                // endif (input->atomType[ia] == nitrogen)
                // end NEW Directional N acceptor
            }
            else if (input->hbond[ia] == A2 && !bondVectors->disorder[ia])
            {
                // A2
                //  ia-th receptor atom = Oxygen
                //  => receptor H-bond acceptor, oxygen.

                // check to see that probe is in front of oxygen, not behind
                double cosTheta = 0;
                for (int i = 0; i < XYZ; i++)
                    cosTheta -= distance[i] * bondVectors->rvector[ia][i];
                // t0 is the angle out of the lone pair plane, calculated
                // as 90 deg - acos (vector to grid point DOT lone pair
                // plane normal)
                double t0 = 0;
                for (int i = 0; i < XYZ; i++)
                    t0 += distance[i] * bondVectors->rvector2[ia][i];
                t0 = (PI / 2) - acos(clamp(t0, -1.0, 1.0));

                // ti is the angle in the lone pair plane, away from the
                // vector between the lone pairs,
                // calculated as (grid vector CROSS lone pair plane normal)
                // DOT C=O vector - 90 deg
                double cross[XYZ];
                cross[0] = distance[1] * bondVectors->rvector2[ia][2] - distance[2] * bondVectors->rvector2[ia][1];
                cross[1] = distance[2] * bondVectors->rvector2[ia][0] - distance[0] * bondVectors->rvector2[ia][2];
                cross[2] = distance[0] * bondVectors->rvector2[ia][1] - distance[1] * bondVectors->rvector2[ia][0];
                double rd2 = hypotenuseSq(cross[0], cross[1], cross[2]);
                if (rd2 < APPROX_ZERO)
                    rd2 = APPROX_ZERO;
                double inv_rd = rsqrt(rd2);
                double ti = 0;
                for (int i = 0; i < XYZ; i++)
                    ti += cross[i] * inv_rd * bondVectors->rvector[ia][i];

                // rdon expressions from Goodford
                rdon = 0;
                if (cosTheta >= 0)
                {
                    ti = fabs(acos(clamp(ti, -1.0, 1.0)) - (PI / 2));
                    // the 2.0*ti can be replaced by (ti + ti) in: rdon = (0.9 + 0.1*sin(2.0*ti))*cos(t0);
                    rdon = (0.9 + 0.1 * sin(ti + ti)) * cos(t0);
                    // 0.34202 = cos (100 deg)
                }
                else if (cosTheta >= -0.34202)
                    rdon = 562.25 * cube(0.116978 - sq(cosTheta)) * cos(t0);

                // endif input->atomType == oxygen, not disordered
            }
            else if (input->hbond[ia] == A2 && bondVectors->disorder[ia])
            {
                // A2
                // cylindrically disordered hydroxyl
                double cosTheta = 0;
                for (int i = 0; i < XYZ; i++)
                    cosTheta -= distance[i] * bondVectors->rvector[ia][i];
                racc = 0;
                rdon = 0;
                double theta = acos(clamp(cosTheta, -1.0, 1.0));
                if (theta <= 1.24791 + (PI / 2))
                {
                    // 1.24791 rad = 180 deg minus C-O-H bond angle, ** 108.5 deg
                    rdon = sq(sq(cos(theta - 1.24791)));    // pow(.., 4)
                    racc = rdon;
                }
            }           // input->atomType test

            // For each probe atom-type,
            // Sum pairwise interactions between each probe
            // at this grid point (gridPos[0:2])
            // and the current receptor atom, ia...
            for (int m = 0; m < gridmaps.getNumAtomMaps(); m++)
            {
                // We do not want to change the current enrg value for any covalent maps, make sure iscovalent is false...
                if (!gridmaps[m].isCovalent)
                {
                    double pwiEnergy = energyLookup(input->atomType[ia], indexR, m);

                    if (gridmaps[m].isHBonder)
                    {
                        // PROBE forms H-bonds...

                        // rsph ramps in angular dependence for distances with negative energy
                        double rsph = clamp(pwiEnergy / 100, 0.0, 1.0);

                        if ((gridmaps[m].hbond == AS || gridmaps[m].hbond == A2)    // AS or A2
                            && (input->hbond[ia] == DS || input->hbond[ia] == D1))
                        {
                            // DS or D1
                            // PROBE can be an H-BOND ACCEPTOR,
                            double f = Hramp * (racc + (1 - racc) * rsph);
                            if (!bondVectors->disorder[ia])
                                gridmaps[m].energies[outputIndex] += f * pwiEnergy;
                            else
                                gridmaps[m].energies[outputIndex] += f * energyLookup(hydrogen, max(0, indexR - 110), m);
                        }
                        else if ((gridmaps[m].hbond == A1)    // A1
                                 && (input->hbond[ia] == DS || input->hbond[ia] == D1))
                        {
                            // DS,D1
                            double hbondEnergy = pwiEnergy * (racc + (1 - racc) * rsph);
                            hbond.min[m] = min(hbond.min[m], hbondEnergy);
                            hbond.max[m] = max(hbond.max[m], hbondEnergy);
                            hbond.flag[m] = true;
                        }
                        else if ((gridmaps[m].hbond == DS || gridmaps[m].hbond == D1) && (input->hbond[ia] > D1))
                        {
                            // DS,D1 vs AS,A1,A2
                            // PROBE is H-BOND DONOR,
                            double hbondEnergy = pwiEnergy * (rdon + (1 - rdon) * rsph);
                            hbond.min[m] = min(hbond.min[m], hbondEnergy);
                            hbond.max[m] = max(hbond.max[m], hbondEnergy);
                            hbond.flag[m] = true;
                        }
                        else
                            // hbonder PROBE-ia cannot form a H-bond...,
                            gridmaps[m].energies[outputIndex] += pwiEnergy;
                    }
                    else
                        // PROBE does not form H-bonds...,
                        gridmaps[m].energies[outputIndex] += pwiEnergy;

                    // add desolvation energy
                    // forcefield desolv coefficient/weight in desolvExpFunc
                    gridmaps[m].energies[outputIndex] += gridmaps[m].solparProbe * input->vol[ia] * desolvExpFunc(indexR) +
                        (input->solpar[ia] + input->solparQ * fabs(input->charge[ia])) * gridmaps[m].volProbe * desolvExpFunc(indexR);
                } // is not covalent
            }
            gridmaps.getDesolvationMap().energies[outputIndex] += input->solparQ * input->vol[ia] * desolvExpFunc(indexR);
        } // ia loop, over all receptor atoms...

        for (int m = 0; m < gridmaps.getNumAtomMaps(); m++)
        {
            double &e = gridmaps[m].energies[outputIndex];
            if (hbond.flag[m])
            {
                e += hbond.min[m];
                e += hbond.max[m];
            }
            e = roundOutput(e);
        }
        gridmaps.getDesolvationMap().energies[outputIndex] = roundOutput(gridmaps.getDesolvationMap().energies[outputIndex]);
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

    beginTimer(1);

    // Calculating the lookup table of the pairwise interaction energies
    PairwiseInteractionEnergies energyLookup;
    energyLookup.calculate(gridmaps, logFile, input->numReceptorTypes, input->receptorTypes, input->rSmooth);

    // Precalculating the exponential function for receptor and ligand desolvation
    DesolvExpFunc desolvExpFunc(parameterLibrary.coeff_desolv);

    // Calculating bond vectors for directional H-bonds
    BondVectors *bondVectors = new BondVectors(&logFile);
    bondVectors->calculate(input, parameterLibrary);

    endTimer(1);

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
            logFile.printFormatted(" %6d   %8.3lf   %5.1lf%%   ", gridCoordZ, input->cgridmin[Z] + gridPosZ, ((z+1) * 100.0) / input->numGridPoints[Z]);
            logFile.printTimeInHMS((gridEndTime - gridStartTime) * (input->numGridPoints[Z] - z));
            logFile.print("  ");
            logFile.printExecutionTimes(gridStartTime, gridEndTime, &timesGridStart, &timerGridEnd);
        }
    */

    beginTimer(2);
    calculateElectrostaticMap(input, gridmaps, parameterLibrary);
    endTimer(2);

    // Calculation of gridmaps
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
