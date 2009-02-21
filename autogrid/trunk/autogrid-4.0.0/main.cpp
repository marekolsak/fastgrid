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
#include <omp.h>

//#define AG_OPENMP

void saveAVSGridmapsFile(const GridMapList &gridmaps, const InputData *input, const ProgramParameters &programParams, LogFile &logFile)
{
    FILE *fldFileAVS;
    if ((fldFileAVS = boincOpenFile(input->fldFilenameAVS, "w")) == 0)
    {
        logFile.printErrorFormatted(ERROR, "can't create grid dimensions data file %s\n", input->fldFilenameAVS);
        logFile.printError(FATAL_ERROR, "Unsuccessful completion.\n\n");
    }
    else
        logFile.printFormatted("\nCreating (AVS-readable) grid maps file : %s\n", input->fldFilenameAVS);

    int numMaps = gridmaps.getNumMapsInclFloatingGrid();
    fprintf(fldFileAVS, "# AVS field file\n#\n");
    fprintf(fldFileAVS, "# AutoDock Atomic Affinity and Electrostatic Grids\n#\n");
    fprintf(fldFileAVS, "# Created by %s.\n#\n", programParams.getProgramName());
    fprintf(fldFileAVS, "#SPACING %.3f\n", float(input->spacing));
    fprintf(fldFileAVS, "#NELEMENTS %d %d %d\n", input->nelements[X], input->nelements[Y], input->nelements[Z]);
    fprintf(fldFileAVS, "#CENTER %.3lf %.3lf %.3lf\n", input->center[X], input->center[Y], input->center[Z]);
    fprintf(fldFileAVS, "#MACROMOLECULE %s\n", input->receptorFilename);
    fprintf(fldFileAVS, "#GRID_PARAMETER_FILE %s\n#\n", programParams.getGridParameterFilename());
    fprintf(fldFileAVS, "ndim=3\t\t\t# number of dimensions in the field\n");
    fprintf(fldFileAVS, "dim1=%d\t\t\t# number of x-elements\n", input->numGridPoints[X]);
    fprintf(fldFileAVS, "dim2=%d\t\t\t# number of y-elements\n", input->numGridPoints[Y]);
    fprintf(fldFileAVS, "dim3=%d\t\t\t# number of z-elements\n", input->numGridPoints[Z]);
    fprintf(fldFileAVS, "nspace=3\t\t# number of physical coordinates per point\n");
    fprintf(fldFileAVS, "veclen=%d\t\t# number of affinity values at each point\n", numMaps);
    fprintf(fldFileAVS, "data=float\t\t# data type (byte, integer, float, double)\n");
    fprintf(fldFileAVS, "field=uniform\t\t# field type (uniform, rectilinear, irregular)\n");
    for (int i = 0; i < XYZ; i++)
        fprintf(fldFileAVS, "coord %d file=%s filetype=ascii offset=%d\n", (i + 1), input->xyzFilename, (i * 2));
    for (int i = 0; i < gridmaps.getNumAtomMaps(); i++)
        fprintf(fldFileAVS, "label=%s-affinity\t# component label for variable %d\n", gridmaps[i].type, (i + 1));
    fprintf(fldFileAVS, "label=Electrostatics\t# component label for variable %d\n", numMaps - 2);
    fprintf(fldFileAVS, "label=Desolvation\t# component label for variable %d\n", numMaps - 1);
    if (gridmaps.containsFloatingGrid())
        fprintf(fldFileAVS, "label=Floating_Grid\t# component label for variable %d\n", numMaps);
    fprintf(fldFileAVS, "#\n# location of affinity grid files and how to read them\n#\n");

    for (int i = 0; i < gridmaps.getNumMaps(); i++)
        fprintf(fldFileAVS, "variable %d file=%s filetype=ascii skip=6\n", (i + 1), gridmaps[i].filename);

    if (gridmaps.containsFloatingGrid())
        fprintf(fldFileAVS, "variable %d file=%s filetype=ascii skip=6\n", numMaps, input->floatingGridFilename);
    fclose(fldFileAVS);
}

void calculateGridmaps(const InputData *input, const GridMapList &gridmaps, const ParameterLibrary &parameterLibrary,
                       const PairwiseInteractionEnergies &energyLookup, const DesolvExpFunc &desolvExpFunc, const BondVectors *bondVectors,
                       LogFile &logFile)
{
    int hydrogen = parameterLibrary.getAtomParameterRecIndex("HD");
    const double logHalf = log(0.5);

    int infoCount = 0, infoCalc = 0;

    logFile.printFormatted("Beginning grid calculations.\n"
                           "\nCalculating %d grids over %d elements, around %d receptor atoms.\n\n"
                           "                    Percent   Estimated Time  Time/this plane\n"
                           "XY-plane  Z-coord   Done      Remaining       Real, User, System\n"
                           "            /Ang              /sec            /sec\n"
                           "________  ________  ________  ______________  __________________________\n\n",
                           gridmaps.getNumMapsInclFloatingGrid(), input->numGridPointsPerMap, input->numReceptorAtoms);

    // Iterate over all grid points, Z(Y (X)) (X is fastest)...
    beginTimer(1);

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

#if defined(AG_OPENMP)
    #pragma omp parallel for schedule(dynamic, 1)
#endif
    // Z axis
    for (int z = 0; z < input->numGridPoints[Z]; z++)
    {
        //  gridPos contains the current grid point.
        int gridCoordZ = z - input->ne[Z];
        double gridPosZ = gridCoordZ * input->spacing;
        int outputIndexZBase = z * input->numGridPoints[X] * input->numGridPoints[Y];

        // Y axis
        for (int y = 0; y < input->numGridPoints[Y]; y++)
        {
            int gridCoordY = y - input->ne[Y];
            double gridPosY = gridCoordY * input->spacing;
            int outputIndexZYBase = outputIndexZBase +  y * input->numGridPoints[X];

            // X axis
            for (int x = 0; x < input->numGridPoints[X]; x++)
            {
                int gridCoordX = x - input->ne[X];
                double gridPosX = gridCoordX * input->spacing;
                double gridPos[XYZ] = {gridPosX, gridPosY, gridPosZ};

                // Calculate outputIndex independently of other grid points. This will come in handy for parallelism.
                int outputIndex = outputIndexZYBase + x;

                for (int j = 0; j < gridmaps.getNumMaps(); j++)
                    if (gridmaps[j].isCovalent)
                    {
                        // Calculate the distance from the current grid point, c, to the covalent attachment point, input->covpos
                        double distance[XYZ];
                        for (int i = 0; i < XYZ; i++)
                            distance[i] = input->covpos[i] - gridPos[i];

                        // Distance squared from current grid point to the covalent attachment point
                        double rcovSq = hypotenuseSq(distance[X], distance[Y], distance[Z]);
                        rcovSq = rcovSq / (input->covHalfWidth * input->covHalfWidth);
                        if (rcovSq < APPROX_ZERO_SQ)
                            rcovSq = APPROX_ZERO_SQ;
                        gridmaps[j].energies[outputIndex] = input->covBarrier * (1 - exp(logHalf * rcovSq));
                    }
                    else // is not covalent
                        gridmaps[j].energies[outputIndex] = 0; // used to initialize to 'constant'for this gridmaps

                // Initialize Min Hbond variables for each new point
                struct HBondInfo
                {
                    double min[MAX_MAPS], max[MAX_MAPS];
                    bool flag[MAX_MAPS];
                } hbond;

                for (int i = 0; i < gridmaps.getNumAtomMaps(); i++)
                {
                    hbond.min[i] = 999999;
                    hbond.max[i] = -999999;
                    hbond.flag[i] = false;
                }

                // NEW2: Find Closest Hbond
                int closestH = 0;
                {
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
                }
                // END NEW2: Find Min Hbond

                double invRMin = 1.0 / BIG;

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

                    if (gridmaps.containsFloatingGrid())
                        // Calculate the so-called "Floating Grid"...
                        invRMin = max(invR, invRMin);

                    double elecEnergy = input->charge[ia] * min(invR, 2.0) * parameterLibrary.coeff_estat;
                    int indexR = -1;
                    //++infoCount;
                    if (input->distDepDiel)
                    {
                        // Distance-dependent dielectric...
                        // apply the estat forcefield coefficient/weight here
                        indexR = min(lookup(1 / invR), MAX_DIST-1); // make sure lookup index is in the table
                        gridmaps.getElectrostaticMap().energies[outputIndex] += elecEnergy * input->epsilon[indexR];
                        //++infoCalc;
                    }
                    else
                        // Constant dielectric...
                        gridmaps.getElectrostaticMap().energies[outputIndex] += elecEnergy * input->invDielCal;

                    // If distance from grid point to atom ia is too large,
                    // or if atom is a disordered hydrogen,
                    //   add nothing to the grid-point's non-bond energy;
                    //   just continue to next atom...

                    // Approximately 3% may pass through this condition
                    if (rSq > sq(NBCUTOFF))
                        continue;   // onto the next atom...
                    if (input->atomType[ia] == hydrogen && bondVectors->disorder[ia])
                        continue;   // onto the next atom...

                    if (indexR == -1)
                        indexR = min(lookup(1 / invR), MAX_DIST-1); // make sure lookup index is in the table

                    for (int i = 0; i < XYZ; i++)
                        distance[i] *= invR;

                    double racc = 1;
                    double rdon = 1;
                    double Hramp = 1;   // Hramp ramps in Hbond acceptor probes

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
                        double rd2 = sq(cross[0]) + sq(cross[1]) + sq(cross[2]);
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
                            //rdon = 562.25 * pow(0.116978 - sq(cosTheta), 3) * cos(t0);
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
                    for (int i = 0; i < gridmaps.getNumAtomMaps(); i++)
                    {
                        // We do not want to change the current enrg value for any covalent maps, make sure iscovalent is false...
                        if (!gridmaps[i].isCovalent)
                        {
                            double energy = energyLookup(input->atomType[ia], indexR, i);

                            if (gridmaps[i].isHBonder)
                            {
                                // PROBE forms H-bonds...

                                // rsph ramps in angular dependence for distances with negative energy
                                double rsph = clamp(energy / 100, 0.0, 1.0);

                                if ((gridmaps[i].hbond == AS || gridmaps[i].hbond == A2)    // AS or A2
                                    && (input->hbond[ia] == DS || input->hbond[ia] == D1))
                                {
                                    // DS or D1
                                    // PROBE can be an H-BOND ACCEPTOR,
                                    double f = Hramp * (racc + (1 - racc) * rsph);
                                    if (!bondVectors->disorder[ia])
                                        gridmaps[i].energies[outputIndex] += f * energy;
                                    else
                                        gridmaps[i].energies[outputIndex] += f * energyLookup(hydrogen, max(0, indexR - 110), i);
                                }
                                else if ((gridmaps[i].hbond == A1)    // A1
                                         && (input->hbond[ia] == DS || input->hbond[ia] == D1))
                                {
                                    // DS,D1
                                    double hbondEnergy = energy * (racc + (1 - racc) * rsph);
                                    hbond.min[i] = min(hbond.min[i], hbondEnergy);
                                    hbond.max[i] = max(hbond.max[i], hbondEnergy);
                                    hbond.flag[i] = true;
                                }
                                else if ((gridmaps[i].hbond == DS || gridmaps[i].hbond == D1) && (input->hbond[ia] > D1))
                                {
                                    // DS,D1 vs AS,A1,A2
                                    // PROBE is H-BOND DONOR,
                                    double hbondEnergy = energy * (rdon + (1 - rdon) * rsph);
                                    hbond.min[i] = min(hbond.min[i], hbondEnergy);
                                    hbond.max[i] = max(hbond.max[i], hbondEnergy);
                                    hbond.flag[i] = true;
                                }
                                else
                                    // hbonder PROBE-ia cannot form a H-bond...,
                                    gridmaps[i].energies[outputIndex] += energy;
                            }
                            else
                                // PROBE does not form H-bonds...,
                                gridmaps[i].energies[outputIndex] += energy;

                            // add desolvation energy
                            // forcefield desolv coefficient/weight in desolvExpFunc
                            gridmaps[i].energies[outputIndex] += gridmaps[i].solparProbe * input->vol[ia] * desolvExpFunc(indexR) +
                                (input->solpar[ia] + input->solparQ * fabs(input->charge[ia])) * gridmaps[i].volProbe * desolvExpFunc(indexR);
                        } // is not covalent
                    }
                    gridmaps.getDesolvationMap().energies[outputIndex] += input->solparQ * input->vol[ia] * desolvExpFunc(indexR);
                } // ia loop, over all receptor atoms...

                for (int i = 0; i < gridmaps.getNumAtomMaps(); i++)
                    if (hbond.flag[i])
                    {
                        gridmaps[i].energies[outputIndex] += hbond.min[i];
                        gridmaps[i].energies[outputIndex] += hbond.max[i];
                    }

                // O U T P U T . . .
                // Now output this grid point's energies to the maps:
                for (int k = 0; k < gridmaps.getNumMaps(); k++)
                {
                    // TODO: writing out energyMin/Max is currently broken, fix it by moving min/max calculations inside logSummary
                    //gridmaps[k].energyMax = max(gridmaps[k].energyMax, gridmaps[k].energies[outputIndex]);
                    //gridmaps[k].energyMin = min(gridmaps[k].energyMin, gridmaps[k].energies[outputIndex]);
                    if (fabs(gridmaps[k].energies[outputIndex]) < PRECISION)
                        gridmaps[k].energies[outputIndex] = 0;
                }
                if (gridmaps.containsFloatingGrid())
                    gridmaps.getFloatingGridMins()[outputIndex] = float(1 / invRMin);
            }
        }
    }
    endTimer(1);

    //fprintf(stderr, "Count: %i, Calc: %i = %i %%\n", infoCount, infoCalc, (infoCalc * 100) / infoCount);
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

    // Calculation of gridmaps
    calculateGridmaps(input, gridmaps, parameterLibrary, energyLookup, desolvExpFunc, bondVectors, logFile);

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
