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
#include <exception>

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

    int numMaps = gridmaps.getNumMaps() + (input->floatingGridFilename[0] != 0);
    fprintf(fldFileAVS, "# AVS field file\n#\n");
    fprintf(fldFileAVS, "# AutoDock Atomic Affinity and Electrostatic Grids\n#\n");
    fprintf(fldFileAVS, "# Created by %s.\n#\n", programParams.getProgramName());
    fprintf(fldFileAVS, "#SPACING %.3f\n", float(input->spacing));
    fprintf(fldFileAVS, "#NELEMENTS %d %d %d\n", input->nelements[X], input->nelements[Y], input->nelements[Z]);
    fprintf(fldFileAVS, "#CENTER %.3lf %.3lf %.3lf\n", input->center[X], input->center[Y], input->center[Z]);
    fprintf(fldFileAVS, "#MACROMOLECULE %s\n", input->receptorFilename);
    fprintf(fldFileAVS, "#GRID_PARAMETER_FILE %s\n#\n", programParams.getGridParameterFilename());
    fprintf(fldFileAVS, "ndim=3\t\t\t# number of dimensions in the field\n");
    fprintf(fldFileAVS, "dim1=%d\t\t\t# number of x-elements\n", input->n1[X]);
    fprintf(fldFileAVS, "dim2=%d\t\t\t# number of y-elements\n", input->n1[Y]);
    fprintf(fldFileAVS, "dim3=%d\t\t\t# number of z-elements\n", input->n1[Z]);
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
    if (input->floatingGridFilename[0])
        fprintf(fldFileAVS, "label=Floating_Grid\t# component label for variable %d\n", numMaps);
    fprintf(fldFileAVS, "#\n# location of affinity grid files and how to read them\n#\n");

    for (int i = 0; i < gridmaps.getNumMaps(); i++)
        fprintf(fldFileAVS, "variable %d file=%s filetype=ascii skip=6\n", (i + 1), gridmaps[i].filename);

    if (input->floatingGridFilename[0])
        fprintf(fldFileAVS, "variable %d file=%s filetype=ascii skip=6\n", numMaps, input->floatingGridFilename);
    fclose(fldFileAVS);
}

void calculateGridmaps(const InputData *input, GridMapList &gridmaps, const ParameterLibrary &parameterLibrary,
                       const PairwiseInteractionEnergies &energyLookup, const DesolvExpFunc &desolvExpFunc, const BondVectors *bondVectors,
                       LogFile &logFile)
{
    int hydrogen = parameterLibrary.getAtomParameterRecIndex("HD");
    char *maptypeptr;           // ptr for current map->type
    double cross[XYZ];
    double c[XYZ];
    int icoord[XYZ] = {0};
    int ctr = 0;
    double PI_halved = PI / 2;
    double rcov = 0.0;          // Distance from current grid point to the covalent attachment point
    double racc, cos_theta, r_min = 0, tmp, rdon, inv_r, inv_rmax, rsph, theta;
    double t0, ti;
    double ln_half = log(0.5);
    double temp_hbond_enrg, hbondmin[MAX_MAPS], hbondmax[MAX_MAPS];
    double rmin, Hramp;
    int nDone = 0;
    bool hbondflag[MAX_MAPS];
    int ii = 0;
    int ic = 0;
    int closestH = 0;
    Clock grd_start;
    Clock grd_end;
    tms tms_grd_start;
    tms tms_grd_end;
    bool warned = false;

    logFile.printFormatted("Beginning grid calculations.\n"
                           "\nCalculating %d grids over %d elements, around %d receptor atoms.\n\n"
                           "                    Percent   Estimated Time  Time/this plane\n"
                           "XY-plane  Z-coord   Done      Remaining       Real, User, System\n"
                           "            /Ang              /sec            /sec\n"
                           "________  ________  ________  ______________  __________________________\n\n",
                           gridmaps.getNumMaps() + (input->floatingGridFilename[0] != 0), input->numGridPointsPerMap, input->numReceptorAtoms);

    // Iterate over all grid points, Z(Y (X)) (X is fastest)...
    for (icoord[Z] = -input->ne[Z]; icoord[Z] <= input->ne[Z]; icoord[Z]++)
    {
        //  c[0:2] contains the current grid point.
        c[Z] = icoord[Z] * input->spacing;
        grd_start = times(&tms_grd_start);

        for (icoord[Y] = -input->ne[Y]; icoord[Y] <= input->ne[Y]; icoord[Y]++)
        {
            c[Y] = icoord[Y] * input->spacing;

            for (icoord[X] = -input->ne[X]; icoord[X] <= input->ne[X]; icoord[X]++)
            {
                c[X] = icoord[X] * input->spacing;

                for (int j = 0; j < gridmaps.getNumMaps(); j++)
                    if (gridmaps[j].isCovalent)
                    {
                        // Calculate the distance from the current grid point, c, to the covalent attachment point, input->covpos
                        double d[XYZ];
                        for (ii = 0; ii < XYZ; ii++)
                            d[ii] = input->covpos[ii] - c[ii];
                        rcov = hypotenuse(d[X], d[Y], d[Z]);
                        rcov = rcov / input->covHalfWidth;
                        if (rcov < APPROX_ZERO)
                            rcov = APPROX_ZERO;
                        gridmaps[j].energy = input->covBarrier * (1 - exp(ln_half * rcov * rcov));
                    }
                    else // is not covalent
                        gridmaps[j].energy = 0; // used to initialize to 'constant'for this gridmaps

                if (input->floatingGridFilename[0])
                    r_min = BIG;

                // Initialize Min Hbond variables for each new point
                for (int mapIndex = 0; mapIndex < gridmaps.getNumAtomMaps(); mapIndex++)
                {
                    hbondmin[mapIndex] = 999999;
                    hbondmax[mapIndex] = -999999;
                    hbondflag[mapIndex] = false;
                }

                // NEW2: Find Closest Hbond
                rmin = 999999;
                closestH = 0;
                for (int ia = 0; ia < input->numReceptorAtoms; ia++)
                    if ((input->hbond[ia] == 1) || (input->hbond[ia] == 2))
                    {
                        // DS or D1
                        double d[XYZ];
                        for (int i = 0; i < XYZ; i++)
                            d[i] = input->coord[ia][i] - c[i];
                        double r = hypotenuse(d[X], d[Y], d[Z]);
                        if (r < rmin)
                        {
                            rmin = r;
                            closestH = ia;
                        }
                    }           // Hydrogen test
                // END NEW2: Find Min Hbond

                //  Do all Receptor (protein, DNA, etc.) atoms...
                for (int ia = 0; ia < input->numReceptorAtoms; ia++)
                {
                    //  Get distance, r, from current grid point, c, to this receptor atom, input->coord,
                    double d[XYZ];
                    for (int i = 0; i < XYZ; i++)
                        d[i] = input->coord[ia][i] - c[i];
                    double r = hypotenuse(d[X], d[Y], d[Z]);
                    if (r < APPROX_ZERO)
                        r = APPROX_ZERO;
                    inv_r = 1 / r;
                    inv_rmax = 1 / max(r, 0.5);

                    for (int i = 0; i < XYZ; i++)
                        d[i] *= inv_r;
                    // make sure lookup index is in the table
                    int indexR = min(lookup(r), MAX_DIST-1);

                    if (input->floatingGridFilename[0])
                        // Calculate the so-called "Floating Grid"...
                        r_min = min(r, r_min);

                    // elecPE is the next-to-last last grid map, i.e. electrostatics
                    if (input->distDepDiel)
                        // Distance-dependent dielectric...
                        // gridmaps.getElectrostaticMap().energy += input->charge[ia] * inv_r * input->epsilon[indexR];
                        // apply the estat forcefield coefficient/weight here
                        gridmaps.getElectrostaticMap().energy += input->charge[ia] * inv_rmax * input->epsilon[indexR] * parameterLibrary.coeff_estat;
                    else
                        // Constant dielectric...
                        // gridmaps.getElectrostaticMap().energy += input->charge[ia] * inv_r * input->invDielCal;
                        gridmaps.getElectrostaticMap().energy += input->charge[ia] * inv_rmax * input->invDielCal * parameterLibrary.coeff_estat;

                    // If distance from grid point to atom ia is too large,
                    // or if atom is a disordered hydrogen,
                    //   add nothing to the grid-point's non-bond energy;
                    //   just continue to next atom...
                    if (r > NBCUTOFF)
                        continue;   // onto the next atom...
                    if (input->atomType[ia] == hydrogen && bondVectors->disorder[ia])
                        continue;   // onto the next atom...

                    // racc = rdon = 1;
                    racc = 1;
                    rdon = 1;
                    // NEW2 Hramp ramps in Hbond acceptor probes
                    Hramp = 1;
                    // END NEW2 Hramp ramps in Hbond acceptor probes

                    if (input->hbond[ia] == 2)
                    {           // D1
                        //  ia-th receptor atom = Hydrogen (4 = H)
                        //  => receptor H-bond donor, OH or NH.
                        //  calculate racc for H-bond ACCEPTOR PROBES at this grid pt.
                        cos_theta = 0;
                        //  d[] = Unit vector from current grid pt to ia_th m/m atom.
                        //  cos_theta = d dot bondVectors->rvector == cos(angle) subtended.
                        for (int i = 0; i < XYZ; i++)
                            cos_theta -= d[i] * bondVectors->rvector[ia][i];

                        if (cos_theta <= 0)
                            //  H->current-grid-pt vector >= 90 degrees from
                            //  N->H or O->H vector,
                            racc = 0;
                        else
                        {
                            //  racc = [cos(theta)]^2.0 for N-H
                            //  racc = [cos(theta)]^4.0 for O-H,
                            switch (bondVectors->rexp[ia])
                            {
                            case 1:
                            default:
                                racc = cos_theta;
                                break;
                            case 2:
                                racc = cos_theta * cos_theta;
                                break;
                            case 4:
                                tmp = cos_theta * cos_theta;
                                racc = tmp * tmp;
                                break;
                            }
                            // racc = pow(cos_theta, bondVectors->rexp[ia]);

                            // NEW2 calculate dot product of bond vector with bond vector of best input->hbond
                            if (ia == closestH)
                                Hramp = 1;
                            else
                            {
                                cos_theta = 0;
                                for (int i = 0; i < XYZ; i++)
                                    cos_theta += bondVectors->rvector[closestH][i] * bondVectors->rvector[ia][i];
                                cos_theta = min(cos_theta, 1.0);
                                cos_theta = max(cos_theta, -1.0);
                                theta = acos(cos_theta);
                                Hramp = 0.5 - 0.5 * cos(theta * 120 / 90);
                            }   // ia test
                            // END NEW2 calculate dot product of bond vector with bond vector of best input->hbond
                        }
                        // endif (input->atomType[ia] == hydrogen)
                        // NEW Directional N acceptor
                    }
                    else if (input->hbond[ia] == 4)
                    {           // A1
                        //  ia-th macromolecule atom = Nitrogen (4 = H)
                        //  calculate rdon for H-bond Donor PROBES at this grid pt.
                        cos_theta = 0;
                        //  d[] = Unit vector from current grid pt to ia_th m/m atom.
                        //  cos_theta = d dot bondVectors->rvector == cos(angle) subtended.
                        for (int i = 0; i < XYZ; i++)
                            cos_theta -= d[i] * bondVectors->rvector[ia][i];

                        if (cos_theta <= 0)
                            //  H->current-grid-pt vector >= 90 degrees from
                            //  X->N vector,
                            rdon = 0;
                        else
                            //  racc = [cos(theta)]^2.0 for H->N
                            rdon = cos_theta * cos_theta;
                        // endif (input->atomType[ia] == nitrogen)
                        // end NEW Directional N acceptor

                    }
                    else if (input->hbond[ia] == 5 && !bondVectors->disorder[ia])
                    {           // A2
                        //  ia-th receptor atom = Oxygen
                        //  => receptor H-bond acceptor, oxygen.
                        rdon = 0;

                        // check to see that probe is in front of oxygen, not behind
                        cos_theta = 0;
                        for (int i = 0; i < XYZ; i++)
                            cos_theta -= d[i] * bondVectors->rvector[ia][i];
                        // t0 is the angle out of the lone pair plane, calculated
                        // as 90 deg - acos (vector to grid point DOT lone pair
                        // plane normal)
                        t0 = 0;
                        for (int i = 0; i < XYZ; i++)
                            t0 += d[i] * bondVectors->rvector2[ia][i];
                        if (t0 > 1)
                        {
                            logFile.printErrorFormatted(WARNING, "I just prevented an attempt to take the arccosine of %f, a value greater than 1.\n", t0);
                            t0 = 1;
                        }
                        else if (t0 < -1)
                        {
                            logFile.printErrorFormatted(WARNING, "I just prevented an attempt to take the arccosine of %f, a value less than -1.\n", t0);
                            t0 = -1;
                        }
                        t0 = PI_halved - acos(t0);

                        // ti is the angle in the lone pair plane, away from the
                        // vector between the lone pairs,
                        // calculated as (grid vector CROSS lone pair plane normal)
                        // DOT C=O vector - 90 deg
                        cross[0] = d[1] * bondVectors->rvector2[ia][2] - d[2] * bondVectors->rvector2[ia][1];
                        cross[1] = d[2] * bondVectors->rvector2[ia][0] - d[0] * bondVectors->rvector2[ia][2];
                        cross[2] = d[0] * bondVectors->rvector2[ia][1] - d[1] * bondVectors->rvector2[ia][0];
                        double rd2 = sq(cross[0]) + sq(cross[1]) + sq(cross[2]);
                        if (rd2 < APPROX_ZERO)
                        {
                            if ((rd2 == 0) && !warned)
                            {
                                logFile.printError(WARNING, "Attempt to divide by zero was just prevented.\n\n");
                                warned = true;
                            }
                            rd2 = APPROX_ZERO;
                        }
                        double inv_rd = 1 / sqrt(rd2);
                        ti = 0;
                        for (int i = 0; i < XYZ; i++)
                            ti += cross[i] * inv_rd * bondVectors->rvector[ia][i];

                        // rdon expressions from Goodford
                        rdon = 0;
                        if (cos_theta >= 0)
                        {
                            if (ti > 1)
                            {
                                logFile.printErrorFormatted(WARNING, "I just prevented an attempt to take the arccosine of %f, a value greater than 1.\n", ti);
                                ti = 1;
                            }
                            else if (ti < -1)
                            {
                                logFile.printErrorFormatted(WARNING, "I just prevented an attempt to take the arccosine of %f, a value less than -1.\n", ti);
                                ti = -1;
                            }
                            ti = acos(ti) - PI_halved;
                            if (ti < 0)
                                ti = -ti;
                            // the 2.0*ti can be replaced by (ti + ti) in: rdon = (0.9 + 0.1*sin(2.0*ti))*cos(t0);
                            rdon = (0.9 + 0.1 * sin(ti + ti)) * cos(t0);
                            // 0.34202 = cos (100 deg)
                        }
                        else if (cos_theta >= -0.34202)
                            rdon = 562.25 * pow(0.116978 - sq(cos_theta), 3) * cos(t0);

                        // endif input->atomType == oxygen, not disordered
                    }
                    else if (input->hbond[ia] == 5 && bondVectors->disorder[ia])
                    {
                        // A2
                        // cylindrically disordered hydroxyl
                        cos_theta = 0;
                        for (int i = 0; i < XYZ; i++)
                            cos_theta -= d[i] * bondVectors->rvector[ia][i];
                        if (cos_theta > 1)
                        {
                            logFile.printErrorFormatted(WARNING, "I just prevented an attempt to take the arccosine of %f, a value greater than 1.\n", cos_theta);
                            cos_theta = 1;
                        }
                        else if (cos_theta < -1)
                        {
                            logFile.printErrorFormatted(WARNING, "I just prevented an attempt to take the arccosine of %f, a value less than -1.\n", cos_theta);
                            cos_theta = -1;
                        }
                        theta = acos(cos_theta);
                        racc = 0;
                        rdon = 0;
                        if (theta <= 1.24791 + PI_halved)
                        {
                            // 1.24791 rad = 180 deg minus C-O-H bond angle, ** 108.5 deg
                            rdon = pow(cos(theta - 1.24791), 4);
                            racc = rdon;
                        }
                    }           // input->atomType test

                    // For each probe atom-type,
                    // Sum pairwise interactions between each probe
                    // at this grid point (c[0:2])
                    // and the current receptor atom, ia...
                    for (int mapIndex = 0; mapIndex < gridmaps.getNumAtomMaps(); mapIndex++)
                    {
                        // We do not want to change the current enrg value for any covalent maps, make sure iscovalent is false...
                        maptypeptr = gridmaps[mapIndex].type;

                        if (!gridmaps[mapIndex].isCovalent)
                        {
                            if (gridmaps[mapIndex].isHBonder)
                            {
                                // PROBE forms H-bonds...

                                // rsph ramps in angular dependence for distances with negative energy
                                rsph = energyLookup(input->atomType[ia], indexR, mapIndex) / 100;
                                rsph = max(rsph, 0);
                                rsph = min(rsph, 1);
                                if ((gridmaps[mapIndex].hbond == 3 || gridmaps[mapIndex].hbond == 5)    // AS or A2
                                    && (input->hbond[ia] == 1 || input->hbond[ia] == 2))
                                {   // DS or D1
                                    // PROBE can be an H-BOND ACCEPTOR,
                                    if (!bondVectors->disorder[ia])
                                        gridmaps[mapIndex].energy += energyLookup(input->atomType[ia], indexR, mapIndex) * Hramp * (racc + (1 - racc) * rsph);
                                    else
                                        gridmaps[mapIndex].energy += energyLookup(hydrogen, max(0, indexR - 110), mapIndex) * Hramp * (racc + (1 - racc) * rsph);
                                }
                                else if ((gridmaps[mapIndex].hbond == 4)    // A1
                                         && (input->hbond[ia] == 1 || input->hbond[ia] == 2))
                                {   // DS,D1
                                    hbondmin[mapIndex] = min(hbondmin[mapIndex], energyLookup(input->atomType[ia], indexR, mapIndex) * (racc + (1 - racc) * rsph));
                                    hbondmax[mapIndex] = max(hbondmax[mapIndex], energyLookup(input->atomType[ia], indexR, mapIndex) * (racc + (1 - racc) * rsph));
                                    hbondflag[mapIndex] = true;
                                }
                                else if ((gridmaps[mapIndex].hbond == 1 || gridmaps[mapIndex].hbond == 2) && (input->hbond[ia] > 2))
                                {   // DS,D1 vs AS,A1,A2
                                    // PROBE is H-BOND DONOR,
                                    temp_hbond_enrg = energyLookup(input->atomType[ia], indexR, mapIndex) * (rdon + (1 - rdon) * rsph);
                                    hbondmin[mapIndex] = min(hbondmin[mapIndex], temp_hbond_enrg);
                                    hbondmax[mapIndex] = max(hbondmax[mapIndex], temp_hbond_enrg);
                                    hbondflag[mapIndex] = true;
                                }
                                else
                                    // hbonder PROBE-ia cannot form a H-bond...,
                                    gridmaps[mapIndex].energy += energyLookup(input->atomType[ia], indexR, mapIndex);
                            }
                            else
                                // PROBE does not form H-bonds...,
                                gridmaps[mapIndex].energy += energyLookup(input->atomType[ia], indexR, mapIndex);

                            // add desolvation energy
                            // forcefield desolv coefficient/weight in desolvExpFunc
                            gridmaps[mapIndex].energy += gridmaps[mapIndex].solparProbe * input->vol[ia] * desolvExpFunc(indexR) +
                                (input->solpar[ia] + input->solparQ * fabs(input->charge[ia])) * gridmaps[mapIndex].volProbe * desolvExpFunc(indexR);
                        }       // is not covalent
                    }           // mapIndex
                    gridmaps.getDesolvationMap().energy += input->solparQ * input->vol[ia] * desolvExpFunc(indexR);
                }               // ia loop, over all receptor atoms...
                for (int mapIndex = 0; mapIndex < gridmaps.getNumAtomMaps(); mapIndex++)
                    if (hbondflag[mapIndex])
                    {
                        gridmaps[mapIndex].energy += hbondmin[mapIndex];
                        gridmaps[mapIndex].energy += hbondmax[mapIndex];
                    }

                // O U T P U T . . .
                // Now output this grid point's energies to the maps:
                // 2 includes new dsolvPE
                for (int k = 0; k < gridmaps.getNumMaps(); k++)
                {
                    gridmaps[k].energyMax = max(gridmaps[k].energyMax, gridmaps[k].energy);
                    gridmaps[k].energyMin = min(gridmaps[k].energyMin, gridmaps[k].energy);

                    float f = float(round3dp(gridmaps[k].energy));
                    fprintf(gridmaps[k].file, "%.3f\n", f);

                }
                if (gridmaps.getFloatingGridFile())
                    fprintf(gridmaps.getFloatingGridFile(), "%.3f\n", float(round3dp(r_min)));
                ctr++;
            }                   // icoord[X] loop
        }                       // icoord[Y] loop

        grd_end = times(&tms_grd_end);
        ++nDone;
        logFile.printFormatted(" %6d   %8.3lf   %5.1lf%%   ", icoord[Z], input->cgridmin[Z] + c[Z], (100 / double(input->n1[Z])) * double(++ic));
        logFile.printTimeInHMS((grd_end - grd_start) * (input->n1[Z] - nDone));
        logFile.print("  ");
        logFile.printExecutionTimes(grd_start, grd_end, &tms_grd_start, &tms_grd_end);
    }                           // icoord[Z] loop
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

    // TODO: consider adding floating grid file handling into the GridMapList class
    gridmaps.setFloatingGridFilename(input->floatingGridFilename);

    // Preparing the files
    gridmaps.prepareFiles(input, programParams.getGridParameterFilename());

    // TODO: add a smarter mechanism of checking for the available disk space, we need to know it as soon as possible.
    // the formerly implemented checks in the middle of calculations were too late

    // Loading the parameter library from the file
    if (input->parameterLibraryFilename[0])
        parameterLibrary.load(input->parameterLibraryFilename);

    // Writing to AVS-readable gridmaps file (fld)
    saveAVSGridmapsFile(gridmaps, input, programParams, logFile);

#if defined(BOINCCOMPOUND)
    boinc_fraction_done(0.1);
#endif

    // Calculating the lookup table of the pairwise interaction energies
    PairwiseInteractionEnergies energyLookup;
    energyLookup.calculate(gridmaps, logFile, input->numReceptorTypes, input->receptorTypes, input->rSmooth);

    // Precalculating the exponential function for receptor and ligand desolvation
    DesolvExpFunc desolvExpFunc(parameterLibrary.coeff_desolv);

    // Calculating bond vectors for directional H-bonds
    BondVectors *bondVectors = new BondVectors(&logFile);
    bondVectors->calculate(input, parameterLibrary);

    // Calculation of gridmaps
    calculateGridmaps(input, gridmaps, parameterLibrary, energyLookup, desolvExpFunc, bondVectors, logFile);

    delete bondVectors;
    delete inputDataLoader;

#if defined(BOINCCOMPOUND)
    boinc_fraction_done(0.9);
#endif

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
