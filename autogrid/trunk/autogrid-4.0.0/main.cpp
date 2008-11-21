/*

   $Id: main.cpp,v 1.58 2007/05/04 07:54:25 garrett Exp $

   AutoGrid

   Copyright (C) 1989-2007, Garrett M. Morris, David S. Goodsell, Ruth Huey, Arthur J. Olson, All Rights Reserved.

   AutoGrid is a Trade Mark of The Scripps Research Institute.

   This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See
   the GNU General Public License for more details.

   You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
   MA 02110-1301, USA.
*/

#pragma region Includes

#include <cmath>
#include <cassert>
#include <search.h>
#include <cstring>
#include <cstdlib>
#include <cctype>   // tolower
#include <cstddef>
#include <iostream>

#if defined(_WIN32)
    #define WIN32_LEAN_AND_MEAN
    #include <windows.h>
#else
    #include <sys/param.h>
#endif

// the BOINC API header file
#if defined(BOINC)
    #include "diagnostics.h"
    #include "boinc_api.h"
    #include "filesys.h"    // boinc_fopen(), etc...
#endif

#include "autogrid.h"
#include "utils.h"
#include "programparameters.h"
#include "read_parameter_library.h"

// classes
#include "exceptions.h" // ExitProgram
#include "logfile.h"    // LogFile
#include "gridmap.h"    // GridMap, GridMapList
#include "pairwiseinteractionenergies.h" // PairwiseInteractionEnergies
#include "inputdata.h"  // InputData

#pragma endregion Includes

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
void appMain(int argc, char **argv)
{
#if defined(_WIN32)
    SetThreadAffinityMask(GetCurrentThread(), 1);
#endif

    // Get the time at the start of the run...
    tms tmsJobStart;
    Clock jobStart = times(&tmsJobStart);

    double versionNumber = 4.00;

    // Initialize BOINC if needed
    initBoinc();

    // Initialize the ProgramParameters object, which parses the command-line arguments
    ProgramParameters programParams(argc, argv);

    // Initialize the log file
    LogFile logFile(versionNumber, programParams.getProgramName(), programParams.getLogFilename());

    // Declaration of gridmaps
    GridMapList gridmaps(&logFile);

    // TODO: put atomParameterManager-related functions into a new class
    // TODO: put LinearFreeEnergyModel, setupParameterLibrary and readParameterLibrary into a new class

    // Initialization of free energy coefficients
    LinearFreeEnergyModel model;

    // Set default values
    // atomParameterManager is also used inside
    setupParameterLibrary(-1, programParams.getProgramName(), programParams.getDebugLevel(), logFile, model);

    // Reading in the grid parameter file
    InputData *input = new InputData();
    input->load(programParams, gridmaps, logFile);  // TODO: shouldn't we put the gridmaps initialization code out of the load function?

    // Load values from the file
    if (input->parameterLibraryFilename[0])
        readParameterLibrary(input->parameterLibraryFilename, -1, programParams.getProgramName(), programParams.getDebugLevel(), logFile, model);

#pragma region Writing to AVS_fld file
{
    int numMaps = gridmaps.getNumMaps() + (input->floatingGridFilename[0] != 0);
    FILE *fldFileAVS;
    if ((fldFileAVS = openFile(input->fldFilenameAVS, "w")) == 0)
    {
        logFile.printErrorFormatted(ERROR, "can't create grid dimensions data file %s\n", input->fldFilenameAVS);
        logFile.printError(FATAL_ERROR, "Unsuccessful completion.\n\n");
    }
    else
        logFile.printFormatted("\nCreating (AVS-readable) grid maps file : %s\n", input->fldFilenameAVS);

    fprintf(fldFileAVS, "# AVS field file\n#\n");
    fprintf(fldFileAVS, "# AutoDock Atomic Affinity and Electrostatic Grids\n#\n");
    fprintf(fldFileAVS, "# Created by %s.\n#\n", programParams.getProgramName());
    fprintf(fldFileAVS, "#SPACING %.3f\n", (float)input->spacing);
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
        fprintf(fldFileAVS, "input->coord %d file=%s filetype=ascii offset=%d\n", (i + 1), input->xyzFilename, (i * 2));
    for (int i = 0; i < gridmaps.getNumAtomMaps(); i++)
        fprintf(fldFileAVS, "label=%s-affinity\t# component label for variable %d\n", gridmaps[i].type, (i + 1));                           // i
    fprintf(fldFileAVS, "label=Electrostatics\t# component label for variable %d\n", numMaps - 2);
    fprintf(fldFileAVS, "label=Desolvation\t# component label for variable %d\n", numMaps - 1);
    if (input->floatingGridFilename[0])
        fprintf(fldFileAVS, "label=Floating_Grid\t# component label for variable %d\n", numMaps);
    fprintf(fldFileAVS, "#\n# location of affinity grid files and how to read them\n#\n");

    for (int i = 0; i < gridmaps.getNumMaps(); i++)
        fprintf(fldFileAVS, "variable %d file=%s filetype=ascii skip=6\n", (i + 1), gridmaps[i].mapFilename);

    if (input->floatingGridFilename[0])
        fprintf(fldFileAVS, "variable %d file=%s filetype=ascii skip=6\n", numMaps, input->floatingGridFilename);
    fclose(fldFileAVS);
}
#pragma endregion Writing to AVS_fld file

#if defined(BOINCCOMPOUND)
    boinc_fraction_done(0.1);
#endif

    PairwiseInteractionEnergies *energyLookup = new PairwiseInteractionEnergies();
    energyLookup->calculate(gridmaps, logFile, input->numReceptorTypes, input->receptorTypes, input->rSmooth);

    double sol_fn[MAX_DIST];

#pragma region Precalculating of the exponential function for receptor and ligand desolvation
{
    double minus_inv_two_sigma_sqd;
    double sigma;

    // exponential function for receptor and ligand desolvation
    // note: the solvation term will not be smoothed
    sigma = 3.6;
    minus_inv_two_sigma_sqd = -1 / (2 * sigma * sigma);
    for (int indx_r = 1; indx_r < MAX_DIST; indx_r++)
    {
        double r = angstrom(indx_r);
        // sol_fn[indx_r] = exp(-sq(r)/(2*sigma*sigma));
        sol_fn[indx_r] = exp(sq(r) * minus_inv_two_sigma_sqd);
        sol_fn[indx_r] *= model.coeff_desolv;
    }
}
#pragma endregion Precalculating of the exponential function for receptor and ligand desolvation

    // canned atom type number
    int hydrogen, carbon, arom_carbon, oxygen, nitrogen;
    int nonHB_hydrogen, nonHB_nitrogen, sulphur, nonHB_sulphur;

    int disorder[AG_MAX_ATOMS];
    int rexp[AG_MAX_ATOMS] = {0};
    double rvector[AG_MAX_ATOMS][XYZ];
    double rvector2[AG_MAX_ATOMS][XYZ];
    double d[XYZ];

    char warned = 'F';

    double inv_rd, rd2; // temporary??

#pragma region Calculating bond vectors for directional H-bonds
{
    double dc[XYZ];
    double rdot;
    int from, to;
    int nbond;
    int i1 = 0, i2 = 0, i3 = 0;

    // Loop over all RECEPTOR atoms to
    // calculate bond vectors for directional H-bonds
    // setup the canned atom types here....
    // at this point set up hydrogen, carbon, oxygen and nitrogen
    hydrogen = getRecIndex("HD");
    nonHB_hydrogen = getRecIndex("H");
    carbon = getRecIndex("C");
    arom_carbon = getRecIndex("A");
    oxygen = getRecIndex("OA");
    nitrogen = getRecIndex("NA");
    nonHB_nitrogen = getRecIndex("N");
    sulphur = getRecIndex("SA");
    nonHB_sulphur = getRecIndex("S");

    // 7:CHANGE HERE: scan the 'mapIndex' from the input
    for (int ia = 0; ia < input->numReceptorAtoms; ia++)
    {                                      //** ia = i_receptor_atom_a **
        disorder[ia] = false;   // initialize disorder flag.
        warned = 'F';

        // Set scan limits looking for bonded atoms
        from = max(ia - 20, 0);
        to = min(ia + 20, input->numReceptorAtoms - 1);

        // If 'ia' is a hydrogen atom, it could be a
        // RECEPTOR hydrogen-BOND DONOR,
        // 8:CHANGE HERE: fix the input->atomType vs atom_types problem in following
        if ((int)input->hbond[ia] == 2) // D1 hydrogen bond donor
        {
            for (int ib = from; ib <= to; ib++)
                if (ib != ia) // ib = i_receptor_atom_b
                {
                    // =>  NH-> or OH->
                    // if ((input->atomType[ib] == nitrogen) || (input->atomType[ib]==nonHB_nitrogen) ||(input->atomType[ib] == oxygen)||(input->atomType[ib] == sulphur)||(input->atomType[ib]==nonHB_sulphur)) {

                    // Calculate the square of the N-H or O-H bond distance, rd2,
                    //                            ib-ia  ib-ia
                    for (int i = 0; i < XYZ; i++)
                        d[i] = input->coord[ia][i] - input->coord[ib][i];
                    rd2 = sq(d[X]) + sq(d[Y]) + sq(d[Z]);
                    // If ia & ib are less than 1.3 A apart -- they are covalently bonded,
                    if (rd2 < 1.90)
                    {           // INCREASED for H-S bonds
                        if (rd2 < APPROX_ZERO)
                        {
                            if (rd2 == 0)
                                logFile.printErrorFormatted(WARNING,
                                    "While calculating an H-O or H-N bond vector...\nAttempt to divide by zero was just prevented.\nAre the coordinates of atoms %d and %d the same?\n\n",
                                    ia + 1, ib + 1);
                            rd2 = APPROX_ZERO;
                        }
                        inv_rd = 1 / sqrt(rd2);

                        // N-H: Set exponent rexp to 2 for m/m H-atom,
                        // if (input->atomType[ib] == nitrogen) rexp[ia] = 2;
                        if ((input->atomType[ib] != oxygen) && (input->atomType[ib] != sulphur))
                            rexp[ia] = 2;

                        // O-H: Set exponent rexp to 4 for m/m H-atom,
                        // and flag disordered hydroxyls
                        if ((input->atomType[ib] == oxygen) || (input->atomType[ib] == sulphur))
                        {
                            rexp[ia] = 4;
                            if (input->disorderH)
                                disorder[ia] = TRUE;
                        }

                        // Normalize the vector from ib to ia, N->H or O->H...
                        for (int i = 0; i < XYZ; i++)
                            rvector[ia][i] = d[i] * inv_rd;

                        // First O-H/N-H H-bond-donor found; Go on to next atom,
                        break;
                    }           // Found covalent bond.
                    // } Found NH or OH in receptor.
                }
            // Finished scanning for the NH or OH in receptor.
            // If 'ia' is an Oxygen atom, it could be a
            // RECEPTOR H_BOND ACCEPTOR,
        }
        else if (input->hbond[ia] == 5)
        {                       // A2
            // Scan from at most, (ia-20)th m/m atom, or ia-th (if ia<20)
            //        to (ia + 5)th m/m-atom
            // determine number of atoms bonded to the oxygen
            nbond = 0;
            int ib = from;
            for (; ib <= to; ib++)
                if (ib != ia)
                {
                    rd2 = 0;

                    for (int i = 0; i < XYZ; i++)
                    {
                        dc[i] = input->coord[ia][i] - input->coord[ib][i];
                        rd2 += sq(dc[i]);
                    }

                    if (((rd2 < 3.61) && ((input->atomType[ib] != hydrogen) && (input->atomType[ib] != nonHB_hydrogen))) ||
                        ((rd2 < 1.69) && ((input->atomType[ib] == hydrogen) || (input->atomType[ib] == nonHB_hydrogen))))
                    {
                        if (nbond == 2)
                            logFile.printErrorFormatted(WARNING, "Found an H-bonding atom with three bonded atoms, atom serial number %d\n", ia + 1);
                        if (nbond == 1)
                        {
                            nbond = 2;
                            i2 = ib;
                        }
                        if (nbond == 0)
                        {
                            nbond = 1;
                            i1 = ib;
                        }
                    }
                }               // (ib != ia)

            // if no bonds, something is wrong
            if (nbond == 0)
                logFile.printErrorFormatted(WARNING, "Oxygen atom found with no bonded atoms, atom serial number %d, input->atomType %d\n", ia + 1, input->atomType[ia]);

            // one bond: Carbonyl Oxygen O=C-X
            if (nbond == 1)
            {
                // calculate normalized carbonyl bond vector rvector[ia][]
                rd2 = 0;
                for (int i = 0; i < XYZ; i++)
                {
                    rvector[ia][i] = input->coord[ia][i] - input->coord[i1][i];
                    rd2 += sq(rvector[ia][i]);
                }
                if (rd2 < APPROX_ZERO)
                {
                    if ((rd2 == 0) && (warned == 'F'))
                    {
                        logFile.printErrorFormatted(WARNING, "Attempt to divide by zero was just prevented.\nAre the coordinates of atoms %d and %d the same?\n\n", ia + 1, i1 + 1);
                        warned = 'T';
                    }
                    rd2 = APPROX_ZERO;
                }
                inv_rd = 1 / sqrt(rd2);
                for (int i = 0; i < XYZ; i++)
                    rvector[ia][i] *= inv_rd;

                // find a second atom (i2) bonded to carbonyl carbon (i1)
                for (int i2 = from; i2 <= to; i2++)
                    if ((i2 != i1) && (i2 != ia))
                    {
                        rd2 = 0;
                        for (int i = 0; i < XYZ; i++)
                        {
                            dc[i] = input->coord[i1][i] - input->coord[i2][i];
                            /*NEW*/ rd2 += sq(dc[i]);
                        }
                        if (((rd2 < 2.89) && (input->atomType[i2] != hydrogen)) || ((rd2 < 1.69) && (input->atomType[i2] == hydrogen)))
                        {

                            // found one
                            // d[i] vector from carbon to second atom
                            rd2 = 0;
                            for (int i = 0; i < XYZ; i++)
                            {
                                d[i] = input->coord[i2][i] - input->coord[i1][i];
                                rd2 += sq(d[i]);
                            }
                            if (rd2 < APPROX_ZERO)
                            {
                                if ((rd2 == 0) && (warned == 'F'))
                                {
                                    logFile.printErrorFormatted(WARNING, "Attempt to divide by zero was just prevented.\nAre the coordinates of atoms %d and %d the same?\n\n", i1 + 1, i2 + 1);
                                    warned = 'T';
                                }
                                rd2 = APPROX_ZERO;
                            }
                            inv_rd = 1 / sqrt(rd2);
                            for (int i = 0; i < XYZ; i++)
                                d[i] *= inv_rd;

                            // C=O cross C-X gives the lone pair plane normal
                            rvector2[ia][0] = rvector[ia][1] * d[2] - rvector[ia][2] * d[1];
                            rvector2[ia][1] = rvector[ia][2] * d[0] - rvector[ia][0] * d[2];
                            rvector2[ia][2] = rvector[ia][0] * d[1] - rvector[ia][1] * d[0];
                            rd2 = 0;
                            for (int i = 0; i < XYZ; i++)
                                rd2 += sq(rvector2[ia][i]);
                            if (rd2 < APPROX_ZERO)
                            {
                                if ((rd2 == 0) && (warned == 'F'))
                                {
                                    logFile.printError(WARNING, "Attempt to divide by zero was just prevented.\n\n");
                                    warned = 'T';
                                }
                                rd2 = APPROX_ZERO;
                            }
                            inv_rd = 1 / sqrt(rd2);
                            for (int i = 0; i < XYZ; i++)
                                rvector2[ia][i] *= inv_rd;
                        }
                    }
            }                   // endif nbond==1

            // two bonds: Hydroxyl or Ether Oxygen X1-O-X2
            if (nbond == 2)
                // disordered hydroxyl
                if ((input->atomType[i1] == hydrogen || input->atomType[i2] == hydrogen) && input->atomType[i1] != input->atomType[i2] && input->disorderH)
                {
                    if ((input->atomType[i1] == carbon) || (input->atomType[i1] == arom_carbon))
                        ib = i1;
                    if ((input->atomType[i2] == carbon) || (input->atomType[i1] == arom_carbon))
                        ib = i2;
                    disorder[ia] = TRUE;
                    rd2 = 0;
                    for (int i = 0; i < XYZ; i++)
                    {
                        rvector[ia][i] = input->coord[ia][i] - input->coord[ib][i];
                        rd2 += sq(rvector[ia][i]);
                    }
                    if (rd2 < APPROX_ZERO)
                    {
                        if ((rd2 == 0) && (warned == 'F'))
                        {
                            logFile.printErrorFormatted(WARNING, "Attempt to divide by zero was just prevented.\nAre the coordinates of atoms %d and %d the same?\n\n", ia + 1, ib + 1);
                            warned = 'T';
                        }
                        rd2 = APPROX_ZERO;
                    }
                    inv_rd = 1 / sqrt(rd2);
                    for (int i = 0; i < XYZ; i++)
                        rvector[ia][i] *= inv_rd;
                }
                else
                {
                    // not a disordered hydroxyl
                    // normalized X1 to X2 vector, defines lone pair plane
                    rd2 = 0;
                    for (int i = 0; i < XYZ; i++)
                    {
                        rvector2[ia][i] = input->coord[i2][i] - input->coord[i1][i];
                        rd2 += sq(rvector2[ia][i]);
                    }
                    if (rd2 < APPROX_ZERO)
                    {
                        if ((rd2 == 0) && (warned == 'F'))
                        {
                            logFile.printErrorFormatted(WARNING, "Attempt to divide by zero was just prevented.\nAre the coordinates of atoms %d and %d the same?\n\n", i1 + 1, i2 + 1);
                            warned = 'T';
                        }
                        rd2 = APPROX_ZERO;
                    }
                    inv_rd = 1 / sqrt(rd2);
                    for (int i = 0; i < XYZ; i++)
                        rvector2[ia][i] *= inv_rd;

                    // vector pointing between the lone pairs:
                    // front of the vector is the oxygen atom,
                    // X1->O vector dotted with normalized X1->X2 vector plus
                    // coords of X1 gives the point on the X1-X2 line for the
                    // back of the vector.
                    rdot = 0;
                    for (int i = 0; i < XYZ; i++)
                        rdot += (input->coord[ia][i] - input->coord[i1][i]) * rvector2[ia][i];
                    rd2 = 0;
                    for (int i = 0; i < XYZ; i++)
                    {
                        rvector[ia][i] = input->coord[ia][i] - ((rdot * rvector2[ia][i]) + input->coord[i1][i]);
                        rd2 += sq(rvector[ia][i]);
                    }
                    if (rd2 < APPROX_ZERO)
                    {
                        if ((rd2 == 0) && (warned == 'F'))
                        {
                            logFile.printError(WARNING, "Attempt to divide by zero was just prevented.\n\n");
                            warned = 'T';
                        }
                        rd2 = APPROX_ZERO;
                    }
                    inv_rd = 1 / sqrt(rd2);
                    for (int i = 0; i < XYZ; i++)
                        rvector[ia][i] *= inv_rd;
                }               // end disordered hydroxyl
        }
        else if (input->hbond[ia] == 4)
        {                       // A1
            // Scan from at most, (ia-20)th m/m atom, or ia-th (if ia<20)
            //        to (ia+5)th m/m-atom
            // determine number of atoms bonded to the oxygen
            nbond = 0;
            int ib = from;
            for (; ib <= to; ib++)
                if (ib != ia)
                {
                    rd2 = 0;
                    for (int i = 0; i < XYZ; i++)
                    {
                        dc[i] = input->coord[ia][i] - input->coord[ib][i];
                        rd2 += sq(dc[i]);
                    }

                    if (((rd2 < 2.89) && ((input->atomType[ib] != hydrogen) && (input->atomType[ib] != nonHB_hydrogen))) ||
                        ((rd2 < 1.69) && ((input->atomType[ib] == hydrogen) || (input->atomType[ib] == nonHB_hydrogen))))
                    {
                        if (nbond == 2)
                        {
                            nbond = 3;
                            i3 = ib;
                        }
                        if (nbond == 1)
                        {
                            nbond = 2;
                            i2 = ib;
                        }
                        if (nbond == 0)
                        {
                            nbond = 1;
                            i1 = ib;
                        }
                    }
                }               // (ib != ia)

            // if no bonds, something is wrong
            if (nbond == 0)
                logFile.printErrorFormatted(WARNING, "Nitrogen atom found with no bonded atoms, atom serial number %d\n", ia);

            // one bond: Azide Nitrogen :N=C-X
            if (nbond == 1)
            {
                // calculate normalized N=C bond vector rvector[ia][]
                rd2 = 0;
                for (int i = 0; i < XYZ; i++)
                {
                    rvector[ia][i] = input->coord[ia][i] - input->coord[i1][i];
                    rd2 += sq(rvector[ia][i]);
                }
                if (rd2 < APPROX_ZERO)
                {
                    if ((rd2 == 0) && (warned == 'F'))
                    {
                        logFile.printErrorFormatted(WARNING, "Attempt to divide by zero was just prevented.\nAre the coordinates of atoms %d and %d the same?\n\n", ia, ib);
                        warned = 'T';
                    }
                    rd2 = APPROX_ZERO;
                }
                inv_rd = 1 / sqrt(rd2);
                for (int i = 0; i < XYZ; i++)
                    rvector[ia][i] *= inv_rd;
            }                   // endif nbond==1

            // two bonds: X1-N=X2
            if (nbond == 2)
            {
                // normalized vector from Nitrogen to midpoint between X1 and X2
                rd2 = 0;
                for (int i = 0; i < XYZ; i++)
                {
                    rvector[ia][i] = input->coord[ia][i] - (input->coord[i2][i] + input->coord[i1][i]) / 2;
                    rd2 += sq(rvector[ia][i]);
                }
                if (rd2 < APPROX_ZERO)
                {
                    if ((rd2 == 0) && (warned == 'F'))
                    {
                        logFile.printErrorFormatted(WARNING, "Attempt to divide by zero was just prevented.\nAre the coordinates of atoms %d and %d the same?\n\n", ia, ib);
                        warned = 'T';
                    }
                    rd2 = APPROX_ZERO;
                }
                inv_rd = 1 / sqrt(rd2);
                for (int i = 0; i < XYZ; i++)
                    rvector[ia][i] *= inv_rd;
            }                   // end two bonds for nitrogen

            // three bonds: X1,X2,X3
            if (nbond == 3)
            {
                // normalized vector from Nitrogen to midpoint between X1, X2, and X3
                rd2 = 0;
                for (int i = 0; i < XYZ; i++)
                {
                    rvector[ia][i] = input->coord[ia][i] - (input->coord[i1][i] + input->coord[i2][i] + input->coord[i3][i]) / 3;
                    rd2 += sq(rvector[ia][i]);
                }
                if (rd2 < APPROX_ZERO)
                {
                    if ((rd2 == 0) && (warned == 'F'))
                    {
                        logFile.printErrorFormatted(WARNING, "Attempt to divide by zero was just prevented.\nAre the coordinates of atoms %d and %d the same?\n\n", ia, ib);
                        warned = 'T';
                    }
                    rd2 = APPROX_ZERO;
                }
                inv_rd = 1 / sqrt(rd2);
                for (int i = 0; i < XYZ; i++)
                    rvector[ia][i] *= inv_rd;

            }                   // end three bonds for Nitrogen
            // endNEW directional N Acceptor
        }                       // end test for atom type
    }                           // Do Next receptor atom...
}
#pragma endregion Calculating bond vectors for directional H-bonds

    // End bond vector loop
    for (int k = 0; k < gridmaps.getNumAtomMaps() + 1; k++)
    {
        gridmaps[k].energyMax = (double)-BIG;
        gridmaps[k].energyMin = (double)BIG;
    }

#pragma region Writing out the correct grid filenames and other parameters
    // Write out the  correct grid_data '.fld' file_name at the  head of each map
    // file, to avoid centering errors in subsequent dockings...
    // AutoDock can then  check to see  if the  input->center of each  map  matches that
    // specified in its parameter file...

    // change numAtomMaps +1 to numAtomMaps + 2 for new dsolvPE map
    for (int k = 0; k < gridmaps.getNumMaps(); k++)
    {
        fprintf(gridmaps[k].file, "GRID_PARAMETER_FILE %s\n", programParams.getGridParameterFilename());
        fprintf(gridmaps[k].file, "GRID_DATA_FILE %s\n", input->fldFilenameAVS);
        fprintf(gridmaps[k].file, "MACROMOLECULE %s\n", input->receptorFilename);
        fprintf(gridmaps[k].file, "SPACING %.3lf\n", input->spacing);
        fprintf(gridmaps[k].file, "NELEMENTS %d %d %d\n", input->nelements[X], input->nelements[Y], input->nelements[Z]);
        fprintf(gridmaps[k].file, "CENTER %.3lf %.3lf %.3lf\n", input->center[X], input->center[Y], input->center[Z]);
    }
    FILE *floatingGridFile = 0;
    if (input->floatingGridFilename[0])
    {
        if ((floatingGridFile = openFile(input->floatingGridFilename, "w")) == 0)
        {
            logFile.printErrorFormatted(ERROR, "can't open grid map \"%s\" for writing.\n", input->floatingGridFilename);
            logFile.printError(FATAL_ERROR, "Unsuccessful completion.\n\n");
        }
        fprintf(floatingGridFile, "GRID_PARAMETER_FILE %s\n", programParams.getGridParameterFilename());
        fprintf(floatingGridFile, "GRID_DATA_FILE %s\n", input->fldFilenameAVS);
        fprintf(floatingGridFile, "MACROMOLECULE %s\n", input->receptorFilename);
        fprintf(floatingGridFile, "SPACING %.3lf\n", input->spacing);
        fprintf(floatingGridFile, "NELEMENTS %d %d %d\n", input->nelements[X], input->nelements[Y], input->nelements[Z]);
        fprintf(floatingGridFile, "CENTER %.3lf %.3lf %.3lf\n", input->center[X], input->center[Y], input->center[Z]);
    }
#pragma endregion Writing out the correct grid filenames and other parameters

#pragma region Calculation of gridmaps
{
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
    int fprintf_retval = 0;
    int nDone = 0;
    bool problemWithWriting = false;
    int hbondflag[MAX_MAPS];
    int ii = 0;
    int ic = 0;
    int closestH = 0;
    Clock grd_start;
    Clock grd_end;
    tms tms_grd_start;
    tms tms_grd_end;

    logFile.printFormatted("Beginning grid calculations.\n"
                           "\nCalculating %d grids over %d elements, around %d receptor atoms.\n\n"
                           "                    Percent   Estimated Time  Time/this plane\n"
                           "XY-plane  Z-input->coord   Done      Remaining       Real, User, System\n"
                           "            /Ang              /sec            /sec\n"
                           "________  ________  ________  ______________  __________________________\n\n",
                           gridmaps.getNumMaps() + (input->floatingGridFilename[0] != 0), input->numGridPointsPerMap, input->numReceptorAtoms);

    // Iterate over all grid points, Z(Y (X)) (X is fastest)...
    for (icoord[Z] = -input->ne[Z]; icoord[Z] <= input->ne[Z]; icoord[Z]++)
    {
        //  c[0:2] contains the current grid point.
        c[Z] = ((double)icoord[Z]) * input->spacing;
        grd_start = times(&tms_grd_start);

        for (icoord[Y] = -input->ne[Y]; icoord[Y] <= input->ne[Y]; icoord[Y]++)
        {
            c[Y] = ((double)icoord[Y]) * input->spacing;

            for (icoord[X] = -input->ne[X]; icoord[X] <= input->ne[X]; icoord[X]++)
            {
                c[X] = ((double)icoord[X]) * input->spacing;

                for (int j = 0; j < gridmaps.getNumMaps(); j++)
                    if (gridmaps[j].isCovalent)
                    {
                        // Calculate the distance from the current grid point, c, to the covalent attachment point, input->covpos
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
                    hbondflag[mapIndex] = FALSE;
                }

                // NEW2: Find Closest Hbond
                rmin = 999999;
                closestH = 0;
                for (int ia = 0; ia < input->numReceptorAtoms; ia++)
                    if ((input->hbond[ia] == 1) || (input->hbond[ia] == 2))
                    {           // DS or D1
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
                    int indx_r = min(lookup(r), MAX_DIST-1);

                    if (input->floatingGridFilename[0])
                        // Calculate the so-called "Floating Grid"...
                        r_min = min(r, r_min);

                    // elecPE is the next-to-last last grid map, i.e. electrostatics
                    if (input->distDepDiel)
                        // Distance-dependent dielectric...
                        // gridmaps.getElectrostaticMap().energy += input->charge[ia] * inv_r * input->epsilon[indx_r];
                        // apply the estat forcefield coefficient/weight here
                        gridmaps.getElectrostaticMap().energy += input->charge[ia] * inv_rmax * input->epsilon[indx_r] * model.coeff_estat;
                    else
                        // Constant dielectric...
                        // gridmaps.getElectrostaticMap().energy += input->charge[ia] * inv_r * input->invDielCal;
                        gridmaps.getElectrostaticMap().energy += input->charge[ia] * inv_rmax * input->invDielCal * model.coeff_estat;

                    // If distance from grid point to atom ia is too large,
                    // or if atom is a disordered hydrogen,
                    //   add nothing to the grid-point's non-bond energy;
                    //   just continue to next atom...
                    if (r > NBCUTOFF)
                        continue;   // onto the next atom...
                    if ((input->atomType[ia] == hydrogen) && (disorder[ia] == TRUE))
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
                        //  cos_theta = d dot rvector == cos(angle) subtended.
                        for (int i = 0; i < XYZ; i++)
                            cos_theta -= d[i] * rvector[ia][i];

                        if (cos_theta <= 0)
                            //  H->current-grid-pt vector >= 90 degrees from
                            //  N->H or O->H vector,
                            racc = 0;
                        else
                        {
                            //  racc = [cos(theta)]^2.0 for N-H
                            //  racc = [cos(theta)]^4.0 for O-H,
                            switch (rexp[ia])
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
                            // racc = pow(cos_theta, (double)rexp[ia]);

                            // NEW2 calculate dot product of bond vector with bond vector of best input->hbond
                            if (ia == closestH)
                                Hramp = 1;
                            else
                            {
                                cos_theta = 0;
                                for (int i = 0; i < XYZ; i++)
                                    cos_theta += rvector[closestH][i] * rvector[ia][i];
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
                        //  cos_theta = d dot rvector == cos(angle) subtended.
                        for (int i = 0; i < XYZ; i++)
                            cos_theta -= d[i] * rvector[ia][i];

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
                    else if ((input->hbond[ia] == 5) && (disorder[ia] == FALSE))
                    {           // A2
                        //  ia-th receptor atom = Oxygen
                        //  => receptor H-bond acceptor, oxygen.
                        rdon = 0;

                        // check to see that probe is in front of oxygen, not behind
                        cos_theta = 0;
                        for (int i = 0; i < XYZ; i++)
                            cos_theta -= d[i] * rvector[ia][i];
                        // t0 is the angle out of the lone pair plane, calculated
                        // as 90 deg - acos (vector to grid point DOT lone pair
                        // plane normal)
                        t0 = 0;
                        for (int i = 0; i < XYZ; i++)
                            t0 += d[i] * rvector2[ia][i];
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
                        cross[0] = d[1] * rvector2[ia][2] - d[2] * rvector2[ia][1];
                        cross[1] = d[2] * rvector2[ia][0] - d[0] * rvector2[ia][2];
                        cross[2] = d[0] * rvector2[ia][1] - d[1] * rvector2[ia][0];
                        rd2 = sq(cross[0]) + sq(cross[1]) + sq(cross[2]);
                        if (rd2 < APPROX_ZERO)
                        {
                            if ((rd2 == 0) && (warned == 'F'))
                            {
                                logFile.printError(WARNING, "Attempt to divide by zero was just prevented.\n\n");
                                warned = 'T';
                            }
                            rd2 = APPROX_ZERO;
                        }
                        inv_rd = 1 / sqrt(rd2);
                        ti = 0;
                        for (int i = 0; i < XYZ; i++)
                            ti += cross[i] * inv_rd * rvector[ia][i];

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
                    else if ((input->hbond[ia] == 5) && (disorder[ia] == TRUE))
                    {           // A2

                        // cylindrically disordered hydroxyl
                        cos_theta = 0;
                        for (int i = 0; i < XYZ; i++)
                            cos_theta -= d[i] * rvector[ia][i];
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
                                rsph = energyLookup->table()[input->atomType[ia]][indx_r][mapIndex] / 100;
                                rsph = max(rsph, 0);
                                rsph = min(rsph, 1);
                                if ((gridmaps[mapIndex].hbond == 3 || gridmaps[mapIndex].hbond == 5)    // AS or A2
                                    && (input->hbond[ia] == 1 || input->hbond[ia] == 2))
                                {   // DS or D1
                                    // PROBE can be an H-BOND ACCEPTOR,
                                    if (disorder[ia] == FALSE)
                                        gridmaps[mapIndex].energy += energyLookup->table()[input->atomType[ia]][indx_r][mapIndex] * Hramp * (racc + (1 - racc) * rsph);
                                    else
                                        gridmaps[mapIndex].energy += energyLookup->table()[hydrogen][max(0, indx_r - 110)][mapIndex] * Hramp * (racc + (1 - racc) * rsph);
                                }
                                else if ((gridmaps[mapIndex].hbond == 4)    // A1
                                         && (input->hbond[ia] == 1 || input->hbond[ia] == 2))
                                {   // DS,D1
                                    hbondmin[mapIndex] = min(hbondmin[mapIndex], energyLookup->table()[input->atomType[ia]][indx_r][mapIndex] * (racc + (1 - racc) * rsph));
                                    hbondmax[mapIndex] = max(hbondmax[mapIndex], energyLookup->table()[input->atomType[ia]][indx_r][mapIndex] * (racc + (1 - racc) * rsph));
                                    hbondflag[mapIndex] = TRUE;
                                }
                                else if ((gridmaps[mapIndex].hbond == 1 || gridmaps[mapIndex].hbond == 2) && (input->hbond[ia] > 2))
                                {   // DS,D1 vs AS,A1,A2
                                    // PROBE is H-BOND DONOR,
                                    temp_hbond_enrg = energyLookup->table()[input->atomType[ia]][indx_r][mapIndex] * (rdon + (1 - rdon) * rsph);
                                    hbondmin[mapIndex] = min(hbondmin[mapIndex], temp_hbond_enrg);
                                    hbondmax[mapIndex] = max(hbondmax[mapIndex], temp_hbond_enrg);
                                    hbondflag[mapIndex] = TRUE;
                                }
                                else
                                    // hbonder PROBE-ia cannot form a H-bond...,
                                    gridmaps[mapIndex].energy += energyLookup->table()[input->atomType[ia]][indx_r][mapIndex];
                            }
                            else
                                // PROBE does not form H-bonds...,
                                gridmaps[mapIndex].energy += energyLookup->table()[input->atomType[ia]][indx_r][mapIndex];

                            // add desolvation energy
                            // forcefield desolv coefficient/weight in sol_fn
                            gridmaps[mapIndex].energy += gridmaps[mapIndex].solparProbe * input->vol[ia] * sol_fn[indx_r] +
                                (input->solpar[ia] + input->solparQ * fabs(input->charge[ia])) * gridmaps[mapIndex].volProbe * sol_fn[indx_r];
                        }       // is not covalent
                    }           // mapIndex
                    gridmaps.getDesolvationMap().energy += input->solparQ * input->vol[ia] * sol_fn[indx_r];
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
                    if (!problemWithWriting)
                    {
                        if (fabs(gridmaps[k].energy) < PRECISION)
                            fprintf_retval = fprintf(gridmaps[k].file, "0.\n");
                        else
                            fprintf_retval = fprintf(gridmaps[k].file, "%.3f\n", (float)round3dp(gridmaps[k].energy));
                        if (fprintf_retval < 0)
                            problemWithWriting = TRUE;
                    }

                    gridmaps[k].energyMax = max(gridmaps[k].energyMax, gridmaps[k].energy);
                    gridmaps[k].energyMin = min(gridmaps[k].energyMin, gridmaps[k].energy);
                }
                if (floatingGridFile)
                    if ((!problemWithWriting) && (fprintf(floatingGridFile, "%.3f\n", (float)round3dp(r_min)) < 0))
                        problemWithWriting = TRUE;
                ctr++;
            }                   // icoord[X] loop
        }                       // icoord[Y] loop

        if (problemWithWriting)
            logFile.printError(WARNING, "Problems writing grid maps - there may not be enough disk space.\n");
        grd_end = times(&tms_grd_end);
        ++nDone;
        logFile.printFormatted(" %6d   %8.3lf   %5.1lf%%   ", icoord[Z], input->cgridmin[Z] + c[Z], (100 / double(input->n1[Z])) * (double)++ic);
        logFile.printTimeInHMS((grd_end - grd_start) * (input->n1[Z] - nDone));
        logFile.print("  ");
        logFile.printExecutionTimes(grd_start, grd_end, &tms_grd_start, &tms_grd_end);
    }                           // icoord[Z] loop
}
#pragma endregion Calculation of gridmaps

    delete energyLookup;

    if (input->floatingGridFilename[0])
        fclose(floatingGridFile);

#if defined(BOINCCOMPOUND)
    boinc_fraction_done(0.9);
#endif

#pragma region Writing out summary
    // Print a summary of extrema-values from the atomic-affinity and
    // electrostatics grid-maps,
    logFile.print("\nGrid\tAtom\tMinimum   \tMaximum\n"
                  "Map \tType\tEnergy    \tEnergy \n"
                  "\t\t(kcal/mol)\t(kcal/mol)\n"
                  "____\t____\t_____________\t_____________\n");

    for (int i = 0; i < gridmaps.getNumAtomMaps(); i++)
        logFile.printFormatted(" %d\t %s\t  %6.2lf\t%9.2le\n", i + 1, gridmaps[i].type, gridmaps[i].energyMin, gridmaps[i].energyMax);

    logFile.printFormatted(" %d\t %c\t  %6.2lf\t%9.2le\tElectrostatic Potential\n"
                           " %d\t %c\t  %6.2lf\t%9.2le\tDesolvation Potential\n"
                           "\n\n * Note:  Every pairwise-atomic interaction was clamped at %.2f\n\n\n",
                           gridmaps.getNumAtomMaps() + 1, 'e', gridmaps.getElectrostaticMap().energyMin, gridmaps[gridmaps.getNumAtomMaps()].energyMax,
                           gridmaps.getNumAtomMaps() + 2, 'd', gridmaps.getDesolvationMap().energyMin, gridmaps[gridmaps.getNumAtomMaps()+1].energyMax,
                           EINTCLAMP);
#pragma endregion Writing out summary

    delete input;

    fprintf(stderr, "\n%s: Successful Completion.\n", programParams.getProgramName());
    logFile.printTitled("Successful Completion.\n");

    tms tmsJobEnd;
    Clock jobEnd = times(&tmsJobEnd);
    logFile.printExecutionTimesInHMS(jobStart, jobEnd, &tmsJobStart, &tmsJobEnd);

#if defined(BOINCCOMPOUND)
    boinc_fraction_done(1);
#endif

#if defined(BOINC)
    boinc_finish(0); // should not return
#endif
}

int main(int argc, char **argv)
{
    try
    {
        appMain(argc, argv);
        return 0;
    }
    catch (ExitProgram &e)
    {
        return e.getExitCode();
    }
}
