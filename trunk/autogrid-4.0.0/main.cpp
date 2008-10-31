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
#include "program_parameters.h"
#include "read_parameter_library.h"
#include "grid_parameter_file.h"

#pragma endregion

//****************************************************************************
// Name: main (executable's name is "autogrid").
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
// Date: 07/07/04
//
// Inputs: Control file, receptor PDBQT file, parameter file
// Returns: Atomic affinity, desolvation and electrostatic grid maps.
// Globals: MAX_DIST, MAX_MAPS
// increased from 8 to 16 6/4/2004
//
// Modification Record
// Date Inits Comments
// 07/06/89 DSG FORTRAN implementation
// 07/05/92 GMM C translation
// 20/09/95 GMM/DSG AutoGrid3
// 07/07/04 DSG/RH AutoGrid4
//****************************************************************************
/* Note: 21/03/03 GMM note: ATOM_MAPS is no longer used here; was used for is_covalent and is_hbonder, but these are now folded into the MapObject and arrayed up to MAX_MAPS (currently).
   MAX_MAPS is always larger than ATOM_MAPS, so this is safe. */
int main(int argc, char **argv)
{
#if defined(_WIN32)
    SetThreadAffinityMask(GetCurrentThread(), 1);
#endif

#pragma region Declarations of variables

    // MAX_DIST is really NBCUTOFF times 100
    double energy_lookup[NUM_RECEPTOR_TYPES][MAX_DIST][MAX_MAPS];

    char *maptypeptr;           // ptr for current map->type
    MapObject *gridmap;         // contains information about gridmaps; gridmap[MAX_MAPS] is dynamically allocated

    int disorder[AG_MAX_ATOMS];
    double rvector[AG_MAX_ATOMS][XYZ];
    double rvector2[AG_MAX_ATOMS][XYZ];

    // XYZ
    double cross[XYZ];
    double c[XYZ];
    double d[XYZ];
    double dc[XYZ];
    int icoord[XYZ] = {0};            // int icenter;

    // LINE_LEN
    char message[LINE_LEN];

    // MAX_DIST
    int MD_1 = MAX_DIST - 1;
    double sol_fn[MAX_DIST];
    double energy_smooth[MAX_DIST];

    int ctr;
    char warned = 'F';

    // for NEW3 desolvation terms
    double dxA;
    double dxB;
    double minus_inv_two_sigma_sqd;
    double PI_halved = PI / 2;
    double rA;
    double rB;                  // double e;
    double rcov = 0.0;          // Distance from current grid point to the covalent attachment point
    double inv_rd, rd2, r;  // re, r2, rd,
    double r_min, inv_r, inv_rmax, racc, rdon, rsph, cos_theta, theta, tmp;
    double rdot;
    double Rij, epsij;
    double t0, ti;
    double ln_half = log(0.5);
    double cA, cB, tmpconst;
    double sigma;
    double version_num = 4.00;

    // are these necessary??
    double temp_hbond_enrg, hbondmin[MAX_MAPS], hbondmax[MAX_MAPS];
    double rmin, Hramp;

    float timeRemaining = 0.;

    int from, to;
    int fprintf_retval = 0;
    int nbond;
    int nDone = 0;
    int problem_wrt = FALSE;
    int xA, xB;
    int hbondflag[MAX_MAPS];

    int i1 = 0, i2 = 0, i3 = 0;
    int closestH = 0;

    Clock job_start;
    Clock job_end;
    struct tms tms_job_start;
    struct tms tms_job_end;

    Clock grd_start;
    Clock grd_end;
    struct tms tms_grd_start;
    struct tms tms_grd_end;
    Real idct = 1/float(CLOCKS_PER_SEC);
    char strtmp[MAX_CHARS];

#pragma endregion

    beginTimer("main() code");
    beginTimer("initialization stuff");
    ag_boinc_init();

    // Get the time at the start of the run...
    job_start = times(&tms_job_start);

    // Parse the command line...
    ProgramParameters programParams;
    process_program_parameters(argc, argv, programParams);

    // Open the log file
    FILE *logFile = stdout;
    if (programParams.logFilename[0])
        if (!(logFile = ag_fopen(programParams.logFilename, "w")))
        {
            fprintf(stderr, "\n%s: Sorry, I can't create the log file \"%s\"\n", programParams.programName, programParams.logFilename);
            fprintf(stderr, "\n%s: Unsuccessful Completion.\n\n", programParams.programName);
            exit(911);
        }

    // Output basic information
    fprint_banner(logFile, version_num);
    fprintf(logFile, "                           $Revision: 1.58 $\n\n\n");
    fprintf(logFile, "\nMaximum number of maps that can be computed = %d (defined by MAX_MAPS in \"autocomm.h\").\n\n\n", MAX_MAPS);
    fprintf(logFile, "This file was created at:\t\t\t");
    fprintf(logFile, getdate(1, strtmp, MAX_CHARS));
    fprintf(logFile, "                   using:\t\t\t\"%s\"\n", ag_gethostname(strtmp, MAX_CHARS));
    fprintf(logFile, "\n\n");

    // Read in default parameters and the grid parameter file
    beginTimer("reading in default parameters and the grid parameter file");
    Linear_FE_Model AD4;
    GridParameterInfo params;
    setup_parameter_library(params.outlev, programParams.programName, programParams.debug, logFile, AD4);
    read_grid_parameter_file(programParams, logFile, AD4, gridmap, params);
    endTimer();

#if defined(BOINCCOMPOUND)
    boinc_fraction_done(0.1);
#endif

    fprintf(logFile, "\n\nCalculating Pairwise Interaction Energies\n");
    fprintf(logFile, "=========================================\n\n");
    endTimer();

    beginTimer("precalculating vdW and H-bond interactions into the energy_lookup array");
    // do the map stuff here:
    // set up xA, xB, npb_r, npb_eps and hbonder
    // before this pt
    for (int ia = 0; ia < params.num_atom_maps; ia++)
        if (gridmap[ia].is_covalent == FALSE)
        {
            // i is the index of the receptor atom type, that the ia type ligand probe will interact with. *//* GPF_MAP
#if defined(DEBUG)
            printf("params.receptor_types_ct=%d\n", params.receptor_types_ct);
#endif
            for (int i = 0; i < params.receptor_types_ct; i++)
            {
                // for each receptor_type
                xA = gridmap[ia].xA[i];
                xB = gridmap[ia].xB[i];
                Rij = gridmap[ia].nbp_r[i];
                epsij = gridmap[ia].nbp_eps[i];
#if defined(DEBUG)
                printf("%d-%d-built xA=%d, xB=%d, npb_r=%6.3lf, nbp_eps=%10.8f for %s\n", ia, i, xA, xB, Rij, epsij, gridmap[ia].type);
#endif
                // for each receptor_type get its parms and fill in tables
                cA = (tmpconst = epsij / (double)(xA - xB)) * pow(Rij, (double)xA) * (double)xB;
                cB = tmpconst * pow(Rij, (double)xB) * (double)xA;
                if (isnan(cA))
                    print_error(programParams.programName, logFile, FATAL_ERROR, "Van der Waals coefficient cA is not a number.  AutoGrid must exit.");
                if (isnan(cB))
                    print_error(programParams.programName, logFile, FATAL_ERROR, "Van der Waals coefficient cB is not a number.  AutoGrid must exit.");
                // printf("tmpconst = %6.4f, cA = %6.4f, cB = %6.4f\n",tmpconst, cA, cB);
                dxA = (double)xA;
                dxB = (double)xB;
                if (xA == 0)
                    print_error(programParams.programName, logFile, FATAL_ERROR, "Van der Waals exponent xA is 0.  AutoGrid must exit.");
                if (xB == 0)
                    print_error(programParams.programName, logFile, FATAL_ERROR, "Van der Waals exponent xB is 0.  AutoGrid must exit.");
                fprintf(logFile, "\n             %9.1lf       %9.1lf \n", cA, cB);
                fprintf(logFile, "    E    =  -----------  -  -----------\n");
                fprintf(logFile, "     %s, %s         %2d              %2d\n", gridmap[ia].type, params.receptor_types[i], xA, xB);
                fprintf(logFile, "                r               r \n\n");
                // loop over distance index, indx_r, from 0 to MAX_DIST *//* GPF_MAP
                fprintf(logFile, "Calculating energies for %s-%s interactions.\n", gridmap[ia].type, params.receptor_types[i]);

                for (int indx_r = 1; indx_r < MAX_DIST; indx_r++)
                {
                    r = angstrom(indx_r);
                    rA = pow(r, dxA);
                    rB = pow(r, dxB);
                    energy_lookup[i][indx_r][ia] = min(EINTCLAMP, (cA / rA - cB / rB));
                }               // for each distance
                energy_lookup[i][0][ia] = EINTCLAMP;
                energy_lookup[i][MD_1][ia] = 0.;

                /* PRINT OUT INITIAL VALUES before smoothing here  fprintf(logFile, "before smoothing\n r "); for (iat = 0; iat < params.receptor_types_ct; iat++) {  fprintf(logFile,
                   " %s ", params.receptor_types[iat]); }  fprintf(logFile, "\n ___"); for (iat = 0; iat < params.receptor_types_ct; iat++) {  fprintf(logFile, " ________"); }
                   fprintf(logFile, "\n"); for (j = 0; j <= 500; j += 10) {  fprintf(logFile, "%4.1lf", angstrom(j)); for (iat = 0; iat < params.receptor_types_ct; iat++) {  fprintf(
                   logFile, (energy_lookup[iat][j][ia]<100000.)?"%9.2lf":"%9.2lg", energy_lookup[iat][j][ia]); }  fprintf(logFile, "\n"); }  fprintf(logFile, "\n"); */

                // smooth with min function *//* GPF_MAP
                if (params.i_smooth > 0)
                {
                    for (int indx_r = 1; indx_r < MAX_DIST; indx_r++)
                    {
                        energy_smooth[indx_r] = 100000.;
                        for (int j = max(0, indx_r - params.i_smooth); j < min(MAX_DIST, indx_r + params.i_smooth); j++)
                            energy_smooth[indx_r] = min(energy_smooth[indx_r], energy_lookup[i][j][ia]);
                    }
                    for (int indx_r = 1; indx_r < MAX_DIST; indx_r++)
                        energy_lookup[i][indx_r][ia] = energy_smooth[indx_r];
                }               // endif smoothing
            }                   // for i in receptor types: build energy table for this map

            // Print out a table, of distance versus energy...
            // GPF_MAP
            fprintf(logFile, "\n\nFinding the lowest pairwise interaction energy within %.1f Angstrom (\"smoothing\").\n\n  r ", params.r_smooth);
            for (int iat = 0; iat < params.receptor_types_ct; iat++)
            {
                fprintf(logFile, "    %s    ", params.receptor_types[iat]);
                //  fprintf(logFile, " %c ", receptor_atom_type_string[iat]);
            }                   // iat
            fprintf(logFile, "\n ___");
            for (int iat = 0; iat < params.receptor_types_ct; iat++)
                fprintf(logFile, " ________");                   // iat
            fprintf(logFile, "\n");
            for (int j = 0; j <= 500; j += 10)
            {
                fprintf(logFile, "%4.1lf", angstrom(j));
                for (int iat = 0; iat < params.receptor_types_ct; iat++)
                    fprintf(logFile, (energy_lookup[iat][j][ia] < 100000.) ? "%9.2lf" : "%9.2lg", energy_lookup[iat][j][ia]);               // iat
                fprintf(logFile, "\n");
            }                   // j
            fprintf(logFile, "\n");
        }
        else
        {
            // parsing for intnbp not needed for covalent maps
            fprintf(logFile, "\nAny internal non-bonded parameters will be ignored for this map, since this is a covalent map.\n");
        }                       // end of else parsing intnbp
    endTimer();

    // exponential function for receptor and ligand desolvation
    // note: the solvation term will not be smoothed
    sigma = 3.6;
    minus_inv_two_sigma_sqd = -1. / (2. * sigma * sigma);
    beginTimer("precalculating the exponential-function values for calculation of desolvation parameters");
    for (int indx_r = 1; indx_r < MAX_DIST; indx_r++)
    {
        r = angstrom(indx_r);
        // sol_fn[indx_r] = exp(-sq(r)/(2.*sigma*sigma));
        sol_fn[indx_r] = exp(sq(r) * minus_inv_two_sigma_sqd);
        sol_fn[indx_r] *= AD4.coeff_desolv;
    }
    endTimer();

    // Loop over all RECEPTOR atoms to
    // calculate bond vectors for directional H-bonds
    // setup the canned atom types here....
    // at this point set up params.hydrogen, params.carbon, params.oxygen and params.nitrogen
    params.hydrogen = get_rec_index("HD");
    params.nonHB_hydrogen = get_rec_index("H");
    params.carbon = get_rec_index("C");
    params.arom_carbon = get_rec_index("A");
    params.oxygen = get_rec_index("OA");
    params.nitrogen = get_rec_index("NA");
    params.nonHB_nitrogen = get_rec_index("N");
    params.sulphur = get_rec_index("SA");
    params.nonHB_sulphur = get_rec_index("S");

    beginTimer("calculate bond vectors for directional H-bonds");
    // 7:CHANGE HERE: scan the 'params.map_index' from the input
    for (int ia = 0; ia < params.num_receptor_atoms; ia++)
    {                                      //** ia = i_receptor_atom_a **
        disorder[ia] = FALSE;   // initialize disorder flag.
        warned = 'F';

        // Set scan limits looking for bonded atoms
        from = max(ia - 20, 0);
        to = min(ia + 20, params.num_receptor_atoms - 1);

        // If 'ia' is a params.hydrogen atom, it could be a
        // RECEPTOR params.hydrogen-BOND DONOR,
        // 8:CHANGE HERE: fix the params.atom_type vs atom_types problem in following
        if ((int)params.hbond[ia] == 2) // D1 params.hydrogen bond donor
        {
            for (int ib = from; ib <= to; ib++)
                if (ib != ia) // ib = i_receptor_atom_b
                {
                    // =>  NH-> or OH->
                    // if ((params.atom_type[ib] == params.nitrogen) || (params.atom_type[ib]==params.nonHB_nitrogen) ||(params.atom_type[ib] == params.oxygen)||(params.atom_type[ib] == params.sulphur)||(params.atom_type[ib]==params.nonHB_sulphur)) {

                    // Calculate the square of the N-H or O-H bond distance, rd2,
                    //                            ib-ia  ib-ia
                    for (int i = 0; i < XYZ; i++)
                        d[i] = params.coord[ia][i] - params.coord[ib][i];
                    rd2 = sq(d[X]) + sq(d[Y]) + sq(d[Z]);
                    // If ia & ib are less than 1.3 A apart -- they are covalently bonded,
                    if (rd2 < 1.90)
                    {           // INCREASED for H-S bonds
                        if (rd2 < APPROX_ZERO)
                        {
                            if (rd2 == 0.)
                            {
                                sprintf(message,
                                              "While calculating an H-O or H-N bond vector...\nAttempt to divide by zero was just prevented.\nAre the coordinates of atoms %d and %d the same?\n\n",
                                              ia + 1, ib + 1);
                                print_error(programParams.programName, logFile, WARNING, message);
                            }
                            rd2 = APPROX_ZERO;
                        }
                        inv_rd = 1. / sqrt(rd2);

                        // N-H: Set exponent params.rexp to 2 for m/m H-atom,
                        // if (params.atom_type[ib] == params.nitrogen) params.rexp[ia] = 2;
                        if ((params.atom_type[ib] != params.oxygen) && (params.atom_type[ib] != params.sulphur))
                            params.rexp[ia] = 2;

                        // O-H: Set exponent params.rexp to 4 for m/m H-atom,
                        // and flag disordered hydroxyls
                        if ((params.atom_type[ib] == params.oxygen) || (params.atom_type[ib] == params.sulphur))
                        {
                            params.rexp[ia] = 4;
                            if (params.disorder_h == TRUE)
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
        else if (params.hbond[ia] == 5)
        {                       // A2
            // Scan from at most, (ia-20)th m/m atom, or ia-th (if ia<20)
            //        to (ia + 5)th m/m-atom
            // determine number of atoms bonded to the params.oxygen
            nbond = 0;
            int ib;
            for (ib = from; ib <= to; ib++)
                if (ib != ia)
                {
                    rd2 = 0.;

                    for (int i = 0; i < XYZ; i++)
                    {
                        dc[i] = params.coord[ia][i] - params.coord[ib][i];
                        rd2 += sq(dc[i]);
                    }

                    // for (int i = 0; i < XYZ; i++) { rd2 += sq(params.coord[ia][i] - params.coord[ib][i]); }
                    if (((rd2 < 3.61) && ((params.atom_type[ib] != params.hydrogen) && (params.atom_type[ib] != params.nonHB_hydrogen))) ||
                        ((rd2 < 1.69) && ((params.atom_type[ib] == params.hydrogen) || (params.atom_type[ib] == params.nonHB_hydrogen))))
                    {
                        if (nbond == 2)
                        {
                            sprintf(message, "Found an H-bonding atom with three bonded atoms, atom serial number %d\n", ia + 1);
                            print_error(programParams.programName, logFile, WARNING, message);
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
            {
                sprintf(message, "Oxygen atom found with no bonded atoms, atom serial number %d, params.atom_type %d\n", ia + 1, params.atom_type[ia]);
                print_error(programParams.programName, logFile, WARNING, message);
            }

            // one bond: Carbonyl Oxygen O=C-X
            if (nbond == 1)
            {
                // calculate normalized carbonyl bond vector rvector[ia][]
                rd2 = 0.;
                for (int i = 0; i < XYZ; i++)
                {
                    rvector[ia][i] = params.coord[ia][i] - params.coord[i1][i];
                    rd2 += sq(rvector[ia][i]);
                }
                if (rd2 < APPROX_ZERO)
                {
                    if ((rd2 == 0.) && (warned == 'F'))
                    {
                        sprintf(message, "Attempt to divide by zero was just prevented.\nAre the coordinates of atoms %d and %d the same?\n\n", ia + 1, i1 + 1);
                        print_error(programParams.programName, logFile, WARNING, message);
                        warned = 'T';
                    }
                    rd2 = APPROX_ZERO;
                }
                inv_rd = 1. / sqrt(rd2);
                for (int i = 0; i < XYZ; i++)
                    rvector[ia][i] *= inv_rd;

                // find a second atom (i2) bonded to carbonyl params.carbon (i1)
                for (i2 = from; i2 <= to; i2++)
                    if ((i2 != i1) && (i2 != ia))
                    {
                        rd2 = 0.;
                        for (int i = 0; i < XYZ; i++)
                        {
                            dc[i] = params.coord[i1][i] - params.coord[i2][i];
                            /*NEW*/ rd2 += sq(dc[i]);
                        }
                        if (((rd2 < 2.89) && (params.atom_type[i2] != params.hydrogen)) || ((rd2 < 1.69) && (params.atom_type[i2] == params.hydrogen)))
                        {

                            // found one
                            // d[i] vector from params.carbon to second atom
                            rd2 = 0.;
                            for (int i = 0; i < XYZ; i++)
                            {
                                d[i] = params.coord[i2][i] - params.coord[i1][i];
                                rd2 += sq(d[i]);
                            }
                            if (rd2 < APPROX_ZERO)
                            {
                                if ((rd2 == 0.) && (warned == 'F'))
                                {
                                    sprintf(message, "Attempt to divide by zero was just prevented.\nAre the coordinates of atoms %d and %d the same?\n\n", i1 + 1, i2 + 1);
                                    print_error(programParams.programName, logFile, WARNING, message);
                                    warned = 'T';
                                }
                                rd2 = APPROX_ZERO;
                            }
                            inv_rd = 1. / sqrt(rd2);
                            for (int i = 0; i < XYZ; i++)
                                d[i] *= inv_rd;

                            // C=O cross C-X gives the lone pair plane normal
                            rvector2[ia][0] = rvector[ia][1] * d[2] - rvector[ia][2] * d[1];
                            rvector2[ia][1] = rvector[ia][2] * d[0] - rvector[ia][0] * d[2];
                            rvector2[ia][2] = rvector[ia][0] * d[1] - rvector[ia][1] * d[0];
                            rd2 = 0.;
                            for (int i = 0; i < XYZ; i++)
                                rd2 += sq(rvector2[ia][i]);
                            if (rd2 < APPROX_ZERO)
                            {
                                if ((rd2 == 0.) && (warned == 'F'))
                                {
                                    sprintf(message, "Attempt to divide by zero was just prevented.\n\n");
                                    print_error(programParams.programName, logFile, WARNING, message);
                                    warned = 'T';
                                }
                                rd2 = APPROX_ZERO;
                            }
                            inv_rd = 1. / sqrt(rd2);
                            for (int i = 0; i < XYZ; i++)
                                rvector2[ia][i] *= inv_rd;
                        }
                    }
            }                   // endif nbond==1

            // two bonds: Hydroxyl or Ether Oxygen X1-O-X2
            if (nbond == 2)
                // disordered hydroxyl
                if (((params.atom_type[i1] == params.hydrogen) || (params.atom_type[i2] == params.hydrogen)) && (params.atom_type[i1] != params.atom_type[i2]) && (params.disorder_h == TRUE))
                {

                    if ((params.atom_type[i1] == params.carbon) || (params.atom_type[i1] == params.arom_carbon))
                        ib = i1;
                    if ((params.atom_type[i2] == params.carbon) || (params.atom_type[i1] == params.arom_carbon))
                        ib = i2;
                    disorder[ia] = TRUE;
                    rd2 = 0.;
                    for (int i = 0; i < XYZ; i++)
                    {
                        rvector[ia][i] = params.coord[ia][i] - params.coord[ib][i];
                        rd2 += sq(rvector[ia][i]);
                    }
                    if (rd2 < APPROX_ZERO)
                    {
                        if ((rd2 == 0.) && (warned == 'F'))
                        {
                            sprintf(message, "Attempt to divide by zero was just prevented.\nAre the coordinates of atoms %d and %d the same?\n\n", ia + 1, ib + 1);
                            print_error(programParams.programName, logFile, WARNING, message);
                            warned = 'T';
                        }
                        rd2 = APPROX_ZERO;
                    }
                    inv_rd = 1. / sqrt(rd2);
                    for (int i = 0; i < XYZ; i++)
                        rvector[ia][i] *= inv_rd;
                }
                else
                {
                    // not a disordered hydroxyl
                    // normalized X1 to X2 vector, defines lone pair plane
                    rd2 = 0.;
                    for (int i = 0; i < XYZ; i++)
                    {
                        rvector2[ia][i] = params.coord[i2][i] - params.coord[i1][i];
                        rd2 += sq(rvector2[ia][i]);
                    }
                    if (rd2 < APPROX_ZERO)
                    {
                        if ((rd2 == 0.) && (warned == 'F'))
                        {
                            sprintf(message, "Attempt to divide by zero was just prevented.\nAre the coordinates of atoms %d and %d the same?\n\n", i1 + 1, i2 + 1);
                            print_error(programParams.programName, logFile, WARNING, message);
                            warned = 'T';
                        }
                        rd2 = APPROX_ZERO;
                    }
                    inv_rd = 1. / sqrt(rd2);
                    for (int i = 0; i < XYZ; i++)
                        rvector2[ia][i] *= inv_rd;

                    // vector pointing between the lone pairs:
                    // front of the vector is the params.oxygen atom,
                    // X1->O vector dotted with normalized X1->X2 vector plus
                    // coords of X1 gives the point on the X1-X2 line for the
                    // back of the vector.
                    rdot = 0.;
                    for (int i = 0; i < XYZ; i++)
                        rdot += (params.coord[ia][i] - params.coord[i1][i]) * rvector2[ia][i];
                    rd2 = 0.;
                    for (int i = 0; i < XYZ; i++)
                    {
                        rvector[ia][i] = params.coord[ia][i] - ((rdot * rvector2[ia][i]) + params.coord[i1][i]);
                        rd2 += sq(rvector[ia][i]);
                    }
                    if (rd2 < APPROX_ZERO)
                    {
                        if ((rd2 == 0.) && (warned == 'F'))
                        {
                            sprintf(message, "Attempt to divide by zero was just prevented.\n\n");
                            print_error(programParams.programName, logFile, WARNING, message);
                            warned = 'T';
                        }
                        rd2 = APPROX_ZERO;
                    }
                    inv_rd = 1. / sqrt(rd2);
                    for (int i = 0; i < XYZ; i++)
                        rvector[ia][i] *= inv_rd;
                }               // end disordered hydroxyl
        }
        else if (params.hbond[ia] == 4)
        {                       // A1
            // Scan from at most, (ia-20)th m/m atom, or ia-th (if ia<20)
            //        to (ia+5)th m/m-atom
            // determine number of atoms bonded to the params.oxygen
            nbond = 0;
            int ib;
            for (ib = from; ib <= to; ib++)
                if (ib != ia)
                {
                    rd2 = 0.;
                    for (int i = 0; i < XYZ; i++)
                    {
                        dc[i] = params.coord[ia][i] - params.coord[ib][i];
                        rd2 += sq(dc[i]);
                    }

                    // for (i = 0; i < XYZ; i++) { rd2 += sq(params.coord[ia][i] - params.coord[ib][i]); }
                    if (((rd2 < 2.89) && ((params.atom_type[ib] != params.hydrogen) && (params.atom_type[ib] != params.nonHB_hydrogen))) ||
                        ((rd2 < 1.69) && ((params.atom_type[ib] == params.hydrogen) || (params.atom_type[ib] == params.nonHB_hydrogen))))
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
            {
                sprintf(message, "Nitrogen atom found with no bonded atoms, atom serial number %d\n", ia);
                print_error(programParams.programName, logFile, WARNING, message);
            }

            // one bond: Azide Nitrogen :N=C-X
            if (nbond == 1)
            {
                // calculate normalized N=C bond vector rvector[ia][]
                rd2 = 0.;
                for (int i = 0; i < XYZ; i++)
                {
                    rvector[ia][i] = params.coord[ia][i] - params.coord[i1][i];
                    rd2 += sq(rvector[ia][i]);
                }
                if (rd2 < APPROX_ZERO)
                {
                    if ((rd2 == 0.) && (warned == 'F'))
                    {
                        sprintf(message, "Attempt to divide by zero was just prevented.\nAre the coordinates of atoms %d and %d the same?\n\n", ia, ib);
                        print_error(programParams.programName, logFile, WARNING, message);
                        warned = 'T';
                    }
                    rd2 = APPROX_ZERO;
                }
                inv_rd = 1. / sqrt(rd2);
                for (int i = 0; i < XYZ; i++)
                    rvector[ia][i] *= inv_rd;
            }                   // endif nbond==1

            // two bonds: X1-N=X2
            if (nbond == 2)
            {
                // normalized vector from Nitrogen to midpoint between X1 and X2
                rd2 = 0.;
                for (int i = 0; i < XYZ; i++)
                {
                    rvector[ia][i] = params.coord[ia][i] - (params.coord[i2][i] + params.coord[i1][i]) / 2.;
                    rd2 += sq(rvector[ia][i]);
                }
                if (rd2 < APPROX_ZERO)
                {
                    if ((rd2 == 0.) && (warned == 'F'))
                    {
                        sprintf(message, "Attempt to divide by zero was just prevented.\nAre the coordinates of atoms %d and %d the same?\n\n", ia, ib);
                        print_error(programParams.programName, logFile, WARNING, message);
                        warned = 'T';
                    }
                    rd2 = APPROX_ZERO;
                }
                inv_rd = 1. / sqrt(rd2);
                for (int i = 0; i < XYZ; i++)
                    rvector[ia][i] *= inv_rd;
            }                   // end two bonds for params.nitrogen

            // three bonds: X1,X2,X3
            if (nbond == 3)
            {
                // normalized vector from Nitrogen to midpoint between X1, X2, and X3
                rd2 = 0.;
                for (int i = 0; i < XYZ; i++)
                {
                    rvector[ia][i] = params.coord[ia][i] - (params.coord[i1][i] + params.coord[i2][i] + params.coord[i3][i]) / 3.;
                    rd2 += sq(rvector[ia][i]);
                }
                if (rd2 < APPROX_ZERO)
                {
                    if ((rd2 == 0.) && (warned == 'F'))
                    {
                        sprintf(message, "Attempt to divide by zero was just prevented.\nAre the coordinates of atoms %d and %d the same?\n\n", ia, ib);
                        print_error(programParams.programName, logFile, WARNING, message);
                        warned = 'T';
                    }
                    rd2 = APPROX_ZERO;
                }
                inv_rd = 1. / sqrt(rd2);
                for (int i = 0; i < XYZ; i++)
                    rvector[ia][i] *= inv_rd;

            }                   // end three bonds for Nitrogen
            // endNEW directional N Acceptor
        }                       // end test for atom type
    }                           // Do Next receptor atom...
    endTimer();

    // End bond vector loop
    for (int k = 0; k < params.num_atom_maps + 1; k++)
    {
        gridmap[k].energy_max = (double)-BIG;
        gridmap[k].energy_min = (double)BIG;
    }

    fprintf(logFile, "Beginning grid calculations.\n");
    fprintf(logFile, "\nCalculating %d grids over %d elements, around %d receptor atoms.\n\n", params.num_maps, params.num_grid_points_per_map, params.num_receptor_atoms);
    fflush(logFile);

    // Write out the  correct grid_data '.fld' file_name at the  head of each map
    // file, to avoid centering errors in subsequent dockings...
    // AutoDock can then  check to see  if the  params.center of each  map  matches that
    // specified in its parameter file...

    // change params.num_atom_maps +1 to params.num_atom_maps + 2 for new params.dsolvPE map
    for (int k = 0; k < params.num_atom_maps + 2; k++)
    {
        fprintf(gridmap[k].map_fileptr, "GRID_PARAMETER_FILE %s\n", programParams.gridParameterFilename);
        fprintf(gridmap[k].map_fileptr, "GRID_DATA_FILE %s\n", params.AVS_fld_filename);
        fprintf(gridmap[k].map_fileptr, "MACROMOLECULE %s\n", params.receptor_filename);
        fprintf(gridmap[k].map_fileptr, "SPACING %.3lf\n", params.spacing);
        fprintf(gridmap[k].map_fileptr, "NELEMENTS %d %d %d\n", params.nelements[X], params.nelements[Y], params.nelements[Z]);
        fprintf(gridmap[k].map_fileptr, "CENTER %.3lf %.3lf %.3lf\n", params.center[X], params.center[Y], params.center[Z]);
    }

    FILE *floating_grid_fileptr;
    if (params.floating_grid)
    {
        if ((floating_grid_fileptr = ag_fopen(params.floating_grid_filename, "w")) == 0)
        {
            sprintf(message, "can't open grid map \"%s\" for writing.\n", params.floating_grid_filename);
            print_error(programParams.programName, logFile, ERROR, message);
            print_error(programParams.programName, logFile, FATAL_ERROR, "Unsuccessful completion.\n\n");
        }

        fprintf(floating_grid_fileptr, "GRID_PARAMETER_FILE %s\n", programParams.gridParameterFilename);
        fprintf(floating_grid_fileptr, "GRID_DATA_FILE %s\n", params.AVS_fld_filename);
        fprintf(floating_grid_fileptr, "MACROMOLECULE %s\n", params.receptor_filename);
        fprintf(floating_grid_fileptr, "SPACING %.3lf\n", params.spacing);
        fprintf(floating_grid_fileptr, "NELEMENTS %d %d %d\n", params.nelements[X], params.nelements[Y], params.nelements[Z]);
        fprintf(floating_grid_fileptr, "CENTER %.3lf %.3lf %.3lf\n", params.center[X], params.center[Y], params.center[Z]);
    }

    fprintf(logFile, "                    Percent   Estimated Time  Time/this plane\n");
    fprintf(logFile, "XY-plane  Z-params.coord   Done      Remaining       Real, User, System\n");
    fprintf(logFile, "            /Ang              /sec            /sec\n");
    fprintf(logFile, "________  ________  ________  ______________  __________________________\n\n");

    // Iterate over all grid points, Z(Y (X)) (X is fastest)...

    beginTimer("calculating gridmap");
    int ic = 0;
    ctr = 0;
    for (icoord[Z] = -params.ne[Z]; icoord[Z] <= params.ne[Z]; icoord[Z]++)
    {
        //  c[0:2] contains the current grid point.
        c[Z] = ((double)icoord[Z]) * params.spacing;
        grd_start = times(&tms_grd_start);

        for (icoord[Y] = -params.ne[Y]; icoord[Y] <= params.ne[Y]; icoord[Y]++)
        {
            c[Y] = ((double)icoord[Y]) * params.spacing;

            for (icoord[X] = -params.ne[X]; icoord[X] <= params.ne[X]; icoord[X]++)
            {
                c[X] = ((double)icoord[X]) * params.spacing;

                for (int j = 0; j < params.num_atom_maps + 2; j++)
                    if (gridmap[j].is_covalent == TRUE)
                    {
                        // Calculate the distance from the current grid point, c, to the covalent attachment point, params.covpos
                        for (int ii = 0; ii < XYZ; ii++)
                            d[ii] = params.covpos[ii] - c[ii];
                        rcov = hypotenuse(d[X], d[Y], d[Z]);
                        rcov = rcov / params.covhalfwidth;
                        if (rcov < APPROX_ZERO)
                            rcov = APPROX_ZERO;
                        gridmap[j].energy = params.covbarrier * (1. - exp(ln_half * rcov * rcov));
                    }
                    else // is not covalent
                        gridmap[j].energy = 0.; // used to initialize to 'constant'for this gridmap

                if (params.floating_grid)
                    r_min = BIG;

                // Initialize Min Hbond variables for each new point
                for (params.map_index = 0; params.map_index < params.num_atom_maps; params.map_index++)
                {
                    hbondmin[params.map_index] = 999999.;
                    hbondmax[params.map_index] = -999999.;
                    hbondflag[params.map_index] = FALSE;
                }

                // NEW2: Find Closest Hbond
                rmin = 999999.;
                closestH = 0;
                for (int ia = 0; ia < params.num_receptor_atoms; ia++)
                    if ((params.hbond[ia] == 1) || (params.hbond[ia] == 2))
                    {           // DS or D1
                        for (int i = 0; i < XYZ; i++)
                            d[i] = params.coord[ia][i] - c[i];
                        r = hypotenuse(d[X], d[Y], d[Z]);
                        if (r < rmin)
                        {
                            rmin = r;
                            closestH = ia;
                        }
                    }           // Hydrogen test
                // END NEW2: Find Min Hbond

                //  Do all Receptor (protein, DNA, etc.) atoms...
                for (int ia = 0; ia < params.num_receptor_atoms; ia++)
                {
                    //  Get distance, r, from current grid point, c, to this receptor atom, params.coord,
                    for (int i = 0; i < XYZ; i++)
                        d[i] = params.coord[ia][i] - c[i];
                    r = hypotenuse(d[X], d[Y], d[Z]);
                    if (r < APPROX_ZERO)
                        r = APPROX_ZERO;
                    inv_r = 1. / r;
                    inv_rmax = 1. / max(r, 0.5);

                    for (int i = 0; i < XYZ; i++)
                        d[i] *= inv_r;
                    // make sure lookup index is in the table
                    int indx_r = min(lookup(r), MD_1);

                    if (params.floating_grid)
                        // Calculate the so-called "Floating Grid"...
                        r_min = min(r, r_min);

                    // params.elecPE is the next-to-last last grid map, i.e. electrostatics
                    if (params.dddiel)
                        // Distance-dependent dielectric...
                        // gridmap[params.elecPE].energy += params.charge[ia] * inv_r * params.epsilon[indx_r];
                        // apply the estat forcefield coefficient/weight here
                        gridmap[params.elecPE].energy += params.charge[ia] * inv_rmax * params.epsilon[indx_r] * AD4.coeff_estat;
                    else
                        // Constant dielectric...
                        // gridmap[params.elecPE].energy += params.charge[ia] * inv_r * params.invdielcal;
                        gridmap[params.elecPE].energy += params.charge[ia] * inv_rmax * params.invdielcal * AD4.coeff_estat;

                    // If distance from grid point to atom ia is too large,
                    // or if atom is a disordered params.hydrogen,
                    //   add nothing to the grid-point's non-bond energy;
                    //   just continue to next atom...
                    if (r > NBCUTOFF)
                        continue;   // onto the next atom...
                    if ((params.atom_type[ia] == params.hydrogen) && (disorder[ia] == TRUE))
                        continue;   // onto the next atom...

                    // racc = rdon = 1.;
                    racc = 1.;
                    rdon = 1.;
                    // NEW2 Hramp ramps in Hbond acceptor probes
                    Hramp = 1.;
                    // END NEW2 Hramp ramps in Hbond acceptor probes

                    if (params.hbond[ia] == 2)
                    {           // D1
                        //  ia-th receptor atom = Hydrogen (4 = H)
                        //  => receptor H-bond donor, OH or NH.
                        //  calculate racc for H-bond ACCEPTOR PROBES at this grid pt.
                        cos_theta = 0.;
                        //  d[] = Unit vector from current grid pt to ia_th m/m atom.
                        //  cos_theta = d dot rvector == cos(angle) subtended.
                        for (int i = 0; i < XYZ; i++)
                            cos_theta -= d[i] * rvector[ia][i];

                        if (cos_theta <= 0.)
                            //  H->current-grid-pt vector >= 90 degrees from
                            //  N->H or O->H vector,
                            racc = 0.;
                        else
                        {
                            //  racc = [cos(theta)]^2.0 for N-H
                            //  racc = [cos(theta)]^4.0 for O-H,
                            switch (params.rexp[ia])
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
                            // racc = pow(cos_theta, (double)params.rexp[ia]);

                            // NEW2 calculate dot product of bond vector with bond vector of best hbond
                            if (ia == closestH)
                                Hramp = 1.;
                            else
                            {
                                cos_theta = 0.;
                                for (int i = 0; i < XYZ; i++)
                                    cos_theta += rvector[closestH][i] * rvector[ia][i];
                                cos_theta = min(cos_theta, 1.0);
                                cos_theta = max(cos_theta, -1.0);
                                theta = acos(cos_theta);
                                Hramp = 0.5 - 0.5 * cos(theta * 120. / 90.);
                            }   // ia test
                            // END NEW2 calculate dot product of bond vector with bond vector of best hbond
                        }
                        // endif (params.atom_type[ia] == params.hydrogen)
                        // NEW Directional N acceptor
                    }
                    else if (params.hbond[ia] == 4)
                    {           // A1
                        //  ia-th macromolecule atom = Nitrogen (4 = H)
                        //  calculate rdon for H-bond Donor PROBES at this grid pt.
                        cos_theta = 0.;
                        //  d[] = Unit vector from current grid pt to ia_th m/m atom.
                        //  cos_theta = d dot rvector == cos(angle) subtended.
                        for (int i = 0; i < XYZ; i++)
                            cos_theta -= d[i] * rvector[ia][i];

                        if (cos_theta <= 0.)
                            //  H->current-grid-pt vector >= 90 degrees from
                            //  X->N vector,
                            rdon = 0.;
                        else
                            //  racc = [cos(theta)]^2.0 for H->N
                            rdon = cos_theta * cos_theta;
                        // endif (params.atom_type[ia] == params.nitrogen)
                        // end NEW Directional N acceptor

                    }
                    else if ((params.hbond[ia] == 5) && (disorder[ia] == FALSE))
                    {           // A2
                        //  ia-th receptor atom = Oxygen
                        //  => receptor H-bond acceptor, params.oxygen.
                        rdon = 0.;

                        // check to see that probe is in front of params.oxygen, not behind
                        cos_theta = 0.;
                        for (int i = 0; i < XYZ; i++)
                            cos_theta -= d[i] * rvector[ia][i];
                        // t0 is the angle out of the lone pair plane, calculated
                        // as 90 deg - acos (vector to grid point DOT lone pair
                        // plane normal)
                        t0 = 0.;
                        for (int i = 0; i < XYZ; i++)
                            t0 += d[i] * rvector2[ia][i];
                        if (t0 > 1.)
                        {
                            t0 = 1.;
                            sprintf(message, "I just prevented an attempt to take the arccosine of %f, a value greater than 1.\n", t0);
                            print_error(programParams.programName, logFile, WARNING, message);
                        }
                        else if (t0 < -1.)
                        {
                            t0 = -1.;
                            sprintf(message, "I just prevented an attempt to take the arccosine of %f, a value less than -1.\n", t0);
                            print_error(programParams.programName, logFile, WARNING, message);
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
                            if ((rd2 == 0.) && (warned == 'F'))
                            {
                                sprintf(message, "Attempt to divide by zero was just prevented.\n\n");
                                print_error(programParams.programName, logFile, WARNING, message);
                                warned = 'T';
                            }
                            rd2 = APPROX_ZERO;
                        }
                        inv_rd = 1. / sqrt(rd2);
                        ti = 0.;
                        for (int i = 0; i < XYZ; i++)
                            ti += cross[i] * inv_rd * rvector[ia][i];

                        // rdon expressions from Goodford
                        rdon = 0.;
                        if (cos_theta >= 0.)
                        {
                            if (ti > 1.)
                            {
                                ti = 1.;
                                sprintf(message, "I just prevented an attempt to take the arccosine of %f, a value greater than 1.\n", ti);
                                print_error(programParams.programName, logFile, WARNING, message);
                            }
                            else if (ti < -1.)
                            {
                                ti = -1.;
                                sprintf(message, "I just prevented an attempt to take the arccosine of %f, a value less than -1.\n", ti);
                                print_error(programParams.programName, logFile, WARNING, message);
                            }
                            ti = acos(ti) - PI_halved;
                            if (ti < 0.)
                                ti = -ti;
                            // the 2.0*ti can be replaced by (ti + ti) in: rdon = (0.9 + 0.1*sin(2.0*ti))*cos(t0);
                            rdon = (0.9 + 0.1 * sin(ti + ti)) * cos(t0);
                            // 0.34202 = cos (100 deg)
                        }
                        else if (cos_theta >= -0.34202)
                            rdon = 562.25 * pow(0.116978 - sq(cos_theta), 3.) * cos(t0);

                        // endif params.atom_type == params.oxygen, not disordered
                    }
                    else if ((params.hbond[ia] == 5) && (disorder[ia] == TRUE))
                    {           // A2

                        // cylindrically disordered hydroxyl
                        cos_theta = 0.;
                        for (int i = 0; i < XYZ; i++)
                            cos_theta -= d[i] * rvector[ia][i];
                        if (cos_theta > 1.)
                        {
                            cos_theta = 1.;
                            sprintf(message, "I just prevented an attempt to take the arccosine of %f, a value greater than 1.\n", cos_theta);
                            print_error(programParams.programName, logFile, WARNING, message);
                        }
                        else if (cos_theta < -1.)
                        {
                            cos_theta = -1.;
                            sprintf(message, "I just prevented an attempt to take the arccosine of %f, a value less than -1.\n", cos_theta);
                            print_error(programParams.programName, logFile, WARNING, message);
                        }
                        theta = acos(cos_theta);
                        racc = 0.;
                        rdon = 0.;
                        if (theta <= 1.24791 + PI_halved)
                        {
                            // 1.24791 rad = 180 deg minus C-O-H bond angle, ** 108.5 deg
                            rdon = pow(cos(theta - 1.24791), 4.);
                            racc = rdon;
                        }
                    }           // params.atom_type test

                    // For each probe atom-type,
                    // Sum pairwise interactions between each probe
                    // at this grid point (c[0:2])
                    // and the current receptor atom, ia...
                    for (params.map_index = 0; params.map_index < params.num_atom_maps; params.map_index++)
                    {
                        // We do not want to change the current enrg value for any covalent maps, make sure iscovalent is false...
                        maptypeptr = gridmap[params.map_index].type;

                        if (gridmap[params.map_index].is_covalent == FALSE)
                        {
                            if (gridmap[params.map_index].is_hbonder == TRUE)
                            {
                                // PROBE forms H-bonds...

                                // rsph ramps in angular dependence for distances with negative energy
                                rsph = energy_lookup[params.atom_type[ia]][indx_r][params.map_index] / 100.;
                                rsph = max(rsph, 0.);
                                rsph = min(rsph, 1.);
                                if ((gridmap[params.map_index].hbond == 3 || gridmap[params.map_index].hbond == 5)    // AS or A2
                                    && (params.hbond[ia] == 1 || params.hbond[ia] == 2))
                                {   // DS or D1
                                    // PROBE can be an H-BOND ACCEPTOR,
                                    if (disorder[ia] == FALSE)
                                        gridmap[params.map_index].energy += energy_lookup[params.atom_type[ia]][indx_r][params.map_index] * Hramp * (racc + (1. - racc) * rsph);
                                    else
                                        gridmap[params.map_index].energy += energy_lookup[params.hydrogen][max(0, indx_r - 110)][params.map_index] * Hramp * (racc + (1. - racc) * rsph);
                                }
                                else if ((gridmap[params.map_index].hbond == 4)    // A1
                                         && (params.hbond[ia] == 1 || params.hbond[ia] == 2))
                                {   // DS,D1
                                    hbondmin[params.map_index] = min(hbondmin[params.map_index], energy_lookup[params.atom_type[ia]][indx_r][params.map_index] * (racc + (1. - racc) * rsph));
                                    hbondmax[params.map_index] = max(hbondmax[params.map_index], energy_lookup[params.atom_type[ia]][indx_r][params.map_index] * (racc + (1. - racc) * rsph));
                                    hbondflag[params.map_index] = TRUE;
                                }
                                else if ((gridmap[params.map_index].hbond == 1 || gridmap[params.map_index].hbond == 2) && (params.hbond[ia] > 2))
                                {   // DS,D1 vs AS,A1,A2
                                    // PROBE is H-BOND DONOR,
                                    temp_hbond_enrg = energy_lookup[params.atom_type[ia]][indx_r][params.map_index] * (rdon + (1. - rdon) * rsph);
                                    hbondmin[params.map_index] = min(hbondmin[params.map_index], temp_hbond_enrg);
                                    hbondmax[params.map_index] = max(hbondmax[params.map_index], temp_hbond_enrg);
                                    hbondflag[params.map_index] = TRUE;
                                }
                                else
                                    // hbonder PROBE-ia cannot form a H-bond...,
                                    gridmap[params.map_index].energy += energy_lookup[params.atom_type[ia]][indx_r][params.map_index];
                            }
                            else
                                // PROBE does not form H-bonds...,
                                gridmap[params.map_index].energy += energy_lookup[params.atom_type[ia]][indx_r][params.map_index];

                            // add desolvation energy
                            // forcefield desolv coefficient/weight in sol_fn
                            gridmap[params.map_index].energy += gridmap[params.map_index].solpar_probe * params.vol[ia] * sol_fn[indx_r] +
                                (params.solpar[ia] + params.solpar_q * fabs(params.charge[ia])) * gridmap[params.map_index].vol_probe * sol_fn[indx_r];
                        }       // is not covalent
                    }           // params.map_index
                    gridmap[params.dsolvPE].energy += params.solpar_q * params.vol[ia] * sol_fn[indx_r];
                }               // ia loop, over all receptor atoms...
                for (params.map_index = 0; params.map_index < params.num_atom_maps; params.map_index++)
                    if (hbondflag[params.map_index])
                    {
                        gridmap[params.map_index].energy += hbondmin[params.map_index];
                        gridmap[params.map_index].energy += hbondmax[params.map_index];
                    }

                // O U T P U T . . .
                // Now output this grid point's energies to the maps:
                // 2 includes new params.dsolvPE
                for (int k = 0; k < params.num_atom_maps + 2; k++)
                {
                    if (!problem_wrt)
                    {
                        if (fabs(gridmap[k].energy) < PRECISION)
                            fprintf_retval = fprintf(gridmap[k].map_fileptr, "0.\n");
                        else
                            fprintf_retval = fprintf(gridmap[k].map_fileptr, "%.3f\n", (float)round3dp(gridmap[k].energy));
                        if (fprintf_retval < 0)
                            problem_wrt = TRUE;
                    }

                    gridmap[k].energy_max = max(gridmap[k].energy_max, gridmap[k].energy);
                    gridmap[k].energy_min = min(gridmap[k].energy_min, gridmap[k].energy);
                }
                if (params.floating_grid)
                    if ((!problem_wrt) && (fprintf(floating_grid_fileptr, "%.3f\n", (float)round3dp(r_min)) < 0))
                        problem_wrt = TRUE;
                ctr++;
            }                   // icoord[X] loop
        }                       // icoord[Y] loop

        if (problem_wrt)
        {
            sprintf(message, "Problems writing grid maps - there may not be enough disk space.\n");
            print_error(programParams.programName, logFile, WARNING, message);
        }
        grd_end = times(&tms_grd_end);
        ++nDone;
        timeRemaining = (float)(grd_end - grd_start) * idct * (float)(params.n1[Z] - nDone);
        fprintf(logFile, " %6d   %8.3lf   %5.1lf%%   ", icoord[Z], params.cgridmin[Z] + c[Z], params.percentdone * (double)++ic);
        prHMSfixed(timeRemaining, logFile);
        fprintf(logFile, "  ");
        timesys(grd_end - grd_start, &tms_grd_start, &tms_grd_end, idct, logFile);
        fflush(logFile);
    }                           // icoord[Z] loop
    endTimer();

#if defined(BOINCCOMPOUND)
    boinc_fraction_done(0.9);
#endif

    // Print a summary of extrema-values from the atomic-affinity and
    // electrostatics grid-maps,
    fprintf(logFile, "\nGrid\tAtom\tMinimum   \tMaximum\n");
    fprintf(logFile, "Map \tType\tEnergy    \tEnergy \n");
    fprintf(logFile, "\t\t(kcal/mol)\t(kcal/mol)\n");
    fprintf(logFile, "____\t____\t_____________\t_____________\n");

    {
        int i;
        for (i = 0; i < params.num_atom_maps; i++)
            fprintf(logFile, " %d\t %s\t  %6.2lf\t%9.2le\n", i + 1, gridmap[i].type, gridmap[i].energy_min, gridmap[i].energy_max);
        fprintf(logFile, " %d\t %c\t  %6.2lf\t%9.2le\tElectrostatic Potential\n", params.num_atom_maps + 1, 'e', gridmap[params.elecPE].energy_min, gridmap[i].energy_max);
        fprintf(logFile, " %d\t %c\t  %6.2lf\t%9.2le\tDesolvation Potential\n", params.num_atom_maps + 2, 'd', gridmap[params.dsolvPE].energy_min, gridmap[i + 1].energy_max);
    }
    fprintf(logFile, "\n\n * Note:  Every pairwise-atomic interaction was clamped at %.2f\n\n", EINTCLAMP);

    // Close all files

    for (int i = 0; i < params.num_atom_maps + 2; i++)
        fclose(gridmap[i].map_fileptr);
    if (params.floating_grid)
        fclose(floating_grid_fileptr);
    // Free up the memory allocated to the gridmap objects...
    free(gridmap);

    fprintf(stderr, "\n%s: Successful Completion.\n", programParams.programName);
    fprintf(logFile, "\n%s: Successful Completion.\n", programParams.programName);

    job_end = times(&tms_job_end);
    timesyshms(job_end - job_start, &tms_job_start, &tms_job_end, idct, logFile);

    fclose(logFile);

#if defined(BOINCCOMPOUND)
    boinc_fraction_done(1.);
#endif

#if defined(BOINC)
    boinc_finish(0);            // should not return
#endif

    endTimer();
    return 0;
}
