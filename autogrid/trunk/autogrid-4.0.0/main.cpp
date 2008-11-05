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

// objects
#include "exceptions.h" // ExitProgram
#include "logfile.h"    // LogFile

#pragma endregion

#pragma region struct MapObject

struct MapObject
{
    int atom_type;          // corresponds to receptor numbers????
    int mapIndex;
    int is_covalent;
    int is_hbonder;
    FILE *map_fileptr;
    char map_filename[MAX_CHARS];
    char type[3];           // eg HD or OA or NA or N
    double constant;        // this will become obsolete
    double energy_max;
    double energy_min;
    double energy;
    double vol_probe;
    double solpar_probe;
    // new 6/28
    double Rij;
    double epsij;
    HBondType hbond;       // hbonding character:
    double RijHB;
    double epsijHB;
    // per receptor type parameters, ordered as in receptor_types
    double nbp_r[NUM_RECEPTOR_TYPES];   // radius of energy-well minimum
    double nbp_eps[NUM_RECEPTOR_TYPES]; // depth of energy-well minimum
    int xA[NUM_RECEPTOR_TYPES]; // generally 12
    int xB[NUM_RECEPTOR_TYPES]; // 6 for non-hbonders 10 for h-bonders
    int hbonder[NUM_RECEPTOR_TYPES];
};

#pragma endregion

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

    double versionNumber = 4.00;

    Linear_FE_Model AD4;
    int outlev = -1;

    // Get the time at the start of the run...
    tms tms_job_start;
    Clock job_start = times(&tms_job_start);

    // Initialization
    initBoinc(); // Initialize BOINC if needed
    ProgramParameters programParams(argc, argv); // Initialize the ProgramParameters object, which parses the command-line arguments
    LogFile logFile(versionNumber, programParams.getProgramName(), programParams.getLogFilename()); // Initialize the log file

    // Read in default parameters
    setupParameterLibrary(outlev, programParams.getProgramName(), programParams.getDebugLevel(), logFile, AD4);

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // Declarations of variables needed below

    MapObject *gridmaps = 0; // array gridmaps[MAX_MAPS] is dynamically alocated

    // variables for RECEPTOR:
    // each type is now at most two characters, eg 'NA\0'
    // NB: these are sparse arrays, some entries are not set
    char receptor_types[NUM_RECEPTOR_TYPES][3] = {0};

    // number of different receptor atom types actually found in receptor PDBQT
    int receptor_types_ct = 0;

    double charge[AG_MAX_ATOMS];
    double vol[AG_MAX_ATOMS];
    double solpar[AG_MAX_ATOMS];
    int atom_type[AG_MAX_ATOMS];
    HBondType hbond[AG_MAX_ATOMS];
    double coord[AG_MAX_ATOMS][XYZ];

    double cgridmin[XYZ];
    double center[XYZ];
    double covpos[XYZ] = {0, 0, 0};         // Cartesian-coordinate of covalent affinity well.
    int ne[XYZ];
    int n1[XYZ];
    int nelements[XYZ];

    char AVS_fld_filename[MAX_CHARS] = {0};
    char floating_grid_filename[MAX_CHARS] = {0};
    char receptor_filename[MAX_CHARS] = {0};
    char xyz_filename[MAX_CHARS] = {0};

    double epsilon[MAX_DIST];

    // for NEW3 desolvation terms
    double solpar_q = .01097;   // unweighted value restored 3:9:05
    double invdielcal = 0;
    double percentdone = 0.0;
    double r_smooth = 0;
    double spacing = 0.375;     // One quarter of a C-C bond length.
    double covhalfwidth = 1.0;
    double covbarrier = 1000.0;

    int num_maps = 0;
    int num_atom_maps = -1;
    int dddiel = FALSE, disorder_h = FALSE;
    int elecPE = 0;
    int dsolvPE = 0;

    int num_grid_points_per_map = INIT_NUM_GRID_PTS;
    int i_smooth = 0;
    int num_receptor_atoms = 0;

#pragma region Reading in the grid parameter file
{
    // LIGAND: maximum is MAX_MAPS
    // each type is now at most two characters plus '\0'
    // currently ligand_atom_types is sparse... some types are not set
    char ligand_types[MAX_MAPS][3] = {0};

    // needed to make regression tests work between platforms
    Real *dummy_map;

    // number of different receptor atom types declared on receptor_types line in GPF
    int receptor_types_gpf_ct = 0;
    int has_receptor_types_in_gpf = 0;

    // array of numbers of each type
    // NB: this is a sparse int array, some entries are 0
    int receptor_atom_type_count[NUM_RECEPTOR_TYPES] = {0};

    double cext[XYZ];
    double cgridmax[XYZ];
    double cmax[XYZ] = {-BIG, -BIG, -BIG};
    double cmin[XYZ] = {BIG, BIG, BIG};
    double csum[XYZ] = {0, 0, 0};
    double cmean[XYZ];

    // LINE_LEN
    char line[LINE_LEN];
    char GPF_line[LINE_LEN];
    int length = LINE_LEN;

    char atom_name[6];
    char record[7];
    char temp_char = ' ';
    char token[5];
    char xyz[5] = "xyz";
    FILE *receptor_fileptr, *xyz_fileptr = 0;

    double q_tot = 0.0;
    double diel;
    double q_max = -BIG, q_min = BIG;
    double ri;
    double temp_vol, temp_solpar;
    double factor = 332.0L;     // Used to convert between calories and SI units

    int GPF_keyword = -1;
    int indcom = 0;
    int infld;
    int mapIndex = -1;

    // Initializes the grid parameter file
    FILE *GPF = stdin;
    if (programParams.getGridParameterFilename()[0])
    {
        GPF = openFile(programParams.getGridParameterFilename(), "r");
        if (!GPF)
        {
            fprintf(stderr, "\n%s: Sorry, I can't find or open Grid Parameter File \"%s\"\n", programParams.getProgramName(), programParams.getGridParameterFilename());
            fprintf(stderr, "\n%s: Unsuccessful Completion.\n\n", programParams.getProgramName());
            throw ExitProgram(911);
        }
    }

    // Read in the grid parameter file...
    ParameterEntry thisparm;
    ParameterEntry *found_parm;
    char FN_parameter_library[MAX_CHARS];   // the AD4 parameters .dat file name
    int parameter_library_found = 0;

    while (fgets(GPF_line, LINE_LEN, GPF) != 0)
    {
        GPF_keyword = parseGPFLine(GPF_line);

        // This first "switch" figures out how to echo the current GPF line.
        switch (GPF_keyword)
        {
        case -1:
            logFile.printFormatted("GPF> %s", GPF_line);
            logFile.printError(WARNING, "Unrecognized keyword in grid parameter file.\n");
            continue;           // while fgets GPF_line...

        case GPF_NULL:
        case GPF_COMMENT:
            logFile.printFormatted("GPF> %s", GPF_line);
            break;

        default:
            logFile.printFormatted("GPF> %s", GPF_line);
            indcom = strIndex(GPF_line, "#");
            if (indcom != -1)
                GPF_line[indcom] = '\0';    // Truncate str. at the comment
            break;

        }                       // first switch

        // This second switch interprets the current GPF line.
        switch (GPF_keyword)
        {
        case GPF_NULL:
        case GPF_COMMENT:
            break;

        case GPF_RECEPTOR:
            {
                // read in the receptor filename
                sscanf(GPF_line, "%*s %s", receptor_filename);
                logFile.printFormatted("\nReceptor Input File :\t%s\n\nReceptor Atom Type Assignments:\n\n", receptor_filename);

                // try to open receptor file
                if ((receptor_fileptr = openFile(receptor_filename, "r")) == 0)
                {
                    logFile.printErrorFormatted(ERROR, "can't find or open receptor PDBQT file \"%s\".\n", receptor_filename);
                    logFile.printError(FATAL_ERROR, "Unsuccessful completion.\n\n");
                }

                // start to read in the lines of the receptor file
                int ia = 0;
                while ((fgets(line, length, receptor_fileptr)) != 0)
                {
                    sscanf(line, "%6s", record);
                    if (equal(record, "ATOM", 4) || // Amino Acid or DNA/RNA atoms
                        equal(record, "HETA", 4) || // Non-standard heteroatoms
                        equal(record, "CHAR", 4))
                    {               // Partial Atomic Charge - not a PDB record

                        strncpy(atom_name, &line[12], 4);
                        /* atom_name is declared as an array of 6 characters, the PDB atom name is 4 characters (C indices 0, 1, 2 and 3) but let's ensure that the fifth character (C index 4)
                           is a null character, which terminates the string. */
                        atom_name[4] = '\0';

                        // Output the serial number of this atom...
                        logFile.printFormatted("Atom no. %2d, \"%s\"", ia + 1, atom_name);

                        // Read in this receptor atom's coordinates,partial charges, and solvation parameters in PDBQS format...
                        sscanf(&line[30], "%lf", &coord[ia][X]);
                        sscanf(&line[38], "%lf", &coord[ia][Y]);
                        sscanf(&line[46], "%lf", &coord[ia][Z]);

                        // Output the coordinates of this atom...
                        logFile.printFormatted(" at (%.3lf, %.3lf, %.3lf), ", coord[ia][X], coord[ia][Y], coord[ia][Z]);

                        // 1:CHANGE HERE: need to set up vol and solpar
                        sscanf(&line[70], "%lf", &charge[ia]);
                        // printf("new type is: %s\n", &line[77]);
                        sscanf(&line[77], "%s", thisparm.autogridType);
                        found_parm = atomParameterManager_find(thisparm.autogridType);
                        if (found_parm != 0)
                        {
                            logFile.printFormatted("DEBUG: found_parm->recIndex = %d", found_parm->recIndex);
                            if (found_parm->recIndex < 0)
                            {
                                strcpy(receptor_types[receptor_types_ct], found_parm->autogridType);
                                found_parm->recIndex = receptor_types_ct++;
                                logFile.printFormatted("DEBUG: found_parm->recIndex => %d", found_parm->recIndex);
                            }
                            atom_type[ia] = found_parm->recIndex;
                            solpar[ia] = found_parm->solpar;
                            vol[ia] = found_parm->vol;
                            hbond[ia] = found_parm->hbond;  // NON=0, DS,D1, AS, A1, A2
                            ++receptor_atom_type_count[found_parm->recIndex];
                        }
                        else
                        {
                            logFile.printFormatted("\n\nreceptor file contains unknown type: '%s'\nadd parameters for it to the parameter library first\n", thisparm.autogridType);
                            exit(-1);
                        }

                        // if from pdbqs: convert cal/molA**3 to kcal/molA**3
                        // solpar[ia] *= 0.001;
                        q_max = max(q_max, charge[ia]);
                        q_min = min(q_min, charge[ia]);

                        if (atom_name[0] == ' ')
                        {
                            // truncate the first character...
                            atom_name[0] = atom_name[1];
                            atom_name[1] = atom_name[2];
                            atom_name[2] = atom_name[3];
                            atom_name[3] = '\0';
                        }
                        else if ((atom_name[0] == '0' || atom_name[0] == '1' || atom_name[0] == '2' ||
                                 atom_name[0] == '3' || atom_name[0] == '4' || atom_name[0] == '5' ||
                                 atom_name[0] == '6' || atom_name[0] == '7' || atom_name[0] == '8' ||
                                 atom_name[0] == '9') && atom_name[1] == 'H')
                        {
                            /* Assume this is the 'mangled' name of a hydrogen atom, after the atom name has been changed from 'HD21' to '1HD2' for example. [0-9]H\(.\)\(.\) 0 1 2 3 : : :
                               : V V V V tmp 0 1 2 tmp : V 0 1 2 3 : : : : V V V V H\(.\)\(.\)[0-9] */
                            temp_char = atom_name[0];
                            atom_name[0] = atom_name[1];
                            atom_name[1] = atom_name[2];
                            atom_name[2] = atom_name[3];
                            atom_name[3] = temp_char;
                        }

                        // Tell the user what you thought this atom was...
                        logFile.printFormatted(" was assigned atom type \"%s\" (recIndex= %d, atom_type= %d).\n", found_parm->autogridType, found_parm->recIndex, atom_type[ia]);

                        // Count the number of each atom type
                        // ++receptor_atom_type_count[ atom_type[ia] ];

                        // Keep track of the extents of the receptor
                        for (int i = 0; i < XYZ; i++)
                        {
                            cmax[i] = max(cmax[i], coord[ia][i]);
                            cmin[i] = min(cmin[i], coord[ia][i]);
                            csum[i] += coord[ia][i];
                        }
                        // Total up the partial charges as we go...
                        q_tot += charge[ia];

                        // Increment the atom counter
                        ia++;

                        // Check that there aren't too many atoms...
                        if (ia > AG_MAX_ATOMS)
                        {
                            logFile.printErrorFormatted(ERROR, "Too many atoms in receptor PDBQT file %s;", receptor_filename);
                            logFile.printErrorFormatted(ERROR, "-- the maximum number of atoms, AG_MAX_ATOMS, allowed is %d.", AG_MAX_ATOMS);
                            logFile.printErrorFormatted(SUGGESTION, "Increase the value in the \"#define AG_MAX_ATOMS %d\" line", AG_MAX_ATOMS);
                            logFile.printError(SUGGESTION, "in the source file \"autogrid.h\", and re-compile AutoGrid.");
                            // FATAL_ERROR will cause AutoGrid to exit...
                            logFile.printError(FATAL_ERROR, "Sorry, AutoGrid cannot continue.");
                        }           // endif
                    }               // endif
                }                   // endwhile
                // Finished reading in the lines of the receptor file
                fclose(receptor_fileptr);
                if (has_receptor_types_in_gpf == 1)
                    // Check that the number of atom types found in the receptor PDBQT
                    // file match the number parsed in by the "receptor_types" command
                    // in the GPF; if they do not match, exit!
                    if (receptor_types_ct != receptor_types_gpf_ct)
                    {
                        logFile.printErrorFormatted(ERROR,
                            "The number of atom types found in the receptor PDBQT (%d) does not match the number specified by the \"receptor_types\" command (%d) in the GPF!\n\n",
                            receptor_types_ct, receptor_types_gpf_ct);
                        // FATAL_ERROR will cause AutoGrid to exit...
                        logFile.printError(FATAL_ERROR, "Sorry, AutoGrid cannot continue.");
                    }

                // Update the total number of atoms in the receptor
                num_receptor_atoms = ia;
                logFile.printFormatted("\nMaximum partial atomic charge found = %+.3lf e\n", q_max);
                logFile.printFormatted("Minimum partial atomic charge found = %+.3lf e\n\n", q_min);
                // Check there are partial charges...
                if (q_max == 0 && q_min == 0)
                {
                    logFile.printErrorFormatted(ERROR, "No partial atomic charges were found in the receptor PDBQT file %s!\n\n", receptor_filename);
                    // FATAL_ERROR will cause AutoGrid to exit...
                    logFile.printError(FATAL_ERROR, "Sorry, AutoGrid cannot continue.");
                }                   // if there are no charges EXIT

                logFile.printFormatted("Atom\tAtom\tNumber of this Type\n");
                logFile.printFormatted("Type\t ID \t in Receptor\n");
                logFile.printFormatted("____\t____\t___________________\n");
                // 2. CHANGE HERE: need to count number of each receptor_type
                for (int ia = 0; ia < receptor_types_ct; ia++)
                    if (receptor_atom_type_count[ia] != 0)
                        logFile.printFormatted(" %d\t %s\t\t%6d\n", (ia), receptor_types[ia], receptor_atom_type_count[ia]);
                logFile.printFormatted("\nTotal number of atoms :\t\t%d atoms \n", num_receptor_atoms);
                logFile.printFormatted("Total charge :\t\t\t%.2lf e\n", q_tot);
                logFile.printFormatted("\n\nReceptor coordinates fit within the following volume:\n\n");
                logFile.printFormatted("                   _______(%.1lf, %.1lf, %.1lf)\n", cmax[X], cmax[Y], cmax[Z]);
                logFile.printFormatted("                  /|     /|\n");
                logFile.printFormatted("                 / |    / |\n");
                logFile.printFormatted("                /______/  |\n");
                logFile.printFormatted("                |  |___|__| Midpoint = (%.1lf, %.1lf, %.1lf)\n", (cmax[X] + cmin[X]) / 2, (cmax[Y] + cmin[Y]) / 2, (cmax[Z] + cmin[Z]) / 2);
                logFile.printFormatted("                |  /   |  /\n");
                logFile.printFormatted("                | /    | /\n");
                logFile.printFormatted("                |/_____|/\n");
                logFile.printFormatted("(%.1lf, %.1lf, %.1lf)      \n", cmin[X], cmin[Y], cmin[Z]);
                logFile.printFormatted("\nMaximum coordinates :\t\t(%.3lf, %.3lf, %.3lf)\n", cmax[X], cmax[Y], cmax[Z]);
                logFile.printFormatted("Minimum coordinates :\t\t(%.3lf, %.3lf, %.3lf)\n\n", cmin[X], cmin[Y], cmin[Z]);
                logFile.printFormatted("\n");
                cmean[0] = csum[0] / (double)num_receptor_atoms;
                cmean[1] = csum[1] / (double)num_receptor_atoms;
                cmean[2] = csum[2] / (double)num_receptor_atoms;
            }
            break;

        case GPF_GRIDFLD:
            sscanf(GPF_line, "%*s %s", AVS_fld_filename);
            infld = strIndex(AVS_fld_filename, ".fld");
            if (infld == -1)
                logFile.printError(FATAL_ERROR, "Grid data file needs the extension \".fld\" for AVS input\n\n");
            else
            {
                infld = strIndex(AVS_fld_filename, "fld");
                strcpy(xyz_filename, AVS_fld_filename);
                xyz_filename[infld] = 'x';
                xyz_filename[infld + 1] = 'y';
                xyz_filename[infld + 2] = 'z';
            }
            if ((xyz_fileptr = openFile(xyz_filename, "w")) == 0)
            {
                logFile.printErrorFormatted(ERROR, "can't create grid extrema data file %s\n", xyz_filename);
                logFile.printError(ERROR, "SORRY!    unable to create the \".xyz\" file.\n\n");
                logFile.printError(FATAL_ERROR, "Unsuccessful completion.\n\n");
            }
            else
                logFile.printFormatted("\nCreating (AVS-readable) grid-coordinates extrema file : %s\n\n", xyz_filename);
            break;

        case GPF_NPTS:
            sscanf(GPF_line, "%*s %d %d %d", &nelements[X], &nelements[Y], &nelements[Z]);
            for (int i = 0; i < XYZ; i++)
            {
                nelements[i] = checkSize(nelements[i], xyz[i], logFile);
                ne[i] = nelements[i] / 2;
                n1[i] = nelements[i] + 1;
            }
            logFile.printFormatted("\n");
            logFile.printFormatted("Number of grid points in x-direction:\t%d\n", n1[X]);
            logFile.printFormatted("Number of grid points in y-direction:\t%d\n", n1[Y]);
            logFile.printFormatted("Number of grid points in z-direction:\t%d\n", n1[Z]);
            logFile.printFormatted("\n");
            num_grid_points_per_map = n1[X] * n1[Y] * n1[Z];
            percentdone = 100 / (double)n1[Z];
            break;

        case GPF_SPACING:
            sscanf(GPF_line, "%*s %lf", &spacing);
            logFile.printFormatted("Grid Spacing :\t\t\t%.3lf Angstrom\n", spacing);
            logFile.printFormatted("\n");
            break;

        case GPF_GRIDCENTER:
            sscanf(GPF_line, "%*s %s", token);
            if (equal(token, "auto", 4))
            {
                for (int i = 0; i < XYZ; i++)
                    center[i] = cmean[i];
                logFile.printFormatted("Grid maps will be centered on the center of mass.\n");
                logFile.printFormatted("Coordinates of center of mass : (%.3lf, %.3lf, %.3lf)\n", center[X], center[Y], center[Z]);
            }
            else
            {
                sscanf(GPF_line, "%*s %lf %lf %lf", &center[X], &center[Y], &center[Z]);
                logFile.printFormatted("\nGrid maps will be centered on user-defined coordinates:\n\n\t\t(%.3lf, %.3lf, %.3lf)\n", center[X], center[Y], center[Z]);
            }
            // centering stuff...
            for (int ia = 0; ia < num_receptor_atoms; ia++)
                for (int i = 0; i < XYZ; i++)
                    coord[ia][i] -= center[i];  // transform to center of gridmaps
            for (int i = 0; i < XYZ; i++)
            {
                cext[i] = spacing * (double)ne[i];
                cgridmax[i] = center[i] + cext[i];
                cgridmin[i] = center[i] - cext[i];
            }
            logFile.printFormatted("\nGrid maps will cover the following volume:\n\n");
            logFile.printFormatted("                   _______(%.1lf, %.1lf, %.1lf)\n", cgridmax[X], cgridmax[Y], cgridmax[Z]);
            logFile.printFormatted("                  /|     /|\n");
            logFile.printFormatted("                 / |    / |\n");
            logFile.printFormatted("                /______/  |\n");
            logFile.printFormatted("                |  |___|__| Midpoint = (%.1lf, %.1lf, %.1lf)\n", center[X], center[Y], center[Z]);
            logFile.printFormatted("                |  /   |  /\n");
            logFile.printFormatted("                | /    | /\n");
            logFile.printFormatted("                |/_____|/\n");
            logFile.printFormatted("(%.1lf, %.1lf, %.1lf)      \n\n", cgridmin[X], cgridmin[Y], cgridmin[Z]);
            for (int i = 0; i < XYZ; i++)
                logFile.printFormatted("Grid map %c-dimension :\t\t%.1lf Angstroms\n", xyz[i], 2 * cext[i]);
            logFile.printFormatted("\nMaximum coordinates :\t\t(%.3lf, %.3lf, %.3lf)\n", cgridmax[X], cgridmax[Y], cgridmax[Z]);
            logFile.printFormatted("Minimum coordinates :\t\t(%.3lf, %.3lf, %.3lf)\n\n", cgridmin[X], cgridmin[Y], cgridmin[Z]);
            for (int i = 0; i < XYZ; i++)
                fprintf(xyz_fileptr, "%.3lf %.3lf\n", cgridmin[i], cgridmax[i]);
            fclose(xyz_fileptr);
            break;

        case GPF_LIGAND_TYPES:
            {
                // Read in the list of atom types in the ligand.
                // GPF_line e.g.: "ligand_types N O A C HH NH"

                // array of ptrs used to parse input line
                char *ligand_atom_types[MAX_MAPS];

                num_atom_maps = parseTypes(GPF_line, ligand_atom_types, MAX_ATOM_TYPES);
                for (int i = 0; i < num_atom_maps; i++)
                    strcpy(ligand_types[i], ligand_atom_types[i]);
                for (int i = 0; i < num_atom_maps; i++)
                {
                    found_parm = atomParameterManager_find(ligand_types[i]);
                    if (found_parm)
                        found_parm->mapIndex = i;
                    else
                    {
                        // return error here
                        logFile.printFormatted("unknown ligand atom type %s\nadd parameters for it to the parameter library first!\n", ligand_atom_types[i]);
                        exit(-1);
                    }
                }

                elecPE = num_atom_maps;
                dsolvPE = elecPE + 1;

                /* num_maps is the number of maps to be created: the number of ligand atom types, plus 1 for the electrostatic map. AutoDock can only read in MAX_MAPS maps, which must include
                   the ligand atom maps and electrostatic map */
                num_maps = num_atom_maps + 2;

                // Check to see if there is enough memory to store these map objects
                gridmaps = new MapObject[num_maps];

                if (gridmaps == 0)
                {
                    logFile.printError(ERROR, "Could not allocate memory to create the MapObject \"gridmaps\".\n");
                    logFile.printError(FATAL_ERROR, "Unsuccessful completion.\n\n");
                }

                // Initialize the gridmaps MapObject
                for (int i = 0; i < num_maps; i++)
                {
                    gridmaps[i].atom_type = 0;   // corresponds to receptor numbers????
                    gridmaps[i].mapIndex = 0;
                    gridmaps[i].is_covalent = 0;
                    gridmaps[i].is_hbonder = 0;
                    gridmaps[i].map_fileptr = (FILE *) 0;
                    strcpy(gridmaps[i].map_filename, "");
                    strcpy(gridmaps[i].type, "");    // eg HD or OA or NA or N
                    gridmaps[i].constant = 0.0L; // this will become obsolete
                    gridmaps[i].energy_max = 0.0L;
                    gridmaps[i].energy_min = 0.0L;
                    gridmaps[i].energy = 0.0L;
                    gridmaps[i].vol_probe = 0.0L;
                    gridmaps[i].solpar_probe = 0.0L;
                    gridmaps[i].Rij = 0.0L;
                    gridmaps[i].epsij = 0.0L;
                    gridmaps[i].hbond = NON; // hbonding character:
                    gridmaps[i].RijHB = 0.0L;
                    gridmaps[i].epsijHB = 0.0L;
                    // per gridmaps[i].receptor type parameters, ordered as in receptor_types
                    for (int j = 0; j < NUM_RECEPTOR_TYPES; j++)
                    {
                        gridmaps[i].nbp_r[j] = 0.0L; // radius of energy-well minimum
                        gridmaps[i].nbp_eps[j] = 0.0L;   // depth of energy-well minimum
                        gridmaps[i].xA[j] = 0;   // generally 12
                        gridmaps[i].xB[j] = 0;   // 6 for non-hbonders 10 for h-bonders
                        gridmaps[i].hbonder[j] = 0;
                    }               // j
                }                   // i

                // Check to see if the number of grid points requested will be feasible; give warning if not enough memory.
                if (num_grid_points_per_map != INIT_NUM_GRID_PTS)
                {
                    dummy_map = (Real *) malloc(sizeof(Real) * (num_maps * num_grid_points_per_map));
                    if (!dummy_map)
                        // Too many maps requested
                        logFile.printErrorFormatted(WARNING,
                            "There will not be enough memory to store these grid maps in AutoDock; \ntry reducing the number of ligand atom types (you have %d including electrostatics) \nor reducing the size of the grid maps (you asked for %d x %d x %d grid points); \n or try running AutoDock on a machine with more RAM than this one.\n",
                            num_maps, n1[X], n1[Y], n1[Z]);
                    else
                        // free up this memory right away; we were just testing to see if we had enough when we try to run AutoDock
                        free(dummy_map);
                }
                else
                {
                    logFile.printError(ERROR, "You need to set the number of grid points using \"npts\" before setting the ligand atom types, using \"ligand_types\".\n");
                    logFile.printError(FATAL_ERROR, "Unsuccessful completion.\n\n");
                }                   // ZZZZZZZZZZZZZZZZZ
                if (!gridmaps)
                {
                    logFile.printErrorFormatted(ERROR, "Too many ligand atom types; there is not enough memory to create these maps.  Try using fewer atom types than %d.\n", num_atom_maps);
                    logFile.printError(FATAL_ERROR, "Unsuccessful completion.\n\n");
                }

                for (int i = 0; i < num_atom_maps; i++)
                {
                    gridmaps[i].is_covalent = FALSE;
                    gridmaps[i].is_hbonder = FALSE;
                    gridmaps[i].mapIndex = i;
                    strcpy(gridmaps[i].type, ligand_types[i]);   // eg HD or OA or NA or N
                    found_parm = atomParameterManager_find(ligand_types[i]);
                    gridmaps[i].atom_type = found_parm->mapIndex;
                    gridmaps[i].solpar_probe = found_parm->solpar;
                    gridmaps[i].vol_probe = found_parm->vol;
                    gridmaps[i].Rij = found_parm->Rij;
                    gridmaps[i].epsij = found_parm->epsij;
                    gridmaps[i].hbond = found_parm->hbond;
                    gridmaps[i].RijHB = found_parm->RijHB;
                    gridmaps[i].epsijHB = found_parm->epsijHB;
                    if (gridmaps[i].hbond > 0)
                        gridmaps[i].is_hbonder = TRUE;

                    for (int j = 0; j < receptor_types_ct; j++)
                    {
                        found_parm = atomParameterManager_find(receptor_types[j]);
                        gridmaps[i].nbp_r[j] = (gridmaps[i].Rij + found_parm->Rij) / 2;
                        gridmaps[i].nbp_eps[j] = sqrt(gridmaps[i].epsij * found_parm->epsij);
                        // apply the vdW forcefield parameter/weight here
                        // This was removed because "setup_p_l" does this for us... gridmaps[i].nbp_eps[j] *= FE_coeff_vdW;
                        gridmaps[i].xA[j] = 12;
                        // setup hbond dependent stuff
                        gridmaps[i].xB[j] = 6;
                        gridmaps[i].hbonder[j] = 0;
                        if ((int)(gridmaps[i].hbond) > 2 && ((int)found_parm->hbond == 1 || (int)found_parm->hbond == 2))
                        {           // AS,A1,A2 map vs DS,D1 probe
                            gridmaps[i].xB[j] = 10;
                            gridmaps[i].hbonder[j] = 1;
                            gridmaps[i].is_hbonder = TRUE;
                            // Rij and epsij for this hb interaction in parm_data.dat file as Rii and epsii for heavy atom hb factors
                            gridmaps[i].nbp_r[j] = gridmaps[i].RijHB;
                            gridmaps[i].nbp_eps[j] = gridmaps[i].epsijHB;

                            // apply the hbond forcefield parameter/weight here
                            // This was removed because "setup_p_l" does this for us... gridmaps[i].nbp_eps[j] *= FE_coeff_hbond;
                        }
                        else if (((int)gridmaps[i].hbond == 1 || (int)gridmaps[i].hbond == 2) && ((int)found_parm->hbond > 2))
                        {           // DS,D1 map vs AS,A1,A2 probe
                            gridmaps[i].xB[j] = 10;
                            gridmaps[i].hbonder[j] = 1;
                            gridmaps[i].is_hbonder = TRUE;
                            // Rij and epsij for this hb interaction in parm_data.dat file as Rii and epsii for heavy atom hb factors
                            gridmaps[i].nbp_r[j] = found_parm->RijHB;
                            gridmaps[i].nbp_eps[j] = found_parm->epsijHB;

                            // apply the hbond forcefield parameter/weight here
                            // This was removed because "setup_p_l" does this for us... gridmaps[i].nbp_eps[j] *= FE_coeff_hbond;
                        }
                    }              // initialize energy parms for each possible receptor type
                }                   // for each map
                logFile.printFormatted("\nAtom type names for ligand atom types 1-%d used for ligand-atom affinity grid maps:\n\n", num_atom_maps);
                for (int i = 0; i < num_atom_maps; i++)
                {
                    logFile.printFormatted("\t\t\tAtom type number %d corresponds to atom type name \"%s\".\n", gridmaps[i].mapIndex, gridmaps[i].type);

                    // FIX THIS!!! Covalent Atom Types are not yet supported with the new AG4/AD4 atom typing mechanism...
                    /* if (gridmaps[i].atom_type == COVALENTTYPE) { gridmaps[i].is_covalent = TRUE;  logFile.printFormatted("\nAtom type number %d will be used to calculate a covalent affinity
                       grid map\n\n", i + 1); } */
                }
                logFile.printFormatted("\n\n");
            }
            break;

        case GPF_RECEPTOR_TYPES:
            {
                // Read in the list of atom types in the receptor.
                // GPF_line e.g.: "receptor_types N O A C HH NH"
                //
                // NOTE: This line is not guaranteed to match the actual
                // atom types present in the receptor PDBQT file
                // specified by the "receptor" command.

                // array of ptrs used to parse input line
                char *receptor_atom_types[NUM_RECEPTOR_TYPES];

                receptor_types_ct = parseTypes(GPF_line, receptor_atom_types, MAX_ATOM_TYPES);
                receptor_types_gpf_ct = receptor_types_ct;
                has_receptor_types_in_gpf = 1;

                for (int i = 0; i < receptor_types_ct; i++)
                    strcpy(receptor_types[i], receptor_atom_types[i]);
                for (int i = 0; i < receptor_types_ct; i++)
                {
                    found_parm = atomParameterManager_find(receptor_atom_types[i]);
                    if (found_parm != 0)
                        found_parm->recIndex = i;
                    else
                    {
                        logFile.printErrorFormatted(ERROR, "Unknown receptor type: \"%s\"\n -- Add parameters for it to the parameter library first!\n", receptor_atom_types[i]);
                        logFile.printError(FATAL_ERROR, "Unsuccessful completion.\n\n");
                    }
                }
            }
            break;

        case GPF_SOL_PAR:      // THIS IS OBSOLETE!!!
            // Read volume and solvation parameter for probe:
            sscanf(GPF_line, "%*s %s %lf %lf", thisparm.autogridType, &temp_vol, &temp_solpar);
            found_parm = atomParameterManager_find(thisparm.autogridType);
            if (found_parm != 0)
            {
                found_parm->vol = temp_vol;
                found_parm->solpar = temp_solpar;
                int i = found_parm->mapIndex;
                if (i >= 0)
                {
                    // DON'T!!!
                    // convert cal/molA^3 to kcal/molA^3
                    // gridmaps[i].solpar_probe = temp_solpar * 0.001;
                    gridmaps[i].solpar_probe = temp_solpar;
                    logFile.printFormatted("\nProbe %s solvation parameters: \n\n\tatomic fragmental volume: %.2f A^3\n\tatomic solvation parameter: %.4f cal/mol A^3\n\n",
                                  found_parm->autogridType, found_parm->vol, found_parm->solpar);
                }
            }
            else
                logFile.printFormatted("%s key not found\n", thisparm.autogridType);
            break;              // end solvation parameter

        case GPF_MAP:
            /* The variable "mapIndex" is the 0-based index of the ligand atom type we are calculating a map for. If the "types" line was CNOSH, there would be 5 ligand atom maps to
               calculate, and since "mapIndex" is initialized to -1, mapIndex will increment each time there is a "map" keyword in the GPF.  The value of mapIndex should therefore go
               from 0 to 4 for each "map" keyword. In this example, num_atom_maps would be 5, and num_atom_maps-1 would be 4, so if mapIndex is > 4, there is something wrong in the
               number of "map" keywords. */
            ++mapIndex;
            if (mapIndex > num_atom_maps - 1)
            {
                logFile.printErrorFormatted(ERROR,
                    "Too many \"map\" keywords (%d);  the \"types\" command declares only %d maps.\nRemove a \"map\" keyword from the GPF.\n", mapIndex + 1,
                    num_atom_maps);
                logFile.printError(FATAL_ERROR, "Unsuccessful completion.\n\n");
            }
            // Read in the filename for this grid map *//* GPF_MAP
            sscanf(GPF_line, "%*s %s", gridmaps[mapIndex].map_filename);
            if ((gridmaps[mapIndex].map_fileptr = openFile(gridmaps[mapIndex].map_filename, "w")) == 0)
            {
                logFile.printErrorFormatted(ERROR, "Cannot open grid map \"%s\" for writing.", gridmaps[mapIndex].map_filename);
                logFile.printError(FATAL_ERROR, "Unsuccessful completion.\n\n");
            }
            logFile.printFormatted("\nOutput Grid Map %d:   %s\n\n", (mapIndex + 1), gridmaps[mapIndex].map_filename);
            break;

        case GPF_ELECMAP:
            sscanf(GPF_line, "%*s %s", gridmaps[elecPE].map_filename);
            if ((gridmaps[elecPE].map_fileptr = openFile(gridmaps[elecPE].map_filename, "w")) == 0)
            {
                logFile.printErrorFormatted(ERROR, "can't open grid map \"%s\" for writing.\n", gridmaps[elecPE].map_filename);
                logFile.printError(FATAL_ERROR, "Unsuccessful completion.\n\n");
            }
            logFile.printFormatted("\nOutput Electrostatic Potential Energy Grid Map: %s\n\n", gridmaps[elecPE].map_filename);
            break;

        case GPF_DSOLVMAP:
            sscanf(GPF_line, "%*s %s", gridmaps[dsolvPE].map_filename);
            if ((gridmaps[dsolvPE].map_fileptr = openFile(gridmaps[dsolvPE].map_filename, "w")) == 0)
            {
                logFile.printErrorFormatted(ERROR, "can't open grid map \"%s\" for writing.\n", gridmaps[dsolvPE].map_filename);
                logFile.printError(FATAL_ERROR, "Unsuccessful completion.\n\n");
            }
            logFile.printFormatted("\nOutput Desolvation Free Energy Grid Map: %s\n\n", gridmaps[dsolvPE].map_filename);
            break;

        case GPF_COVALENTMAP:
            sscanf(GPF_line, "%*s %lf %lf %lf %lf %lf", &covhalfwidth, &covbarrier, &(covpos[X]), &(covpos[Y]), &(covpos[Z]));
            logFile.printFormatted("\ncovalentmap <half-width in Angstroms> <barrier> <x> <y> <z>\n");
            logFile.printFormatted("\nCovalent well's half-width in Angstroms:         %8.3f\n", covhalfwidth);
            logFile.printFormatted("\nCovalent barrier energy in kcal/mol:             %8.3f\n", covbarrier);
            logFile.printFormatted("\nCovalent attachment point will be positioned at: (%8.3f, %8.3f, %8.3f)\n\n", covpos[X], covpos[Y], covpos[Z]);
            for (int i = 0; i < XYZ; i++)
                // center covpos in the grid maps frame of reference,
                covpos[i] -= center[i];
            break;

        case GPF_DISORDER:
            disorder_h = TRUE;
            logFile.printFormatted("\nHydroxyls will be disordered \n\n");
            break;

        case GPF_SMOOTH:
            sscanf(GPF_line, "%*s %lf", &r_smooth);
            logFile.printFormatted("\nPotentials will be smoothed by: %.3lf Angstrom\n\n", r_smooth);
            // Angstrom is divided by A_DIVISOR in look-up table.
            // Typical value of r_smooth is 0.5 Angstroms
            // so i_smooth = 0.5 * 100 / 2 = 25
            i_smooth = (int)(r_smooth * A_DIVISOR / 2);
            break;

        case GPF_QASP:
            sscanf(GPF_line, "%*s %lf", &solpar_q);
            logFile.printFormatted("\nCharge component of the atomic solvation parameter: %.3lf\n\n", solpar_q);
            // Typical value of solpar_q is 0.001118
            break;

        case GPF_DIEL:
            sscanf(GPF_line, "%*s %lf", &diel);
            if (diel < 0)
            {
                // negative...
                dddiel = TRUE;
                // calculate ddd of Mehler & Solmajer
                logFile.printFormatted("\nUsing *distance-dependent* dielectric function of Mehler and Solmajer, Prot.Eng.4, 903-910.\n\n");
                epsilon[0] = 1.0;
                for (int indx_r = 1; indx_r < MAX_DIST; indx_r++)
                    epsilon[indx_r] = calculateDDDMehlerSolmajer(angstrom(indx_r), APPROX_ZERO);
                logFile.printFormatted("  d   Dielectric\n ___  __________\n");
                for (int i = 0; i <= 500; i += 10)
                {
                    ri = angstrom(i);
                    logFile.printFormatted("%4.1lf%9.2lf\n", ri, epsilon[i]);
                }
                logFile.printFormatted("\n");
                // convert epsilon to 1 / epsilon
                for (int i = 1; i < MAX_DIST; i++)
                    epsilon[i] = factor / epsilon[i];
            }
            else
            {
                // positive or zero...
                dddiel = FALSE;
                if (diel <= APPROX_ZERO)
                    diel = 40;
                logFile.printFormatted("Using a *constant* dielectric of:  %.2f\n", diel);
                invdielcal = factor / diel;
            }
            break;

        case GPF_FMAP:
            sscanf(GPF_line, "%*s %s", floating_grid_filename);
            logFile.printFormatted("\nFloating Grid file name = %s\n", floating_grid_filename);
            ++num_maps;
            break;

        case GPF_PARAM_FILE:
            // open and read the AD4 parameters .dat file
            parameter_library_found = sscanf(GPF_line, "%*s %s ", FN_parameter_library);
            readParameterLibrary(FN_parameter_library, outlev, programParams.getProgramName(), programParams.getDebugLevel(), logFile, AD4);
            break;
        }                       // second switch
    }                           // while

    logFile.printFormatted("\n>>> Closing the grid parameter file (GPF)... <<<\n\n");
    logFile.printFormatted(UnderLine);
    fclose(GPF);

    if (!floating_grid_filename[0])
        logFile.printFormatted("\n\nNo Floating Grid was requested.\n");
}
#pragma endregion

#pragma region Writing to AVS_fld file
{
    FILE *AVS_fld_fileptr;

    if ((AVS_fld_fileptr = openFile(AVS_fld_filename, "w")) == 0)
    {
        logFile.printErrorFormatted(ERROR, "can't create grid dimensions data file %s\n", AVS_fld_filename);
        logFile.printError(FATAL_ERROR, "Unsuccessful completion.\n\n");
    }
    else
        logFile.printFormatted("\nCreating (AVS-readable) grid maps file : %s\n", AVS_fld_filename);
    fprintf(AVS_fld_fileptr, "# AVS field file\n#\n");
    fprintf(AVS_fld_fileptr, "# AutoDock Atomic Affinity and Electrostatic Grids\n#\n");
    fprintf(AVS_fld_fileptr, "# Created by %s.\n#\n", programParams.getProgramName());
    fprintf(AVS_fld_fileptr, "#SPACING %.3f\n", (float)spacing);
    fprintf(AVS_fld_fileptr, "#NELEMENTS %d %d %d\n", nelements[X], nelements[Y], nelements[Z]);
    fprintf(AVS_fld_fileptr, "#CENTER %.3lf %.3lf %.3lf\n", center[X], center[Y], center[Z]);
    fprintf(AVS_fld_fileptr, "#MACROMOLECULE %s\n", receptor_filename);
    fprintf(AVS_fld_fileptr, "#GRID_PARAMETER_FILE %s\n#\n", programParams.getGridParameterFilename());
    fprintf(AVS_fld_fileptr, "ndim=3\t\t\t# number of dimensions in the field\n");
    fprintf(AVS_fld_fileptr, "dim1=%d\t\t\t# number of x-elements\n", n1[X]);
    fprintf(AVS_fld_fileptr, "dim2=%d\t\t\t# number of y-elements\n", n1[Y]);
    fprintf(AVS_fld_fileptr, "dim3=%d\t\t\t# number of z-elements\n", n1[Z]);
    fprintf(AVS_fld_fileptr, "nspace=3\t\t# number of physical coordinates per point\n");
    fprintf(AVS_fld_fileptr, "veclen=%d\t\t# number of affinity values at each point\n", num_maps);
    fprintf(AVS_fld_fileptr, "data=float\t\t# data type (byte, integer, float, double)\n");
    fprintf(AVS_fld_fileptr, "field=uniform\t\t# field type (uniform, rectilinear, irregular)\n");
    for (int i = 0; i < XYZ; i++)
        fprintf(AVS_fld_fileptr, "coord %d file=%s filetype=ascii offset=%d\n", (i + 1), xyz_filename, (i * 2));
    for (int i = 0; i < num_atom_maps; i++)
        fprintf(AVS_fld_fileptr, "label=%s-affinity\t# component label for variable %d\n", gridmaps[i].type, (i + 1));                           // i
    fprintf(AVS_fld_fileptr, "label=Electrostatics\t# component label for variable %d\n", num_maps - 2);
    fprintf(AVS_fld_fileptr, "label=Desolvation\t# component label for variable %d\n", num_maps - 1);
    if (floating_grid_filename[0])
        fprintf(AVS_fld_fileptr, "label=Floating_Grid\t# component label for variable %d\n", num_maps);
    fprintf(AVS_fld_fileptr, "#\n# location of affinity grid files and how to read them\n#\n");
    for (int i = 0; i < num_atom_maps; i++)
        fprintf(AVS_fld_fileptr, "variable %d file=%s filetype=ascii skip=6\n", (i + 1), gridmaps[i].map_filename);
    fprintf(AVS_fld_fileptr, "variable %d file=%s filetype=ascii skip=6\n", num_atom_maps + 1, gridmaps[elecPE].map_filename);
    fprintf(AVS_fld_fileptr, "variable %d file=%s filetype=ascii skip=6\n", num_atom_maps + 2, gridmaps[dsolvPE].map_filename);
    if (floating_grid_filename[0])
        fprintf(AVS_fld_fileptr, "variable %d file=%s filetype=ascii skip=6\n", num_maps, floating_grid_filename);
    fclose(AVS_fld_fileptr);
}
#pragma endregion

#if defined(BOINCCOMPOUND)
    boinc_fraction_done(0.1);
#endif

    // MAX_DIST is not NBCUTOFF times 100 as it should be, it's a power of two for a little faster memory access
    // values beyond this threshold are unused
    double energy_lookup[NUM_RECEPTOR_TYPES][MAX_DIST][MAX_MAPS];
    double r;   // temporary??

#pragma region Calculating pairwise interaction energies
{
    double energy_smooth[MAX_DIST];
    double dxA;
    double dxB;
    double rA;
    double rB;
    double Rij, epsij;
    double cA, cB, tmpconst;
    int xA, xB;

    logFile.printFormatted("\n\nCalculating Pairwise Interaction Energies\n");
    logFile.printFormatted("=========================================\n\n");

    // do the map stuff here:
    // set up xA, xB, npb_r, npb_eps and hbonder
    // before this pt
    for (int ia = 0; ia < num_atom_maps; ia++)
        if (gridmaps[ia].is_covalent == FALSE)
        {
            // i is the index of the receptor atom type, that the ia type ligand probe will interact with. *//* GPF_MAP
            for (int i = 0; i < receptor_types_ct; i++)
            {
                // for each receptor_type
                xA = gridmaps[ia].xA[i];
                xB = gridmaps[ia].xB[i];
                Rij = gridmaps[ia].nbp_r[i];
                epsij = gridmaps[ia].nbp_eps[i];

                // for each receptor_type get its parms and fill in tables
                cA = (tmpconst = epsij / (double)(xA - xB)) * pow(Rij, (double)xA) * (double)xB;
                cB = tmpconst * pow(Rij, (double)xB) * (double)xA;
                if (isnan(cA))
                    logFile.printError(FATAL_ERROR, "Van der Waals coefficient cA is not a number.  AutoGrid must exit.");
                if (isnan(cB))
                    logFile.printError(FATAL_ERROR, "Van der Waals coefficient cB is not a number.  AutoGrid must exit.");
                dxA = (double)xA;
                dxB = (double)xB;
                if (xA == 0)
                    logFile.printError(FATAL_ERROR, "Van der Waals exponent xA is 0.  AutoGrid must exit.");
                if (xB == 0)
                    logFile.printError(FATAL_ERROR, "Van der Waals exponent xB is 0.  AutoGrid must exit.");
                logFile.printFormatted("\n             %9.1lf       %9.1lf \n", cA, cB);
                logFile.printFormatted("    E    =  -----------  -  -----------\n");
                logFile.printFormatted("     %s, %s         %2d              %2d\n", gridmaps[ia].type, receptor_types[i], xA, xB);
                logFile.printFormatted("                r               r \n\n");
                // loop over distance index, indx_r, from 0 to MAX_DIST *//* GPF_MAP
                logFile.printFormatted("Calculating energies for %s-%s interactions.\n", gridmaps[ia].type, receptor_types[i]);

                for (int indx_r = 1; indx_r < MAX_DIST; indx_r++)
                {
                    r = angstrom(indx_r);
                    rA = pow(r, dxA);
                    rB = pow(r, dxB);
                    energy_lookup[i][indx_r][ia] = min(EINTCLAMP, (cA / rA - cB / rB));
                }               // for each distance
                energy_lookup[i][0][ia] = EINTCLAMP;
                energy_lookup[i][MAX_DIST-1][ia] = 0;

                // smooth with min function *//* GPF_MAP
                if (i_smooth > 0)
                {
                    for (int indx_r = 1; indx_r < MAX_DIST; indx_r++)
                    {
                        energy_smooth[indx_r] = 100000;
                        for (int j = max(0, indx_r - i_smooth); j < min(MAX_DIST, indx_r + i_smooth); j++)
                            energy_smooth[indx_r] = min(energy_smooth[indx_r], energy_lookup[i][j][ia]);
                    }
                    for (int indx_r = 1; indx_r < MAX_DIST; indx_r++)
                        energy_lookup[i][indx_r][ia] = energy_smooth[indx_r];
                }               // endif smoothing
            }                   // for i in receptor types: build energy table for this map

            // Print out a table, of distance versus energy...
            // GPF_MAP
            logFile.printFormatted("\n\nFinding the lowest pairwise interaction energy within %.1f Angstrom (\"smoothing\").\n\n  r ", r_smooth);
            for (int iat = 0; iat < receptor_types_ct; iat++)
                logFile.printFormatted("    %s    ", receptor_types[iat]);
            logFile.printFormatted("\n ___");
            for (int iat = 0; iat < receptor_types_ct; iat++)
                logFile.printFormatted(" ________");                   // iat
            logFile.printFormatted("\n");
            for (int j = 0; j <= 500; j += 10)
            {
                logFile.printFormatted("%4.1lf", angstrom(j));
                for (int iat = 0; iat < receptor_types_ct; iat++)
                    logFile.printFormatted((energy_lookup[iat][j][ia] < 100000) ? "%9.2lf" : "%9.2lg", energy_lookup[iat][j][ia]);               // iat
                logFile.printFormatted("\n");
            }                   // j
            logFile.printFormatted("\n");
        }
        else
        {
            // parsing for intnbp not needed for covalent maps
            logFile.printFormatted("\nAny internal non-bonded parameters will be ignored for this map, since this is a covalent map.\n");
        }                       // end of else parsing intnbp
}
#pragma endregion

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
        r = angstrom(indx_r);
        // sol_fn[indx_r] = exp(-sq(r)/(2*sigma*sigma));
        sol_fn[indx_r] = exp(sq(r) * minus_inv_two_sigma_sqd);
        sol_fn[indx_r] *= AD4.coeff_desolv;
    }
}
#pragma endregion

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
    for (int ia = 0; ia < num_receptor_atoms; ia++)
    {                                      //** ia = i_receptor_atom_a **
        disorder[ia] = FALSE;   // initialize disorder flag.
        warned = 'F';

        // Set scan limits looking for bonded atoms
        from = max(ia - 20, 0);
        to = min(ia + 20, num_receptor_atoms - 1);

        // If 'ia' is a hydrogen atom, it could be a
        // RECEPTOR hydrogen-BOND DONOR,
        // 8:CHANGE HERE: fix the atom_type vs atom_types problem in following
        if ((int)hbond[ia] == 2) // D1 hydrogen bond donor
        {
            for (int ib = from; ib <= to; ib++)
                if (ib != ia) // ib = i_receptor_atom_b
                {
                    // =>  NH-> or OH->
                    // if ((atom_type[ib] == nitrogen) || (atom_type[ib]==nonHB_nitrogen) ||(atom_type[ib] == oxygen)||(atom_type[ib] == sulphur)||(atom_type[ib]==nonHB_sulphur)) {

                    // Calculate the square of the N-H or O-H bond distance, rd2,
                    //                            ib-ia  ib-ia
                    for (int i = 0; i < XYZ; i++)
                        d[i] = coord[ia][i] - coord[ib][i];
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
                        // if (atom_type[ib] == nitrogen) rexp[ia] = 2;
                        if ((atom_type[ib] != oxygen) && (atom_type[ib] != sulphur))
                            rexp[ia] = 2;

                        // O-H: Set exponent rexp to 4 for m/m H-atom,
                        // and flag disordered hydroxyls
                        if ((atom_type[ib] == oxygen) || (atom_type[ib] == sulphur))
                        {
                            rexp[ia] = 4;
                            if (disorder_h == TRUE)
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
        else if (hbond[ia] == 5)
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
                        dc[i] = coord[ia][i] - coord[ib][i];
                        rd2 += sq(dc[i]);
                    }

                    if (((rd2 < 3.61) && ((atom_type[ib] != hydrogen) && (atom_type[ib] != nonHB_hydrogen))) ||
                        ((rd2 < 1.69) && ((atom_type[ib] == hydrogen) || (atom_type[ib] == nonHB_hydrogen))))
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
                logFile.printErrorFormatted(WARNING, "Oxygen atom found with no bonded atoms, atom serial number %d, atom_type %d\n", ia + 1, atom_type[ia]);

            // one bond: Carbonyl Oxygen O=C-X
            if (nbond == 1)
            {
                // calculate normalized carbonyl bond vector rvector[ia][]
                rd2 = 0;
                for (int i = 0; i < XYZ; i++)
                {
                    rvector[ia][i] = coord[ia][i] - coord[i1][i];
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
                            dc[i] = coord[i1][i] - coord[i2][i];
                            /*NEW*/ rd2 += sq(dc[i]);
                        }
                        if (((rd2 < 2.89) && (atom_type[i2] != hydrogen)) || ((rd2 < 1.69) && (atom_type[i2] == hydrogen)))
                        {

                            // found one
                            // d[i] vector from carbon to second atom
                            rd2 = 0;
                            for (int i = 0; i < XYZ; i++)
                            {
                                d[i] = coord[i2][i] - coord[i1][i];
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
                if (((atom_type[i1] == hydrogen) || (atom_type[i2] == hydrogen)) && (atom_type[i1] != atom_type[i2]) && (disorder_h == TRUE))
                {
                    if ((atom_type[i1] == carbon) || (atom_type[i1] == arom_carbon))
                        ib = i1;
                    if ((atom_type[i2] == carbon) || (atom_type[i1] == arom_carbon))
                        ib = i2;
                    disorder[ia] = TRUE;
                    rd2 = 0;
                    for (int i = 0; i < XYZ; i++)
                    {
                        rvector[ia][i] = coord[ia][i] - coord[ib][i];
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
                        rvector2[ia][i] = coord[i2][i] - coord[i1][i];
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
                        rdot += (coord[ia][i] - coord[i1][i]) * rvector2[ia][i];
                    rd2 = 0;
                    for (int i = 0; i < XYZ; i++)
                    {
                        rvector[ia][i] = coord[ia][i] - ((rdot * rvector2[ia][i]) + coord[i1][i]);
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
        else if (hbond[ia] == 4)
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
                        dc[i] = coord[ia][i] - coord[ib][i];
                        rd2 += sq(dc[i]);
                    }

                    if (((rd2 < 2.89) && ((atom_type[ib] != hydrogen) && (atom_type[ib] != nonHB_hydrogen))) ||
                        ((rd2 < 1.69) && ((atom_type[ib] == hydrogen) || (atom_type[ib] == nonHB_hydrogen))))
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
                    rvector[ia][i] = coord[ia][i] - coord[i1][i];
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
                    rvector[ia][i] = coord[ia][i] - (coord[i2][i] + coord[i1][i]) / 2;
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
                    rvector[ia][i] = coord[ia][i] - (coord[i1][i] + coord[i2][i] + coord[i3][i]) / 3;
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
#pragma endregion

    // End bond vector loop
    for (int k = 0; k < num_atom_maps + 1; k++)
    {
        gridmaps[k].energy_max = (double)-BIG;
        gridmaps[k].energy_min = (double)BIG;
    }

#pragma region Writing out the correct grid filenames and other parameters
    // Write out the  correct grid_data '.fld' file_name at the  head of each map
    // file, to avoid centering errors in subsequent dockings...
    // AutoDock can then  check to see  if the  center of each  map  matches that
    // specified in its parameter file...

    // change num_atom_maps +1 to num_atom_maps + 2 for new dsolvPE map
    for (int k = 0; k < num_atom_maps + 2; k++)
    {
        fprintf(gridmaps[k].map_fileptr, "GRID_PARAMETER_FILE %s\n", programParams.getGridParameterFilename());
        fprintf(gridmaps[k].map_fileptr, "GRID_DATA_FILE %s\n", AVS_fld_filename);
        fprintf(gridmaps[k].map_fileptr, "MACROMOLECULE %s\n", receptor_filename);
        fprintf(gridmaps[k].map_fileptr, "SPACING %.3lf\n", spacing);
        fprintf(gridmaps[k].map_fileptr, "NELEMENTS %d %d %d\n", nelements[X], nelements[Y], nelements[Z]);
        fprintf(gridmaps[k].map_fileptr, "CENTER %.3lf %.3lf %.3lf\n", center[X], center[Y], center[Z]);
    }
    FILE *floating_grid_fileptr = 0;
    if (floating_grid_filename[0])
    {
        if ((floating_grid_fileptr = openFile(floating_grid_filename, "w")) == 0)
        {
            logFile.printErrorFormatted(ERROR, "can't open grid map \"%s\" for writing.\n", floating_grid_filename);
            logFile.printError(FATAL_ERROR, "Unsuccessful completion.\n\n");
        }
        fprintf(floating_grid_fileptr, "GRID_PARAMETER_FILE %s\n", programParams.getGridParameterFilename());
        fprintf(floating_grid_fileptr, "GRID_DATA_FILE %s\n", AVS_fld_filename);
        fprintf(floating_grid_fileptr, "MACROMOLECULE %s\n", receptor_filename);
        fprintf(floating_grid_fileptr, "SPACING %.3lf\n", spacing);
        fprintf(floating_grid_fileptr, "NELEMENTS %d %d %d\n", nelements[X], nelements[Y], nelements[Z]);
        fprintf(floating_grid_fileptr, "CENTER %.3lf %.3lf %.3lf\n", center[X], center[Y], center[Z]);
    }
#pragma endregion

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
    int problem_wrt = FALSE;
    int hbondflag[MAX_MAPS];
    int ii = 0;
    int ic = 0;
    int closestH = 0;
    Clock grd_start;
    Clock grd_end;
    tms tms_grd_start;
    tms tms_grd_end;

    logFile.printFormatted("Beginning grid calculations.\n");
    logFile.printFormatted("\nCalculating %d grids over %d elements, around %d receptor atoms.\n\n", num_maps, num_grid_points_per_map, num_receptor_atoms);
    logFile.printFormatted("                    Percent   Estimated Time  Time/this plane\n");
    logFile.printFormatted("XY-plane  Z-coord   Done      Remaining       Real, User, System\n");
    logFile.printFormatted("            /Ang              /sec            /sec\n");
    logFile.printFormatted("________  ________  ________  ______________  __________________________\n\n");

    // Iterate over all grid points, Z(Y (X)) (X is fastest)...
    for (icoord[Z] = -ne[Z]; icoord[Z] <= ne[Z]; icoord[Z]++)
    {
        //  c[0:2] contains the current grid point.
        c[Z] = ((double)icoord[Z]) * spacing;
        grd_start = times(&tms_grd_start);

        for (icoord[Y] = -ne[Y]; icoord[Y] <= ne[Y]; icoord[Y]++)
        {
            c[Y] = ((double)icoord[Y]) * spacing;

            for (icoord[X] = -ne[X]; icoord[X] <= ne[X]; icoord[X]++)
            {
                c[X] = ((double)icoord[X]) * spacing;

                for (int j = 0; j < num_atom_maps + 2; j++)
                    if (gridmaps[j].is_covalent == TRUE)
                    {
                        // Calculate the distance from the current grid point, c, to the covalent attachment point, covpos
                        for (ii = 0; ii < XYZ; ii++)
                            d[ii] = covpos[ii] - c[ii];
                        rcov = hypotenuse(d[X], d[Y], d[Z]);
                        rcov = rcov / covhalfwidth;
                        if (rcov < APPROX_ZERO)
                            rcov = APPROX_ZERO;
                        gridmaps[j].energy = covbarrier * (1 - exp(ln_half * rcov * rcov));
                    }
                    else // is not covalent
                        gridmaps[j].energy = 0; // used to initialize to 'constant'for this gridmaps

                if (floating_grid_filename[0])
                    r_min = BIG;

                // Initialize Min Hbond variables for each new point
                for (int mapIndex = 0; mapIndex < num_atom_maps; mapIndex++)
                {
                    hbondmin[mapIndex] = 999999;
                    hbondmax[mapIndex] = -999999;
                    hbondflag[mapIndex] = FALSE;
                }

                // NEW2: Find Closest Hbond
                rmin = 999999;
                closestH = 0;
                for (int ia = 0; ia < num_receptor_atoms; ia++)
                    if ((hbond[ia] == 1) || (hbond[ia] == 2))
                    {           // DS or D1
                        for (int i = 0; i < XYZ; i++)
                            d[i] = coord[ia][i] - c[i];
                        r = hypotenuse(d[X], d[Y], d[Z]);
                        if (r < rmin)
                        {
                            rmin = r;
                            closestH = ia;
                        }
                    }           // Hydrogen test
                // END NEW2: Find Min Hbond

                //  Do all Receptor (protein, DNA, etc.) atoms...
                for (int ia = 0; ia < num_receptor_atoms; ia++)
                {
                    //  Get distance, r, from current grid point, c, to this receptor atom, coord,
                    for (int i = 0; i < XYZ; i++)
                        d[i] = coord[ia][i] - c[i];
                    r = hypotenuse(d[X], d[Y], d[Z]);
                    if (r < APPROX_ZERO)
                        r = APPROX_ZERO;
                    inv_r = 1 / r;
                    inv_rmax = 1 / max(r, 0.5);

                    for (int i = 0; i < XYZ; i++)
                        d[i] *= inv_r;
                    // make sure lookup index is in the table
                    int indx_r = min(lookup(r), MAX_DIST-1);

                    if (floating_grid_filename[0])
                        // Calculate the so-called "Floating Grid"...
                        r_min = min(r, r_min);

                    // elecPE is the next-to-last last grid map, i.e. electrostatics
                    if (dddiel)
                        // Distance-dependent dielectric...
                        // gridmaps[elecPE].energy += charge[ia] * inv_r * epsilon[indx_r];
                        // apply the estat forcefield coefficient/weight here
                        gridmaps[elecPE].energy += charge[ia] * inv_rmax * epsilon[indx_r] * AD4.coeff_estat;
                    else
                        // Constant dielectric...
                        // gridmaps[elecPE].energy += charge[ia] * inv_r * invdielcal;
                        gridmaps[elecPE].energy += charge[ia] * inv_rmax * invdielcal * AD4.coeff_estat;

                    // If distance from grid point to atom ia is too large,
                    // or if atom is a disordered hydrogen,
                    //   add nothing to the grid-point's non-bond energy;
                    //   just continue to next atom...
                    if (r > NBCUTOFF)
                        continue;   // onto the next atom...
                    if ((atom_type[ia] == hydrogen) && (disorder[ia] == TRUE))
                        continue;   // onto the next atom...

                    // racc = rdon = 1;
                    racc = 1;
                    rdon = 1;
                    // NEW2 Hramp ramps in Hbond acceptor probes
                    Hramp = 1;
                    // END NEW2 Hramp ramps in Hbond acceptor probes

                    if (hbond[ia] == 2)
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

                            // NEW2 calculate dot product of bond vector with bond vector of best hbond
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
                            // END NEW2 calculate dot product of bond vector with bond vector of best hbond
                        }
                        // endif (atom_type[ia] == hydrogen)
                        // NEW Directional N acceptor
                    }
                    else if (hbond[ia] == 4)
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
                        // endif (atom_type[ia] == nitrogen)
                        // end NEW Directional N acceptor

                    }
                    else if ((hbond[ia] == 5) && (disorder[ia] == FALSE))
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

                        // endif atom_type == oxygen, not disordered
                    }
                    else if ((hbond[ia] == 5) && (disorder[ia] == TRUE))
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
                    }           // atom_type test

                    // For each probe atom-type,
                    // Sum pairwise interactions between each probe
                    // at this grid point (c[0:2])
                    // and the current receptor atom, ia...
                    for (int mapIndex = 0; mapIndex < num_atom_maps; mapIndex++)
                    {
                        // We do not want to change the current enrg value for any covalent maps, make sure iscovalent is false...
                        maptypeptr = gridmaps[mapIndex].type;

                        if (gridmaps[mapIndex].is_covalent == FALSE)
                        {
                            if (gridmaps[mapIndex].is_hbonder == TRUE)
                            {
                                // PROBE forms H-bonds...

                                // rsph ramps in angular dependence for distances with negative energy
                                rsph = energy_lookup[atom_type[ia]][indx_r][mapIndex] / 100;
                                rsph = max(rsph, 0);
                                rsph = min(rsph, 1);
                                if ((gridmaps[mapIndex].hbond == 3 || gridmaps[mapIndex].hbond == 5)    // AS or A2
                                    && (hbond[ia] == 1 || hbond[ia] == 2))
                                {   // DS or D1
                                    // PROBE can be an H-BOND ACCEPTOR,
                                    if (disorder[ia] == FALSE)
                                        gridmaps[mapIndex].energy += energy_lookup[atom_type[ia]][indx_r][mapIndex] * Hramp * (racc + (1 - racc) * rsph);
                                    else
                                        gridmaps[mapIndex].energy += energy_lookup[hydrogen][max(0, indx_r - 110)][mapIndex] * Hramp * (racc + (1 - racc) * rsph);
                                }
                                else if ((gridmaps[mapIndex].hbond == 4)    // A1
                                         && (hbond[ia] == 1 || hbond[ia] == 2))
                                {   // DS,D1
                                    hbondmin[mapIndex] = min(hbondmin[mapIndex], energy_lookup[atom_type[ia]][indx_r][mapIndex] * (racc + (1 - racc) * rsph));
                                    hbondmax[mapIndex] = max(hbondmax[mapIndex], energy_lookup[atom_type[ia]][indx_r][mapIndex] * (racc + (1 - racc) * rsph));
                                    hbondflag[mapIndex] = TRUE;
                                }
                                else if ((gridmaps[mapIndex].hbond == 1 || gridmaps[mapIndex].hbond == 2) && (hbond[ia] > 2))
                                {   // DS,D1 vs AS,A1,A2
                                    // PROBE is H-BOND DONOR,
                                    temp_hbond_enrg = energy_lookup[atom_type[ia]][indx_r][mapIndex] * (rdon + (1 - rdon) * rsph);
                                    hbondmin[mapIndex] = min(hbondmin[mapIndex], temp_hbond_enrg);
                                    hbondmax[mapIndex] = max(hbondmax[mapIndex], temp_hbond_enrg);
                                    hbondflag[mapIndex] = TRUE;
                                }
                                else
                                    // hbonder PROBE-ia cannot form a H-bond...,
                                    gridmaps[mapIndex].energy += energy_lookup[atom_type[ia]][indx_r][mapIndex];
                            }
                            else
                                // PROBE does not form H-bonds...,
                                gridmaps[mapIndex].energy += energy_lookup[atom_type[ia]][indx_r][mapIndex];

                            // add desolvation energy
                            // forcefield desolv coefficient/weight in sol_fn
                            gridmaps[mapIndex].energy += gridmaps[mapIndex].solpar_probe * vol[ia] * sol_fn[indx_r] +
                                (solpar[ia] + solpar_q * fabs(charge[ia])) * gridmaps[mapIndex].vol_probe * sol_fn[indx_r];
                        }       // is not covalent
                    }           // mapIndex
                    gridmaps[dsolvPE].energy += solpar_q * vol[ia] * sol_fn[indx_r];
                }               // ia loop, over all receptor atoms...
                for (int mapIndex = 0; mapIndex < num_atom_maps; mapIndex++)
                    if (hbondflag[mapIndex])
                    {
                        gridmaps[mapIndex].energy += hbondmin[mapIndex];
                        gridmaps[mapIndex].energy += hbondmax[mapIndex];
                    }

                // O U T P U T . . .
                // Now output this grid point's energies to the maps:
                // 2 includes new dsolvPE
                for (int k = 0; k < num_atom_maps + 2; k++)
                {
                    if (!problem_wrt)
                    {
                        if (fabs(gridmaps[k].energy) < PRECISION)
                            fprintf_retval = fprintf(gridmaps[k].map_fileptr, "0.\n");
                        else
                            fprintf_retval = fprintf(gridmaps[k].map_fileptr, "%.3f\n", (float)round3dp(gridmaps[k].energy));
                        if (fprintf_retval < 0)
                            problem_wrt = TRUE;
                    }

                    gridmaps[k].energy_max = max(gridmaps[k].energy_max, gridmaps[k].energy);
                    gridmaps[k].energy_min = min(gridmaps[k].energy_min, gridmaps[k].energy);
                }
                if (floating_grid_fileptr)
                    if ((!problem_wrt) && (fprintf(floating_grid_fileptr, "%.3f\n", (float)round3dp(r_min)) < 0))
                        problem_wrt = TRUE;
                ctr++;
            }                   // icoord[X] loop
        }                       // icoord[Y] loop

        if (problem_wrt)
            logFile.printError(WARNING, "Problems writing grid maps - there may not be enough disk space.\n");
        grd_end = times(&tms_grd_end);
        ++nDone;
        logFile.printFormatted(" %6d   %8.3lf   %5.1lf%%   ", icoord[Z], cgridmin[Z] + c[Z], percentdone * (double)++ic);
        logFile.printTimeInHMS((grd_end - grd_start) * (n1[Z] - nDone));
        logFile.printFormatted("  ");
        logFile.printExecutionTimes(grd_start, grd_end, &tms_grd_start, &tms_grd_end);
    }                           // icoord[Z] loop
}
#pragma endregion

    if (floating_grid_filename[0])
        fclose(floating_grid_fileptr);

#if defined(BOINCCOMPOUND)
    boinc_fraction_done(0.9);
#endif

#pragma region Writing out summary
    // Print a summary of extrema-values from the atomic-affinity and
    // electrostatics grid-maps,
    logFile.printFormatted("\nGrid\tAtom\tMinimum   \tMaximum\n");
    logFile.printFormatted("Map \tType\tEnergy    \tEnergy \n");
    logFile.printFormatted("\t\t(kcal/mol)\t(kcal/mol)\n");
    logFile.printFormatted("____\t____\t_____________\t_____________\n");

    {
        int i = 0;
        for (; i < num_atom_maps; i++)
            logFile.printFormatted(" %d\t %s\t  %6.2lf\t%9.2le\n", i + 1, gridmaps[i].type, gridmaps[i].energy_min, gridmaps[i].energy_max);
        logFile.printFormatted(" %d\t %c\t  %6.2lf\t%9.2le\tElectrostatic Potential\n", num_atom_maps + 1, 'e', gridmaps[elecPE].energy_min, gridmaps[i].energy_max);
        logFile.printFormatted(" %d\t %c\t  %6.2lf\t%9.2le\tDesolvation Potential\n", num_atom_maps + 2, 'd', gridmaps[dsolvPE].energy_min, gridmaps[i + 1].energy_max);
    }
    logFile.printFormatted("\n\n * Note:  Every pairwise-atomic interaction was clamped at %.2f\n\n\n", EINTCLAMP);
#pragma endregion

#pragma region Closing all files
    // Close all files
    for (int i = 0; i < num_atom_maps + 2; i++)
        fclose(gridmaps[i].map_fileptr);

    // Free up the memory allocated to the gridmaps objects...
    delete [] gridmaps;

    fprintf(stderr, "\n%s: Successful Completion.\n", programParams.getProgramName());
    logFile.printTitled("Successful Completion.\n");

    tms tms_job_end;
    Clock job_end = times(&tms_job_end);
    logFile.printExecutionTimesInHMS(job_start, job_end, &tms_job_start, &tms_job_end);

#pragma endregion

#if defined(BOINCCOMPOUND)
    boinc_fraction_done(1);
#endif

#if defined(BOINC)
    boinc_finish(0);            // should not return
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
