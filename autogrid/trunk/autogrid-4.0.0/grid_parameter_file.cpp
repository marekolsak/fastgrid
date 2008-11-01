#include "grid_parameter_file.h"
#include "autogrid.h"
#include "utils.h"
#include "read_parameter_library.h"
#include <cstdlib>
#include <cstring>
#include <cmath>

GridParameterInfo::GridParameterInfo()
{
    memset(receptor_types, 0, sizeof(receptor_types));
    receptor_types_ct = 0;
    memset(covpos, 0, sizeof(covpos));

    solpar_q = .01097;   // unweighted value restored 3:9:05
    percentdone = 0.0;
    r_smooth = 0.;
    spacing = 0.375;     // One quarter of a C-C bond length.
    covhalfwidth = 1.0;
    covbarrier = 1000.0;

    factor = 332.0L;

    num_maps = 0;
    num_atom_maps = -1;
    floating_grid = FALSE, dddiel = FALSE, disorder_h = FALSE;
    elecPE = 0;
    dsolvPE = 0;

    outlev = -1;

    num_grid_points_per_map = INIT_NUM_GRID_PTS;
    i_smooth = 0;

    map_index = -1;
}

void read_grid_parameter_file(ProgramParameters &programParams, FILE *logFile, Linear_FE_Model &AD4, MapObject *(&gridmap), GridParameterInfo &params)
{
    // Open the grid parameter file
    FILE *GPF = stdin;
    if (programParams.gridParameterFilename[0])
        if (!(GPF = ag_fopen(programParams.gridParameterFilename, "r")))
        {
            fprintf(stderr, "\n%s: Sorry, I can't find or open Grid Parameter File \"%s\"\n", programParams.programName, programParams.gridParameterFilename);
            fprintf(stderr, "\n%s: Unsuccessful Completion.\n\n", programParams.programName);
            exit(911);
        }

    // Variables
    char GPF_line[LINE_LEN];
    int GPF_keyword = -1;
    int indcom = 0;
    FILE *receptor_fileptr;
    char message[LINE_LEN];
    char line[LINE_LEN];
    int length = LINE_LEN;
    char record[7];
    char atom_name[6];
    char temp_char = ' ';
    double q_max = -BIG, q_min = BIG;
    double q_tot = 0.0;
    int has_receptor_types_in_gpf = 0;
    int infld;
    FILE *xyz_fileptr;
    char xyz[5] = "xyz";
    char token[5];
    double temp_vol, temp_solpar;
    double diel;
    double ri;
    FILE *AVS_fld_fileptr;
    char xyz_filename[MAX_CHARS];

    // needed to make regression tests work between platforms
    Real *dummy_map;

    // array of ptrs used to parse input line
    char *receptor_atom_types[NUM_RECEPTOR_TYPES];

    // LIGAND: maximum is MAX_MAPS
    // each type is now at most two characters plus '\0'
    // currently ligand_atom_types is sparse... some types are not set
    char ligand_types[MAX_MAPS][3] = {0};

    // array of ptrs used to parse input line
    char *ligand_atom_types[MAX_MAPS];

    // number of different receptor atom types declared on params.receptor_types line in GPF
    int receptor_types_gpf_ct = 0;

    // XYZ
    double cext[XYZ];
    double cgridmax[XYZ];
    double cmax[XYZ] = {-BIG, -BIG, -BIG};
    double cmin[XYZ] = {BIG, BIG, BIG};
    double csum[XYZ] = {0, 0, 0};
    double cmean[XYZ];

    // array of numbers of each type
    // NB: this is a sparse int array, some entries are 0
    int receptor_atom_type_count[NUM_RECEPTOR_TYPES] = {0};

    // Read in the grid parameter file...
    ParameterEntry thisparm;
    ParameterEntry *found_parm;
    char FN_parameter_library[MAX_CHARS];   // the AD4 parameters .dat file name
    int parameter_library_found = 0;

    while (fgets(GPF_line, LINE_LEN, GPF) != 0)
    {
        GPF_keyword = gpfparser(GPF_line);

        // This first "switch" figures out how to echo the current GPF line.
        switch (GPF_keyword)
        {
        case -1:
            fprintf(logFile, "GPF> %s", GPF_line);
            print_error(programParams.programName, logFile, WARNING, "Unrecognized keyword in grid parameter file.\n");
            continue;           // while fgets GPF_line...

        case GPF_NULL:
        case GPF_COMMENT:
            fprintf(logFile, "GPF> %s", GPF_line);
            fflush(logFile);
            break;

        default:
            fprintf(logFile, "GPF> %s", GPF_line);
            indcom = strindex(GPF_line, "#");
            if (indcom != -1)
                GPF_line[indcom] = '\0';    // Truncate str. at the comment
            fflush(logFile);
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
                sscanf(GPF_line, "%*s %s", params.receptor_filename);
                fprintf(logFile, "\nReceptor Input File :\t%s\n\nReceptor Atom Type Assignments:\n\n", params.receptor_filename);

                // try to open receptor file
                if ((receptor_fileptr = ag_fopen(params.receptor_filename, "r")) == 0)
                {
                    sprintf(message, "can't find or open receptor PDBQT file \"%s\".\n", params.receptor_filename);
                    print_error(programParams.programName, logFile, ERROR, message);
                    print_error(programParams.programName, logFile, FATAL_ERROR, "Unsuccessful completion.\n\n");
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
                        fprintf(logFile, "Atom no. %2d, \"%s\"", ia + 1, atom_name);
                        fflush(logFile);

                        // Read in this receptor atom's coordinates,partial charges, and solvation parameters in PDBQS format...
                        sscanf(&line[30], "%lf", &params.coord[ia][X]);
                        sscanf(&line[38], "%lf", &params.coord[ia][Y]);
                        sscanf(&line[46], "%lf", &params.coord[ia][Z]);

                        // Output the coordinates of this atom...
                        fprintf(logFile, " at (%.3lf, %.3lf, %.3lf), ", params.coord[ia][X], params.coord[ia][Y], params.coord[ia][Z]);
                        fflush(logFile);

                        // 1:CHANGE HERE: need to set up params.vol and params.solpar
                        sscanf(&line[70], "%lf", &params.charge[ia]);
                        // printf("new type is: %s\n", &line[77]);
                        sscanf(&line[77], "%s", thisparm.autogrid_type);
                        found_parm = apm_find(thisparm.autogrid_type);
                        if (found_parm != 0)
                        {
                            fprintf(logFile, "DEBUG: found_parm->rec_index = %d", found_parm->rec_index);
                            if (found_parm->rec_index < 0)
                            {
                                strcpy(params.receptor_types[params.receptor_types_ct], found_parm->autogrid_type);
                                found_parm->rec_index = params.receptor_types_ct++;
                                fprintf(logFile, "DEBUG: found_parm->rec_index => %d", found_parm->rec_index);
                            }
                            params.atom_type[ia] = found_parm->rec_index;
                            params.solpar[ia] = found_parm->solpar;
                            params.vol[ia] = found_parm->vol;
                            params.hbond[ia] = found_parm->hbond;  // NON=0, DS,D1, AS, A1, A2
#if defined(DEBUG)
                            printf("%d:key=%s, type=%d,params.solpar=%f,params.vol=%f\n", ia, thisparm.autogrid_type, params.atom_type[ia], params.solpar[ia], params.vol[ia]);
#endif
                            ++receptor_atom_type_count[found_parm->rec_index];
                        }
                        else
                        {
                            fprintf(logFile, "\n\nreceptor file contains unknown type: '%s'\nadd parameters for it to the parameter library first\n", thisparm.autogrid_type);
                            exit(-1);
                        }

                        // if from pdbqs: convert cal/molA**3 to kcal/molA**3
                        // params.solpar[ia] *= 0.001;
                        q_max = max(q_max, params.charge[ia]);
                        q_min = min(q_min, params.charge[ia]);

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
                            /* Assume this is the 'mangled' name of a params.hydrogen atom, after the atom name has been changed from 'HD21' to '1HD2' for example. [0-9]H\(.\)\(.\) 0 1 2 3 : : :
                               : V V V V tmp 0 1 2 tmp : V 0 1 2 3 : : : : V V V V H\(.\)\(.\)[0-9] */
                            temp_char = atom_name[0];
                            atom_name[0] = atom_name[1];
                            atom_name[1] = atom_name[2];
                            atom_name[2] = atom_name[3];
                            atom_name[3] = temp_char;
                        }

                        // Tell the user what you thought this atom was...
                        fprintf(logFile, " was assigned atom type \"%s\" (rec_index= %d, params.atom_type= %d).\n", found_parm->autogrid_type, found_parm->rec_index, params.atom_type[ia]);
                        fflush(logFile);

                        // Count the number of each atom type
                        // ++receptor_atom_type_count[ params.atom_type[ia] ];

                        // Keep track of the extents of the receptor
                        for (int i = 0; i < XYZ; i++)
                        {
                            cmax[i] = max(cmax[i], params.coord[ia][i]);
                            cmin[i] = min(cmin[i], params.coord[ia][i]);
                            csum[i] += params.coord[ia][i];
                        }
                        // Total up the partial charges as we go...
                        q_tot += params.charge[ia];

                        // Increment the atom counter
                        ia++;

                        // Check that there aren't too many atoms...
                        if (ia > AG_MAX_ATOMS)
                        {
                            sprintf(message, "Too many atoms in receptor PDBQT file %s;", params.receptor_filename);
                            print_error(programParams.programName, logFile, ERROR, message);
                            sprintf(message, "-- the maximum number of atoms, AG_MAX_ATOMS, allowed is %d.", AG_MAX_ATOMS);
                            print_error(programParams.programName, logFile, ERROR, message);
                            sprintf(message, "Increase the value in the \"#define AG_MAX_ATOMS %d\" line", AG_MAX_ATOMS);
                            print_error(programParams.programName, logFile, SUGGESTION, message);
                            print_error(programParams.programName, logFile, SUGGESTION, "in the source file \"autogrid.h\", and re-compile AutoGrid.");
                            fflush(logFile);
                            // FATAL_ERROR will cause AutoGrid to exit...
                            print_error(programParams.programName, logFile, FATAL_ERROR, "Sorry, AutoGrid cannot continue.");
                        }           // endif
                    }               // endif
                }                   // endwhile
                // Finished reading in the lines of the receptor file
                fclose(receptor_fileptr);
                if (has_receptor_types_in_gpf == 1)
                    // Check that the number of atom types found in the receptor PDBQT
                    // file match the number parsed in by the "params.receptor_types" command
                    // in the GPF; if they do not match, exit!
                    if (params.receptor_types_ct != receptor_types_gpf_ct)
                    {
                        sprintf(message,
                                      "The number of atom types found in the receptor PDBQT (%d) does not match the number specified by the \"params.receptor_types\" command (%d) in the GPF!\n\n",
                                      params.receptor_types_ct, receptor_types_gpf_ct);
                        print_error(programParams.programName, logFile, ERROR, message);
                        // FATAL_ERROR will cause AutoGrid to exit...
                        print_error(programParams.programName, logFile, FATAL_ERROR, "Sorry, AutoGrid cannot continue.");
                    }

                // Update the total number of atoms in the receptor
                params.num_receptor_atoms = ia;
                fprintf(logFile, "\nMaximum partial atomic params.charge found = %+.3lf e\n", q_max);
                fprintf(logFile, "Minimum partial atomic params.charge found = %+.3lf e\n\n", q_min);
                fflush(logFile);
                // Check there are partial charges...
                if (q_max == 0. && q_min == 0.)
                {
                    sprintf(message, "No partial atomic charges were found in the receptor PDBQT file %s!\n\n", params.receptor_filename);
                    print_error(programParams.programName, logFile, ERROR, message);
                    // FATAL_ERROR will cause AutoGrid to exit...
                    print_error(programParams.programName, logFile, FATAL_ERROR, "Sorry, AutoGrid cannot continue.");
                }                   // if there are no charges EXIT

                for (int ia = 0; ia < params.num_receptor_atoms; ia++)
                    params.rexp[ia] = 0;
                fprintf(logFile, "Atom\tAtom\tNumber of this Type\n");
                fprintf(logFile, "Type\t ID \t in Receptor\n");
                fprintf(logFile, "____\t____\t___________________\n");
                fflush(logFile);
                // 2. CHANGE HERE: need to count number of each receptor_type
                for (int ia = 0; ia < params.receptor_types_ct; ia++)
                    if (receptor_atom_type_count[ia] != 0)
                        fprintf(logFile, " %d\t %s\t\t%6d\n", (ia), params.receptor_types[ia], receptor_atom_type_count[ia]);
                fprintf(logFile, "\nTotal number of atoms :\t\t%d atoms \n", params.num_receptor_atoms);
                fflush(logFile);
                fprintf(logFile, "Total params.charge :\t\t\t%.2lf e\n", q_tot);
                fflush(logFile);
                fprintf(logFile, "\n\nReceptor coordinates fit within the following volume:\n\n");
                fflush(logFile);
                fprintf(logFile, "                   _______(%.1lf, %.1lf, %.1lf)\n", cmax[X], cmax[Y], cmax[Z]);
                fprintf(logFile, "                  /|     /|\n");
                fprintf(logFile, "                 / |    / |\n");
                fprintf(logFile, "                /______/  |\n");
                fprintf(logFile, "                |  |___|__| Midpoint = (%.1lf, %.1lf, %.1lf)\n", (cmax[X] + cmin[X]) / 2., (cmax[Y] + cmin[Y]) / 2., (cmax[Z] + cmin[Z]) / 2.);
                fprintf(logFile, "                |  /   |  /\n");
                fprintf(logFile, "                | /    | /\n");
                fprintf(logFile, "                |/_____|/\n");
                fprintf(logFile, "(%.1lf, %.1lf, %.1lf)      \n", cmin[X], cmin[Y], cmin[Z]);
                fprintf(logFile, "\nMaximum coordinates :\t\t(%.3lf, %.3lf, %.3lf)\n", cmax[X], cmax[Y], cmax[Z]);
                fprintf(logFile, "Minimum coordinates :\t\t(%.3lf, %.3lf, %.3lf)\n\n", cmin[X], cmin[Y], cmin[Z]);
                fprintf(logFile, "\n");
                cmean[0] = csum[0] / (double)params.num_receptor_atoms;
                cmean[1] = csum[1] / (double)params.num_receptor_atoms;
                cmean[2] = csum[2] / (double)params.num_receptor_atoms;
                fflush(logFile);
            }
            break;

        case GPF_GRIDFLD:
            sscanf(GPF_line, "%*s %s", params.AVS_fld_filename);
            infld = strindex(params.AVS_fld_filename, ".fld");
            if (infld == -1)
                print_error(programParams.programName, logFile, FATAL_ERROR, "Grid data file needs the extension \".fld\" for AVS input\n\n");
            else
            {
                infld = strindex(params.AVS_fld_filename, "fld");
                strcpy(xyz_filename, params.AVS_fld_filename);
                xyz_filename[infld] = 'x';
                xyz_filename[infld + 1] = 'y';
                xyz_filename[infld + 2] = 'z';
            }
            if ((AVS_fld_fileptr = ag_fopen(params.AVS_fld_filename, "w")) == 0)
            {
                sprintf(message, "can't create grid dimensions data file %s\n", params.AVS_fld_filename);
                print_error(programParams.programName, logFile, ERROR, message);
                print_error(programParams.programName, logFile, FATAL_ERROR, "Unsuccessful completion.\n\n");
            }
            else
                fprintf(logFile, "\nCreating (AVS-readable) grid maps file : %s\n", params.AVS_fld_filename);
            if ((xyz_fileptr = ag_fopen(xyz_filename, "w")) == 0)
            {
                sprintf(message, "can't create grid extrema data file %s\n", xyz_filename);
                print_error(programParams.programName, logFile, ERROR, message);
                sprintf(message, "SORRY!    unable to create the \".xyz\" file.\n\n");
                print_error(programParams.programName, logFile, ERROR, message);
                print_error(programParams.programName, logFile, FATAL_ERROR, "Unsuccessful completion.\n\n");
            }
            else
                fprintf(logFile, "\nCreating (AVS-readable) grid-coordinates extrema file : %s\n\n", xyz_filename);
            fflush(logFile);
            break;

        case GPF_NPTS:
            sscanf(GPF_line, "%*s %d %d %d", &params.nelements[X], &params.nelements[Y], &params.nelements[Z]);
            for (int i = 0; i < XYZ; i++)
            {
                params.nelements[i] = check_size(params.nelements[i], xyz[i], programParams.programName, logFile);
                params.ne[i] = params.nelements[i] / 2;
                params.n1[i] = params.nelements[i] + 1;
            }
            fprintf(logFile, "\n");
            fprintf(logFile, "Number of grid points in x-direction:\t%d\n", params.n1[X]);
            fprintf(logFile, "Number of grid points in y-direction:\t%d\n", params.n1[Y]);
            fprintf(logFile, "Number of grid points in z-direction:\t%d\n", params.n1[Z]);
            fprintf(logFile, "\n");
            params.num_grid_points_per_map = params.n1[X] * params.n1[Y] * params.n1[Z];
            params.percentdone = 100. / (double)params.n1[Z];
            fflush(logFile);
            break;

        case GPF_SPACING:
            sscanf(GPF_line, "%*s %lf", &params.spacing);
            fprintf(logFile, "Grid Spacing :\t\t\t%.3lf Angstrom\n", params.spacing);
            fprintf(logFile, "\n");
            fflush(logFile);
            break;

        case GPF_GRIDCENTER:
            sscanf(GPF_line, "%*s %s", token);
            if (equal(token, "auto", 4))
            {
                for (int i = 0; i < XYZ; i++)
                    params.center[i] = cmean[i];
                fprintf(logFile, "Grid maps will be centered on the params.center of mass.\n");
                fprintf(logFile, "Coordinates of params.center of mass : (%.3lf, %.3lf, %.3lf)\n", params.center[X], params.center[Y], params.center[Z]);
            }
            else
            {
                sscanf(GPF_line, "%*s %lf %lf %lf", &params.center[X], &params.center[Y], &params.center[Z]);
                fprintf(logFile, "\nGrid maps will be centered on user-defined coordinates:\n\n\t\t(%.3lf, %.3lf, %.3lf)\n", params.center[X], params.center[Y], params.center[Z]);
            }
            // centering stuff...
            for (int ia = 0; ia < params.num_receptor_atoms; ia++)
                for (int i = 0; i < XYZ; i++)
                    params.coord[ia][i] -= params.center[i];  // transform to params.center of gridmaps
            for (int i = 0; i < XYZ; i++)
            {
                cext[i] = params.spacing * (double)params.ne[i];
                cgridmax[i] = params.center[i] + cext[i];
                params.cgridmin[i] = params.center[i] - cext[i];
            }
            fprintf(logFile, "\nGrid maps will cover the following volume:\n\n");
            fprintf(logFile, "                   _______(%.1lf, %.1lf, %.1lf)\n", cgridmax[X], cgridmax[Y], cgridmax[Z]);
            fprintf(logFile, "                  /|     /|\n");
            fprintf(logFile, "                 / |    / |\n");
            fprintf(logFile, "                /______/  |\n");
            fprintf(logFile, "                |  |___|__| Midpoint = (%.1lf, %.1lf, %.1lf)\n", params.center[X], params.center[Y], params.center[Z]);
            fprintf(logFile, "                |  /   |  /\n");
            fprintf(logFile, "                | /    | /\n");
            fprintf(logFile, "                |/_____|/\n");
            fprintf(logFile, "(%.1lf, %.1lf, %.1lf)      \n\n", params.cgridmin[X], params.cgridmin[Y], params.cgridmin[Z]);
            for (int i = 0; i < XYZ; i++)
                fprintf(logFile, "Grid map %c-dimension :\t\t%.1lf Angstroms\n", xyz[i], 2. * cext[i]);
            fprintf(logFile, "\nMaximum coordinates :\t\t(%.3lf, %.3lf, %.3lf)\n", cgridmax[X], cgridmax[Y], cgridmax[Z]);
            fprintf(logFile, "Minimum coordinates :\t\t(%.3lf, %.3lf, %.3lf)\n\n", params.cgridmin[X], params.cgridmin[Y], params.cgridmin[Z]);
            for (int i = 0; i < XYZ; i++)
                fprintf(xyz_fileptr, "%.3lf %.3lf\n", params.cgridmin[i], cgridmax[i]);
            fclose(xyz_fileptr);
            fflush(logFile);
            break;

        case GPF_LIGAND_TYPES:
            // Read in the list of atom types in the ligand.
            // GPF_line e.g.: "ligand_types N O A C HH NH"
            params.num_atom_maps = parsetypes(GPF_line, ligand_atom_types, MAX_ATOM_TYPES);
            for (int i = 0; i < params.num_atom_maps; i++)
            {
                strcpy(ligand_types[i], ligand_atom_types[i]);
#if defined(DEBUG)
                fprintf(logFile, "%d %s ->%s\n", i, ligand_atom_types[i], ligand_types[i]);
#endif
            }
            for (int i = 0; i < params.num_atom_maps; i++)
            {
                found_parm = apm_find(ligand_types[i]);
                if (found_parm != 0)
                {
                    found_parm->map_index = i;
#if defined(DEBUG)
                    fprintf(logFile, "found ligand type: %-6s%2d\n", found_parm->autogrid_type, found_parm->map_index);
#endif
                }
                else
                {
                    // return error here
                    fprintf(logFile, "unknown ligand atom type %s\nadd parameters for it to the parameter library first!\n", ligand_atom_types[i]);
                    exit(-1);
                }
            }

            params.elecPE = params.num_atom_maps;
            params.dsolvPE = params.elecPE + 1;

            /* params.num_maps is the number of maps to be created: the number of ligand atom types, plus 1 for the electrostatic map. AutoDock can only read in MAX_MAPS maps, which must include
               the ligand atom maps and electrostatic map */
            params.num_maps = params.num_atom_maps + 2;

            // Check to see if there is enough memory to store these map objects
            gridmap = (MapObject *) malloc(sizeof(MapObject) * params.num_maps);

            if (gridmap == 0)
            {
                print_error(programParams.programName, logFile, ERROR, "Could not allocate memory to create the MapObject \"gridmap\".\n");
                print_error(programParams.programName, logFile, FATAL_ERROR, "Unsuccessful completion.\n\n");
            }

            // Initialize the gridmap MapObject
            for (int i = 0; i < params.num_maps; i++)
            {
                gridmap[i].atom_type = 0;   // corresponds to receptor numbers????
                gridmap[i].map_index = 0;
                gridmap[i].is_covalent = 0;
                gridmap[i].is_hbonder = 0;
                gridmap[i].map_fileptr = (FILE *) 0;
                strcpy(gridmap[i].map_filename, "");
                strcpy(gridmap[i].type, "");    // eg HD or OA or NA or N
                gridmap[i].constant = 0.0L; // this will become obsolete
                gridmap[i].energy_max = 0.0L;
                gridmap[i].energy_min = 0.0L;
                gridmap[i].energy = 0.0L;
                gridmap[i].vol_probe = 0.0L;
                gridmap[i].solpar_probe = 0.0L;
                gridmap[i].Rij = 0.0L;
                gridmap[i].epsij = 0.0L;
                gridmap[i].hbond = NON; // hbonding character:
                gridmap[i].Rij_hb = 0.0L;
                gridmap[i].epsij_hb = 0.0L;
                // per gridmap[i].receptor type parameters, ordered as in params.receptor_types
                for (int j = 0; j < NUM_RECEPTOR_TYPES; j++)
                {
                    gridmap[i].nbp_r[j] = 0.0L; // radius of energy-well minimum
                    gridmap[i].nbp_eps[j] = 0.0L;   // depth of energy-well minimum
                    gridmap[i].xA[j] = 0;   // generally 12
                    gridmap[i].xB[j] = 0;   // 6 for non-hbonders 10 for h-bonders
                    gridmap[i].hbonder[j] = 0;
                }               // j
            }                   // i

            // Check to see if the number of grid points requested will be feasible; give warning if not enough memory.
            if (params.num_grid_points_per_map != INIT_NUM_GRID_PTS)
            {
                dummy_map = (Real *) malloc(sizeof(Real) * (params.num_maps * params.num_grid_points_per_map));
                if (!dummy_map)
                {
                    // Too many maps requested
                    sprintf(message,
                                  "There will not be enough memory to store these grid maps in AutoDock; \ntry reducing the number of ligand atom types (you have %d including electrostatics) \nor reducing the size of the grid maps (you asked for %d x %d x %d grid points); \n or try running AutoDock on a machine with more RAM than this one.\n",
                                  params.num_maps, params.n1[X], params.n1[Y], params.n1[Z]);
                    print_error(programParams.programName, logFile, WARNING, message);
                }
                else
                    // free up this memory right away; we were just testing to see if we had enough when we try to run AutoDock
                    free(dummy_map);
            }
            else
            {
                print_error(programParams.programName, logFile, ERROR, "You need to set the number of grid points using \"npts\" before setting the ligand atom types, using \"ligand_types\".\n");
                print_error(programParams.programName, logFile, FATAL_ERROR, "Unsuccessful completion.\n\n");
            }                   // ZZZZZZZZZZZZZZZZZ
            if (!gridmap)
            {
                sprintf(message, "Too many ligand atom types; there is not enough memory to create these maps.  Try using fewer atom types than %d.\n", params.num_atom_maps);
                print_error(programParams.programName, logFile, ERROR, message);
                print_error(programParams.programName, logFile, FATAL_ERROR, "Unsuccessful completion.\n\n");
            }

            for (int i = 0; i < params.num_atom_maps; i++)
            {
                gridmap[i].is_covalent = FALSE;
                gridmap[i].is_hbonder = FALSE;
                gridmap[i].map_index = i;
                strcpy(gridmap[i].type, ligand_types[i]);   // eg HD or OA or NA or N
                found_parm = apm_find(ligand_types[i]);
                gridmap[i].atom_type = found_parm->map_index;
                gridmap[i].solpar_probe = found_parm->solpar;
                gridmap[i].vol_probe = found_parm->vol;
                gridmap[i].Rij = found_parm->Rij;
                gridmap[i].epsij = found_parm->epsij;
                gridmap[i].hbond = found_parm->hbond;
                gridmap[i].Rij_hb = found_parm->Rij_hb;
                gridmap[i].epsij_hb = found_parm->epsij_hb;
                if (gridmap[i].hbond > 0)
                    gridmap[i].is_hbonder = TRUE;

#if defined(DEBUG)
                fprintf(logFile, " setting ij parms for map %d \n", i);
                fprintf(logFile, "for gridmap[%d], type->%s,Rij->%6.4f, epsij->%6.4f, hbond->%d\n", i, found_parm->autogrid_type, gridmap[i].Rij, gridmap[i].epsij, gridmap[i].hbond);
#endif
                for (int j = 0; j < params.receptor_types_ct; j++)
                {
                    found_parm = apm_find(params.receptor_types[j]);
                    gridmap[i].nbp_r[j] = (gridmap[i].Rij + found_parm->Rij) / 2.;
                    gridmap[i].nbp_eps[j] = sqrt(gridmap[i].epsij * found_parm->epsij);
                    // apply the vdW forcefield parameter/weight here
                    // This was removed because "setup_p_l" does this for us... gridmap[i].nbp_eps[j] *= FE_coeff_vdW;
                    gridmap[i].xA[j] = 12;
                    // setup hbond dependent stuff
                    gridmap[i].xB[j] = 6;
                    gridmap[i].hbonder[j] = 0;
                    if ((int)(gridmap[i].hbond) > 2 && ((int)found_parm->hbond == 1 || (int)found_parm->hbond == 2))
                    {           // AS,A1,A2 map vs DS,D1 probe
                        gridmap[i].xB[j] = 10;
                        gridmap[i].hbonder[j] = 1;
                        gridmap[i].is_hbonder = TRUE;
                        // Rij and epsij for this hb interaction in parm_data.dat file as Rii and epsii for heavy atom hb factors
                        gridmap[i].nbp_r[j] = gridmap[i].Rij_hb;
                        gridmap[i].nbp_eps[j] = gridmap[i].epsij_hb;

                        // apply the hbond forcefield parameter/weight here
                        // This was removed because "setup_p_l" does this for us... gridmap[i].nbp_eps[j] *= FE_coeff_hbond;
#if defined(DEBUG)
                        fprintf(logFile, "set %d-%d hb eps to %6.4f*%6.4f=%6.4f\n", i, j, gridmap[i].epsij_hb, found_parm->epsij_hb, gridmap[i].nbp_eps[j]);
#endif
                    }
                    else if (((int)gridmap[i].hbond == 1 || (int)gridmap[i].hbond == 2) && ((int)found_parm->hbond > 2))
                    {           // DS,D1 map vs AS,A1,A2 probe
                        gridmap[i].xB[j] = 10;
                        gridmap[i].hbonder[j] = 1;
                        gridmap[i].is_hbonder = TRUE;
                        // Rij and epsij for this hb interaction in parm_data.dat file as Rii and epsii for heavy atom hb factors
                        gridmap[i].nbp_r[j] = found_parm->Rij_hb;
                        gridmap[i].nbp_eps[j] = found_parm->epsij_hb;

                        // apply the hbond forcefield parameter/weight here
                        // This was removed because "setup_p_l" does this for us... gridmap[i].nbp_eps[j] *= FE_coeff_hbond;
#if defined(DEBUG)
                        fprintf(logFile, "2: set %d-%d hb eps to %6.4f*%6.4f=%6.4f\n", i, j, gridmap[i].epsij_hb, found_parm->epsij_hb, gridmap[i].nbp_eps[j]);
#endif
                    }
#if defined(DEBUG)
                    fprintf(logFile, "vs receptor_type[%d]:type->%s, hbond->%d ", j, found_parm->autogrid_type, (int)found_parm->hbond);
                    fprintf(logFile, "nbp_r->%6.4f, nbp_eps->%6.4f,xB=%d,hbonder=%d\n", gridmap[i].nbp_r[j], gridmap[i].nbp_eps[j], gridmap[i].xB[j], gridmap[i].hbonder[j]);
#endif
                }              // initialize energy parms for each possible receptor type
            }                   // for each map
            fprintf(logFile, "\nAtom type names for ligand atom types 1-%d used for ligand-atom affinity grid maps:\n\n", params.num_atom_maps);
            for (int i = 0; i < params.num_atom_maps; i++)
            {
                fprintf(logFile, "\t\t\tAtom type number %d corresponds to atom type name \"%s\".\n", gridmap[i].map_index, gridmap[i].type);

                // FIX THIS!!! Covalent Atom Types are not yet supported with the new AG4/AD4 atom typing mechanism...
                /* if (gridmap[i].atom_type == COVALENTTYPE) { gridmap[i].is_covalent = TRUE;  fprintf(logFile, "\nAtom type number %d will be used to calculate a covalent affinity
                   grid map\n\n", i + 1); } */
            }
            fprintf(logFile, "\n\n");
            fflush(logFile);
            break;

        case GPF_RECEPTOR_TYPES:
            // Read in the list of atom types in the receptor.
            // GPF_line e.g.: "params.receptor_types N O A C HH NH"
            //
            // NOTE: This line is not guaranteed to match the actual
            // atom types present in the receptor PDBQT file
            // specified by the "receptor" command.
            params.receptor_types_ct = parsetypes(GPF_line, receptor_atom_types, MAX_ATOM_TYPES);
            receptor_types_gpf_ct = params.receptor_types_ct;
            has_receptor_types_in_gpf = 1;
#if defined(DEBUG)
            printf("receptor_types_gpf_ct=%d\n", receptor_types_gpf_ct);
            printf("params.receptor_types_ct=%d\n", params.receptor_types_ct);
#endif
            for (int i = 0; i < params.receptor_types_ct; i++)
            {
                strcpy(params.receptor_types[i], receptor_atom_types[i]);
#if defined(DEBUG)
                printf("%d %s  ->%s\n", i, receptor_atom_types[i], params.receptor_types[i]);
#endif
            }
            for (int i = 0; i < params.receptor_types_ct; i++)
            {
                found_parm = apm_find(receptor_atom_types[i]);
                if (found_parm != 0)
                    found_parm->rec_index = i;
                else
                {
                    sprintf(message, "Unknown receptor type: \"%s\"\n -- Add parameters for it to the parameter library first!\n", receptor_atom_types[i]);
                    print_error(programParams.programName, logFile, ERROR, message);
                    print_error(programParams.programName, logFile, FATAL_ERROR, "Unsuccessful completion.\n\n");
                }
            }
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
#if defined(DEBUG)
            printf("assigned receptor types:params.arom_carbon->%d, params.hydrogen->%d,params.nonHB_hydrogen->%d, params.carbon->%d, params.oxygen->%d, params.nitrogen->%d\n, params.nonHB_nitrogen->%d, params.sulphur->%d, params.nonHB_sulphur->%d\n",
                   params.arom_carbon, params.hydrogen, params.nonHB_hydrogen, params.carbon, params.oxygen, params.nitrogen, params.nonHB_nitrogen, params.sulphur, params.nonHB_sulphur);
#endif
            fflush(logFile);
            break;

        case GPF_SOL_PAR:      // THIS IS OBSOLETE!!!
            // Read volume and solvation parameter for probe:
            sscanf(GPF_line, "%*s %s %lf %lf", thisparm.autogrid_type, &temp_vol, &temp_solpar);
            found_parm = apm_find(thisparm.autogrid_type);
            if (found_parm != 0)
            {
                found_parm->vol = temp_vol;
                found_parm->solpar = temp_solpar;
                int i = found_parm->map_index;
                if (i >= 0)
                {
                    // DON'T!!!
                    // convert cal/molA^3 to kcal/molA^3
                    // gridmap[i].solpar_probe = temp_solpar * 0.001;
                    gridmap[i].solpar_probe = temp_solpar;
                    fprintf(logFile, "\nProbe %s solvation parameters: \n\n\tatomic fragmental volume: %.2f A^3\n\tatomic solvation parameter: %.4f cal/mol A^3\n\n",
                                  found_parm->autogrid_type, found_parm->vol, found_parm->solpar);
                    fflush(logFile);
                }
            }
            else
                fprintf(logFile, "%s key not found\n", thisparm.autogrid_type);
            fflush(logFile);
            break;              // end solvation parameter

        /*case GPF_CONSTANT:
            break;*/

        case GPF_MAP:
            /* The variable "params.map_index" is the 0-based index of the ligand atom type we are calculating a map for. If the "types" line was CNOSH, there would be 5 ligand atom maps to
               calculate, and since "params.map_index" is initialized to -1, params.map_index will increment each time there is a "map" keyword in the GPF.  The value of params.map_index should therefore go
               from 0 to 4 for each "map" keyword. In this example, params.num_atom_maps would be 5, and params.num_atom_maps-1 would be 4, so if params.map_index is > 4, there is something wrong in the
               number of "map" keywords. */
            ++params.map_index;
            if (params.map_index > params.num_atom_maps - 1)
            {
                sprintf(message, "Too many \"map\" keywords (%d);  the \"types\" command declares only %d maps.\nRemove a \"map\" keyword from the GPF.\n", params.map_index + 1,
                              params.num_atom_maps);
                print_error(programParams.programName, logFile, ERROR, message);
                print_error(programParams.programName, logFile, FATAL_ERROR, "Unsuccessful completion.\n\n");
            }
            // Read in the filename for this grid map *//* GPF_MAP
            sscanf(GPF_line, "%*s %s", gridmap[params.map_index].map_filename);
            if ((gridmap[params.map_index].map_fileptr = ag_fopen(gridmap[params.map_index].map_filename, "w")) == 0)
            {
                sprintf(message, "Cannot open grid map \"%s\" for writing.", gridmap[params.map_index].map_filename);
                print_error(programParams.programName, logFile, ERROR, message);
                print_error(programParams.programName, logFile, FATAL_ERROR, "Unsuccessful completion.\n\n");
            }
            fprintf(logFile, "\nOutput Grid Map %d:   %s\n\n", (params.map_index + 1), gridmap[params.map_index].map_filename);
            fflush(logFile);

            break;

        case GPF_ELECMAP:
            sscanf(GPF_line, "%*s %s", gridmap[params.elecPE].map_filename);
            if ((gridmap[params.elecPE].map_fileptr = ag_fopen(gridmap[params.elecPE].map_filename, "w")) == 0)
            {
                sprintf(message, "can't open grid map \"%s\" for writing.\n", gridmap[params.elecPE].map_filename);
                print_error(programParams.programName, logFile, ERROR, message);
                print_error(programParams.programName, logFile, FATAL_ERROR, "Unsuccessful completion.\n\n");
            }
            fprintf(logFile, "\nOutput Electrostatic Potential Energy Grid Map: %s\n\n", gridmap[params.elecPE].map_filename);
            break;

        case GPF_DSOLVMAP:
            sscanf(GPF_line, "%*s %s", gridmap[params.dsolvPE].map_filename);
            if ((gridmap[params.dsolvPE].map_fileptr = ag_fopen(gridmap[params.dsolvPE].map_filename, "w")) == 0)
            {
                sprintf(message, "can't open grid map \"%s\" for writing.\n", gridmap[params.dsolvPE].map_filename);
                print_error(programParams.programName, logFile, ERROR, message);
                print_error(programParams.programName, logFile, FATAL_ERROR, "Unsuccessful completion.\n\n");
            }
            fprintf(logFile, "\nOutput Desolvation Free Energy Grid Map: %s\n\n", gridmap[params.dsolvPE].map_filename);
            break;

        case GPF_COVALENTMAP:
            sscanf(GPF_line, "%*s %lf %lf %lf %lf %lf", &params.covhalfwidth, &params.covbarrier, &(params.covpos[X]), &(params.covpos[Y]), &(params.covpos[Z]));
            fprintf(logFile, "\ncovalentmap <half-width in Angstroms> <barrier> <x> <y> <z>\n");
            fprintf(logFile, "\nCovalent well's half-width in Angstroms:         %8.3f\n", params.covhalfwidth);
            fprintf(logFile, "\nCovalent barrier energy in kcal/mol:             %8.3f\n", params.covbarrier);
            fprintf(logFile, "\nCovalent attachment point will be positioned at: (%8.3f, %8.3f, %8.3f)\n\n", params.covpos[X], params.covpos[Y], params.covpos[Z]);
            for (int i = 0; i < XYZ; i++)
                // params.center params.covpos in the grid maps frame of reference,
                params.covpos[i] -= params.center[i];
            break;

        case GPF_DISORDER:
            params.disorder_h = TRUE;
            fprintf(logFile, "\nHydroxyls will be disordered \n\n");
            break;

        case GPF_SMOOTH:
            sscanf(GPF_line, "%*s %lf", &params.r_smooth);
            fprintf(logFile, "\nPotentials will be smoothed by: %.3lf Angstrom\n\n", params.r_smooth);
            // Angstrom is divided by A_DIVISOR in look-up table.
            // Typical value of params.r_smooth is 0.5 Angstroms
            // so params.i_smooth = 0.5 * 100. / 2 = 25
            params.i_smooth = (int)(params.r_smooth * A_DIVISOR / 2.);
            break;

        case GPF_QASP:
            sscanf(GPF_line, "%*s %lf", &params.solpar_q);
            fprintf(logFile, "\nCharge component of the atomic solvation parameter: %.3lf\n\n", params.solpar_q);
            // Typical value of params.solpar_q is 0.001118
            break;

        case GPF_DIEL:
            sscanf(GPF_line, "%*s %lf", &diel);
            if (diel < 0.)
            {
                // negative...
                params.dddiel = TRUE;
                // calculate ddd of Mehler & Solmajer
                fprintf(logFile, "\nUsing *distance-dependent* dielectric function of Mehler and Solmajer, Prot.Eng.4, 903-910.\n\n");
                params.epsilon[0] = 1.0;
                for (int indx_r = 1; indx_r < MAX_DIST; indx_r++)
                    params.epsilon[indx_r] = calc_ddd_Mehler_Solmajer(angstrom(indx_r), APPROX_ZERO);
                fprintf(logFile, "  d   Dielectric\n ___  __________\n");
                for (int i = 0; i <= 500; i += 10)
                {
                    ri = angstrom(i);
                    fprintf(logFile, "%4.1lf%9.2lf\n", ri, params.epsilon[i]);
                }
                fprintf(logFile, "\n");
                // convert params.epsilon to 1 / params.epsilon
                for (int i = 1; i < MAX_DIST; i++)
                    params.epsilon[i] = params.factor / params.epsilon[i];
            }
            else
            {
                // positive or zero...
                params.dddiel = FALSE;
                if (diel <= APPROX_ZERO)
                    diel = 40.;
                fprintf(logFile, "Using a *constant* dielectric of:  %.2f\n", diel);
                params.invdielcal = params.factor / diel;
            }
            break;

        case GPF_FMAP:
            sscanf(GPF_line, "%*s %s", params.floating_grid_filename);
            fprintf(logFile, "\nFloating Grid file name = %s\n", params.floating_grid_filename);
            ++params.num_maps;
            params.floating_grid = TRUE;
            break;

        case GPF_PARAM_FILE:
            // open and read the AD4 parameters .dat file
            parameter_library_found = sscanf(GPF_line, "%*s %s ", FN_parameter_library);
            read_parameter_library(FN_parameter_library, params.outlev, programParams.programName, programParams.debug, logFile, AD4);
            break;
        }                       // second switch
    }                           // while

    fprintf(logFile, "\n>>> Closing the grid parameter file (GPF)... <<<\n\n");
    fprintf(logFile, UnderLine);
    fclose(GPF);

    if (!params.floating_grid)
        fprintf(logFile, "\n\nNo Floating Grid was requested.\n");

    fprintf(AVS_fld_fileptr, "# AVS field file\n#\n");
    fprintf(AVS_fld_fileptr, "# AutoDock Atomic Affinity and Electrostatic Grids\n#\n");
    fprintf(AVS_fld_fileptr, "# Created by %s.\n#\n", programParams.programName);
    fprintf(AVS_fld_fileptr, "#SPACING %.3f\n", (float)params.spacing);
    fprintf(AVS_fld_fileptr, "#NELEMENTS %d %d %d\n", params.nelements[X], params.nelements[Y], params.nelements[Z]);
    fprintf(AVS_fld_fileptr, "#CENTER %.3lf %.3lf %.3lf\n", params.center[X], params.center[Y], params.center[Z]);
    fprintf(AVS_fld_fileptr, "#MACROMOLECULE %s\n", params.receptor_filename);
    fprintf(AVS_fld_fileptr, "#GRID_PARAMETER_FILE %s\n#\n", programParams.gridParameterFilename);
    fprintf(AVS_fld_fileptr, "ndim=3\t\t\t# number of dimensions in the field\n");
    fprintf(AVS_fld_fileptr, "dim1=%d\t\t\t# number of x-elements\n", params.n1[X]);
    fprintf(AVS_fld_fileptr, "dim2=%d\t\t\t# number of y-elements\n", params.n1[Y]);
    fprintf(AVS_fld_fileptr, "dim3=%d\t\t\t# number of z-elements\n", params.n1[Z]);
    fprintf(AVS_fld_fileptr, "nspace=3\t\t# number of physical coordinates per point\n");
    fprintf(AVS_fld_fileptr, "veclen=%d\t\t# number of affinity values at each point\n", params.num_maps);
    fprintf(AVS_fld_fileptr, "data=float\t\t# data type (byte, integer, float, double)\n");
    fprintf(AVS_fld_fileptr, "field=uniform\t\t# field type (uniform, rectilinear, irregular)\n");
    for (int i = 0; i < XYZ; i++)
        fprintf(AVS_fld_fileptr, "params.coord %d file=%s filetype=ascii offset=%d\n", (i + 1), xyz_filename, (i * 2));
    for (int i = 0; i < params.num_atom_maps; i++)
        fprintf(AVS_fld_fileptr, "label=%s-affinity\t# component label for variable %d\n", gridmap[i].type, (i + 1));                           // i
    fprintf(AVS_fld_fileptr, "label=Electrostatics\t# component label for variable %d\n", params.num_maps - 2);
    fprintf(AVS_fld_fileptr, "label=Desolvation\t# component label for variable %d\n", params.num_maps - 1);
    if (params.floating_grid)
        fprintf(AVS_fld_fileptr, "label=Floating_Grid\t# component label for variable %d\n", params.num_maps);
    fprintf(AVS_fld_fileptr, "#\n# location of affinity grid files and how to read them\n#\n");
    for (int i = 0; i < params.num_atom_maps; i++)
        fprintf(AVS_fld_fileptr, "variable %d file=%s filetype=ascii skip=6\n", (i + 1), gridmap[i].map_filename);
    fprintf(AVS_fld_fileptr, "variable %d file=%s filetype=ascii skip=6\n", params.num_atom_maps + 1, gridmap[params.elecPE].map_filename);
    fprintf(AVS_fld_fileptr, "variable %d file=%s filetype=ascii skip=6\n", params.num_atom_maps + 2, gridmap[params.dsolvPE].map_filename);
    if (params.floating_grid)
        fprintf(AVS_fld_fileptr, "variable %d file=%s filetype=ascii skip=6\n", params.num_maps, params.floating_grid_filename);
    fclose(AVS_fld_fileptr);
}
