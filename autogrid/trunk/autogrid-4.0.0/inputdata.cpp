#include "inputdata.h"
#include "utils.h"
#include "exceptions.h"
#include "constants.h"
#include <cstring>
#include <cmath>

InputData::InputData()
{
    fldFilenameAVS[0] = 0;
    floatingGridFilename[0] = 0;
    receptorFilename[0] = 0;
    xyzFilename[0] = 0;
    parameterLibraryFilename[0] = 0;

    memset(receptorTypes, 0, sizeof(receptorTypes));
    numReceptorTypes = 0;
    numGridPointsPerMap = INIT_NUM_GRID_PTS;
    numReceptorAtoms = 0;
    memset(covpos, 0, sizeof(covpos));

    // for NEW3 desolvation terms
    solparQ = .01097;   // unweighted value restored 3:9:05
    invDielCal = 0;
    rSmooth = 0;
    spacing = 0.375;     // One quarter of a C-C bond length.
    covHalfWidth = 1.0;
    covBarrier = 1000.0;

    distDepDiel = false;
    disorderH = false;
}

void InputData::load(const char *gridParameterFilename, GridMapList &gridmaps, ParameterLibrary &parameterLibrary, LogFile &logFile)
{
    // LIGAND: maximum is MAX_MAPS
    // each type is now at most two characters plus '\0'
    // currently ligandAtomTypes is sparse... some types are not set
    char ligandTypes[MAX_MAPS][3] = {0};

    // number of different receptor atom types declared on receptorTypes line in GPF
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
    char GPFLine[LINE_LEN];
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
    if (gridParameterFilename[0])
    {
        GPF = openFile(gridParameterFilename, "r");
        if (!GPF)
        {
            logFile.printErrorFormatted(ERROR, "Sorry, I can't find or open Grid Parameter File \"%s\"", gridParameterFilename);
            logFile.printErrorFormatted(ERROR, "Unsuccessful Completion.\n");
            throw ExitProgram(911);
        }
    }

    // Read in the grid parameter file...
    ParameterEntry thisparm;
    ParameterEntry *foundParam;

    while (fgets(GPFLine, LINE_LEN, GPF) != 0)
    {
        GPF_keyword = parseGPFLine(GPFLine);

        // This first "switch" figures out how to echo the current GPF line.
        switch (GPF_keyword)
        {
        case -1:
            logFile.printFormatted("GPF> %s", GPFLine);
            logFile.printError(WARNING, "Unrecognized keyword in grid parameter file.\n");
            continue;           // while fgets GPFLine...

        case GPF_NULL:
        case GPF_COMMENT:
            logFile.printFormatted("GPF> %s", GPFLine);
            break;

        default:
            logFile.printFormatted("GPF> %s", GPFLine);
            indcom = strIndex(GPFLine, "#");
            if (indcom != -1)
                GPFLine[indcom] = '\0';    // Truncate str. at the comment
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
                // read in the receptor gridParameterFilename
                sscanf(GPFLine, "%*s %s", receptorFilename);
                logFile.printFormatted("\nReceptor Input File :\t%s\n\nReceptor Atom Type Assignments:\n\n", receptorFilename);

                // try to open receptor file
                if ((receptor_fileptr = openFile(receptorFilename, "r")) == 0)
                {
                    logFile.printErrorFormatted(ERROR, "can't find or open receptor PDBQT file \"%s\".\n", receptorFilename);
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
                        sscanf(&line[77], "%s", thisparm.autogridType);
                        foundParam = parameterLibrary.findAtomParameter(thisparm.autogridType);
                        if (foundParam != 0)
                        {
                            logFile.printFormatted("DEBUG: foundParam->recIndex = %d", foundParam->recIndex);
                            if (foundParam->recIndex < 0)
                            {
                                strcpy(receptorTypes[numReceptorTypes], foundParam->autogridType);
                                foundParam->recIndex = numReceptorTypes++;
                                logFile.printFormatted("DEBUG: foundParam->recIndex => %d", foundParam->recIndex);
                            }
                            atomType[ia] = foundParam->recIndex;
                            solpar[ia] = foundParam->solpar;
                            vol[ia] = foundParam->vol;
                            hbond[ia] = foundParam->hbond;  // NON=0, DS,D1, AS, A1, A2
                            ++receptor_atom_type_count[foundParam->recIndex];
                        }
                        else
                        {
                            logFile.printFormatted("\n\nreceptor file contains unknown type: '%s'\nadd parameters for it to the parameter library first\n", thisparm.autogridType);
                            throw ExitProgram(-1);
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
                        logFile.printFormatted(" was assigned atom type \"%s\" (recIndex= %d, atomType= %d).\n", foundParam->autogridType, foundParam->recIndex, atomType[ia]);

                        // Count the number of each atom type
                        // ++receptor_atom_type_count[ atomType[ia] ];

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
                            logFile.printErrorFormatted(ERROR, "Too many atoms in receptor PDBQT file %s;", receptorFilename);
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
                    // file match the number parsed in by the "receptorTypes" command
                    // in the GPF; if they do not match, exit!
                    if (numReceptorTypes != receptor_types_gpf_ct)
                    {
                        logFile.printErrorFormatted(ERROR,
                            "The number of atom types found in the receptor PDBQT (%d) does not match the number specified by the \"receptor_types\" command (%d) in the GPF!\n\n",
                            numReceptorTypes, receptor_types_gpf_ct);
                        // FATAL_ERROR will cause AutoGrid to exit...
                        logFile.printError(FATAL_ERROR, "Sorry, AutoGrid cannot continue.");
                    }

                // Update the total number of atoms in the receptor
                numReceptorAtoms = ia;
                logFile.printFormatted("\nMaximum partial atomic charge found = %+.3lf e\n", q_max);
                logFile.printFormatted("Minimum partial atomic charge found = %+.3lf e\n\n", q_min);

                // Check there are partial charges...
                if (q_max == 0 && q_min == 0)
                {
                    logFile.printErrorFormatted(ERROR, "No partial atomic charges were found in the receptor PDBQT file %s!\n\n", receptorFilename);
                    // FATAL_ERROR will cause AutoGrid to exit...
                    logFile.printError(FATAL_ERROR, "Sorry, AutoGrid cannot continue.");
                }                   // if there are no charges EXIT

                logFile.print("Atom\tAtom\tNumber of this Type\n"
                              "Type\t ID \t in Receptor\n"
                              "____\t____\t___________________\n");

                // 2. CHANGE HERE: need to count number of each receptor_type
                for (int ia = 0; ia < numReceptorTypes; ia++)
                    if (receptor_atom_type_count[ia] != 0)
                        logFile.printFormatted(" %d\t %s\t\t%6d\n", (ia), receptorTypes[ia], receptor_atom_type_count[ia]);

                logFile.printFormatted("\nTotal number of atoms :\t\t%d atoms \n"
                                       "Total charge :\t\t\t%.2lf e\n"
                                       "\n\nReceptor coordinates fit within the following volume:\n\n"
                                       "                   _______(%.1lf, %.1lf, %.1lf)\n"
                                       "                  /|     /|\n"
                                       "                 / |    / |\n"
                                       "                /______/  |\n"
                                       "                |  |___|__| Midpoint = (%.1lf, %.1lf, %.1lf)\n"
                                       "                |  /   |  /\n"
                                       "                | /    | /\n"
                                       "                |/_____|/\n"
                                       "(%.1lf, %.1lf, %.1lf)      \n"
                                       "\nMaximum coordinates :\t\t(%.3lf, %.3lf, %.3lf)\n"
                                       "Minimum coordinates :\t\t(%.3lf, %.3lf, %.3lf)\n\n"
                                       "\n",
                                       numReceptorAtoms, q_tot, cmax[X], cmax[Y], cmax[Z],
                                       (cmax[X] + cmin[X]) / 2, (cmax[Y] + cmin[Y]) / 2, (cmax[Z] + cmin[Z]) / 2,
                                       cmin[X], cmin[Y], cmin[Z], cmax[X], cmax[Y], cmax[Z], cmin[X], cmin[Y], cmin[Z]);

                cmean[0] = csum[0] / (double)numReceptorAtoms;
                cmean[1] = csum[1] / (double)numReceptorAtoms;
                cmean[2] = csum[2] / (double)numReceptorAtoms;
            }
            break;

        case GPF_GRIDFLD:
            sscanf(GPFLine, "%*s %s", fldFilenameAVS);
            infld = strIndex(fldFilenameAVS, ".fld");
            if (infld == -1)
                logFile.printError(FATAL_ERROR, "Grid data file needs the extension \".fld\" for AVS input\n\n");
            else
            {
                infld = strIndex(fldFilenameAVS, "fld");
                strcpy(xyzFilename, fldFilenameAVS);
                xyzFilename[infld] = 'x';
                xyzFilename[infld + 1] = 'y';
                xyzFilename[infld + 2] = 'z';
            }
            if ((xyz_fileptr = openFile(xyzFilename, "w")) == 0)
            {
                logFile.printErrorFormatted(ERROR, "can't create grid extrema data file %s\n", xyzFilename);
                logFile.printError(ERROR, "SORRY!    unable to create the \".xyz\" file.\n\n");
                logFile.printError(FATAL_ERROR, "Unsuccessful completion.\n\n");
            }
            else
                logFile.printFormatted("\nCreating (AVS-readable) grid-coordinates extrema file : %s\n\n", xyzFilename);
            break;

        case GPF_NPTS:
            sscanf(GPFLine, "%*s %d %d %d", &nelements[X], &nelements[Y], &nelements[Z]);
            for (int i = 0; i < XYZ; i++)
            {
                nelements[i] = checkSize(nelements[i], xyz[i], logFile);
                ne[i] = nelements[i] / 2;
                n1[i] = nelements[i] + 1;
            }
            logFile.printFormatted("\nNumber of grid points in x-direction:\t%d\n"
                                   "Number of grid points in y-direction:\t%d\n"
                                   "Number of grid points in z-direction:\t%d\n\n",
                                   n1[X], n1[Y], n1[Z]);

            numGridPointsPerMap = n1[X] * n1[Y] * n1[Z];
            break;

        case GPF_SPACING:
            sscanf(GPFLine, "%*s %lf", &spacing);
            logFile.printFormatted("Grid Spacing :\t\t\t%.3lf Angstrom\n\n", spacing);
            break;

        case GPF_GRIDCENTER:
            sscanf(GPFLine, "%*s %s", token);
            if (equal(token, "auto", 4))
            {
                for (int i = 0; i < XYZ; i++)
                    center[i] = cmean[i];
                logFile.printFormatted("Grid maps will be centered on the center of mass.\n"
                                       "Coordinates of center of mass : (%.3lf, %.3lf, %.3lf)\n", center[X], center[Y], center[Z]);
            }
            else
            {
                sscanf(GPFLine, "%*s %lf %lf %lf", &center[X], &center[Y], &center[Z]);
                logFile.printFormatted("\nGrid maps will be centered on user-defined coordinates:\n\n\t\t(%.3lf, %.3lf, %.3lf)\n", center[X], center[Y], center[Z]);
            }
            // centering stuff...
            for (int ia = 0; ia < numReceptorAtoms; ia++)
                for (int i = 0; i < XYZ; i++)
                    coord[ia][i] -= center[i];  // transform to center of gridmaps
            for (int i = 0; i < XYZ; i++)
            {
                cext[i] = spacing * (double)ne[i];
                cgridmax[i] = center[i] + cext[i];
                cgridmin[i] = center[i] - cext[i];
            }

            logFile.printFormatted("\nGrid maps will cover the following volume:\n\n"
                                   "                   _______(%.1lf, %.1lf, %.1lf)\n"
                                   "                  /|     /|\n"
                                   "                 / |    / |\n"
                                   "                /______/  |\n"
                                   "                |  |___|__| Midpoint = (%.1lf, %.1lf, %.1lf)\n"
                                   "                |  /   |  /\n"
                                   "                | /    | /\n"
                                   "                |/_____|/\n"
                                   "(%.1lf, %.1lf, %.1lf)      \n\n", cgridmax[X], cgridmax[Y], cgridmax[Z],
                                   center[X], center[Y], center[Z], cgridmin[X], cgridmin[Y], cgridmin[Z]);

            for (int i = 0; i < XYZ; i++)
                logFile.printFormatted("Grid map %c-dimension :\t\t%.1lf Angstroms\n", xyz[i], 2 * cext[i]);

            logFile.printFormatted("\nMaximum coordinates :\t\t(%.3lf, %.3lf, %.3lf)\n"
                                   "Minimum coordinates :\t\t(%.3lf, %.3lf, %.3lf)\n\n", cgridmax[X], cgridmax[Y], cgridmax[Z], cgridmin[X], cgridmin[Y], cgridmin[Z]);
            for (int i = 0; i < XYZ; i++)
                fprintf(xyz_fileptr, "%.3lf %.3lf\n", cgridmin[i], cgridmax[i]);
            fclose(xyz_fileptr);
            break;

        case GPF_LIGAND_TYPES:
            {
                // Read in the list of atom types in the ligand.
                // GPFLine e.g.: "ligand_types N O A C HH NH"

                // array of ptrs used to parse input line
                char *ligandAtomTypes[MAX_MAPS];

                int numAtomMaps = parseTypes(GPFLine, ligandAtomTypes, MAX_ATOM_TYPES);
                for (int i = 0; i < numAtomMaps; i++)
                    strcpy(ligandTypes[i], ligandAtomTypes[i]);
                for (int i = 0; i < numAtomMaps; i++)
                {
                    foundParam = parameterLibrary.findAtomParameter(ligandTypes[i]);
                    if (foundParam)
                        foundParam->mapIndex = i;
                    else
                    {
                        // return error here
                        logFile.printFormatted("unknown ligand atom type %s\nadd parameters for it to the parameter library first!\n", ligandAtomTypes[i]);
                        throw ExitProgram(-1);
                    }
                }

                // Check to see if there is enough memory to store these map objects
                gridmaps.setNumMaps(numAtomMaps+2);

                if (numGridPointsPerMap == INIT_NUM_GRID_PTS)
                {
                    logFile.printError(ERROR, "You need to set the number of grid points using \"npts\" before setting the ligand atom types, using \"ligand_types\".\n");
                    logFile.printError(FATAL_ERROR, "Unsuccessful completion.\n\n");
                }

                for (int i = 0; i < numAtomMaps; i++)
                {
                    gridmaps[i].mapIndex = i;
                    strcpy(gridmaps[i].type, ligandTypes[i]);   // eg HD or OA or NA or N
                    foundParam = parameterLibrary.findAtomParameter(ligandTypes[i]);
                    gridmaps[i].atomType = foundParam->mapIndex;
                    gridmaps[i].solparProbe = foundParam->solpar;
                    gridmaps[i].volProbe = foundParam->vol;
                    gridmaps[i].Rij = foundParam->Rij;
                    gridmaps[i].epsij = foundParam->epsij;
                    gridmaps[i].hbond = foundParam->hbond;
                    gridmaps[i].RijHB = foundParam->RijHB;
                    gridmaps[i].epsijHB = foundParam->epsijHB;
                    if (gridmaps[i].hbond > 0)
                        gridmaps[i].isHBonder = true;

                    for (int j = 0; j < numReceptorTypes; j++)
                    {
                        foundParam = parameterLibrary.findAtomParameter(receptorTypes[j]);
                        gridmaps[i].nbpR[j] = (gridmaps[i].Rij + foundParam->Rij) / 2;
                        gridmaps[i].nbpEps[j] = sqrt(gridmaps[i].epsij * foundParam->epsij);
                        // apply the vdW forcefield parameter/weight here
                        // This was removed because "setup_p_l" does this for us... gridmaps[i].nbpEps[j] *= FE_coeff_vdW;
                        gridmaps[i].xA[j] = 12;
                        // setup hbond dependent stuff
                        gridmaps[i].xB[j] = 6;
                        gridmaps[i].hbonder[j] = 0;
                        if ((int)(gridmaps[i].hbond) > 2 && ((int)foundParam->hbond == 1 || (int)foundParam->hbond == 2))
                        {           // AS,A1,A2 map vs DS,D1 probe
                            gridmaps[i].xB[j] = 10;
                            gridmaps[i].hbonder[j] = 1;
                            gridmaps[i].isHBonder = true;
                            // Rij and epsij for this hb interaction in parm_data.dat file as Rii and epsii for heavy atom hb factors
                            gridmaps[i].nbpR[j] = gridmaps[i].RijHB;
                            gridmaps[i].nbpEps[j] = gridmaps[i].epsijHB;

                            // apply the hbond forcefield parameter/weight here
                            // This was removed because "setup_p_l" does this for us... gridmaps[i].nbpEps[j] *= FE_coeff_hbond;
                        }
                        else if (((int)gridmaps[i].hbond == 1 || (int)gridmaps[i].hbond == 2) && ((int)foundParam->hbond > 2))
                        {           // DS,D1 map vs AS,A1,A2 probe
                            gridmaps[i].xB[j] = 10;
                            gridmaps[i].hbonder[j] = 1;
                            gridmaps[i].isHBonder = true;
                            // Rij and epsij for this hb interaction in parm_data.dat file as Rii and epsii for heavy atom hb factors
                            gridmaps[i].nbpR[j] = foundParam->RijHB;
                            gridmaps[i].nbpEps[j] = foundParam->epsijHB;

                            // apply the hbond forcefield parameter/weight here
                            // This was removed because "setup_p_l" does this for us... gridmaps[i].nbpEps[j] *= FE_coeff_hbond;
                        }
                    }              // initialize energy parms for each possible receptor type
                }                   // for each map

                logFile.printFormatted("\nAtom type names for ligand atom types 1-%d used for ligand-atom affinity grid maps:\n\n", numAtomMaps);
                for (int i = 0; i < numAtomMaps; i++)
                {
                    logFile.printFormatted("\t\t\tAtom type number %d corresponds to atom type name \"%s\".\n", gridmaps[i].mapIndex, gridmaps[i].type);

                    // FIX THIS!!! Covalent Atom Types are not yet supported with the new AG4/AD4 atom typing mechanism...
                    /* if (gridmaps[i].atomType == COVALENTTYPE) { gridmaps[i].isCovalent = true;  logFile.printFormatted("\nAtom type number %d will be used to calculate a covalent affinity
                       grid map\n\n", i + 1); } */
                }
                logFile.print("\n\n");
            }
            break;

        case GPF_RECEPTOR_TYPES:
            {
                // Read in the list of atom types in the receptor.
                // GPFLine e.g.: "receptorTypes N O A C HH NH"
                //
                // NOTE: This line is not guaranteed to match the actual
                // atom types present in the receptor PDBQT file
                // specified by the "receptor" command.

                // array of ptrs used to parse input line
                char *receptor_atom_types[NUM_RECEPTOR_TYPES];

                numReceptorTypes = parseTypes(GPFLine, receptor_atom_types, MAX_ATOM_TYPES);
                receptor_types_gpf_ct = numReceptorTypes;
                has_receptor_types_in_gpf = 1;

                for (int i = 0; i < numReceptorTypes; i++)
                    strcpy(receptorTypes[i], receptor_atom_types[i]);
                for (int i = 0; i < numReceptorTypes; i++)
                {
                    foundParam = parameterLibrary.findAtomParameter(receptor_atom_types[i]);
                    if (foundParam != 0)
                        foundParam->recIndex = i;
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
            sscanf(GPFLine, "%*s %s %lf %lf", thisparm.autogridType, &temp_vol, &temp_solpar);
            foundParam = parameterLibrary.findAtomParameter(thisparm.autogridType);
            if (foundParam != 0)
            {
                foundParam->vol = temp_vol;
                foundParam->solpar = temp_solpar;
                int i = foundParam->mapIndex;
                if (i >= 0)
                {
                    // DON'T!!!
                    // convert cal/molA^3 to kcal/molA^3
                    // gridmaps[i].solparProbe = temp_solpar * 0.001;
                    gridmaps[i].solparProbe = temp_solpar;
                    logFile.printFormatted("\nProbe %s solvation parameters: \n\n\tatomic fragmental volume: %.2f A^3\n\tatomic solvation parameter: %.4f cal/mol A^3\n\n",
                                  foundParam->autogridType, foundParam->vol, foundParam->solpar);
                }
            }
            else
                logFile.printFormatted("%s key not found\n", thisparm.autogridType);
            break;              // end solvation parameter

        case GPF_MAP:
            /* The variable "mapIndex" is the 0-based index of the ligand atom type we are calculating a map for. If the "types" line was CNOSH, there would be 5 ligand atom maps to
               calculate, and since "mapIndex" is initialized to -1, mapIndex will increment each time there is a "map" keyword in the GPF.  The value of mapIndex should therefore go
               from 0 to 4 for each "map" keyword. In this example, numAtomMaps would be 5, and numAtomMaps-1 would be 4, so if mapIndex is > 4, there is something wrong in the
               number of "map" keywords. */
            ++mapIndex;
            if (mapIndex >= gridmaps.getNumAtomMaps())
            {
                logFile.printErrorFormatted(ERROR,
                    "Too many \"map\" keywords (%d);  the \"types\" command declares only %d maps.\nRemove a \"map\" keyword from the GPF.\n",
                    mapIndex + 1, gridmaps.getNumAtomMaps());
                logFile.printError(FATAL_ERROR, "Unsuccessful completion.\n\n");
            }
            // Read in the gridParameterFilename for this grid map *//* GPF_MAP
            sscanf(GPFLine, "%*s %s", gridmaps[mapIndex].mapFilename);
            if ((gridmaps[mapIndex].file = openFile(gridmaps[mapIndex].mapFilename, "w")) == 0)
            {
                logFile.printErrorFormatted(ERROR, "Cannot open grid map \"%s\" for writing.", gridmaps[mapIndex].mapFilename);
                logFile.printError(FATAL_ERROR, "Unsuccessful completion.\n\n");
            }
            logFile.printFormatted("\nOutput Grid Map %d:   %s\n\n", (mapIndex + 1), gridmaps[mapIndex].mapFilename);
            break;

        case GPF_ELECMAP:
            sscanf(GPFLine, "%*s %s", gridmaps.getElectrostaticMap().mapFilename);
            if ((gridmaps.getElectrostaticMap().file = openFile(gridmaps.getElectrostaticMap().mapFilename, "w")) == 0)
            {
                logFile.printErrorFormatted(ERROR, "can't open grid map \"%s\" for writing.\n", gridmaps.getElectrostaticMap().mapFilename);
                logFile.printError(FATAL_ERROR, "Unsuccessful completion.\n\n");
            }
            logFile.printFormatted("\nOutput Electrostatic Potential Energy Grid Map: %s\n\n", gridmaps.getElectrostaticMap().mapFilename);
            break;

        case GPF_DSOLVMAP:
            sscanf(GPFLine, "%*s %s", gridmaps.getDesolvationMap().mapFilename);
            if ((gridmaps.getDesolvationMap().file = openFile(gridmaps.getDesolvationMap().mapFilename, "w")) == 0)
            {
                logFile.printErrorFormatted(ERROR, "can't open grid map \"%s\" for writing.\n", gridmaps.getDesolvationMap().mapFilename);
                logFile.printError(FATAL_ERROR, "Unsuccessful completion.\n\n");
            }
            logFile.printFormatted("\nOutput Desolvation Free Energy Grid Map: %s\n\n", gridmaps.getDesolvationMap().mapFilename);
            break;

        case GPF_COVALENTMAP:
            sscanf(GPFLine, "%*s %lf %lf %lf %lf %lf", &covHalfWidth, &covBarrier, &(covpos[X]), &(covpos[Y]), &(covpos[Z]));
            logFile.printFormatted("\ncovalentmap <half-width in Angstroms> <barrier> <x> <y> <z>\n"
                                   "\nCovalent well's half-width in Angstroms:         %8.3f\n"
                                   "\nCovalent barrier energy in kcal/mol:             %8.3f\n"
                                   "\nCovalent attachment point will be positioned at: (%8.3f, %8.3f, %8.3f)\n\n",
                                   covHalfWidth, covBarrier, covpos[X], covpos[Y], covpos[Z]);

            // center covpos in the grid maps frame of reference,
            for (int i = 0; i < XYZ; i++)
                covpos[i] -= center[i];
            break;

        case GPF_DISORDER:
            disorderH = true;
            logFile.print("\nHydroxyls will be disordered \n\n");
            break;

        case GPF_SMOOTH:
            sscanf(GPFLine, "%*s %lf", &rSmooth);
            logFile.printFormatted("\nPotentials will be smoothed by: %.3lf Angstrom\n\n", rSmooth);
            break;

        case GPF_QASP:
            sscanf(GPFLine, "%*s %lf", &solparQ);
            logFile.printFormatted("\nCharge component of the atomic solvation parameter: %.3lf\n\n", solparQ);
            // Typical value of solparQ is 0.001118
            break;

        case GPF_DIEL:
            sscanf(GPFLine, "%*s %lf", &diel);
            if (diel < 0)
            {
                // negative...
                distDepDiel = true;
                // calculate ddd of Mehler & Solmajer
                epsilon[0] = 1.0;
                for (int indx_r = 1; indx_r < MAX_DIST; indx_r++)
                    epsilon[indx_r] = calculateDDDMehlerSolmajer(angstrom(indx_r), APPROX_ZERO);
                logFile.print("\nUsing *distance-dependent* dielectric function of Mehler and Solmajer, Prot.Eng.4, 903-910.\n\n"
                              "  d   Dielectric\n ___  __________\n");
                for (int i = 0; i <= 500; i += 10)
                {
                    ri = angstrom(i);
                    logFile.printFormatted("%4.1lf%9.2lf\n", ri, epsilon[i]);
                }
                logFile.print("\n");
                // convert epsilon to 1 / epsilon
                for (int i = 1; i < MAX_DIST; i++)
                    epsilon[i] = factor / epsilon[i];
            }
            else
            {
                // positive or zero...
                distDepDiel = false;
                if (diel <= APPROX_ZERO)
                    diel = 40;
                logFile.printFormatted("Using a *constant* dielectric of:  %.2f\n", diel);
                invDielCal = factor / diel;
            }
            break;

        case GPF_FMAP:
            sscanf(GPFLine, "%*s %s", floatingGridFilename);
            logFile.printFormatted("\nFloating Grid file name = %s\n", floatingGridFilename);
            break;

        case GPF_PARAM_FILE:
            // open and read the AD4 parameters .dat file
            sscanf(GPFLine, "%*s %s ", parameterLibraryFilename);
            break;
        }                       // second switch
    }                           // while

    logFile.print("\n>>> Closing the grid parameter file (GPF)... <<<\n\n" UnderLine);
    fclose(GPF);

    if (!floatingGridFilename[0])
        logFile.print("\n\nNo Floating Grid was requested.\n");
}
