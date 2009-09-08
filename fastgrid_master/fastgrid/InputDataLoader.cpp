/*
    FastGrid (formerly AutoGrid)

    Copyright (C) 2009 The Scripps Research Institute. All rights reserved.
    Copyright (C) 2009 Masaryk University. All rights reserved.

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

#include <cstring>
#include <cctype>
#include <climits>
#include "InputDataLoader.h"
#include "Utils.h"
#include "Exceptions.h"

// GPF tokens
enum GPFTokens
{
    GPF_NULL = 0,
    GPF_COMMENT,
    GPF_RECEPTOR,
    GPF_GRIDFLD,
    GPF_NPTS,
    GPF_SPACING,
    GPF_GRIDCENTER,
    GPF_LIGAND_TYPES,
    GPF_MAP,
    GPF_NBP_COEFFS,
    GPF_NBP_R_EPS,
    GPF_ELECMAP,
    GPF_DIEL,
    GPF_FMAP,
    GPF_DISORDER,
    GPF_SMOOTH,
    GPF_SOL_PAR,
    GPF_CONSTANT,
    GPF_COVALENTMAP,
    GPF_RECEPTOR_TYPES,
    GPF_PARAM_FILE,
    GPF_DSOLVMAP,
    GPF_QASP
};

InputDataLoader::InputDataLoader(LogFile *logFile): logFile(logFile)
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
    covalentPoint = 0;

    // for NEW3 desolvation terms
    solparQ = .01097;   // unweighted value restored 3:9:05
    invDielCal = 0;
    rSmooth = 0;
    gridSpacing = 0.375;     // One quarter of a C-C bond length.
    covHalfWidthSquaredInv = 1.0;
    covBarrier = 1000.0;

    distDepDiel = false;
    disorderH = false;

    receptorAtom = 0;
}

InputDataLoader::~InputDataLoader()
{
    if (receptorAtom)
    {
        delete [] charge;
        delete [] vol;
        delete [] solpar;
        delete [] atomType;
        delete [] hbond;
        delete [] receptorAtom;
    }
}

void InputDataLoader::load(const char *gridParameterFilename, GridMapList &gridmaps, ParameterLibrary &parameterLibrary)
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

    Vec3d gridExtent;
    Vec3d gridCornerMax;
    Vec3d cmax = -BIG;
    Vec3d cmin = BIG;
    Vec3d csum = 0;
    Vec3d cmean;

    // LINE_LEN
    char line[LINE_LEN];
    char GPFLine[LINE_LEN];
    int length = LINE_LEN;

    char atom_name[6];
    char record[LINE_LEN];
    char temp_char = ' ';
    char token[LINE_LEN];
    static const char xyz[] = "xyz"; // used to print headings
    FILE *receptorFile, *xyz_fileptr = 0;

    double q_tot = 0.0;
    double diel;
    double q_max = -BIG, q_min = BIG;
    double ri;
    double temp_vol, temp_solpar;

    int GPF_keyword = -1;
    int indcom = 0;
    int infld;
    int mapIndex = -1;

    // Initializes the grid parameter file
    FILE *GPF = stdin;
    if (gridParameterFilename[0])
    {
        GPF = boincOpenFile(gridParameterFilename, "r");
        if (!GPF)
        {
            logFile->printErrorFormatted(ERROR, "Sorry, I can't find or open Grid Parameter File \"%s\"", gridParameterFilename);
            logFile->printErrorFormatted(ERROR, "Unsuccessful Completion.\n");
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
            logFile->printFormatted("GPF> %s", GPFLine);
            logFile->printError(WARNING, "Unrecognized keyword in grid parameter file.\n");
            continue;           // while fgets GPFLine...

        case GPF_NULL:
        case GPF_COMMENT:
            logFile->printFormatted("GPF> %s", GPFLine);
            break;

        default:
            logFile->printFormatted("GPF> %s", GPFLine);
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
                logFile->printFormatted("\nReceptor Input File :\t%s\n\nReceptor Atom Type Assignments:\n\n", receptorFilename);

                // try to open receptor file
                if ((receptorFile = boincOpenFile(receptorFilename, "r")) == 0)
                {
                    logFile->printErrorFormatted(ERROR, "can't find or open receptor PDBQT file \"%s\".\n", receptorFilename);
                    logFile->printError(FATAL_ERROR, "Unsuccessful completion.\n\n");
                }

                // get the upper bound of total number of atoms
                fseek(receptorFile, 0, SEEK_END);
                int sizeOfFile = ftell(receptorFile);
                fseek(receptorFile, 0, SEEK_SET);
                char *fileContent = new char[sizeOfFile];
                if (!fread(fileContent, sizeOfFile, 1, receptorFile))
                    logFile->printError(FATAL_ERROR, "Cannot read the receptor file.");
                fseek(receptorFile, 0, SEEK_SET);

                int numAtomsMax = 1;
                for (int i = 0; i < sizeOfFile; i++)
                    if (fileContent[i] == '\n')
                        ++numAtomsMax;
                
                delete [] fileContent; 
                
                // reserve space for atoms
                charge = new double[numAtomsMax];
                vol = new double[numAtomsMax];
                solpar = new double[numAtomsMax];
                atomType = new int[numAtomsMax];
                hbond = new HBondType[numAtomsMax];
                receptorAtom = new Vec4d[numAtomsMax];

                // start to read in the lines of the receptor file
                int ia = 0;
                while ((fgets(line, length, receptorFile)) != 0)
                {
                    sscanf(line, "%6s", record);
                    if (strncmp(record, "ATOM", 4) == 0 || // Amino Acid or DNA/RNA atoms
                        strncmp(record, "HETA", 4) == 0 || // Non-standard heteroatoms
                        strncmp(record, "CHAR", 4) == 0)
                    {
                        // Partial Atomic Charge - not a PDB record

                        strncpy(atom_name, &line[12], 4);
                        /* atom_name is declared as an array of 6 characters, the PDB atom name is 4 characters (C indices 0, 1, 2 and 3) but let's ensure that the fifth character (C index 4)
                           is a null character, which terminates the string. */
                        atom_name[4] = '\0';

                        // Output the serial number of this atom...
                        logFile->printFormatted("Atom no. %2d, \"%s\"", ia + 1, atom_name);

                        // Read in this receptor atom's coordinates,partial charges, and solvation parameters in PDBQS format...
                        sscanf(&line[30], "%lf", &receptorAtom[ia].x);
                        sscanf(&line[38], "%lf", &receptorAtom[ia].y);
                        sscanf(&line[46], "%lf", &receptorAtom[ia].z);

                        // Output the coordinates of this atom...
                        logFile->printFormatted(" at (%.3lf, %.3lf, %.3lf), ", receptorAtom[ia].x, receptorAtom[ia].y, receptorAtom[ia].z);

                        // 1:CHANGE HERE: need to set up vol and solpar
                        sscanf(&line[70], "%lf", &charge[ia]);
                        sscanf(&line[77], "%s", thisparm.autogridType);
                        foundParam = parameterLibrary.findAtomParameter(thisparm.autogridType);
                        if (foundParam != 0)
                        {
                            logFile->printFormatted("DEBUG: foundParam->recIndex = %d", foundParam->recIndex);
                            if (foundParam->recIndex < 0)
                            {
                                strncpy(receptorTypes[numReceptorTypes], foundParam->autogridType, 3);
                                foundParam->recIndex = numReceptorTypes++;
                                logFile->printFormatted("DEBUG: foundParam->recIndex => %d", foundParam->recIndex);
                            }
                            atomType[ia] = foundParam->recIndex;
                            solpar[ia] = foundParam->solpar;
                            vol[ia] = foundParam->vol;
                            hbond[ia] = foundParam->hbond;  // NON=0, DS,D1, AS, A1, A2
                            ++receptor_atom_type_count[foundParam->recIndex];
                        }
                        else
                        {
                            logFile->printFormatted("\n\nreceptor file contains unknown type: '%s'\nadd parameters for it to the parameter library first\n", thisparm.autogridType);
                            throw ExitProgram(-1);
                        }

                        // if from pdbqs: convert cal/molA**3 to kcal/molA**3
                        // solpar[ia] *= 0.001;
                        q_max = Mathd::Max(q_max, charge[ia]);
                        q_min = Mathd::Min(q_min, charge[ia]);

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
                        logFile->printFormatted(" was assigned atom type \"%s\" (recIndex= %d, atomType= %d).\n", foundParam->autogridType, foundParam->recIndex, atomType[ia]);

                        // Count the number of each atom type
                        // ++receptor_atom_type_count[ atomType[ia] ];

                        // Keep track of the extents of the receptor
                        for (int i = 0; i < 3; i++)
                        {
                            cmax[i] = Mathd::Max(cmax[i], receptorAtom[ia][i]);
                            cmin[i] = Mathd::Min(cmin[i], receptorAtom[ia][i]);
                        }
                        csum += Vec3d(receptorAtom[ia]);
                        // Total up the partial charges as we go...
                        q_tot += charge[ia];

                        // Increment the atom counter
                        ia++;

                        // Check that there aren't too many atoms...
                        if (ia > AG_MAX_ATOMS)
                        {
                            logFile->printErrorFormatted(ERROR, "Too many atoms in receptor PDBQT file %s;", receptorFilename);
                            logFile->printErrorFormatted(ERROR, "-- the maximum number of atoms, AG_MAX_ATOMS, allowed is %ul.", AG_MAX_ATOMS);
                            logFile->printErrorFormatted(SUGGESTION, "Increase the value in the \"#define AG_MAX_ATOMS %ul\" line", AG_MAX_ATOMS);
                            logFile->printError(SUGGESTION, "in the source file \"autogrid.h\", and re-compile " APPNAME ".");
                            // FATAL_ERROR will cause this program to exit...
                            logFile->printError(FATAL_ERROR, "Sorry, " APPNAME " cannot continue.");
                        }           // endif
                    }               // endif
                }                   // endwhile
                // Finished reading in the lines of the receptor file
                fclose(receptorFile);
                if (has_receptor_types_in_gpf == 1)
                    // Check that the number of atom types found in the receptor PDBQT
                    // file match the number parsed in by the "receptorTypes" command
                    // in the GPF; if they do not match, exit!
                    if (numReceptorTypes != receptor_types_gpf_ct)
                    {
                        logFile->printErrorFormatted(ERROR,
                            "The number of atom types found in the receptor PDBQT (%d) does not match the number specified by the \"receptor_types\" command (%d) in the GPF!\n\n",
                            numReceptorTypes, receptor_types_gpf_ct);
                        // FATAL_ERROR will cause this program to exit...
                        logFile->printError(FATAL_ERROR, "Sorry, " APPNAME " cannot continue.");
                    }

                // Update the total number of atoms in the receptor
                numReceptorAtoms = ia;
                logFile->printFormatted("\nMaximum partial atomic charge found = %+.3lf e\n", q_max);
                logFile->printFormatted("Minimum partial atomic charge found = %+.3lf e\n\n", q_min);

                // Check there are partial charges...
                if (q_max == 0 && q_min == 0)
                {
                    logFile->printErrorFormatted(ERROR, "No partial atomic charges were found in the receptor PDBQT file %s!\n\n", receptorFilename);
                    // FATAL_ERROR will cause this program to exit...
                    logFile->printError(FATAL_ERROR, "Sorry, " APPNAME " cannot continue.");
                }                   // if there are no charges EXIT

                logFile->print("Atom\tAtom\tNumber of this Type\n"
                              "Type\t ID \t in Receptor\n"
                              "____\t____\t___________________\n");

                // 2. CHANGE HERE: need to count number of each receptor_type
                for (int ia = 0; ia < numReceptorTypes; ia++)
                    if (receptor_atom_type_count[ia] != 0)
                        logFile->printFormatted(" %d\t %s\t\t%6d\n", (ia), receptorTypes[ia], receptor_atom_type_count[ia]);

                logFile->printFormatted("\nTotal number of atoms :\t\t%d atoms \n"
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
                                       numReceptorAtoms, q_tot, cmax.x, cmax.y, cmax.z,
                                       (cmax.x + cmin.x) / 2, (cmax.y + cmin.y) / 2, (cmax.z + cmin.z) / 2,
                                       cmin.x, cmin.y, cmin.z, cmax.x, cmax.y, cmax.z, cmin.x, cmin.y, cmin.z);

                cmean[0] = csum[0] / numReceptorAtoms;
                cmean[1] = csum[1] / numReceptorAtoms;
                cmean[2] = csum[2] / numReceptorAtoms;
            }
            break;

        case GPF_GRIDFLD:
            sscanf(GPFLine, "%*s %s", fldFilenameAVS);
            infld = strIndex(fldFilenameAVS, ".fld");
            if (infld == -1)
                logFile->printError(FATAL_ERROR, "Grid data file needs the extension \".fld\" for AVS input\n\n");
            else
            {
                infld = strIndex(fldFilenameAVS, "fld");
                strncpy(xyzFilename, fldFilenameAVS, MAX_CHARS);
                xyzFilename[infld] = 'x';
                xyzFilename[infld + 1] = 'y';
                xyzFilename[infld + 2] = 'z';
            }
            if ((xyz_fileptr = boincOpenFile(xyzFilename, "w")) == 0)
            {
                logFile->printErrorFormatted(ERROR, "can't create grid extrema data file %s\n", xyzFilename);
                logFile->printError(ERROR, "SORRY!    unable to create the \".xyz\" file.\n\n");
                logFile->printError(FATAL_ERROR, "Unsuccessful completion.\n\n");
            }
            else
                logFile->printFormatted("\nCreating (AVS-readable) grid-coordinates extrema file : %s\n\n", xyzFilename);
            break;

        case GPF_NPTS:
            {
                Vec3i numGridPointsMinusOne;
                sscanf(GPFLine, "%*s %d %d %d", &numGridPointsMinusOne.x, &numGridPointsMinusOne.y, &numGridPointsMinusOne.z);

                // numGridPointsMinusOne mustn't be negative, shouldn't be zero or larger than MAX_GRID_PTS and should be even
                for (int i = 0; i < 3; i++)
                    numGridPointsMinusOne[i] = checkSize(numGridPointsMinusOne[i], xyz[i]);

                numGridPointsDiv2 = numGridPointsMinusOne / 2;
                numGridPoints = numGridPointsMinusOne + 1;
                logFile->printFormatted("\nNumber of grid points in x-direction:\t%d\n"
                                       "Number of grid points in y-direction:\t%d\n"
                                       "Number of grid points in z-direction:\t%d\n\n",
                                       numGridPoints.x, numGridPoints.y, numGridPoints.z);

                numGridPointsPerMap = numGridPoints.Cube();
            }
            break;

        case GPF_SPACING:
            sscanf(GPFLine, "%*s %lf", &gridSpacing);
            logFile->printFormatted("Grid Spacing :\t\t\t%.3lf Angstrom\n\n", gridSpacing);
            break;

        case GPF_GRIDCENTER:
            sscanf(GPFLine, "%*s %s", token);
            if (strncmp(token, "auto", 4) == 0)
            {
                gridCenter = cmean;
                logFile->printFormatted("Grid maps will be centered on the center of mass.\n"
                                       "Coordinates of center of mass : (%.3lf, %.3lf, %.3lf)\n", gridCenter.x, gridCenter.y, gridCenter.z);
            }
            else
            {
                sscanf(GPFLine, "%*s %lf %lf %lf", &gridCenter.x, &gridCenter.y, &gridCenter.z);
                logFile->printFormatted("\nGrid maps will be centered on user-defined coordinates:\n\n\t\t(%.3lf, %.3lf, %.3lf)\n", gridCenter.x, gridCenter.y, gridCenter.z);
            }
            // centering stuff...
            for (int ia = 0; ia < numReceptorAtoms; ia++)
                *((Vec3d*)&receptorAtom[ia]) -= gridCenter;  // transform to center of gridmaps
            gridExtent = Vec3d(numGridPointsDiv2) * gridSpacing;
            gridCornerMax = gridCenter + gridExtent;
            gridCornerMin = gridCenter - gridExtent;

            logFile->printFormatted("\nGrid maps will cover the following volume:\n\n"
                                   "                   _______(%.1lf, %.1lf, %.1lf)\n"
                                   "                  /|     /|\n"
                                   "                 / |    / |\n"
                                   "                /______/  |\n"
                                   "                |  |___|__| Midpoint = (%.1lf, %.1lf, %.1lf)\n"
                                   "                |  /   |  /\n"
                                   "                | /    | /\n"
                                   "                |/_____|/\n"
                                   "(%.1lf, %.1lf, %.1lf)      \n\n", gridCornerMax.x, gridCornerMax.y, gridCornerMax.z,
                                   gridCenter.x, gridCenter.y, gridCenter.z, gridCornerMin.x, gridCornerMin.y, gridCornerMin.z);

            for (int i = 0; i < 3; i++)
                logFile->printFormatted("Grid map %c-dimension :\t\t%.1lf Angstroms\n", xyz[i], 2 * gridExtent[i]);

            logFile->printFormatted("\nMaximum coordinates :\t\t(%.3lf, %.3lf, %.3lf)\n"
                                   "Minimum coordinates :\t\t(%.3lf, %.3lf, %.3lf)\n\n", gridCornerMax.x, gridCornerMax.y, gridCornerMax.z, gridCornerMin.x, gridCornerMin.y, gridCornerMin.z);
            for (int i = 0; i < 3; i++)
                fprintf(xyz_fileptr, "%.3lf %.3lf\n", gridCornerMin[i], gridCornerMax[i]);
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
                    strncpy(ligandTypes[i], ligandAtomTypes[i], 3);
                for (int i = 0; i < numAtomMaps; i++)
                {
                    foundParam = parameterLibrary.findAtomParameter(ligandTypes[i]);
                    if (foundParam)
                        foundParam->mapIndex = i;
                    else
                    {
                        // return error here
                        logFile->printFormatted("unknown ligand atom type %s\nadd parameters for it to the parameter library first!\n", ligandAtomTypes[i]);
                        throw ExitProgram(-1);
                    }
                }

                // Check to see if there is enough memory to store these map objects
                gridmaps.setNumMaps(numAtomMaps+2);

                if (numGridPointsPerMap == INIT_NUM_GRID_PTS)
                {
                    logFile->printError(ERROR, "You need to set the number of grid points using \"npts\" before setting the ligand atom types, using \"ligand_types\".\n");
                    logFile->printError(FATAL_ERROR, "Unsuccessful completion.\n\n");
                }

                for (int i = 0; i < numAtomMaps; i++)
                {
                    strncpy(gridmaps[i].type, ligandTypes[i], 3);   // eg HD or OA or NA or N
                    foundParam = parameterLibrary.findAtomParameter(ligandTypes[i]);
                    if (strcmp(ligandTypes[i], "Z") == 0)
                    {
                        logFile->print("Found covalent map atomtype\n");
                        gridmaps[i].isCovalent = true;
                    }
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
                        if (gridmaps[i].hbond > D1 && (foundParam->hbond == DS || foundParam->hbond == D1))
                        {
                            // AS,A1,A2 map vs DS,D1 probe
                            gridmaps[i].xB[j] = 10;
                            gridmaps[i].hbonder[j] = 1;
                            gridmaps[i].isHBonder = true;
                            // Rij and epsij for this hb interaction in parm_data.dat file as Rii and epsii for heavy atom hb factors
                            gridmaps[i].nbpR[j] = gridmaps[i].RijHB;
                            gridmaps[i].nbpEps[j] = gridmaps[i].epsijHB;

                            // apply the hbond forcefield parameter/weight here
                            // This was removed because "setup_p_l" does this for us... gridmaps[i].nbpEps[j] *= FE_coeff_hbond;
                        }
                        else if ((gridmaps[i].hbond == DS || gridmaps[i].hbond == D1) && foundParam->hbond > D1)
                        {
                            // DS,D1 map vs AS,A1,A2 probe
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

                logFile->printFormatted("\nAtom type names for ligand atom types 1-%d used for ligand-atom affinity grid maps:\n\n", numAtomMaps);
                for (int i = 0; i < numAtomMaps; i++)
                {
                    logFile->printFormatted("\t\t\tAtom type number %d corresponds to atom type name \"%s\".\n", i, gridmaps[i].type);

                    if (gridmaps[i].isCovalent)
                        logFile->printFormatted("\nAtom type number %d will be used to calculate a covalent affinity grid map\n\n", i + 1);
                }
                logFile->print("\n\n");
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
                    strncpy(receptorTypes[i], receptor_atom_types[i], 3);
                for (int i = 0; i < numReceptorTypes; i++)
                {
                    foundParam = parameterLibrary.findAtomParameter(receptor_atom_types[i]);
                    if (foundParam != 0)
                        foundParam->recIndex = i;
                    else
                    {
                        logFile->printErrorFormatted(ERROR, "Unknown receptor type: \"%s\"\n -- Add parameters for it to the parameter library first!\n", receptor_atom_types[i]);
                        logFile->printError(FATAL_ERROR, "Unsuccessful completion.\n\n");
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
                    logFile->printFormatted("\nProbe %s solvation parameters: \n\n\tatomic fragmental volume: %.2f A^3\n\tatomic solvation parameter: %.4f cal/mol A^3\n\n",
                                  foundParam->autogridType, foundParam->vol, foundParam->solpar);
                }
            }
            else
                logFile->printFormatted("%s key not found\n", thisparm.autogridType);
            break;              // end solvation parameter

        case GPF_MAP:
            /* The variable "mapIndex" is the 0-based index of the ligand atom type we are calculating a map for. If the "types" line was CNOSH, there would be 5 ligand atom maps to
               calculate, and since "mapIndex" is initialized to -1, mapIndex will increment each time there is a "map" keyword in the GPF.  The value of mapIndex should therefore go
               from 0 to 4 for each "map" keyword. In this example, numAtomMaps would be 5, and numAtomMaps-1 would be 4, so if mapIndex is > 4, there is something wrong in the
               number of "map" keywords. */
            ++mapIndex;
            if (mapIndex >= gridmaps.getNumAtomMaps())
            {
                logFile->printErrorFormatted(ERROR,
                    "Too many \"map\" keywords (%d);  the \"ligand_types\" command declares only %d atom types.\nRemove a \"map\" keyword from the GPF.\n",
                    mapIndex + 1, gridmaps.getNumAtomMaps());
                logFile->printError(FATAL_ERROR, "Unsuccessful completion.\n\n");
            }
            // Read in the gridParameterFilename for this grid map
            sscanf(GPFLine, "%*s %s", gridmaps[mapIndex].filename);
            logFile->printFormatted("\nOutput Grid Map %d:   %s\n\n", (mapIndex + 1), gridmaps[mapIndex].filename);
            break;

        case GPF_ELECMAP:
            sscanf(GPFLine, "%*s %s", gridmaps.getElectrostaticMap().filename);
            logFile->printFormatted("\nOutput Electrostatic Potential Energy Grid Map: %s\n\n", gridmaps.getElectrostaticMap().filename);
            break;

        case GPF_DSOLVMAP:
            sscanf(GPFLine, "%*s %s", gridmaps.getDesolvationMap().filename);
            logFile->printFormatted("\nOutput Desolvation Free Energy Grid Map: %s\n\n", gridmaps.getDesolvationMap().filename);
            break;

        case GPF_COVALENTMAP:
            {
                double covHalfWidth;
                sscanf(GPFLine, "%*s %lf %lf %lf %lf %lf", &covHalfWidth, &covBarrier, &(covalentPoint.x), &(covalentPoint.y), &(covalentPoint.z));
                covHalfWidthSquaredInv = 1 / Mathd::Sqr(covHalfWidth);

                logFile->printFormatted("\ncovalentmap <half-width in Angstroms> <barrier> <x> <y> <z>\n"
                                       "\nCovalent well's half-width in Angstroms:         %8.3f\n"
                                       "\nCovalent barrier energy in kcal/mol:             %8.3f\n"
                                       "\nCovalent attachment point will be positioned at: (%8.3f, %8.3f, %8.3f)\n\n",
                                       covHalfWidth, covBarrier, covalentPoint.x, covalentPoint.y, covalentPoint.z);

                // center covalentPoint in the grid maps frame of reference,
                covalentPoint -= gridCenter;
            }
            break;

        case GPF_DISORDER:
            disorderH = true;
            logFile->print("\nHydroxyls will be disordered \n\n");
            break;

        case GPF_SMOOTH:
            sscanf(GPFLine, "%*s %lf", &rSmooth);
            logFile->printFormatted("\nPotentials will be smoothed by: %.3lf Angstrom\n\n", rSmooth);
            break;

        case GPF_QASP:
            sscanf(GPFLine, "%*s %lf", &solparQ);
            logFile->printFormatted("\nCharge component of the atomic solvation parameter: %.3lf\n\n", solparQ);
            // Typical value of solparQ is 0.001118
            break;

        case GPF_DIEL:
            sscanf(GPFLine, "%*s %lf", &diel);
            if (diel < 0)
            {
                // negative...
                distDepDiel = true;
                // calculate inverted ddd of Mehler & Solmajer
                for (int indexR = 0; indexR < MAX_DIST; indexR++)
                    epsilon[indexR] = calculateDistDepDielInv(indexToAngstrom<double>(indexR));
                logFile->print("\nUsing *distance-dependent* dielectric function of Mehler and Solmajer, Prot.Eng.4, 903-910.\n\n"
                              "  d   Dielectric\n ___  __________\n");
                for (int i = 0; i <= 500; i += 10)
                {
                    ri = indexToAngstrom<double>(i);
                    logFile->printFormatted("%4.1lf%9.2lf\n", ri, DDD_FACTOR / epsilon[i]);
                }
                logFile->print("\n");
            }
            else
            {
                // positive or zero...
                distDepDiel = false;
                if (diel <= APPROX_ZERO)
                    diel = 40;
                logFile->printFormatted("Using a *constant* dielectric of:  %.2f\n", diel);
                invDielCal = DDD_FACTOR / diel;
            }
            break;

        case GPF_FMAP:
            sscanf(GPFLine, "%*s %s", floatingGridFilename);
            logFile->printFormatted("\nFloating Grid file name = %s\n", floatingGridFilename);
            break;

        case GPF_PARAM_FILE:
            // open and read the AD4 parameters .dat file
            sscanf(GPFLine, "%*s %s ", parameterLibraryFilename);
            break;
        }                       // second switch
    }                           // while: finished reading gpf

    // Map files checkpoint  (number of maps, desolv and elec maps) SF

    // Number of maps defined for atom types
    if (mapIndex < gridmaps.getNumAtomMaps() - 1)
    {
        logFile->printFormatted("Too few \"map\" keywords (%d);  the \"ligand_types\" command declares %d atom types.\nAdd a \"map\" keyword from the GPF.\n", mapIndex + 1, gridmaps.getNumAtomMaps());
        logFile->printError(ERROR, "Not enough map keywords found.\n");
        logFile->printError(FATAL_ERROR, "Unsuccessful completion.\n\n");
	}

    // Desolvation map
    if (!gridmaps.getDesolvationMap().filename[0])
    {
        logFile->print("The desolvation map file is not defined in the GPF.\n");
        logFile->printError(ERROR, "No desolvation map file defined.\n");
        logFile->printError(FATAL_ERROR, "Unsuccessful completion.\n\n");
	}

    // Electrostatic map
    if (!gridmaps.getElectrostaticMap().filename[0])
    {
        logFile->print("The electrostatic map file is not defined in the GPF.\n");
        logFile->printError(ERROR, "No electrostatic map file defined.\n");
        logFile->printError(FATAL_ERROR, "Unsuccessful completion.\n\n");
	}
    // End of map files checkpoint SF

    logFile->print("\n>>> Closing the grid parameter file (GPF)... <<<\n\n" UnderLine);
    fclose(GPF);

    if (!floatingGridFilename[0])
        logFile->print("\n\nNo Floating Grid was requested.\n");

    // apply the estat forcefield coefficient/weight here
    if (distDepDiel)
        for (int ia = 0; ia < numReceptorAtoms; ia++)
            receptorAtom[ia].w = charge[ia] * parameterLibrary.coeff_estat;
    else
        // apply the constant dielectric
        for (int ia = 0; ia < numReceptorAtoms; ia++)
            receptorAtom[ia].w = charge[ia] * parameterLibrary.coeff_estat * invDielCal;
}

int InputDataLoader::parseGPFLine(const char *line)
{
    size_t l;
	int i, token = -1; // return -1 if nothing is recognized
    char c[LINE_LEN];

    l = strIndex(line, " ");
    if (l == -1)
    {
        l = strIndex(line, "\t");
        if (l == -1)
            l = strlen(line);
    }
    for(i=0; i<l; i++)
        c[i] = char(tolower(line[i]));

    if ((c[0]=='\n')||(c[0]=='\0'))
        token = GPF_NULL;
    else if (c[0]=='#')
        token = GPF_COMMENT;
    else if (strncmp(c, "receptor_types", 14) == 0)
        token = GPF_RECEPTOR_TYPES;
    else if (strncmp(c, "receptor", 8) == 0)
        token = GPF_RECEPTOR;
    else if (strncmp(c, "gridfld", 7) == 0)
        token = GPF_GRIDFLD;
    else if (strncmp(c, "npts", 4) == 0)
        token = GPF_NPTS;
    else if (strncmp(c, "spacing", 7) == 0)
        token = GPF_SPACING;
    else if (strncmp(c, "gridcenter", 10) == 0)
        token = GPF_GRIDCENTER;
    else if (strncmp(c, "types", 5) == 0)
        token = GPF_LIGAND_TYPES;
    else if (strncmp(c, "ligand_types", 12) == 0)
        token = GPF_LIGAND_TYPES;
    else if (strncmp(c, "map", 3) == 0)
        token = GPF_MAP;
    else if (strncmp(c, "elecmap", 7) == 0)
        token = GPF_ELECMAP;
    else if (strncmp(c, "dsolvmap", 8) == 0)
        token = GPF_DSOLVMAP;
    else if (strncmp(c, "covalentmap", 11) == 0)
        token = GPF_COVALENTMAP;
    else if (strncmp(c, "nbp_coeffs", 10) == 0)
        token = GPF_NBP_COEFFS;
    else if (strncmp(c, "nbp_r_eps", 9) == 0)
        token = GPF_NBP_R_EPS;
    else if (strncmp(c, "dielectric", 10) == 0)
        token = GPF_DIEL;
    else if (strncmp(c, "qasp", 4) == 0)
        token = GPF_QASP;
    else if (strncmp(c, "fmap", 4) == 0)
        token = GPF_FMAP;
    else if (strncmp(c, "disorder_h", 10) == 0)
        token = GPF_DISORDER;
    else if (strncmp(c, "smooth", 6) == 0)
        token = GPF_SMOOTH;
    else if (strncmp(c, "sol_par", 7) == 0)
        token = GPF_SOL_PAR;
    else if (strncmp(c, "constant", 8) == 0)
        token = GPF_CONSTANT;
    else if (strncmp(c, "parameter_file", 14) == 0)
        token = GPF_PARAM_FILE;

    return token;
}

// checks that number of grid elements is valid
int InputDataLoader::checkSize(int numGridPointsMinusOne, char axischar)
{
    // numGridPointsMinusOne mustn't be negative, shouldn't be zero or larger than MAX_GRID_PTS and should be even
    if (numGridPointsMinusOne < 0)
        logFile->printErrorFormatted(FATAL_ERROR, "Negative number of %c-grid elements!  Aborting.\n\n", axischar);
    else if (numGridPointsMinusOne == 0)
        logFile->printErrorFormatted(WARNING, "0 %c-grid elements!\n\n", axischar);
    else if (numGridPointsMinusOne > MAX_GRID_PTS)
    {
        logFile->printErrorFormatted(WARNING, "Maximum number of %c-grid elements allowed is %d. Using this value.\n", axischar, MAX_GRID_PTS);
        numGridPointsMinusOne = MAX_GRID_PTS;
    }
    else if (numGridPointsMinusOne % 2 == 1)
    {
        logFile->printTitledFormatted("Number of grid elements must be even; %c-elements changed to: %d\n", axischar, numGridPointsMinusOne);
        numGridPointsMinusOne -= 1;
    }

    return numGridPointsMinusOne;
}

// utility func for parsing types
int InputDataLoader::parseTypes(char *line, char **words, int maxwords)
{
    char *char_ptr = line;
    int num_types = 0;
    // flag for first word which is always a keyword
    int found_keyword = 0;
    int index = 0;

    for (;;)
    {
        // skip spaces
        while (isspace(*char_ptr))
        {
            char_ptr++;
            index++;
        }

        // done parsing when get eol 'null' character
        // could get null after a space
        if (*char_ptr == '\0')
            return num_types; // return number of 'types' found

        // the first word is the keyword not a type
        if (!found_keyword)
            found_keyword++;
        else
            words[num_types++] = char_ptr; // words is a list of indicies of beginning of 1 or 2 char types

        // once in a type, skip possible 2nd characters up to a space or null character
        while (!isspace(*char_ptr) && *char_ptr != '\0')
        {
            char_ptr++;
            index++;
        }

        // done parsing when get eol 'null' character
        // could get null after a character
        if (*char_ptr == '\0')
            return num_types;

        // make each 'type' a null terminated string
        *char_ptr++ = '\0';
        index++;

        //if there are too many types, return
        if (num_types >= maxwords)
            return num_types;
    }
}

// returns index of t in s, -1 if none.
int InputDataLoader::strIndex(const char *s, const char *t)
{
    const char *r = strstr(s, t);
    return r? int(r-s) : -1;
}
