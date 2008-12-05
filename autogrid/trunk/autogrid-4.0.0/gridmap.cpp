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

#include "gridmap.h"
#include "utils.h"
#include <new>
#include <cstring>

GridMap::GridMap()
{
    memset(this, 0, sizeof(*this));
    energyMax = -BIG;
    energyMin = BIG;
    hbond = NON;
}

///////////////////////////////////////////////////////////////////////////////

GridMapList::GridMapList(LogFile *logFile): numMaps(0), numAtomMaps(0), elecIndex(0), desolvIndex(0), gridmaps(0), logFile(logFile), floatingGridFile(0)
{
    floatingGridFilename[0] = 0;
}

GridMapList::~GridMapList()
{
    for (int i = 0; i < numMaps; i++)
        if (gridmaps[i].file)
            fclose(gridmaps[i].file);

    if (floatingGridFile)
        fclose(floatingGridFile);

    if (gridmaps)
        delete [] gridmaps;
}

void GridMapList::prepareFiles(const InputData *input, const char *gridParameterFilename)
{
    // Write out the  correct grid_data '.fld' file_name at the head of each map
    // file, to avoid centering errors in subsequent dockings...
    // AutoDock can then check to see if the center of each map matches that
    // specified in its parameter file...

    char fileHeader[1<<14];
    snprintf(fileHeader, 1<<14,
        "GRID_PARAMETER_FILE %s\n"
        "GRID_DATA_FILE %s\n"
        "MACROMOLECULE %s\n"
        "SPACING %.3lf\n"
        "NELEMENTS %d %d %d\n"
        "CENTER %.3lf %.3lf %.3lf\n",
        gridParameterFilename, input->fldFilenameAVS, input->receptorFilename, input->spacing,
        input->nelements[X], input->nelements[Y], input->nelements[Z],
        input->center[X], input->center[Y], input->center[Z]);

    int fileHeaderLength = strlen(fileHeader);

    // Open files
    for (int i = 0; i < numMaps; i++)
        if (gridmaps[i].filename[0])
        {
            if ((gridmaps[i].file = boincOpenFile(gridmaps[i].filename, "w")) == 0)
            {
                logFile->printErrorFormatted(ERROR, "Cannot open grid map \"%s\" for writing.", gridmaps[i].filename);
                logFile->printError(FATAL_ERROR, "Unsuccessful completion.\n\n");
            }

            fwrite(fileHeader, fileHeaderLength, 1, gridmaps[i].file);
        }

    if (input->floatingGridFilename[0])
    {
        if ((floatingGridFile = boincOpenFile(input->floatingGridFilename, "w")) == 0)
        {
            logFile->printErrorFormatted(ERROR, "can't open grid map \"%s\" for writing.\n", input->floatingGridFilename);
            logFile->printError(FATAL_ERROR, "Unsuccessful completion.\n\n");
        }
        fwrite(fileHeader, fileHeaderLength, 1, floatingGridFile);
    }
}

void GridMapList::setNumMaps(int num)
{
    if (gridmaps)
        delete [] gridmaps;

    numAtomMaps = num - 2;
    elecIndex = num - 2;
    desolvIndex = num - 1;
    numMaps = num;

    gridmaps = new(std::nothrow) GridMap[num];
    if (!gridmaps)
    {
        logFile->printError(ERROR, "Could not allocate memory to create the GridMap \"gridmaps\".\n");
        logFile->printError(FATAL_ERROR, "Unsuccessful completion.\n\n");
    }
}

void GridMapList::logSummary()
{
    // Print a summary of extrema-values from the atomic-affinity and
    // electrostatics grid-maps,
    logFile->print("\nGrid\tAtom\tMinimum   \tMaximum\n"
                  "Map \tType\tEnergy    \tEnergy \n"
                  "\t\t(kcal/mol)\t(kcal/mol)\n"
                  "____\t____\t_____________\t_____________\n");

    for (int i = 0; i < numAtomMaps; i++)
        logFile->printFormatted(" %d\t %s\t  %6.2lf\t%9.2le\n", i + 1, gridmaps[i].type, gridmaps[i].energyMin, gridmaps[i].energyMax);

    logFile->printFormatted(" %d\t %c\t  %6.2lf\t%9.2le\tElectrostatic Potential\n"
                           " %d\t %c\t  %6.2lf\t%9.2le\tDesolvation Potential\n"
                           "\n\n * Note:  Every pairwise-atomic interaction was clamped at %.2f\n\n\n",
                           getElectrostaticMapIndex() + 1, 'e', getElectrostaticMap().energyMin, getElectrostaticMap().energyMax,
                           getDesolvationMapIndex() + 1, 'd', getDesolvationMap().energyMin, getDesolvationMap().energyMax,
                           EINTCLAMP);

    fprintf(stderr, "\n%s: Successful Completion.\n", logFile->getProgramName());
    logFile->printTitled("Successful Completion.\n");
}

void GridMapList::setFloatingGridFilename(const char *filename)
{
    strncpy(floatingGridFilename, filename, MAX_CHARS);
}
