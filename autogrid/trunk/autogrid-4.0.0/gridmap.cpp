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
#include <new>

GridMap::GridMap()
{
    atomType = 0;   // corresponds to receptor numbers????
    mapIndex = 0;
    isCovalent = false;
    isHBonder = false;
    file = 0;
    mapFilename[0] = 0;
    type[0] = 0;    // eg HD or OA or NA or N
    constant = 0; // this will become obsolete
    energyMax = -BIG;
    energyMin = BIG;
    energy = 0;
    volProbe = 0;
    solparProbe = 0;
    Rij = 0;
    epsij = 0;
    hbond = NON; // hbonding character:
    RijHB = 0;
    epsijHB = 0;

    // per gridmaps[i].receptor type parameters, ordered as in receptorTypes
    for (int j = 0; j < NUM_RECEPTOR_TYPES; j++)
    {
        nbpR[j] = 0; // radius of energy-well minimum
        nbpEps[j] = 0;   // depth of energy-well minimum
        xA[j] = 0;   // generally 12
        xB[j] = 0;   // 6 for non-hbonders 10 for h-bonders
        hbonder[j] = 0;
    }
}

GridMap::~GridMap()
{
    if (file)
        fclose(file);
}

///////////////////////////////////////////////////////////////////////////////

GridMapList::GridMapList(LogFile *logFile): numMaps(0), numAtomMaps(0), elecIndex(0), desolvIndex(0), gridmaps(0), logFile(logFile)
{}

GridMapList::~GridMapList()
{
    if (gridmaps)
        delete [] gridmaps;
}

void GridMapList::initFileHeader(const InputData *input, const char *gridParameterFilename)
{
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

    fileHeaderLength = strlen(fileHeader);
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
