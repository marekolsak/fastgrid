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

#include <new>
#include <cstring>
#include "GridMap.h"
#include "Utils.h"

GridMap::GridMap()
{
    memset(this, 0, sizeof(*this));
}

///////////////////////////////////////////////////////////////////////////////

GridMapList::GridMapList(LogFile *logFile): numMaps(0), numAtomMaps(0), elecIndex(0), desolvIndex(0), gridmaps(0),
    logFile(logFile), floatingGridMins(0), useFloatingGrid(false)
{}

GridMapList::~GridMapList()
{
    for (int i = 0; i < numMaps; i++)
        if (gridmaps[i].filename[0])
            delete [] gridmaps[i].energies;

    if (useFloatingGrid)
        delete [] floatingGridMins;

    if (gridmaps)
        delete [] gridmaps;
}

void GridMapList::setNumMaps(int numMaps)
{
    if (gridmaps)
        delete [] gridmaps;

    numAtomMaps = numMaps - 2;
    elecIndex = numMaps - 2;
    desolvIndex = numMaps - 1;
    this->numMaps = numMaps;

    gridmaps = new(std::nothrow) GridMap[numMaps];
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

    for (int i = 0; i < numMaps; i++)
    {
        double energyMin, energyMax;
        calculateEnergyMinMax(i, energyMin, energyMax);
        logFile->printFormatted(" %d\t %s\t  %6.2lf\t%9.2le\n", i + 1, gridmaps[i].type, energyMin, energyMax);
    }

    double energyMinElec, energyMaxElec;
    double energyMinDesolv, energyMaxDesolv;
    calculateEnergyMinMax(elecIndex, energyMinElec, energyMaxElec);
    calculateEnergyMinMax(desolvIndex, energyMinDesolv, energyMaxDesolv);

    logFile->printFormatted(" %d\t %c\t  %6.2lf\t%9.2le\tElectrostatic Potential\n"
                           " %d\t %c\t  %6.2lf\t%9.2le\tDesolvation Potential\n"
                           "\n\n * Note:  Every pairwise-atomic interaction was clamped at %.2f\n\n\n",
                           getElectrostaticMapIndex() + 1, 'e', energyMinElec, energyMaxElec,
                           getDesolvationMapIndex() + 1, 'd', energyMinDesolv, energyMaxDesolv,
                           EINTCLAMP);

    fprintf(stderr, "\n%s: Successful Completion.\n", logFile->getProgramName());
    logFile->printTitled("Successful Completion.\n");
}

void GridMapList::enableFloatingGrid()
{
    useFloatingGrid = true;
}

int GridMapList::getNumMapsInclFloatingGrid() const
{
    return numMaps + (useFloatingGrid ? 1 : 0);
}

void GridMapList::prepareGridmaps(int numGridPointsPerMap)
{
    this->numGridPointsPerMap = numGridPointsPerMap;

    for (int i = 0; i < numMaps; i++)
        if (gridmaps[i].filename[0])
        {
            gridmaps[i].energies = new double[numGridPointsPerMap];
            memset(gridmaps[i].energies, 0, numGridPointsPerMap * sizeof(double));
        }
    if (useFloatingGrid)
        floatingGridMins = new float[numGridPointsPerMap];
}

void GridMapList::initFileHeader(const InputData *input, const char *gridParameterFilename)
{
    // Write out the correct grid data '.fld' filename at the head of each map file,
    // to avoid centering errors in subsequent dockings...
    // AutoDock can then check to see if the center of each map matches that
    // specified in its parameter file...

    // The header of all files
    fileHeaderLength = snprintf(fileHeader, 1<<14,
        "GRID_PARAMETER_FILE %s\n"
        "GRID_DATA_FILE %s\n"
        "MACROMOLECULE %s\n"
        "SPACING %.3lf\n"
        "NELEMENTS %d %d %d\n"
        "CENTER %.3lf %.3lf %.3lf\n",
        gridParameterFilename, input->fldFilenameAVS, input->receptorFilename, input->gridSpacing,
        input->numGridPoints.x-1, input->numGridPoints.y-1, input->numGridPoints.z-1,
        input->gridCenter.x, input->gridCenter.y, input->gridCenter.z);
}

void GridMapList::saveElectrostaticMap() const
{
    saveGridMap(elecIndex);
}

void GridMapList::saveAtomMapsAndDesolvMap() const
{
    for (int i = 0; i < numAtomMaps; i++)
        saveGridMap(i);
    
    saveGridMap(desolvIndex);
}

void GridMapList::saveFloatingGrid(const char *floatingGridFilename) const
{
    // Floating grid
    if (useFloatingGrid)
    {
        // Open the file
        FILE *file = 0;
        if ((file = boincOpenFile(floatingGridFilename, "w")) == 0)
        {
            logFile->printErrorFormatted(ERROR, "can't open grid map \"%s\" for writing.\n", floatingGridFilename);
            logFile->printError(FATAL_ERROR, "Unsuccessful completion.\n\n");
        }
        if (!fwrite(fileHeader, fileHeaderLength, 1, file))
            logFile->printError(FATAL_ERROR, "Not enough disk space.");

        // Save the floating grid
        for (int j = 0; j < numGridPointsPerMap; j++)
            fprintf(file, "%.3f\n", floatingGridMins[j]);

        fclose(file);
    }
}

void GridMapList::saveGridMap(int index) const
{
    if (gridmaps[index].filename[0])
    {
        // Open the file
        FILE *file;
        if ((file = boincOpenFile(gridmaps[index].filename, "w")) == 0)
        {
            logFile->printErrorFormatted(ERROR, "Cannot open grid map \"%s\" for writing.", gridmaps[index].filename);
            logFile->printError(FATAL_ERROR, "Unsuccessful completion.\n\n");
        }
        if (!fwrite(fileHeader, fileHeaderLength, 1, file))
            logFile->printError(FATAL_ERROR, "Not enough disk space.");

        // Save energies
        for (int j = 0; j < numGridPointsPerMap; j++)
        {
            double f = gridmaps[index].energies[j];
            if (f == 0)
            {
                if (!fwrite("0\n", 2, 1, file))
                    logFile->printError(FATAL_ERROR, "Not enough disk space.");
            }
            else
                fprintf(file, "%.3f\n", f);
        }

        fclose(file);
    }
}

void GridMapList::calculateEnergyMinMax(int map, double &energyMin, double &energyMax)
{
    energyMax = -BIG;
    energyMin = BIG;

    for (int j = 0; j < numGridPointsPerMap; j++)
    {
        energyMax = Mathd::Max(energyMax, gridmaps[map].energies[j]);
        energyMin = Mathd::Min(energyMin, gridmaps[map].energies[j]);
    }
}
