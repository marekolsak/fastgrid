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

#pragma once
#include "InputData.h"

struct GridMap
{
    int atomType;          // corresponds to receptor numbers????
    bool isCovalent;
    int isHBonder;
    char filename[MAX_CHARS];
    char type[3];           // eg HD or OA or NA or N
    double constant;        // this will become obsolete
    double volProbe;
    double solparProbe;

    // new 6/28
    double Rij;
    double epsij;
    HBondType hbond;       // hbonding character
    double RijHB;
    double epsijHB;

    // per receptor type parameters, ordered as in receptorTypes
    double nbpR[NUM_RECEPTOR_TYPES];   // radius of energy-well minimum
    double nbpEps[NUM_RECEPTOR_TYPES]; // depth of energy-well minimum
    int xA[NUM_RECEPTOR_TYPES]; // generally 12
    int xB[NUM_RECEPTOR_TYPES]; // 6 for non-hbonders 10 for h-bonders
    int hbonder[NUM_RECEPTOR_TYPES];

    double *energies;  // output energies, this array will be saved to file

    GridMap();
};

class GridMapList
{
public:
    GridMapList(LogFile *logFile);
    ~GridMapList();

    // Allocates gridmaps.
    // "num" is the number of maps to be created: the number of ligand atom types, plus 1 for the electrostatic map,
    // plus 1 for the desolvation map. (the floating grid should not be included here, see enableFloatingGrid)
    // Keep in mind that AutoDock can only read in MAX_MAPS maps.
    void setNumMaps(int numMaps);

    // Enables the floating grid. By default, the floating grid is disabled.
    void enableFloatingGrid();

    // Allocates memory for output energies.
    void prepareGridmaps(int numGridPointsPerMap);

    // TODO: unify functions: setNumMaps, enableFloatingGrid, prepareGridmaps, to:
    // void initialize(int numMaps, bool enableFloatingGrid, int numGridPointsPerMap);

    // Saves all energies and possibly the floating grid to files
    void initFileHeader(const InputData *input, const char *gridParameterFilename);
    void saveElectrostaticMap() const;
    void saveAtomMapsAndDesolvMap() const;
    void saveFloatingGrid(const char *floatingGridFilename) const;

    // Writes out summary
    void logSummary();

    // Read/write access to maps
    GridMap &operator [](int i)             { return gridmaps[i]; }
    GridMap &getElectrostaticMap()          { return gridmaps[elecIndex]; }
    GridMap &getDesolvationMap()            { return gridmaps[desolvIndex]; }
    float *getFloatingGridMins() const      { return floatingGridMins; }

    // Read-only access to maps
    const GridMap &operator [](int i) const { return gridmaps[i]; }
    const GridMap &getElectrostaticMap()const{ return gridmaps[elecIndex]; }
    const GridMap &getDesolvationMap() const{ return gridmaps[desolvIndex]; }

    int getNumAtomMaps() const              { return numAtomMaps; }
    int getNumMaps() const                  { return numMaps; }
    int getNumMapsInclFloatingGrid() const;
    int getElectrostaticMapIndex() const    { return elecIndex; }
    int getDesolvationMapIndex() const      { return desolvIndex; }
    bool containsFloatingGrid() const       { return useFloatingGrid; }

private:
    LogFile *logFile;

    int numMaps, numAtomMaps, elecIndex, desolvIndex;
    GridMap *gridmaps;

    bool useFloatingGrid;
    float *floatingGridMins;
    int numGridPointsPerMap;
    char fileHeader[1<<14];
    int fileHeaderLength;

    void calculateEnergyMinMax(int map, double &energyMin, double &energyMax);
    void saveGridMap(int index) const;
};
