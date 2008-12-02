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

#pragma once
#include "inputdata.h"

struct GridMap
{
    int atomType;          // corresponds to receptor numbers????
    int mapIndex;
    bool isCovalent;
    int isHBonder;
    FILE *file;
    char mapFilename[MAX_CHARS];
    char type[3];           // eg HD or OA or NA or N
    double constant;        // this will become obsolete
    double energyMax;
    double energyMin;
    double energy;
    double volProbe;
    double solparProbe;

    // new 6/28
    double Rij;
    double epsij;
    HBondType hbond;       // hbonding character:
    double RijHB;
    double epsijHB;

    // per receptor type parameters, ordered as in receptorTypes
    double nbpR[NUM_RECEPTOR_TYPES];   // radius of energy-well minimum
    double nbpEps[NUM_RECEPTOR_TYPES]; // depth of energy-well minimum
    int xA[NUM_RECEPTOR_TYPES]; // generally 12
    int xB[NUM_RECEPTOR_TYPES]; // 6 for non-hbonders 10 for h-bonders
    int hbonder[NUM_RECEPTOR_TYPES];

    GridMap();
    ~GridMap();
};

class GridMapList
{
public:
    GridMapList(LogFile *logFile);
    ~GridMapList();

    void initFileHeader(const InputData *input, const char *gridParameterFilename);
    const char *getFileHeader() const       { return fileHeader; }
    int getFileHeaderLength() const         { return fileHeaderLength; }

    // Allocates gridmaps.
    // "num" is the number of maps to be created: the number of ligand atom types, plus 1 for the electrostatic map,
    // plus 1 for the desolvation map.
    // Keep in mind that AutoDock can only read in MAX_MAPS maps.
    void setNumMaps(int num);

    // Writes out summary
    void logSummary();

    // Read/write access
    GridMap &operator [](int i)             { return gridmaps[i]; }
    GridMap &getElectrostaticMap()          { return gridmaps[elecIndex]; }
    GridMap &getDesolvationMap()            { return gridmaps[desolvIndex]; }

    // Read-only access
    const GridMap &operator [](int i) const { return gridmaps[i]; }
    const GridMap &getElectrostaticMap()const{ return gridmaps[elecIndex]; }
    const GridMap &getDesolvationMap() const{ return gridmaps[desolvIndex]; }

    int getNumAtomMaps() const              { return numAtomMaps; }
    int getNumMaps() const                  { return numMaps; }
    int getElectrostaticMapIndex() const    { return elecIndex; }
    int getDesolvationMapIndex() const      { return desolvIndex; }

private:
    int numMaps, numAtomMaps, elecIndex, desolvIndex;
    GridMap *gridmaps;
    LogFile *logFile;
    char fileHeader[1<<14];
    int fileHeaderLength;
};
