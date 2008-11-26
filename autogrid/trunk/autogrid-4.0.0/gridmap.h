#pragma once
#include "autogrid.h"
#include "linearfreeenergymodel.h"

class GridMap
{
public:
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

    // Allocates gridmaps.
    // "num" is the number of maps to be created: the number of ligand atom types, plus 1 for the electrostatic map,
    // plus 1 for the desolvation map.
    // Keep in mind that AutoDock can only read in MAX_MAPS maps.
    void setNumMaps(int num);

    // Read/write access
    GridMap &operator [](int i)             { return gridmaps[i]; }
    GridMap &getElectrostaticMap()          { return gridmaps[elecIndex]; }
    GridMap &getDesolvationMap()            { return gridmaps[desolvIndex]; }

    // Read-only access
    const GridMap &operator [](int i) const { return gridmaps[i]; }
    const GridMap &getElectrostaticMap() const { return gridmaps[elecIndex]; }
    const GridMap &getDesolvationMap() const { return gridmaps[desolvIndex]; }

    int getNumAtomMaps() const              { return numAtomMaps; }
    int getNumMaps() const                  { return numMaps; }
    int getElectrostaticMapIndex() const    { return elecIndex; }
    int getDesolvationMapIndex() const      { return desolvIndex; }

private:
    int numMaps, numAtomMaps, elecIndex, desolvIndex;
    GridMap *gridmaps;
    LogFile *logFile;
};
