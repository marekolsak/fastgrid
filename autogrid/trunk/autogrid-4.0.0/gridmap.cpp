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
    energyMax = 0;
    energyMin = 0;
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
