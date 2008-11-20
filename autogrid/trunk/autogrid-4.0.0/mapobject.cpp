#include "mapobject.h"

MapObject::MapObject()
{
    atomType = 0;   // corresponds to receptor numbers????
    mapIndex = 0;
    isCovalent = 0;
    isHBonder = 0;
    mapFile = 0;
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
