#pragma once
#include <cstdio>
#include "parameters.h"
#include "autogrid.h"

class MapObject
{
public:
    int atomType;          // corresponds to receptor numbers????
    int mapIndex;
    int isCovalent;
    int isHBonder;
    FILE *mapFile;
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

    MapObject();
};
