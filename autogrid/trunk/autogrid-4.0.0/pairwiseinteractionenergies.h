#pragma once
#include "gridmap.h"

class PairwiseInteractionEnergies
{
public:
    // MAX_DIST is not NBCUTOFF times 100 as it should be, it's a power of two for a little faster memory access,
    // values beyond this threshold are unused
    typedef double LookupTable[NUM_RECEPTOR_TYPES][MAX_DIST][MAX_MAPS];

    // calculates the table of pair-wise interaction energies
    void calculate(const GridMapList &gridmaps, LogFile &logFile,
                   // TODO: unify these arguments
                   int numReceptorTypes, const char (&receptorTypes)[NUM_RECEPTOR_TYPES][3], double rSmooth);

    // returns the calculated lookup table, read-only access
    const LookupTable &table() const { return energyLookup; }

private:
    LookupTable energyLookup;
};
