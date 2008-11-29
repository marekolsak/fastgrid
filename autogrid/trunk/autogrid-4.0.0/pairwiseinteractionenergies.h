#pragma once
#include "gridmap.h"

// Calculates the lookup table of the pairwise interaction energies
class PairwiseInteractionEnergies
{
public:
    struct LookupTable
    {
        // MAX_DIST is not NBCUTOFF times 100 as it should be, it's a power of two for a little faster memory access,
        // values beyond this threshold are unused
        double table[NUM_RECEPTOR_TYPES][MAX_DIST][MAX_MAPS];
    };

    PairwiseInteractionEnergies();
    ~PairwiseInteractionEnergies();

    // calculates the table of pair-wise interaction energies
    void calculate(const GridMapList &gridmaps, LogFile &logFile,
                   int numReceptorTypes, const char (&receptorTypes)[NUM_RECEPTOR_TYPES][3], double rSmooth);

    // returns the precalculated value
    double operator ()(int atomType, int indexR, int mapIndex) const
    {
        return energyLookup->table[atomType][indexR][mapIndex];
    }

private:
    LookupTable *energyLookup;
};
