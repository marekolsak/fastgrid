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
#include "GridMap.h"

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
