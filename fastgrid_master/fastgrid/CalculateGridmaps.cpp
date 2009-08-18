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

#include "CalculateGridmaps.h"
#include "SpatialGrid.h"
#include "NearestNeighborSearch3.h"
#include "Utils.h"

#define USE_SPATIAL_GRID
#define USE_NNS

struct HBondInfo
{
    struct Info
    {
        double min, max;
        bool flag;
    } info[MAX_MAPS];

    inline HBondInfo(int numAtomMaps)
    {
        Info def;
        def.min = 999999;
        def.max = -999999;
        def.flag = false;

        for (int m = 0; m < numAtomMaps; m++)
            info[m] = def;
    }

    inline Info &operator [](int i)
    {
        return info[i];
    }

    inline void insert(int mapIndex, double energy)
    {
        Info &i = info[mapIndex];
        i.min = Mathd::Min(i.min, energy);
        i.max = Mathd::Max(i.max, energy);
        i.flag = true;
    }
};

static inline int findClosestHBond(const InputData *input, const Vec3d &gridPos)
{
    int closestH = 0;
    double rminSq = Mathd::Sqr(999999.0);
    for (int ia = 0; ia < input->numReceptorAtoms; ia++)
        if (input->hbond[ia] == DS || input->hbond[ia] == D1)
        {
            // DS or D1
            Vec3d distance = Vec3d(input->receptorAtom[ia]) - gridPos;
            double rSq = distance.MagnitudeSqr();
            if (rSq < rminSq)
            {
                rminSq = rSq;
                closestH = ia;
            }
        }           // Hydrogen test
    return closestH;
}

static inline void getHBondAngularFunction(const InputData *input, const BondVectors *bondVectors, int ia, int closestH, const Vec3d &distance,
                             double &racc, double &rdon, double &Hramp) // output
{
    racc = 1;
    rdon = 1;
    Hramp = 1;   // Hramp ramps in Hbond acceptor probes

    //  cosTheta = distance dot bondVectors->rvector == cos(angle) subtended.
    double cosTheta = -Vec3d::Dot(distance, bondVectors->rvector[ia]);

    // Calculate racc/rdon/Hramp
    if (input->hbond[ia] == D1)
    {
        // D1
        //  ia-th receptor atom = Hydrogen (4 = H)
        //  => receptor H-bond donor, OH or NH.
        //  calculate racc for H-bond ACCEPTOR PROBES at this grid pt.

        racc = 0;

        //  H->current-grid-pt vector < 90 degrees from
        //  N->H or O->H vector,
        if (cosTheta > 0)
        {
            //  racc = [cos(theta)]^2.0 for N-H
            //  racc = [cos(theta)]^4.0 for O-H
            racc = cosTheta;

            switch (bondVectors->rexp[ia]) // racc = pow(cosTheta, bondVectors->rexp[ia]);
            {
            case 4:
                racc = Mathd::Sqr(racc);
            case 2:
                racc = Mathd::Sqr(racc);
            }

            // NEW2 calculate dot product of bond vector with bond vector of best input->hbond
            if (ia != closestH)
            {
                double theta = Vec3d::AngleUnnorm(bondVectors->rvector[closestH], bondVectors->rvector[ia]);
                Hramp = 0.5 - 0.5 * cos(theta * (120.0 / 90.0));
            }   // ia test
            // END NEW2 calculate dot product of bond vector with bond vector of best input->hbond
        }
        // endif (input->atomType[ia] == hydrogen)
    }
    else if (input->hbond[ia] == A1)
    {
        // NEW Directional N acceptor
        //  ia-th macromolecule atom = Nitrogen (4 = H)
        //  calculate rdon for H-bond Donor PROBES at this grid pt.

        //  H->current-grid-pt vector < 90 degrees from X->N vector
        rdon = 0;
        if (cosTheta > 0)
            rdon = Mathd::Sqr(cosTheta); // for H->N
        // endif (input->atomType[ia] == nitrogen)
        // end NEW Directional N acceptor
    }
    else if (input->hbond[ia] == A2)
    {
        rdon = 0;

        if (bondVectors->disorder[ia])
        {
            // A2
            // cylindrically disordered hydroxyl

            racc = 0;
            double theta = Mathd::Acos(cosTheta);
            if (theta <= 1.24791 + (PI / 2))
            {
                // 1.24791 rad = 180 deg minus C-O-H bond angle, ** 108.5 deg
                rdon = Mathd::Sqr(Mathd::Sqr(cos(theta - 1.24791)));    // pow(.., 4)
                racc = rdon;
            }
        }
        else
        {
            // A2
            //  ia-th receptor atom = Oxygen
            //  => receptor H-bond acceptor, oxygen.

            // rdon expressions from Goodford
            if (cosTheta >= 0)
            {
                // ti is the angle in the lone pair plane, away from the
                // vector between the lone pairs,
                // calculated as (grid vector CROSS lone pair plane normal)
                // DOT C=O vector - 90 deg
                Vec3d cross = Vec3d::Cross(distance, bondVectors->rvector2[ia]);
                double rd2 = cross.MagnitudeSqr();
                if (rd2 < APPROX_ZERO)
                    rd2 = APPROX_ZERO;
                double inv_rd = Mathd::Rsqrt(rd2);

                double ti = fabs(Mathd::Acos(inv_rd * Vec3d::Dot(cross, bondVectors->rvector[ia])) - (PI / 2));
                // the 2.0*ti can be replaced by (ti + ti) in: rdon = (0.9 + 0.1*sin(2.0*ti))*cos(t0);
                rdon = 0.9 + 0.1 * sin(ti + ti);
                // 0.34202 = cos (100 deg)
            }
            else if (cosTheta >= -0.34202)
                rdon = 562.25 * Mathd::Cube(0.116978 - Mathd::Sqr(cosTheta));

            // t0 is the angle out of the lone pair plane, calculated
            // as 90 deg - acos (vector to grid point DOT lone pair
            // plane normal)
            double t0 = (PI / 2) - Vec3d::AngleUnnorm(distance, bondVectors->rvector2[ia]);
            rdon *= cos(t0);

            // endif input->atomType == oxygen, not disordered
        }
    }
}

static inline void sumPairwiseInteractions(const InputData *input, const GridMapList &gridmaps, const PairwiseInteractionEnergies &energyLookup,
                             const DesolvExpFunc &desolvExpFunc, const BondVectors *bondVectors, HBondInfo &hbond,
                             int outputIndex, int m, int ia, int indexR, int hydrogen, double racc, double rdon, double Hramp)
{
    double &e = gridmaps[m].energies[outputIndex];
    double pwiEnergy = energyLookup(input->atomType[ia], indexR, m);

    if (gridmaps[m].isHBonder)
    {
        // PROBE forms H-bonds...

        // rsph ramps in angular dependence for distances with negative energy
        double rsph = Mathd::Saturate(pwiEnergy / 100);

        if ((gridmaps[m].hbond == AS || gridmaps[m].hbond == A2) && (input->hbond[ia] == DS || input->hbond[ia] == D1))
            // PROBE can be an H-BOND ACCEPTOR,
            e += (bondVectors->disorder[ia] ? energyLookup(hydrogen, Mathi::Max(0, indexR - 110), m) : pwiEnergy) *
                 Hramp * (racc + (1 - racc) * rsph);
        else if (gridmaps[m].hbond == A1 && (input->hbond[ia] == DS || input->hbond[ia] == D1))
            // A1 vs DS, D1
            hbond.insert(m, pwiEnergy * (racc + (1 - racc) * rsph));
        else if ((gridmaps[m].hbond == DS || gridmaps[m].hbond == D1) && input->hbond[ia] >= AS)
            // DS,D1 vs AS,A1,A2
            // PROBE is H-BOND DONOR,
            hbond.insert(m, pwiEnergy * (rdon + (1 - rdon) * rsph));
        else
            // hbonder PROBE-ia cannot form a H-bond...,
            e += pwiEnergy;
    }
    else
        // PROBE does not form H-bonds...,
        e += pwiEnergy;

    // add desolvation energy
    // forcefield desolv coefficient/weight in desolvExpFunc
    e += gridmaps[m].solparProbe * input->vol[ia] * desolvExpFunc(indexR) +
        (input->solpar[ia] + input->solparQ * fabs(input->charge[ia])) * gridmaps[m].volProbe * desolvExpFunc(indexR);
}

template<int UseNNS, int UseCutoffGrid>
static inline void calculateGridPoint(const InputData *input, GridMapList &gridmaps, int hydrogen,
                       const PairwiseInteractionEnergies &energyLookup, const DesolvExpFunc &desolvExpFunc, const BondVectors *bondVectors,
                       const SpatialGrid &cutoffGrid, const int *indicesHtoA, const NearestNeighborSearch3d &hsearch,
                       const Vec3d &gridPos, int outputIndex)
{
    HBondInfo hbond(gridmaps.getNumAtomMaps());

    int closestH;
    if (UseNNS)
        closestH = indicesHtoA[hsearch.searchNearest(gridPos)];
    else
        closestH = findClosestHBond(input, gridPos);

    //  Do all Receptor (protein, DNA, etc.) atoms...
    int num;
    SpatialCell cell = 0;
    if (UseCutoffGrid)
    {
        cell = cutoffGrid.getCellByPos(gridPos);
        num = cutoffGrid.getNumElements(cell);
    }
    else
        num = input->numReceptorAtoms;

    for (int index = 0; index < num; index++)
    {
        int ia = UseCutoffGrid ? cutoffGrid.getElement(cell, index) : index;

        // If distance from grid point to atom ia is too large,
        // or if atom is a disordered hydrogen,
        //   add nothing to the grid-point's non-bond energy;
        //   just continue to next atom...

        //  Get distance from current grid point to this receptor atom
        Vec3d distance = Vec3d(input->receptorAtom[ia]) - gridPos;

        // rSq = |distance|^2
        double rSq = distance.MagnitudeSqr();

        if (rSq > Mathd::Sqr(NBCUTOFF))
            continue;   // onto the next atom...
        if (input->atomType[ia] == hydrogen && bondVectors->disorder[ia])
            continue;   // onto the next atom...

        // Normalize the distance vector
        if (rSq < Mathd::Sqr(APPROX_ZERO*APPROX_ZERO))
            rSq = Mathd::Sqr(APPROX_ZERO*APPROX_ZERO);
        double invR = Mathd::Rsqrt(rSq);
        distance *= invR;

        int indexR = angstromToIndex<int>(1 / invR);

        double racc, rdon, Hramp;
        getHBondAngularFunction(input, bondVectors, ia, closestH, distance, racc, rdon, Hramp);

        // For each probe atom-type,
        // Sum pairwise interactions between each probe
        // at this grid point (gridPos[0:2])
        // and the current receptor atom, ia...
        for (int m = 0; m < gridmaps.getNumAtomMaps(); m++)
            if (!gridmaps[m].isCovalent)
                sumPairwiseInteractions(input, gridmaps, energyLookup, desolvExpFunc, bondVectors, hbond, outputIndex, m, ia, indexR, hydrogen, racc, rdon, Hramp);

        gridmaps.getDesolvationMap().energies[outputIndex] += input->solparQ * input->vol[ia] * desolvExpFunc(indexR);
    } // ia loop, over all receptor atoms...

    // Adjust maps of hydrogen-bonding atoms by adding largest and
    // smallest interaction of all 'pair-wise' interactions with receptor atoms
    for (int m = 0; m < gridmaps.getNumAtomMaps(); m++)
    {
        double &e = gridmaps[m].energies[outputIndex];
        if (hbond[m].flag)
            e += hbond[m].min + hbond[m].max;
        e = roundOutput(e);
    }
    double &e = gridmaps.getDesolvationMap().energies[outputIndex];
    e = roundOutput(e);
}

template<int UseNNS, int UseCutoffGrid>
static void LoopOverGrid(const InputData *input, GridMapList &gridmaps, int hydrogen,
                       const PairwiseInteractionEnergies &energyLookup, const DesolvExpFunc &desolvExpFunc, const BondVectors *bondVectors,
                       const SpatialGrid &cutoffGrid, const int *indicesHtoA, const NearestNeighborSearch3d &hsearch)
{
#if defined(AG_OPENMP)
    #pragma omp parallel for schedule(dynamic, 1)
#endif
    FOR_EACH_GRID_POINT(gridPos, outputIndex)
        calculateGridPoint<UseNNS, UseCutoffGrid>(input, gridmaps, hydrogen, energyLookup, desolvExpFunc, bondVectors,
                                                  cutoffGrid, indicesHtoA, hsearch, gridPos, outputIndex);
    END_FOR()
}

static void initCutoffGrid(const InputData *input, const ProgramParameters &programParams, SpatialGrid &cutoffGrid)
{
    // Find distance^2 between two closest atoms
#if defined(AG_OPENMP)
    int subsets = omp_get_max_threads()*4;
#else
    int subsets = 1;
#endif
    int countPerSubset = Mathi::Max(1, int(Mathd::Ceil(double(input->numReceptorAtoms) / subsets)));
    double *minDistances = new double[subsets];
    for (int s = 0; s < subsets; s++)
        minDistances[s] = 10000000;

    NearestNeighborSearch3d search;
    search.create(input->receptorAtom, input->numReceptorAtoms);

#if defined(AG_OPENMP)
    #pragma omp parallel for schedule(dynamic, 1)
#endif
    for (int s = 0; s < subsets; s++)
    {
        int start = s * countPerSubset;
        int end = Mathi::Min((s+1) * countPerSubset - 1, input->numReceptorAtoms - 1);

        for (int i = start; i < end; i++)
        {
            double dist = search.getDistanceSqrOfNearest2(input->receptorAtom[i]);
            if (dist < minDistances[s])
                minDistances[s] = dist;
        }
    }

    double minDistanceSq = 10000000;
    for (int s = 0; s < subsets; s++)
        if (minDistances[s] < minDistanceSq)
            minDistanceSq = minDistances[s];

    delete [] minDistances;

    // Estimate a lower bound of atom volume
    double minAtomRadius = sqrt(minDistanceSq)/2;
    double minAtomVolume = (4.0/3.0) * 3.14159265358979323846 * Mathd::Cube(minAtomRadius);
    double invMinAtomVolume = 1 / minAtomVolume;
    Vec3d gridSize = Vec3d(input->numGridPoints) * input->gridSpacing;

    // Calculate the cell size and the max atoms per cell according to memory usage
    double cellSize = NBCUTOFF/2;
    int maxAtomsPerCell = 0;
    for (int i = 0; i < 40; i++, cellSize += 0.5)
    {
        maxAtomsPerCell = int(Mathd::Cube(cellSize + NBCUTOFF*2) * invMinAtomVolume + 1); // + 1 because of rounding
        if ((cutoffGrid.estimateMemorySize(gridSize, cellSize, maxAtomsPerCell) >> 20) <= size_t(programParams.getCutoffGridMemoryLimit()))
            break;
    }

    if (programParams.benchmarkEnabled())
    {
        //fprintf(stderr, "Cutoff Grid: Min atom radius: %f\n", minAtomRadius);
        //fprintf(stderr, "Cutoff Grid: Max atoms per cell: %i\n", maxAtomsPerCell);
        //fprintf(stderr, "Cutoff Grid: Cell size: %f\n", cellSize);
        //fprintf(stderr, "Cutoff Grid: Allocating %" SIZE_T_FLAG "u MiB\n", cutoffGrid.estimateMemorySize(gridSize, cellSize, maxAtomsPerCell) >> 20);
    }

    cutoffGrid.create(gridSize, cellSize, 0, maxAtomsPerCell);
    cutoffGrid.insertSpheres(input->numReceptorAtoms, &input->receptorAtom[0], NBCUTOFF);
}

static void initHSearch(const InputData *input, NearestNeighborSearch3d &hsearch, int *indicesHtoA)
{
    int numH = 0;
    Vec3d *hcoord = new Vec3d[input->numReceptorAtoms];

    for (int ia = 0; ia < input->numReceptorAtoms; ia++)
        if (input->hbond[ia] == DS || input->hbond[ia] == D1)
        {
            hcoord[numH] = Vec3d(input->receptorAtom[ia]);
            indicesHtoA[numH] = ia;
            ++numH;
        }

    hsearch.create(hcoord, numH, true);
}

void calculateGridmaps(const InputData *input, const ProgramParameters &programParams, GridMapList &gridmaps, const ParameterLibrary &parameterLibrary,
                       const PairwiseInteractionEnergies &energyLookup, const DesolvExpFunc &desolvExpFunc, const BondVectors *bondVectors)
{
    SpatialGrid cutoffGrid;
    NearestNeighborSearch3d hsearch;
    int *indicesHtoA = new int[input->numReceptorAtoms]; // table for translating a H index into a receptor atom index

    if (programParams.useCutoffGrid())
    {
        // Create a grid for a cutoff of distant receptor atoms
        Timer *t0 = 0;
        if (programParams.benchmarkEnabled())
            t0 = Timer::startNew("Cutoff grid            ");
        initCutoffGrid(input, programParams, cutoffGrid);
        if (programParams.benchmarkEnabled())
            t0->stopAndLog(stderr);
    }

    if (programParams.useNNS())
    {
        // Create a tree for finding the closest H
        Timer *t1 = 0;
        if (programParams.benchmarkEnabled())
            t1 = Timer::startNew("KD-tree of hydrogens   ");
        initHSearch(input, hsearch, indicesHtoA);
        if (programParams.benchmarkEnabled())
            t1->stopAndLog(stderr);
    }

    int hydrogen = parameterLibrary.getAtomParameterRecIndex("HD");

    Timer *t2 = 0;
    if (programParams.benchmarkEnabled())
        t2 = Timer::startNew("Atom & desolvation maps");

    // Determine which template parameter to use based on program parameters
    // This enabled/disables the nearest neighbor search and the cutoff grid
    if (programParams.useNNS())
        if (programParams.useCutoffGrid())
            LoopOverGrid<1, 1>(input, gridmaps, hydrogen, energyLookup, desolvExpFunc, bondVectors, cutoffGrid, indicesHtoA, hsearch);
        else
            LoopOverGrid<1, 0>(input, gridmaps, hydrogen, energyLookup, desolvExpFunc, bondVectors, cutoffGrid, indicesHtoA, hsearch);
    else
        if (programParams.useCutoffGrid())
            LoopOverGrid<0, 1>(input, gridmaps, hydrogen, energyLookup, desolvExpFunc, bondVectors, cutoffGrid, indicesHtoA, hsearch);
        else
            LoopOverGrid<0, 0>(input, gridmaps, hydrogen, energyLookup, desolvExpFunc, bondVectors, cutoffGrid, indicesHtoA, hsearch);

    if (programParams.benchmarkEnabled())
        t2->stopAndLog(stderr);

    delete [] indicesHtoA;
}
