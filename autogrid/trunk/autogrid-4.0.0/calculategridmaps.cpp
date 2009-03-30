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

#include "calculategridmaps.h"
#include "spatialgrid.h"
#include "nearestneighborsearch3.h"
#include "utils.h"

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
            Vec3d distance = input->receptorAtomCoord[ia].xyz - gridPos;
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
        // NEW Directional N acceptor
    }
    else if (input->hbond[ia] == A1)
    {
        // A1
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

static inline void calculateGridPoint(const InputData *input, GridMapList &gridmaps, int hydrogen,
                       const PairwiseInteractionEnergies &energyLookup, const DesolvExpFunc &desolvExpFunc, const BondVectors *bondVectors,
                       const SpatialGrid &cutoffGrid, const int *indicesHtoA, const NearestNeighborSearch3d &hsearch,
                       const Vec3d &gridPos, int outputIndex)
{
    HBondInfo hbond(gridmaps.getNumAtomMaps());

#if defined(USE_NNS)
    int closestH = indicesHtoA[hsearch.searchNearest(gridPos)];
#else
    int closestH = findClosestHBond(input, gridPos);
#endif

    //  Do all Receptor (protein, DNA, etc.) atoms...
#if defined(USE_SPATIAL_GRID)
    SpatialCell cell = cutoffGrid.getCellByPos(gridPos);
    int num = cutoffGrid.getNumElements(cell);

    for (int index = 0; index < num; index++)
    {
        int ia = cutoffGrid.getElement(cell, index);
#else
    int num = input->numReceptorAtoms;
    for (int ia = 0; ia < num; ia++)
    {
#endif
        // If distance from grid point to atom ia is too large,
        // or if atom is a disordered hydrogen,
        //   add nothing to the grid-point's non-bond energy;
        //   just continue to next atom...

        //  Get distance from current grid point to this receptor atom
        Vec3d distance = input->receptorAtomCoord[ia].xyz - gridPos;

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

        int indexR = lookup(1 / invR);

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

static void initCutoffGrid(const InputData *input, SpatialGrid &cutoffGrid)
{
    Vec3d gridSize = Vec3d(input->numGridPoints) * input->gridSpacing;
    double cellSize = NBCUTOFF / 4.0;

    // TODO: reduce the bucket size (maxElementsInCell) according to a density of atoms
    cutoffGrid.create(gridSize, cellSize, 0, input->numReceptorAtoms);
    cutoffGrid.insertSpheres(input->numReceptorAtoms, &input->receptorAtomCoord[0], NBCUTOFF);
}

static void initHSearch(const InputData *input, NearestNeighborSearch3d &hsearch, int *indicesHtoA)
{
    int numH = 0;
    Vec3d *hcoord = new Vec3d[input->numReceptorAtoms];

    for (int ia = 0; ia < input->numReceptorAtoms; ia++)
        if (input->hbond[ia] == DS || input->hbond[ia] == D1)
        {
            hcoord[numH] = input->receptorAtomCoord[ia].xyz;
            indicesHtoA[numH] = ia;
            ++numH;
        }

    hsearch.create(hcoord, numH, true);
}

void calculateGridmaps(const InputData *input, GridMapList &gridmaps, const ParameterLibrary &parameterLibrary,
                       const PairwiseInteractionEnergies &energyLookup, const DesolvExpFunc &desolvExpFunc, const BondVectors *bondVectors)
{
    SpatialGrid cutoffGrid;
    NearestNeighborSearch3d hsearch;
    int *indicesHtoA = new int[input->numReceptorAtoms]; // table for translating a H index into a receptor atom index

#if defined(USE_SPATIAL_GRID)
    // Create a grid for a cutoff of distant receptor atoms
    Timer *t0 = Timer::startNew("SGRIDGEN");
    initCutoffGrid(input, cutoffGrid);
    t0->stopAndLog(stderr);
#endif

#if defined(USE_NNS)
    // Create a tree for finding the closest H
    Timer *t1 = Timer::startNew("HTREEGEN");
    initHSearch(input, hsearch, indicesHtoA);
    t1->stopAndLog(stderr);
#endif

    int hydrogen = parameterLibrary.getAtomParameterRecIndex("HD");

    Timer *t2 = Timer::startNew("ATOMMAPS");
#if defined(AG_OPENMP)
    #pragma AG_OPENMP_PARALLEL_FOR
#endif
    FOR_EACH_GRID_POINT(gridPos, outputIndex)
    {
        calculateGridPoint(input, gridmaps, hydrogen, energyLookup, desolvExpFunc, bondVectors,
                           cutoffGrid, indicesHtoA, hsearch, gridPos, outputIndex);
    }
    END_FOR();
    t2->stopAndLog(stderr);

    delete [] indicesHtoA;
}
