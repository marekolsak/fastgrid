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
#include "autogrid.h"

typedef void *SpatialCell;

class SpatialGrid
{
public:
    typedef int32 T;

    SpatialGrid(): grid(0) {}
    ~SpatialGrid() { if (grid) delete [] grid; }

    // Creates the grid with the given size of the grid and the maximal expected size of each cell
    void create(const Vec3d &gridSize, double maxCellSize, const Vec3d &centerPos, int maxElementsInCell)
    {
        Vec3d minNumCells = gridSize / maxCellSize;
        Vec3d numCells = Vec3d(Mathd::Ceil(minNumCells.x),
                               Mathd::Ceil(minNumCells.y),
                               Mathd::Ceil(minNumCells.z));
        Vec3d cellSize = gridSize / numCells;

        create(Vec3i(numCells + 0.5), cellSize, centerPos, maxElementsInCell);
    }

    // Creates the grid with the given number of cells and the size of each cell
    void create(const Vec3i &numCells, const Vec3d &cellSize, const Vec3d &centerPos, int maxElementsInCell)
    {
        this->numCells = numCells;
        for (int i = 0; i < 3; i++)
            this->numCells[i] = Mathi::Max(1, this->numCells[i]);
        this->numCellsMinus1 = this->numCells - 1;

        this->cellSize = cellSize;
        cellSizeHalf = cellSize * 0.5;
        cellDiagonalHalf = cellSizeHalf.Magnitude();
        cellSizeInv = cellSize;
        cellSizeInv.Inverse();
        this->maxElementsInCell = maxElementsInCell;

        // The actual number of elements is at the beginning of every cell
        sizeofCell = 1 + maxElementsInCell;
        // Align to a multiple of 8 bytes
        sizeofCell = (((sizeofCell * sizeof(T) - 1) / 8 + 1) * 8) / sizeof(T);
        // Allocate the grid
        int numAllCells = this->numCells.x * this->numCells.y * this->numCells.z;
        grid = new T[numAllCells * sizeofCell];

        // Zeroize the count of elements in each grid
        T *end = grid + numAllCells * sizeofCell;
        for (T *p = grid; p != end; p += sizeofCell)
            *p = 0;

        // extent is a vector from the corner to the center
        Vec3d extent = Vec3d(this->numCells) * cellSizeHalf;

        // cornerPosMin is a position of the corner which has minimal coordinates
        cornerPosMin = centerPos - extent;

        // cornerCellPosMin is a position at the center of the corner cell
        cornerCellPosMin = cornerPosMin + cellSizeHalf;
    }

    // Set indices to be in the valid range
    void clampIndices(Vec3i &in) const
    {
        for (int i = 0; i < 3; i++)
            in[i] = Mathi::Clamp(in[i], 0, numCellsMinus1[i]);
    }

    // Returns an internal cell ID, which is actually a pointer to the grid array
    SpatialCell getCellByIndices(const Vec3i &indices) const
    {
        Vec3i in = indices;
        clampIndices(in);
        return getCellByClampedIndices(in);
    }

    // Returns an internal cell ID, which is actually a pointer to the grid array
    SpatialCell getCellByClampedIndices(const Vec3i &in) const
    {
        int index = getScalarIndexByIndices(in) * sizeofCell;
        return grid + index;
    }

    // Calculates a scalar index from 3D indices
    int getScalarIndexByIndices(const Vec3i &in) const
    {
        return numCells.x * (numCells.y * in.z + in.y) + in.x;
    }

    // Returns indices of the cell at the given position
    void getIndicesByPos(const Vec3d &pos, Vec3i &indices) const
    {
        indices = Vec3i((pos - cornerPosMin) * cellSizeInv);
    }

    // Returns the center of the cell at the given indices
    void getCellCenterPosByIndices(const Vec3i &indices, Vec3d &pos) const
    {
        pos = Vec3d(indices) * cellSize + cornerCellPosMin;
    }

    // Returns the corner of the cell at the given indices
    void getCellCornerPosByIndices(const Vec3i &indices, Vec3d &pos) const
    {
        pos = Vec3d(indices) * cellSize + cornerPosMin;
    }

    // Returns the internal cell ID at the given position
    SpatialCell getCellByPos(const Vec3d &pos) const
    {
        Vec3i ipos;
        getIndicesByPos(pos, ipos);
        return getCellByIndices(ipos);
    }

    // Returns the number of elements in the cell
    T getNumElements(const SpatialCell cell) const
    {
        return *(T*)cell;
    }

    // Returns the index-th element in the given cell
    T getElement(const SpatialCell cell, int index) const
    {
        return *((T*)cell + 1 + index);
    }

    // Inserts ID of your object to the cell, which occupies the given cell
    void insertIntoCell(SpatialCell cell, T id)
    {
        // The first index represents the number of elements
        T &num = *(T*)cell;

        if (num < maxElementsInCell)
        {
            // Elements begin at the second index
            T *elements = (T*)cell+1;
            elements[num] = id;
            ++num;
        }
    }

    // Inserts ID of your object to the cell at the given indices
    void insertAtIndices(const Vec3i &indices, T id)
    {
        insertIntoCell(getCellByIndices(indices), id);
    }

    // Inserts ID of your object to the cell at the given indices
    void insertAtClampedIndices(const Vec3i &indices, T id)
    {
        insertIntoCell(getCellByClampedIndices(indices), id);
    }

    // Inserts ID to the cell at the given position
    void insertPos(const Vec3d &pos, T id)
    {
        insertIntoCell(getCellByPos(pos), id);
    }

    // Returns the numCells vector of the grid
    const Vec3i &getNumCells() const
    {
        return numCells;
    }

    // For each cell: if the cell is in range of the sphere, insert ID into the cell.
    void insertSphere(const Sphere3d &s, T id)
    {
        Vec3i indicesMin, indicesMax;
        getIndicesByPos(s.pos - s.radius, indicesMin);
        getIndicesByPos(s.pos + s.radius, indicesMax);
        clampIndices(indicesMin);
        clampIndices(indicesMax);

        // Too simple to be parallelized
        for (int x = indicesMin.x; x <= indicesMax.x; x++)
            for (int y = indicesMin.y; y <= indicesMax.y; y++)
                for (int z = indicesMin.z; z <= indicesMax.z; z++)
                {
                    Vec3i indices(x, y, z);
                    Vec3d cellPos;
                    getCellCornerPosByIndices(indices, cellPos);

                    if (Intersect(AxisAlignedBox3d(cellPos, cellPos + cellSize), s))
                        insertAtClampedIndices(indices, id);
                }
    }

    // For each sphere, for each cell: if the cell is in range of the sphere, insert its index into the cell.
    template<typename Vec> // Restriction: Typecasting Vec -> Vec3d must be possible (e.g. allows using a Vec4d array)
    void insertSpheres(int num, const Vec *pos, double radius)
    {
        // For each cell, precalculate a box representing the cell in Euclidean space
        AxisAlignedBox3d *boxes = new AxisAlignedBox3d[numCells.Cube()];
#if defined(AG_OPENMP)
    #pragma omp parallel for schedule(dynamic, 4)
#endif
        for (int x = 0; x < numCells.x; x++)
            for (int y = 0; y < numCells.y; y++)
                for (int z = 0; z < numCells.z; z++)
                {
                    Vec3i indices(x, y, z);
                    Vec3d cellPos;
                    getCellCornerPosByIndices(indices, cellPos);
                    boxes[getScalarIndexByIndices(indices)] = AxisAlignedBox3d(cellPos, cellPos + cellSize);
                }

        // For each sphere, precalculate a range of possible indices of cells the sphere might intersect with
        Vec3i *sphereIndicesRange = new Vec3i[num*2];
#if defined(AG_OPENMP)
    #pragma omp parallel for schedule(dynamic, 1024)
#endif
        for (int s = 0; s < num; s++)
        {
            int sindex = s*2;
            getIndicesByPos(Vec3d(pos[s]) - radius, sphereIndicesRange[sindex]);
            getIndicesByPos(Vec3d(pos[s]) + radius, sphereIndicesRange[sindex+1]);
            clampIndices(sphereIndicesRange[sindex]);
            clampIndices(sphereIndicesRange[sindex+1]);
        }

        // Divide the X axis up to 4*CPUs segments to balance fluctuating density of spheres
#if defined(AG_OPENMP)
        int segments = omp_get_max_threads() * 4;
#else
        int segments = 1;
#endif
        int numCellsXPerSegment = Mathi::Max(1, int(Mathd::Ceil(double(numCells.x) / segments)));
        // The reason we do this is that iterating through all X elements concurrently and then through all the spheres
        // would be waste of computational time since not all spheres intersect the particular X coordinate, so we divide
        // the X axis to a reasonable number of segments and iterate only the ones which might really intersect the sphere.

#if defined(AG_OPENMP)
    #pragma omp parallel for schedule(dynamic, 1)
#endif
        // For each segment and each sphere
        for (int i = 0; i < segments; i++)
            for (int s = 0; s < num; s++)
            {
                // Get the range of indices for this sphere
                int sindex = s*2;
                Vec3i indicesMin = sphereIndicesRange[sindex];
                Vec3i indicesMax = sphereIndicesRange[sindex+1];

                // Adjust indices in the X axis for them to be in this segment
                indicesMin.x = Mathi::Max(i * numCellsXPerSegment, indicesMin.x);
                indicesMax.x = Mathi::Min((i+1) * numCellsXPerSegment - 1, indicesMax.x);

                // For each cell, insert the sphere if intersecting
                for (int x = indicesMin.x; x <= indicesMax.x; x++)
                    for (int y = indicesMin.y; y <= indicesMax.y; y++)
                        for (int z = indicesMin.z; z <= indicesMax.z; z++)
                        {
                            int index = getScalarIndexByIndices(Vec3i(x, y, z));
                            if (Intersect(boxes[index], Sphere3d(Vec3d(pos[s]), radius)))
                                insertIntoCell(grid + (index * sizeofCell), T(s));
                        }
            }

        delete [] sphereIndicesRange;
        delete [] boxes;
    }

private:
    T *grid;
    Vec3i numCells, numCellsMinus1;
    T maxElementsInCell, sizeofCell;
    Vec3d cellSize, cellSizeHalf, cellSizeInv;
    Vec3d cornerPosMin, cornerCellPosMin;
    double cellDiagonalHalf;
};
