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
#include <cstring>
#include "autogrid.h"
#include "Exceptions.h"

typedef void *SpatialCell;

class SpatialGrid
{
public:
    typedef int32 T;

    SpatialGrid(): grid(0) {}

	~SpatialGrid()
	{	
		if (grid)
		{
			for (int i = 0; i < numAllCells; i++)
				if (grid[i])
					delete [] grid[i];
			delete [] grid;
		}
	}

    // Creates the grid with the given size of the grid and the maximal expected size of each cell
    void create(const Vec3d &gridSize, double maxCellSize, const Vec3d &centerPos)
    {
        Vec3d minNumCells = gridSize / maxCellSize;
        Vec3d numCells = Vec3d(Mathd::Ceil(minNumCells.x),
                               Mathd::Ceil(minNumCells.y),
                               Mathd::Ceil(minNumCells.z));
        Vec3d cellSize = gridSize / numCells;

        create(Vec3i(numCells + 0.5), cellSize, centerPos);
    }

	// Returns the size of memory
	size_t getSizeOfMemory() const
	{
		if (!grid)
			return 0;

		size_t size = numAllCells * sizeof(T*);
		for (int i = 0; i < numAllCells; i++)
			size += (grid[i][0]+2) * sizeof(T);
		return size;
	}

    // Creates the grid with the given number of cells and the size of each cell
    void create(const Vec3i &numCells, const Vec3d &cellSize, const Vec3d &centerPos)
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

        // Allocate the grid
        numAllCells = this->numCells.Cube();
        grid = new T*[numAllCells];
		memset(grid, 0, sizeof(T*) * numAllCells);

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
    SpatialCell getCellFromIndices(const Vec3i &indices) const
    {
        Vec3i in = indices;
        clampIndices(in);
        return grid[getScalarIndexFromIndices(in)];
    }

    // Calculates a scalar index from 3D indices
    int getScalarIndexFromIndices(const Vec3i &in) const
    {
        return numCells.x * (numCells.y * in.z + in.y) + in.x;
    }

    // Returns indices of the cell at the given position
    void getIndicesFromPos(const Vec3d &pos, Vec3i &indices) const
    {
        indices = Vec3i((pos - cornerPosMin) * cellSizeInv);
    }

    // Returns the corner of the cell at the given indices
    void getCellCornerPosFromIndices(const Vec3i &indices, Vec3d &pos) const
    {
        pos = Vec3d(indices) * cellSize + cornerPosMin;
    }

    // Returns the internal cell ID at the given position
    SpatialCell getCellFromPos(const Vec3d &pos) const
    {
        Vec3i ipos;
        getIndicesFromPos(pos, ipos);
        return getCellFromIndices(ipos);
    }

    // Returns the number of elements in the cell
    T getNumElements(const SpatialCell cell) const
    {
        return *((T*)cell + 1);
    }

    // Returns the index-th element in the given cell
    T getElement(const SpatialCell cell, int index) const
    {
        return *((T*)cell + 2 + index);
    }

    // Inserts ID of your object to the cell, which occupies the given cell
    void insertIntoCell(SpatialCell cell, T id)
    {
		T *c = (T*)cell;

		T maxNum = c[0]; // The maximum number of elements
        T &num = c[1]; // The actual number of elements
		T *elements = c+2;

        if (num < maxNum)
        {
            elements[num] = id;
            ++num;
        }
        else
        {
            fprintf(stderr, "SpatialGrid::insertIntoCell: cell overflow!\n");
            throw ExitProgram(1);
        }
    }

    // Inserts ID of your object to the cell at the given indices
    void insertAtIndices(const Vec3i &indices, T id)
    {
        insertIntoCell(getCellFromIndices(indices), id);
    }

    // Inserts ID of your object to the cell at the given indices
    void insertAtClampedIndices(const Vec3i &indices, T id)
    {
        insertIntoCell(grid[getScalarIndexFromIndices(indices)], id);
    }

    // Inserts ID to the cell at the given position
    void insertPos(const Vec3d &pos, T id)
    {
        insertIntoCell(getCellFromPos(pos), id);
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
        getIndicesFromPos(s.pos - s.radius, indicesMin);
        getIndicesFromPos(s.pos + s.radius, indicesMax);
        clampIndices(indicesMin);
        clampIndices(indicesMax);

        // Too simple to be parallelized
        for (int x = indicesMin.x; x <= indicesMax.x; x++)
            for (int y = indicesMin.y; y <= indicesMax.y; y++)
                for (int z = indicesMin.z; z <= indicesMax.z; z++)
                {
                    Vec3i indices(x, y, z);
                    Vec3d cellPos;
                    getCellCornerPosFromIndices(indices, cellPos);

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
                    getCellCornerPosFromIndices(indices, cellPos);
                    boxes[getScalarIndexFromIndices(indices)] = AxisAlignedBox3d(cellPos, cellPos + cellSize);
                }

        // For each sphere, precalculate a range of possible indices of cells the sphere might intersect with
        Vec3i *sphereIndicesRange = new Vec3i[num*2];
#if defined(AG_OPENMP)
    #pragma omp parallel for schedule(dynamic, 1024)
#endif
        for (int s = 0; s < num; s++)
        {
            int sindex = s*2;
            getIndicesFromPos(Vec3d(pos[s]) - radius, sphereIndicesRange[sindex]);
            getIndicesFromPos(Vec3d(pos[s]) + radius, sphereIndicesRange[sindex+1]);
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

		// FIRST DETERMINE THE SIZE OF EACH CELL
		assert(sizeof(size_t) == sizeof(T*));
		size_t *counts = (size_t*)grid;
		memset(counts, 0, sizeof(size_t) * numAllCells);

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
                            int index = getScalarIndexFromIndices(Vec3i(x, y, z));
                            if (Intersect(boxes[index], Sphere3d(Vec3d(pos[s]), radius)))
                                ++counts[index];
                        }
            }

		// ALLOCATE CELLS
		for (int i = 0; i < numAllCells; i++)
		{
			T *cell = new T[2 + counts[i]];
			cell[0] = T(counts[i]);
			cell[1] = 0;
			grid[i] = cell;
		}

		// THEN, INSERT ATOMS INTO THE CELLS
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
                            int index = getScalarIndexFromIndices(Vec3i(x, y, z));
                            if (Intersect(boxes[index], Sphere3d(Vec3d(pos[s]), radius)))
                                insertIntoCell(grid[index], T(s));
                        }
            }

        delete [] sphereIndicesRange;
        delete [] boxes;
    }

private:
    T **grid;
    Vec3i numCells, numCellsMinus1;
	int numAllCells;
    Vec3d cellSize, cellSizeHalf, cellSizeInv;
    Vec3d cornerPosMin, cornerCellPosMin;
    double cellDiagonalHalf;
};
