#pragma once
#include "autogrid.h"
#include "math.h"

typedef void *SpatialCell;

template<typename T>    // T should be an integer (short, int, ...)
class SpatialGrid
{
public:
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
    void clampIndices(Vec3i &in)
    {
        for (int i = 0; i < 3; i++)
            in[i] = Mathi::Clamp(in[i], 0, numCellsMinus1[i]);
    }

    // Returns the internal cell ID, which is actually the pointer to the grid array
    SpatialCell getCellByIndices(const Vec3i &indices)
    {
        Vec3i in = indices;
        clampIndices(in);

        int index = (numCells.x * (numCells.y * in.z + in.y) + in.x) * sizeofCell;
        return grid + index;
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
    SpatialCell getCellByPos(const Vec3d &pos)
    {
        Vec3i ipos;
        getIndicesByPos(pos, ipos);
        return getCellByIndices(ipos);
    }

    // Returns the number of elements in the cell
    int getNumElements(const SpatialCell cell) const
    {
        return *(T*)cell;
    }

    // Returns the index-th element in the given cell
    T getElement(const SpatialCell cell, int index) const
    {
        return *((T*)cell + 1 + index);
    }

    // Inserts ID of your object to the cell, which occupies the given cell
    void insertInCell(SpatialCell cell, T id)
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
        insertInCell(getCellByIndices(indices), id);
    }

    // Inserts ID to the cell at the given position
    void insertPos(const Vec3d &pos, T id)
    {
        insertInCell(getCellByPos(pos), id);
    }

    // Returns the numCells vector of the grid
    const Vec3i &getNumCells() const
    {
        return numCells;
    }

    // Inserts ID to all cells which are in the range of the sphere
    void insertSphere(const Sphere3d &s, T id)
    {
        Vec3i indicesMin, indicesMax;
        getIndicesByPos(s.pos - s.radius, indicesMin);
        getIndicesByPos(s.pos + s.radius, indicesMax);
        clampIndices(indicesMin);
        clampIndices(indicesMax);

#if defined(AG_OPENMP)
    #pragma AG_OPENMP_PARALLEL_FOR
#endif
        for (int x = indicesMin.x; x <= indicesMax.x; x++)
            for (int y = indicesMin.y; y <= indicesMax.y; y++)
                for (int z = indicesMin.z; z <= indicesMax.z; z++)
                {
                    Vec3i indices(x, y, z);
                    Vec3d cellPos;
                    getCellCenterPosByIndices(indices, cellPos);

                    bool test = Intersect(s, cellPos);
                    if (!test)
                        test = Intersect(AxisAlignedBox3d(cellPos - cellSizeHalf, cellPos + cellSizeHalf), s);
                    if (test)
                            insertAtIndices(indices, id);
                }
    }

private:
    T *grid;
    Vec3i numCells, numCellsMinus1;
    int maxElementsInCell, sizeofCell;
    Vec3d cellSize, cellSizeHalf, cellSizeInv;
    Vec3d cornerPosMin, cornerCellPosMin;
    double cellDiagonalHalf;
};
