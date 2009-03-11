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

    // Creates the grid
    void create(const Vec3i &size, double cellSize, const Vec3d &centerPos, int maxElementsInCell)
    {
        this->size = size;
        for (int i = 0; i < 3; i++)
            this->size[i] = Mathi::Max(1, this->size[i]);
        this->sizeMinus1 = this->size - 1;

        this->cellSize = cellSize;
        cellSizeHalf = cellSize * 0.5;
        cellSizeInv = 1 / cellSize;
        this->maxElementsInCell = maxElementsInCell;

        // The actual number of elements will be at the beginning of every cell
        sizeofCell = 1 + maxElementsInCell;
        sizeofCell = (((sizeofCell * sizeof(T) - 1) / 8 + 1) * 8) / sizeof(T); // align to a multiple of 8 bytes
        int numCells = this->size.x * this->size.y * this->size.z;
        grid = new T[numCells * sizeofCell];

        // Zeroize counts
        T *end = grid + numCells * sizeofCell;
        for (T *p = grid; p != end; p += sizeofCell)
            *p = 0;

        // extent is a vector from the corner to the center
        Vec3d extent = Vec3d(this->size) * cellSizeHalf;

        // cornerPosMin is a position of the corner which has minimal coordinates
        cornerPosMin = centerPos - extent;

        // cornerCellPosMin is a position at the center of the corner cell
        cornerCellPosMin = cornerPosMin + cellSizeHalf;
    }

    void clampIndices(Vec3i &in)
    {
        for (int i = 0; i < 3; i++)
            in[i] = Mathi::Clamp(in[i], 0, sizeMinus1[i]);
    }

    // Returns the internal cell ID, which is actually the pointer to the grid array
    SpatialCell getCellByIndices(const Vec3i &indices)
    {
        Vec3i in = indices;
        clampIndices(in);

        int index = (size.x * (size.y * in.z + in.y) + in.x) * sizeofCell;
        return grid + index;
    }

    // Returns indices of the cell at the given position
    void getIndicesByPos(const Vec3d &pos, Vec3i &indices) const
    {
        indices = Vec3i((pos - cornerPosMin) * cellSizeInv);
    }

    // Returns the center of the cell at the given indices
    void getCellPosByIndices(const Vec3i &indices, Vec3d &pos) const
    {
        pos = Vec3d(indices) * cellSize + cornerCellPosMin;
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

    // Returns the size[XYZ] of the grid
    const Vec3i &getSize() const
    {
        return size;
    }

    // Inserts ID to all cells which are in the range of the sphere
    void insertSphere(const Sphere3d &s, T id)
    {
        Vec3i indicesMin, indicesMax;
        getIndicesByPos(s.pos - s.radius, indicesMin);
        getIndicesByPos(s.pos + s.radius, indicesMax);
        clampIndices(indicesMin);
        clampIndices(indicesMax);

        for (int x = indicesMin.x; x <= indicesMax.x; x++)
            for (int y = indicesMin.y; y <= indicesMax.y; y++)
                for (int z = indicesMin.z; z <= indicesMax.z; z++)
                {
                    Vec3i indices(x, y, z);
                    Vec3d cellPos;
                    getCellPosByIndices(indices, cellPos);
                    AxisAlignedBox3d b(cellPos - cellSizeHalf, cellPos + cellSizeHalf);

                    if (Intersect(b, s))
                        insertAtIndices(indices, id);
                }
    }

private:
    T *grid;
    Vec3i size, sizeMinus1;
    int maxElementsInCell, sizeofCell;
    double cellSize, cellSizeHalf, cellSizeInv;
    Vec3d cornerPosMin, cornerCellPosMin;
};
