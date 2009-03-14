#pragma once
#include "ann/ann.h"

// Given a set of points, this class can find the nearest neighbor point from the given position in space.
template<typename T>
class NearestNeighborSearch3
{
public:
    NearestNeighborSearch3(): releaseCoords(false), coords(0), tree(0), pointarray(0) {}

    ~NearestNeighborSearch3()
    {
        if (tree)
            delete tree;
        if (pointarray)
            delete [] pointarray;
        if (releaseCoords && coords)
            delete [] coords;
    }

    void create(const Vec3<T> *points, int num, bool releaseCoordsOnDestroy = false)
    {
        releaseCoords = releaseCoordsOnDestroy;
        coords = points;
        pointarray = new ANNpoint[num];
        for (int i = 0; i < num; i++)
            pointarray[i] = const_cast<double*>(static_cast<const double*>(coords[i]));

        tree = new ANNkd_tree(pointarray, num, 3);
    }

    int searchNearest(const Vec3<T> &point)
    {
        int result;
        double d;
        tree->annkSearch(const_cast<double*>(static_cast<const double*>(point)), 1, &result, &d, 0);
        return result;
    }

private:
    bool releaseCoords;
    const Vec3<T> *coords;
    ANNkd_tree *tree;
    ANNpointArray pointarray;
};

typedef NearestNeighborSearch3<int32> NearestNeighborSearch3i;
typedef NearestNeighborSearch3<float> NearestNeighborSearch3f;
typedef NearestNeighborSearch3<double> NearestNeighborSearch3d;
