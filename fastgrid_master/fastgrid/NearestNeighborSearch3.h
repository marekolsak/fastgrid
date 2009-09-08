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
#include "ann/ANN.h"

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
            switch (stride)
            {
            case 3:
                delete [] (Vec3<T>*)coords;
                break;

            case 4:
                delete [] (Vec4<T>*)coords;
                break;
            }
    }

    void create(const Vec3<T> *points, int num, bool releaseCoordsOnDestroy = false)
    {
        create(&points->x, 3, num, releaseCoordsOnDestroy);
    }

    void create(const Vec4<T> *points, int num, bool releaseCoordsOnDestroy = false)
    {
        create(&points->x, 4, num, releaseCoordsOnDestroy);
    }

    int searchNearest(const Vec3<T> &point) const
    {
        int result;
        double d;
        const_cast<ANNkd_tree*>(tree)->annkSearch(const_cast<double*>(&point.x), 1, &result, &d, 0);
        return result;
    }

    T getDistanceSqrOfNearest2(const Vec3<T> &point) const
    {
        int result[2];
        double d[2];
        const_cast<ANNkd_tree*>(tree)->annkSearch(const_cast<double*>(&point.x), 2, result, d, 0);
        return Vec3<T>::DistanceSqr(*(Vec3<T>*)&coords[result[0]*stride], *(Vec3<T>*)&coords[result[1]*stride]);
    }

private:
    bool releaseCoords;
    int stride;
    const T *coords;
    ANNkd_tree *tree;
    ANNpointArray pointarray;

    void create(const T *points, int stride, int num, bool releaseCoordsOnDestroy)
    {
        this->stride = stride;
        releaseCoords = releaseCoordsOnDestroy;
        coords = points;
        pointarray = new ANNpoint[num];
        for (int i = 0; i < num; i++)
            pointarray[i] = const_cast<double*>(&coords[i*stride]);

        tree = new ANNkd_tree(pointarray, num, 3);
    }
};

typedef NearestNeighborSearch3<int32> NearestNeighborSearch3i;
typedef NearestNeighborSearch3<float> NearestNeighborSearch3f;
typedef NearestNeighborSearch3<double> NearestNeighborSearch3d;
