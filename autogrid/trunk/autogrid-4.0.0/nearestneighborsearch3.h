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
            pointarray[i] = const_cast<double*>(&coords[i].x);

        tree = new ANNkd_tree(pointarray, num, 3);
    }

    int searchNearest(const Vec3<T> &point) const
    {
        int result;
        double d;
        const_cast<ANNkd_tree*>(tree)->annkSearch(const_cast<double*>(&point.x), 1, &result, &d, 0);
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
