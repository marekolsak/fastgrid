/*
    Linear Algebra / Math library

    Copyright (C) 2003-2009, Marek Olsak (maraeo@gmail.com), All Rights Reserved.
    Copyright (C) 2003-2005, Tomas Pastorek (tomas@tomaspastorek.cz), All Rights Reserved.

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

#include "../All.h"

namespace Rune
{
    /**
        Spocita bounding sphere, pouziva 3 ruzne zpusoby, vypocet neni uplne presny
    **************************************************************************************************/
    template<typename T>
    void Sphere3<T>::ApproximateBestOf3(const Vec3<T> *vertices, int count)
    {
        Sphere3<T> aabb, av, md;
        aabb.ApproximateFromAABB(vertices, count);
        av.ApproximateFromAverage(vertices, count);
        md.ApproximateFromMostDistantVerts(vertices, count);

        // Nastavi tu nejmensi
        if (aabb.radius < av.radius && aabb.radius < md.radius)
            *this = aabb;
        else if (av.radius < md.radius)
            *this = av;
        else
            *this = md;
    }

    template<typename T>
    void Sphere3<T>::ApproximateFromAABB(const Vec3<T> *vertices, int count)
    {
        pos = FindAABBCenter(vertices, count);
        radius = CalculateRadius(pos, vertices, count);
    }

    template<typename T>
    void Sphere3<T>::ApproximateFromAABB(const Sphere3<T> *spheres, int count)
    {
        // TODO: tato aproximace je hodne nepresna, pro presnejsi aproximaci by bylo lepsi vyuzit knihovnu CGAL (je pod LGPL)
        Vec3<T> *vertices = new Vec3<T>[count];
        T maxR = 0;
        for (int i = 0; i < count; i++)
        {
            vertices[i] = spheres[i].pos;
            if (maxR < spheres[i].radius)
                maxR = spheres[i].radius;
        }

        ApproximateFromAABB(vertices, count);
        delete [] vertices;
        radius += maxR;
    }

    template<typename T>
    void Sphere3<T>::ApproximateFromAverage(const Vec3<T> *vertices, int count)
    {
        pos = FindAverageCenter(vertices, count);
        radius = CalculateRadius(pos, vertices, count);
    }

    template<typename T>
    void Sphere3<T>::ApproximateFromMostDistantVerts(const Vec3<T> *vertices, int count)
    {
        pos = FindCenterFromMostDistantVerts(vertices, count);
        radius = CalculateRadius(pos, vertices, count);
    }

    template<typename T>
    Vec3<T> Sphere3<T>::FindAABBCenter(const Vec3<T> *vertices, int count)
    {
        AxisAlignedBox3<T> box;
        box.Approximate(vertices, count);
        return Vec3<T>::Center(box.min, box.max);
    }

    template<typename T>
    Vec3<T> Sphere3<T>::FindAverageCenter(const Vec3<T> *vertices, int count)
    {
        Vec3<T> result = *vertices;
        const Vec3<T> *it, *last = vertices+count;
        for (it = vertices+1; it != last; ++it)
            result += *it;
        result *= 1 / T(count);
        return result;
    }

    template<typename T>
    Vec3<T> Sphere3<T>::FindCenterFromMostDistantVerts(const Vec3<T> *vertices, int count)
    {
        T maxDist = 0;
        Vec3<T> result;
        const Vec3<T> *it, *it2, *last = vertices+count;
        for (it = vertices; it != last; ++it)
            for (it2 = it; it2 != last; ++it2)
            {
                T dist = Vec3<T>::DistanceSqr(*it, *it2);
                if (dist > maxDist)
                {
                    maxDist = dist;
                    result = Vec3<T>::Center(*it, *it2);
                }
            }
        return result;
    }

    template<typename T>
    T Sphere3<T>::CalculateRadius(const Vec3<T> &center, const Vec3<T> *vertices, int count)
    {
        T maxRadius = 0;
        const Vec3<T> *it, *last = vertices+count;
        for (it = vertices; it != last; ++it)
        {
            T radius = Vec3<T>::DistanceSqr(*it, center);
            if (radius > maxRadius) maxRadius = radius;
        }
        return Math<T>::Sqrt(maxRadius);
    }

    template
    class RUNEMATH_API Sphere3<float>;
    template
    class RUNEMATH_API Sphere3<double>;
}
