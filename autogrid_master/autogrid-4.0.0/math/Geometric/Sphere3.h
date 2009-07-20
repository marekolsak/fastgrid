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

#pragma once

namespace Rune
{
    /**
        Trida koule ve 3D
    **************************************************************************************************/
    template<typename T>
    class RUNEMATH_API Sphere3
    {
    public:
        Vec3<T> pos;
        T radius;

        Sphere3() {}
        Sphere3(const Vec3<T> &position, T sphereradius): pos(position), radius(sphereradius) {}

        void ApproximateBestOf3(const Vec3<T> *vertices, int count);
        void ApproximateFromAABB(const Vec3<T> *vertices, int count);
        void ApproximateFromAABB(const Sphere3<T> *spheres, int count);
        void ApproximateFromAverage(const Vec3<T> *vertices, int count);
        void ApproximateFromMostDistantVerts(const Vec3<T> *vertices, int count);

    private:
        static Vec3<T> FindAABBCenter(const Vec3<T> *vertices, int count);
        static Vec3<T> FindAverageCenter(const Vec3<T> *vertices, int count);
        static Vec3<T> FindCenterFromMostDistantVerts(const Vec3<T> *vertices, int count);
        static T CalculateRadius(const Vec3<T> &center, const Vec3<T> *vertices, int count);
    };

    typedef Sphere3<float> Sphere3f;
    typedef Sphere3<double> Sphere3d;
}
