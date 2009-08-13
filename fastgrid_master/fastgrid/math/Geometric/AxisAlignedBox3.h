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
        Trida kvadru ve 3D, ktery je rovnobezny s osami (AABB = axis-aligned bounding box)
    **************************************************************************************************/
    template<typename T>
    class AxisAlignedBox3
    {
    public:
        Vec3<T> min, max;

        AxisAlignedBox3() {}
        AxisAlignedBox3(const Vec3<T> &Min, const Vec3<T> &Max): min(Min), max(Max) {}

        Vec3<T> GetCenter() const   { return Vec3<T>::Center(min, max); }
        Vec3<T> GetExtents() const  { return (max - min) * T(0.5); }

        RUNEMATH_API void GetVertices(Vec3<T> *vertices) const;
        RUNEMATH_API void SetCenterAndExtents(const Vec3<T> &center, const Vec3<T> &extents);
        RUNEMATH_API void Approximate(const Vec3<T> *vertices, int count);
        RUNEMATH_API void Approximate(const AxisAlignedBox3 *boxes, int count);
    };

    typedef AxisAlignedBox3<float> AxisAlignedBox3f;
    typedef AxisAlignedBox3<double> AxisAlignedBox3d;
}
