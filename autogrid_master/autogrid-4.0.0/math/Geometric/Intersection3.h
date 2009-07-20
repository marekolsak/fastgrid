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
    enum Interior
    {
        SOLID = 0,
        HOLLOW
    };

    template<typename T>
    RUNEMATH_API bool Intersect(const AxisAlignedBox3<T> &b, const Vec3<T> &p);
    template<typename T>
    RUNEMATH_API bool Intersect(const AxisAlignedBox3<T> &b1, const AxisAlignedBox3<T> &b2);
    template<typename T>
    RUNEMATH_API bool Intersect(const AxisAlignedBox3<T> &b, const Sphere3<T> &s);
    template<typename T>
    RUNEMATH_API bool Intersect(const AxisAlignedBox3<T> &b, Interior bi, const Sphere3<T> &s, Interior si);
    template<typename T>
    RUNEMATH_API bool Intersect(const Sphere3<T> s1, const Sphere3<T> &s2);
    template<typename T>
    RUNEMATH_API bool Intersect(const Sphere3<T> &s1, const Sphere3<T> &s2, Vec3<T> &point);
    template<typename T>
    RUNEMATH_API bool Intersect(const Sphere3<T> &s, const Vec3<T> &point);
}
