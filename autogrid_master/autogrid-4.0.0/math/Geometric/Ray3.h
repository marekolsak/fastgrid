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
        Trida primky ve 3D
    **************************************************************************************************/
    template<typename T>
    class Ray3
    {
    public:
        Vec3<T> origin;
        Vec3<T> direction;

        Ray3() {}
        Ray3(const Vec3<T> &Origin, const Vec3<T> &Direction): origin(Origin), direction(Direction) {}

        RUNEMATH_API void SetSecondPoint(Vec3<T> &p);
    };

    template<typename T>
    RUNEMATH_API Ray3<T> operator *(const Matrix3<T> &m, const Ray3<T> &r);
    template<typename T>
    RUNEMATH_API Ray3<T> operator *(const Ray3<T> &r, const Matrix3<T> &m);
    template<typename T>
    RUNEMATH_API Ray3<T> operator *(const Matrix4<T> &m, const Ray3<T> &r);
    template<typename T>
    RUNEMATH_API Ray3<T> operator *(const Ray3<T> &r, const Matrix4<T> &m);

    typedef Ray3<float> Ray3f;
    typedef Ray3<double> Ray3d;
}
