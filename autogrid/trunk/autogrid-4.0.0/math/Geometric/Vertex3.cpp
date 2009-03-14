/*
    Auxiliary Math library

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
        Transformuje vertex
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API Vertex3<T> operator *(const Matrix3<T> &m, const Vertex3<T> &v)
    {
        return m * Vec3<T>(v);
    }

    template RUNEMATH_API Vertex3<float> operator *(const Matrix3<float> &m, const Vertex3<float> &v);
    template RUNEMATH_API Vertex3<double> operator *(const Matrix3<double> &m, const Vertex3<double> &v);

    /**
        Transformuje vertex transponovanou matici
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API Vertex3<T> operator *(const Vertex3<T> &v, const Matrix3<T> &m)
    {
        return Vec3<T>(v) * m;
    }

    template RUNEMATH_API Vertex3<float> operator *(const Vertex3<float> &v, const Matrix3<float> &m);
    template RUNEMATH_API Vertex3<double> operator *(const Vertex3<double> &v, const Matrix3<double> &m);

    /**
        Transformuje vertex (je ekvivalentni transformaci vektoru (x, y, z, 1))
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API Vertex3<T> operator *(const Matrix4<T> &m, const Vertex3<T> &v)
    {
        return m * Vec3<T>(v) + Vec3<T>(m[12], m[13], m[14]);
    }

    template RUNEMATH_API Vertex3<float> operator *(const Matrix4<float> &m, const Vertex3<float> &v);
    template RUNEMATH_API Vertex3<double> operator *(const Matrix4<double> &m, const Vertex3<double> &v);

    /**
        Transformuje vertex transponovanou matici (je ekvivalentni transformaci vektoru (x, y, z, 1))
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API Vertex3<T> operator *(const Vertex3<T> &v, const Matrix4<T> &m)
    {
        return Vec3<T>(v) * m + Vec3<T>(m[3], m[7], m[11]);
    }

    template RUNEMATH_API Vertex3<float> operator *(const Vertex3<float> &v, const Matrix4<float> &m);
    template RUNEMATH_API Vertex3<double> operator *(const Vertex3<double> &v, const Matrix4<double> &m);
}
