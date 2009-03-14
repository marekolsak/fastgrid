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

#pragma once

namespace Rune
{
    /**
        Bod v prostoru, obalujici trida pro typove odliseni nasobeni vektoru a vertexu matici
    **************************************************************************************************/
    template<typename T>
    class Vertex3 : public Vec3<T>
    {
    public:
        typedef Rune::Vec3<T> Vec3;

        Vertex3() {}
        Vertex3(const Vec3 &V): Vec3(V) {}
    };

    template<typename T>
    RUNEMATH_API Vertex3<T> operator *(const Matrix3<T> &m, const Vertex3<T> &v);
    template<typename T>
    RUNEMATH_API Vertex3<T> operator *(const Vertex3<T> &v, const Matrix3<T> &m);
    template<typename T>
    RUNEMATH_API Vertex3<T> operator *(const Matrix4<T> &m, const Vertex3<T> &v);
    template<typename T>
    RUNEMATH_API Vertex3<T> operator *(const Vertex3<T> &v, const Matrix4<T> &m);

    typedef Vertex3<float> Vertex3f;
    typedef Vertex3<double> Vertex3d;
    typedef Vertex3<int32> Vertex3i;
}
