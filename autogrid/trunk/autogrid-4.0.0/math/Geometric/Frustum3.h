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
        Trida komoleho jehlanu
    **************************************************************************************************/
    template<typename T>
    class Frustum3
    {
    public:
        Plane3<T> planes[6];

        /**
            Nastavi 6 ploch definujici komoly jehlan pohledu kamery za predpokladu, ze "m"
            obsahuje projekcni matici, pripadne muze taky obsahovat soucin projekcni a modelovaci matice
        **************************************************************************************************/
        RUNEMATH_API void Calculate(const Matrix4<T> &m);
    };

    template<typename T>
    RUNEMATH_API Frustum3<T> operator *(const Matrix3<T> &m, const Frustum3<T> &f);
    template<typename T>
    RUNEMATH_API Frustum3<T> operator *(const Frustum3<T> &f, const Matrix3<T> &m);
    template<typename T>
    RUNEMATH_API Frustum3<T> operator *(const Matrix4<T> &m, const Frustum3<T> &f);
    template<typename T>
    RUNEMATH_API Frustum3<T> operator *(const Frustum3<T> &f, const Matrix4<T> &m);

    typedef Frustum3<float> Frustum3f;
    typedef Frustum3<double> Frustum3d;
}
