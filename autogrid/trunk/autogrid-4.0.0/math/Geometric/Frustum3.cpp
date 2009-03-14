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
        Nastavi 6 ploch definujici komoly jehlan pohledu kamery za predpokladu, ze "m"
        obsahuje projekcni matici, pripadne muze taky obsahovat soucin projekcni a modelovaci matice
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API void Frustum3<T>::Calculate(const Matrix4<T> &m)
    {
        planes[0].Set(m[3] - m[0], m[7] - m[4], m[11] - m[8], m[15] - m[12]);
        planes[0].Normalize();
        planes[1].Set(m[3] + m[0], m[7] + m[4], m[11] + m[8], m[15] + m[12]);
        planes[1].Normalize();
        planes[2].Set(m[3] + m[1], m[7] + m[5], m[11] + m[9], m[15] + m[13]);
        planes[2].Normalize();
        planes[3].Set(m[3] - m[1], m[7] - m[5], m[11] - m[9], m[15] - m[13]);
        planes[3].Normalize();
        planes[4].Set(m[3] - m[2], m[7] - m[6], m[11] - m[10], m[15] - m[14]);
        planes[4].Normalize();
        planes[5].Set(m[3] - m[2], m[7] - m[6], m[11] - m[10], m[15] - m[14]);
        planes[5].Normalize();
    }

    template RUNEMATH_API void Frustum3<float>::Calculate(const Matrix4<float> &m);
    template RUNEMATH_API void Frustum3<double>::Calculate(const Matrix4<double> &m);

    /**
        Transformuje komoly jehlan
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API Frustum3<T> operator *(const Matrix3<T> &m, const Frustum3<T> &f)
    {
        Frustum3<T> result;
        for (int i = 0; i < 6; i++)
            result.planes[i] = m * f.planes[i];
        return result;
    }

    template RUNEMATH_API Frustum3<float> operator *(const Matrix3<float> &m, const Frustum3<float> &f);
    template RUNEMATH_API Frustum3<double> operator *(const Matrix3<double> &m, const Frustum3<double> &f);

    /**
        Transformuje komoly jehlan transponovanou matici
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API Frustum3<T> operator *(const Frustum3<T> &f, const Matrix3<T> &m)
    {
        Frustum3<T> result;
        for (int i = 0; i < 6; i++)
            result.planes[i] = f.planes[i] * m;
        return result;
    }

    template RUNEMATH_API Frustum3<float> operator *(const Frustum3<float> &f, const Matrix3<float> &m);
    template RUNEMATH_API Frustum3<double> operator *(const Frustum3<double> &f, const Matrix3<double> &m);

    /**
        Transformuje komoly jehlan
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API Frustum3<T> operator *(const Matrix4<T> &m, const Frustum3<T> &f)
    {
        Frustum3<T> result;
        for (int i = 0; i < 6; i++)
            result.planes[i] = m * f.planes[i];
        return result;
    }

    template RUNEMATH_API Frustum3<float> operator *(const Matrix4<float> &m, const Frustum3<float> &f);
    template RUNEMATH_API Frustum3<double> operator *(const Matrix4<double> &m, const Frustum3<double> &f);

    /**
        Transformuje komoly jehlan transponovanou matici
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API Frustum3<T> operator *(const Frustum3<T> &f, const Matrix4<T> &m)
    {
        Frustum3<T> result;
        for (int i = 0; i < 6; i++)
            result.planes[i] = f.planes[i] * m;
        return result;
    }

    template RUNEMATH_API Frustum3<float> operator *(const Frustum3<float> &f, const Matrix4<float> &m);
    template RUNEMATH_API Frustum3<double> operator *(const Frustum3<double> &f, const Matrix4<double> &m);
}
