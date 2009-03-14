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
    template<typename T>
    RUNEMATH_API void Ray3<T>::SetSecondPoint(Vec3<T> &p)
    {
        direction = (p - origin).GetNormalized();
    }

    template RUNEMATH_API void Ray3<float>::SetSecondPoint(Vec3<float> &p);
    template RUNEMATH_API void Ray3<double>::SetSecondPoint(Vec3<double> &p);

    /**
        Transformuje primku v prostoru
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API Ray3<T> operator *(const Matrix3<T> &m, const Ray3<T> &r)
    {
        return Ray3<T>(m * r.origin, m * r.direction);
    }

    template RUNEMATH_API Ray3<float> operator *(const Matrix3<float> &m, const Ray3<float> &r);
    template RUNEMATH_API Ray3<double> operator *(const Matrix3<double> &m, const Ray3<double> &r);

    /**
        Transformuje primku v prostoru transponovanou matici
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API Ray3<T> operator *(const Ray3<T> &r, const Matrix3<T> &m)
    {
        return Ray3<T>(r.origin * m, r.direction * m);
    }

    template RUNEMATH_API Ray3<float> operator *(const Ray3<float> &r, const Matrix3<float> &m);
    template RUNEMATH_API Ray3<double> operator *(const Ray3<double> &r, const Matrix3<double> &m);

    /**
        Transformuje primku v prostoru
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API Ray3<T> operator *(const Matrix4<T> &m, const Ray3<T> &r)
    {
        return Ray3<T>(m * Vertex3<T>(r.origin), m * r.direction);
    }

    template RUNEMATH_API Ray3<float> operator *(const Matrix4<float> &m, const Ray3<float> &r);
    template RUNEMATH_API Ray3<double> operator *(const Matrix4<double> &m, const Ray3<double> &r);

    /**
        Transformuje primku v prostoru transponovanou matici
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API Ray3<T> operator *(const Ray3<T> &r, const Matrix4<T> &m)
    {
        return Ray3<T>(Vertex3<T>(r.origin) * m, r.direction * m);
    }

    template RUNEMATH_API Ray3<float> operator *(const Ray3<float> &r, const Matrix4<float> &m);
    template RUNEMATH_API Ray3<double> operator *(const Ray3<double> &r, const Matrix4<double> &m);
}
