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
        Pro bod - vraci -1, kdyz je za plochou, 1 kdyz pred plochou a 0, kdyz lezi na plose
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API int Plane3<T>::GetSide(const Vec3<T> &point) const
    {
        T f = GetDistance(point);
        if (f > Math<T>::Epsilon()) return 1;
        if (f < -Math<T>::Epsilon()) return -1;
        return 0;
    }

    template RUNEMATH_API int Plane3<float>::GetSide(const Vec3<float> &point) const;
    template RUNEMATH_API int Plane3<double>::GetSide(const Vec3<double> &point) const;

    /**
        Pro kouli - vraci -1, kdyz je za plochou, 1 kdyz pred plochou a 0, kdyz se dotyka plochy
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API int Plane3<T>::GetSide(const Sphere3<T> &sphere) const
    {
        T f = GetDistance(sphere.pos);
        if (Math<T>::Abs(f) < sphere.radius) return 0;

        if (f > 0) f -= sphere.radius;
        else f += sphere.radius;

        if (f > Math<T>::Epsilon()) return 1;
        if (f < -Math<T>::Epsilon()) return -1;
        return 0;
    }

    template RUNEMATH_API int Plane3<float>::GetSide(const Sphere3<float> &sphere) const;
    template RUNEMATH_API int Plane3<double>::GetSide(const Sphere3<double> &sphere) const;

    /**
        Pro kvadr rovnobezny s osama - vraci -1, kdyz je za plochou, 1 kdyz pred plochou a 0, kdyz se dotyka plochy
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API int Plane3<T>::GetSide(const AxisAlignedBox3<T> &box) const
    {
        int sides[8];
        sides[0] = GetSide(box.min);
        sides[1] = GetSide(Vec3<T>(box.max.x, box.min.y, box.min.z));
        sides[2] = GetSide(Vec3<T>(box.min.x, box.max.y, box.min.z));
        sides[3] = GetSide(Vec3<T>(box.min.x, box.min.y, box.max.z));
        sides[4] = GetSide(Vec3<T>(box.max.x, box.max.y, box.min.z));
        sides[5] = GetSide(Vec3<T>(box.min.x, box.max.y, box.max.z));
        sides[6] = GetSide(Vec3<T>(box.max.x, box.min.y, box.max.z));
        sides[7] = GetSide(box.max);

        if (sides[0] == sides[1] == sides[2] == sides[3] == sides[4] == sides[5] == sides[6] == sides[7]) return *sides;
        return 0;
    }

    template RUNEMATH_API int Plane3<float>::GetSide(const AxisAlignedBox3<float> &box) const;
    template RUNEMATH_API int Plane3<double>::GetSide(const AxisAlignedBox3<double> &box) const;

    /**
        Transformuje rovinu v prostoru
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API Plane3<T> operator *(const Matrix3<T> &m, const Plane3<T> &p)
    {
        return Plane3<T>(m * p.normal, m * p.GetPoint()).GetNormalized();
    }

    template RUNEMATH_API Plane3<float> operator *(const Matrix3<float> &m, const Plane3<float> &p);
    template RUNEMATH_API Plane3<double> operator *(const Matrix3<double> &m, const Plane3<double> &p);

    /**
        Transformuje rovinu v prostoru transponovanou matici
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API Plane3<T> operator *(const Plane3<T> &p, const Matrix3<T> &m)
    {
        return Plane3<T>(p.normal * m, p.GetPoint() * m).GetNormalized();
    }

    template RUNEMATH_API Plane3<float> operator *(const Plane3<float> &p, const Matrix3<float> &m);
    template RUNEMATH_API Plane3<double> operator *(const Plane3<double> &p, const Matrix3<double> &m);

    /**
        Transformuje rovinu v prostoru
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API Plane3<T> operator *(const Matrix4<T> &m, const Plane3<T> &p)
    {
        return Plane3<T>(m * p.normal, m * Vertex3<T>(p.GetPoint())).GetNormalized();
    }

    template RUNEMATH_API Plane3<float> operator *(const Matrix4<float> &m, const Plane3<float> &p);
    template RUNEMATH_API Plane3<double> operator *(const Matrix4<double> &m, const Plane3<double> &p);

    /**
        Transformuje rovinu v prostoru transponovanou matici
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API Plane3<T> operator *(const Plane3<T> &p, const Matrix4<T> &m)
    {
        return Plane3<T>(p.normal * m, Vertex3<T>(p.GetPoint()) * m).GetNormalized();
    }

    template RUNEMATH_API Plane3<float> operator *(const Plane3<float> &p, const Matrix4<float> &m);
    template RUNEMATH_API Plane3<double> operator *(const Plane3<double> &p, const Matrix4<double> &m);
}
