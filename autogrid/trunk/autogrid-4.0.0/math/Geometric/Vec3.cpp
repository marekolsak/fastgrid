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
    RUNEMATH_API bool Vec3<T>::SafeIsEqual(const Vec3<T> &v) const
    {
        return Math<T>::SafeIsEqual(x, v.x) && Math<T>::SafeIsEqual(y, v.y) && Math<T>::SafeIsEqual(z, v.z);
    }

    template RUNEMATH_API bool Vec3<float>::SafeIsEqual(const Vec3<float> &v) const;
    template RUNEMATH_API bool Vec3<double>::SafeIsEqual(const Vec3<double> &v) const;

    /**
        Rotuje vektor kolem rotacni osy (x,y,z) uhlem angle
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API void Vec3<T>::Rotate(T x, T y, T z, T angle)
    {
        T s = Math<T>::SinDeg(angle), c = Math<T>::CosDeg(angle);
        T cr = 1-c, xx = x*x, xy = x*y, xz = x*z, yy = y*y, yz = y*z, zz = z*z, sx = s*x, sy = s*y, sz = s*z;
        Vec3<T> r1(xx+c*(1-xx), xy*cr-sz,    xz*cr+sy);
        Vec3<T> r2(xy*cr+sz,    yy+c*(1-yy), yz*cr-sx);
        Vec3<T> r3(xz*cr-sy,    yz*cr+sx,    zz+c*(1-zz));
        Set(Dot(*this, r1), Dot(*this, r2), Dot(*this, r3));
    }

    template RUNEMATH_API void Vec3<float>::Rotate(float x, float y, float z, float angle);
    template RUNEMATH_API void Vec3<double>::Rotate(double x, double y, double z, double angle);

    /**
        Rotuje vektor kolem osy X uhlem angle
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API void Vec3<T>::RotateX(T angle)
    {
        T s = Math<T>::SinDeg(angle), c = Math<T>::CosDeg(angle);
        Set(x, y*c-z*s, y*s+z*c);
    }

    template RUNEMATH_API void Vec3<float>::RotateX(float angle);
    template RUNEMATH_API void Vec3<double>::RotateX(double angle);

    /**
        Rotuje vektor kolem osy Y uhlem angle
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API void Vec3<T>::RotateY(T angle)
    {
        T s = Math<T>::SinDeg(angle), c = Math<T>::CosDeg(angle);
        Set(x*c+z*s, y, -x*s+z*c);
    }

    template RUNEMATH_API void Vec3<float>::RotateY(float angle);
    template RUNEMATH_API void Vec3<double>::RotateY(double angle);

    /**
        Rotuje vektor kolem osy Z uhlem angle
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API void Vec3<T>::RotateZ(T angle)
    {
        T s = Math<T>::SinDeg(angle), c = Math<T>::CosDeg(angle);
        Set(x*c-y*s, x*s+y*c, z);
    }

    template RUNEMATH_API void Vec3<float>::RotateZ(float angle);
    template RUNEMATH_API void Vec3<double>::RotateZ(double angle);

    template<typename T>
    RUNEMATH_API T Vec3<T>::Dot(const Vec3<T> &v1, const Vec3<T> &v2)
    {
        return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
    }

    template RUNEMATH_API float Vec3<float>::Dot(const Vec3<float> &v1, const Vec3<float> &v2);
    template RUNEMATH_API double Vec3<double>::Dot(const Vec3<double> &v1, const Vec3<double> &v2);

    template<typename T>
    RUNEMATH_API Vec3<T> Vec3<T>::Cross(const Vec3<T> &v1, const Vec3<T> &v2)
    {
        return Vec3<T>(v1.y*v2.z - v1.z*v2.y, v1.z*v2.x - v1.x*v2.z, v1.x*v2.y - v1.y*v2.x);
    }

    template RUNEMATH_API Vec3<float> Vec3<float>::Cross(const Vec3<float> &v1, const Vec3<float> &v2);
    template RUNEMATH_API Vec3<double> Vec3<double>::Cross(const Vec3<double> &v1, const Vec3<double> &v2);

    template<typename T>
    RUNEMATH_API T Vec3<T>::GetDistanceSqr(const Vec3<T> &v1, const Vec3<T> &v2)
    {
        return (v1 - v2).MagnitudeSqr();
    }

    template RUNEMATH_API float Vec3<float>::GetDistanceSqr(const Vec3<float> &v1, const Vec3<float> &v2);
    template RUNEMATH_API double Vec3<double>::GetDistanceSqr(const Vec3<double> &v1, const Vec3<double> &v2);

    template<typename T>
    RUNEMATH_API T Vec3<T>::GetDistance(const Vec3<T> &v1, const Vec3<T> &v2)
    {
        return Math<T>::Sqrt(GetDistanceSqr(v1, v2));
    }

    template RUNEMATH_API float Vec3<float>::GetDistance(const Vec3<float> &v1, const Vec3<float> &v2);
    template RUNEMATH_API double Vec3<double>::GetDistance(const Vec3<double> &v1, const Vec3<double> &v2);

    template<typename T>
    RUNEMATH_API Vec3<T> Vec3<T>::GetNormal(const Vec3<T> &v0, const Vec3<T> &v1, const Vec3<T> &v2)
    {
        return GetNormalUN(v0, v1, v2).GetNormalized();
    }

    template RUNEMATH_API Vec3<float> Vec3<float>::GetNormal(const Vec3<float> &v0, const Vec3<float> &v1, const Vec3<float> &v2);
    template RUNEMATH_API Vec3<double> Vec3<double>::GetNormal(const Vec3<double> &v0, const Vec3<double> &v1, const Vec3<double> &v2);

    template<typename T>
    RUNEMATH_API Vec3<T> Vec3<T>::GetNormalUN(const Vec3<T> &v0, const Vec3<T> &v1, const Vec3<T> &v2)
    {
        return Cross(v2-v1, v0-v2);
    }

    template RUNEMATH_API Vec3<float> Vec3<float>::GetNormalUN(const Vec3<float> &v0, const Vec3<float> &v1, const Vec3<float> &v2);
    template RUNEMATH_API Vec3<double> Vec3<double>::GetNormalUN(const Vec3<double> &v0, const Vec3<double> &v1, const Vec3<double> &v2);

    template<typename T>
    RUNEMATH_API Vec3<T> Vec3<T>::GetCenter(const Vec3<T> &min, const Vec3<T> &max)
    {
        return (min + max) * T(0.5);
    }

    template RUNEMATH_API Vec3<float> Vec3<float>::GetCenter(const Vec3<float> &min, const Vec3<float> &max);
    template RUNEMATH_API Vec3<double> Vec3<double>::GetCenter(const Vec3<double> &min, const Vec3<double> &max);

    template<typename T>
    RUNEMATH_API T Vec3<T>::GetAngle(const Vec3<T> &a, const Vec3<T> &b)
    {
        return Math<T>::Acos(Dot(a, b)*a.RMagnitude()*b.RMagnitude());
    }

    template RUNEMATH_API float Vec3<float>::GetAngle(const Vec3<float> &a, const Vec3<float> &b);
    template RUNEMATH_API double Vec3<double>::GetAngle(const Vec3<double> &a, const Vec3<double> &b);

    /**
        Vraci tangent vektor (vektor mapovani textury pro normalmapping), vstupem jsou vertexy
        trojuhelniku a jeho texturove koordinaty
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API Vec3<T> Vec3<T>::GetTangent(const Vec3<T> &v0, const Vec3<T> &v1, const Vec3<T> &v2, const Vec2 &t0, const Vec2 &t1, const Vec2 &t2)
    {
        return ((v1 - v0) * (t2.y - t0.y) + (v2 - v0) * (t0.y - t1.y)).GetNormalized();
    }

    template RUNEMATH_API Vec3<float> Vec3<float>::GetTangent(const Vec3<float> &v0, const Vec3<float> &v1, const Vec3<float> &v2, const Vec2 &t0, const Vec2 &t1, const Vec2 &t2);
    template RUNEMATH_API Vec3<double> Vec3<double>::GetTangent(const Vec3<double> &v0, const Vec3<double> &v1, const Vec3<double> &v2, const Vec2 &t0, const Vec2 &t1, const Vec2 &t2);

    /**
        Vraci bitangent vektor (vektor mapovani textury pro normalmapping), vstupem jsou vertexy
        trojuhelniku a jeho texturove koordinaty
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API Vec3<T> Vec3<T>::GetBitangent(const Vec3<T> &v0, const Vec3<T> &v1, const Vec3<T> &v2, const Vec2 &t0, const Vec2 &t1, const Vec2 &t2)
    {
        return ((v1 - v0) * (t2.x - t0.x) + (v2 - v0) * (t0.x - t1.x)).GetNormalized();
    }

    template RUNEMATH_API Vec3<float> Vec3<float>::GetBitangent(const Vec3<float> &v0, const Vec3<float> &v1, const Vec3<float> &v2, const Vec2 &t0, const Vec2 &t1, const Vec2 &t2);
    template RUNEMATH_API Vec3<double> Vec3<double>::GetBitangent(const Vec3<double> &v0, const Vec3<double> &v1, const Vec3<double> &v2, const Vec2 &t0, const Vec2 &t1, const Vec2 &t2);

    /**
        Vraci obsah trojuhelniku
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API T Vec3<T>::GetArea(const Vec3<T> &a, const Vec3<T> &b, const Vec3<T> &c)
    {
        return Cross(b - a, c - a).Magnitude() * T(0.5);
    }

    template RUNEMATH_API float Vec3<float>::GetArea(const Vec3<float> &a, const Vec3<float> &b, const Vec3<float> &c);
    template RUNEMATH_API double Vec3<double>::GetArea(const Vec3<double> &a, const Vec3<double> &b, const Vec3<double> &c);
}
