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
        Ctyrrozmerny vektor
    **************************************************************************************************/
    template<typename T>
    class Vec4
    {
    public:
        typedef Rune::Vec2<T> Vec2;
        typedef Rune::Vec3<T> Vec3;

        T x, y, z, w;

        operator Vec3() const                           { return Vec3(x, y, z); } // This implicit conversion should not be allowed!

        Vec4()                                          {}
        Vec4(T f)                                       { x = f; y = f; z = f; w = f; }
        Vec4(T X, T Y, T Z, T W)                        { x = X; y = Y; z = Z; w = W; }
        Vec4(T X, T Y, const Vec2 &v)                   { x = X; y = Y; z = v.x; w = v.y; }
        Vec4(T X, const Vec2 &v, T W)                   { x = X; y = v.x; z = v.y; w = W; }
        Vec4(const Vec2 &v, T Z, T W)                   { x = v.x; y = v.y; z = Z; w = W; }
        Vec4(const Vec2 &v1, const Vec2 &v2)            { x = v1.x; y = v1.y; z = v2.x; w = v2.y; }
        Vec4(const Vec3 &v, T W)                        { x = v.x; y = v.y; z = v.z; w = W; }
        Vec4(T X, const Vec3 &v)                        { x = X; y = v.x; z = v.y; w = v.z; }

        template<typename U>
        explicit Vec4(const Vec4<U> &v)                 { x = T(v.x); y = T(v.y); z = T(v.z); w = T(v.w); }

        bool operator == (const Vec4 &v) const          { return x == v.x && y == v.y && z == v.z && w == v.w; }
        bool operator != (const Vec4 &v) const          { return x != v.x && y != v.y && z != v.z && w != v.w; }
        bool operator <(const Vec4 &v) const            { return x < v.x && y < v.y && z < v.z && w < v.w; }
        bool operator >(const Vec4 &v) const            { return x > v.x && y > v.y && z > v.z && w > v.w; }
        bool operator <=(const Vec4 &v) const           { return x <= v.x && y <= v.y && z <= v.z && w <= v.w; }
        bool operator >=(const Vec4 &v) const           { return x >= v.x && y >= v.y && z >= v.z && w >= v.w; }

        bool operator ==(T f) const                     { return x == f && y == f && z == f && w == f; }
        bool operator !=(T f) const                     { return x != f && y != f && z != f && w != f; }
        bool operator <(T f) const                      { return x < f && y < f && z < f && w < f; }
        bool operator >(T f) const                      { return x > f && y > f && z > f && w > f; }
        bool operator <=(T f) const                     { return x <= f && y <= f && z <= f && w <= f; }
        bool operator >=(T f) const                     { return x >= f && y >= f && z >= f && w >= f; }

        void operator = (const Vec4 &v)                 { x = v.x; y = v.y; z = v.z; w = v.w; }
        void operator +=(const Vec4 &v)                 { x += v.x; y += v.y; z += v.z; w += v.w; }
        void operator +=(T f)                           { x += f; y += f; z += f; w += f; }
        void operator -=(const Vec4 &v)                 { x -= v.x; y -= v.y; z -= v.z; w -= v.w; }
        void operator -=(T f)                           { x -= f; y -= f; z -= f; w -= f; }
        void operator *=(const Vec4 &v)                 { x *= v.x; y *= v.y; z *= v.z; w *= v.w; }
        void operator *=(T f)                           { x *= f; y *= f; z *= f; w *= f; }
        void operator /=(T f)                           { x /= f; y /= f; z /= f; w /= f; }

        Vec4 operator +(const Vec4 &v) const            { return Vec4(x+v.x, y+v.y, z+v.z, w+v.w); }
        Vec4 operator +(T f) const                      { return Vec4(x+f, y+f, z+f, w+f); }
        Vec4 operator -(const Vec4 &v) const            { return Vec4(x-v.x, y-v.y, z-v.z, w-v.w); }
        Vec4 operator -(T f) const                      { return Vec4(x-f, y-f, z-f, w-f); }
        Vec4 operator *(const Vec4 &v) const            { return Vec4(x*v.x, y*v.y, z*v.z, w*v.w); }
        Vec4 operator *(T f) const                      { return Vec4(x*f, y*f, z*f, w*f); }
        Vec4 operator /(T f) const                      { return Vec4(x/f, y/f, z/f, w/f); }
        Vec4 operator -() const                         { return Vec4(-x, -y, -z, -w); }
        T& operator[] (int i)                           { return (&x)[i]; }
        T operator[] (int i) const                      { return (&x)[i]; }

        Vec4 GetAbs() const                             { return Vec4(Math<T>::Abs(x), Math<T>::Abs(y), Math<T>::Abs(z), Math<T>::Abs(w)); }
        Vec4 GetNormalized() const                      { T f = MagnitudeInv(); return Vec4(x*f, y*f, z*f, w*f); }
        T Magnitude() const                             { return Math<T>::Sqrt(MagnitudeSqr()); }
        T MagnitudeInv() const                          { return Math<T>::Rsqrt(MagnitudeSqr()); }
        T MagnitudeSqr() const                          { return x*x + y*y + z*z + w*w; }
        void Normalize()                                { operator *=(MagnitudeInv()); }
        bool IsEmpty() const                            { return x == 0 && y == 0 && z == 0 && w == 0; }
        void PackTo01()                                 { operator *=(0.5f); operator +=(0.5f); }
        void PackBack()                                 { operator -=(0.5f); operator *=(2); }
        void Set(const T X, const T Y, const T Z, const T W) { x = X; y = Y; z = Z; w = W; }
        void Inverse()                                  { x = 1/x; y = 1/y; z = 1/z; w = 1/w; }
        Vec4 Frac() const                               { return Vec4(Math<T>::Frac(x), Math<T>::Frac(y), Math<T>::Frac(z), Math<T>::Frac(w)); }

        RUNEMATH_API bool SafeIsEqual(const Vec4 &v) const;
        RUNEMATH_API static T Dot(const Vec4 &v1, const Vec4 &v2);
        RUNEMATH_API static Vec4 Cross(const Vec4 &a, const Vec4 &b, const Vec4 &c);
    };

    template<typename T>
    std::ostream& operator <<(std::ostream &stream, const Vec4<T> &v)
    {
        return stream << ("(" + tostrf(v.x) + ", " + tostrf(v.y) + ", " + tostrf(v.z) + ", " + tostrf(v.w) + ")");
    }

    template<typename T>
    std::istream& operator >>(std::istream &stream, Vec4<T> &v)
    {
        char tmp;
        return stream >> tmp >> v.x >> tmp >> v.y >> tmp >> v.z >> tmp >> v.w >> tmp;
    }

    typedef Vec4<float> Vec4f;
    typedef Vec4<double> Vec4d;
    typedef Vec4<int32> Vec4i;
    typedef Vec4<bool> Vec4b;
}
