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
        Trirozmerny vektor
    **************************************************************************************************/
    template<typename T>
    class Vec3
    {
    public:
        T x,y,z;

        typedef Rune::Vec2<T> Vec2;

        Vec3() {}
        Vec3(T f): x(f), y(f), z(f) {}
        Vec3(const Vec2 &v, T Z): x(v.x), y(v.y), z(Z) {}
        Vec3(T X, const Vec2 &v): x(X), y(v.x), z(v.y) {}
        Vec3(T X, T Y, T Z): x(X), y(Y), z(Z) {}

        template<typename U>
        explicit Vec3(const Vec3<U> &v): x(T(v.x)), y(T(v.y)), z(T(v.z)) {}

        bool operator ==(const Vec3 &v) const           { return x == v.x && y == v.y && z == v.z; }
        bool operator !=(const Vec3 &v) const           { return x != v.x && y != v.y && z != v.z; }
        bool operator <(const Vec3 &v) const            { return x < v.x && y < v.y && z < v.z; }
        bool operator >(const Vec3 &v) const            { return x > v.x && y > v.y && z > v.z; }
        bool operator <=(const Vec3 &v) const           { return x <= v.x && y <= v.y && z <= v.z; }
        bool operator >=(const Vec3 &v) const           { return x >= v.x && y >= v.y && z >= v.z; }

        bool operator ==(T f) const                     { return x == f && y == f && z == f; }
        bool operator !=(T f) const                     { return x != f && y != f && z != f; }
        bool operator <(T f) const                      { return x < f && y < f && z < f; }
        bool operator >(T f) const                      { return x > f && y > f && z > f; }
        bool operator <=(T f) const                     { return x <= f && y <= f && z <= f; }
        bool operator >=(T f) const                     { return x >= f && y >= f && z >= f; }

        void operator =(const Vec3 &v)                  { x = v.x; y = v.y; z = v.z; }
        void operator +=(const Vec3 &v)                 { x += v.x; y += v.y; z += v.z; }
        void operator +=(T f)                           { x += f; y += f; z += f; }
        void operator -=(const Vec3 &v)                 { x -= v.x; y -= v.y; z -= v.z; }
        void operator -=(T f)                           { x -= f; y -= f; z -= f; }
        void operator *=(T f)                           { x *= f; y *= f; z *= f; }
        void operator /=(T f)                           { f = 1/f; x *= f; y *= f; z *= f; }

        Vec3 operator +(const Vec3 &v) const            { return Vec3(x+v.x, y+v.y, z+v.z); }
        Vec3 operator +(T f) const                      { return Vec3(x+f, y+f, z+f); }
        Vec3 operator -(const Vec3 &v) const            { return Vec3(x-v.x, y-v.y, z-v.z); }
        Vec3 operator -(T f) const                      { return Vec3(x-f, y-f, z-f); }
        Vec3 operator *(const Vec3 &v) const            { return Vec3(x*v.x, y*v.y, z*v.z); }
        Vec3 operator *(T f) const                      { return Vec3(x*f, y*f, z*f); }
        Vec3 operator /(const Vec3 &v) const            { return Vec3(x/v.x, y/v.y, z/v.z); }
        Vec3 operator /(T f) const                      { f = 1/f; return Vec3(x*f, y*f, z*f); }
        Vec3 operator -() const                         { return Vec3(-x, -y, -z); }
        T& operator [](int i)                           { return (&x)[i]; }
        T operator [](int i) const                      { return (&x)[i]; }

        operator T* ()                                  { return &x; }
        operator const T* () const                      { return &x; }

        Vec2 GetVec2() const                            { return Vec2(x, y); }
        T GetMax() const                                { return Math<T>::Max(Math<T>::Max(x, y), z); }
        Vec3 GetAbs() const                             { return Vec3(Math<T>::Abs(x), Math<T>::Abs(y), Math<T>::Abs(z)); }
        Vec3 GetNormalized() const                      { T f = RMagnitude(); return Vec3(x*f, y*f, z*f); }
        T Magnitude() const                             { return Math<T>::Sqrt(MagnitudeSqr()); }
        T RMagnitude() const                            { return Math<T>::Rsqrt(MagnitudeSqr()); }
        T MagnitudeSqr() const                          { return x*x + y*y + z*z; }
        Vec3 Reflect(const Vec3 &normal) const          { return *this - 2*Dot(normal, *this)*normal; }
        Vec3 Project(const Vec3 &v) const               { return *this * Dot(v, *this)/MagnitudeSqr(); }
        Vec3 Orthogonalize(const Vec3 &v) const         { return v - Project(v); }
        void Normalize()                                { operator *=(RMagnitude()); }
        void PackTo01()                                 { operator *=(T(0.5)); operator +=(T(0.5)); }
        void PackBack()                                 { operator *=(2); operator -=(1); }
        void Set(T X, T Y, T Z)                         { x = X; y = Y; z = Z; }
        void Rotate(const Vec3 &v, const T angle)       { Rotate(v.x, v.y, v.z, angle); }
        void Inverse()                                  { x = 1/x; y = 1/y; z = 1/z; }
        Vec3 Frac() const                               { return Vec3(Math<T>::Frac(x), Math<T>::Frac(y), Math<T>::Frac(z)); }

        RUNEMATH_API bool SafeIsEqual(const Vec3 &v) const;
        RUNEMATH_API void Rotate(T x, T y, T z, T angle);
        RUNEMATH_API void RotateX(T angle);
        RUNEMATH_API void RotateY(T angle);
        RUNEMATH_API void RotateZ(T angle);
        RUNEMATH_API static T Dot(const Vec3 &v1, const Vec3 &v2);
        RUNEMATH_API static Vec3 Cross(const Vec3 &v1, const Vec3 &v2);
        RUNEMATH_API static T GetDistanceSqr(const Vec3 &v1, const Vec3 &v2);
        RUNEMATH_API static T GetDistance(const Vec3 &v1, const Vec3 &v2);
        RUNEMATH_API static Vec3 GetNormal(const Vec3 &v0, const Vec3 &v1, const Vec3 &v2);
        RUNEMATH_API static Vec3 GetNormalUN(const Vec3 &v0, const Vec3 &v1, const Vec3 &v2);
        RUNEMATH_API static Vec3 GetCenter(const Vec3 &min, const Vec3 &max);
        RUNEMATH_API static T GetAngle(const Vec3 &a, const Vec3 &b);
        RUNEMATH_API static Vec3 GetTangent(const Vec3 &v0, const Vec3 &v1, const Vec3 &v2, const Vec2 &t0, const Vec2 &t1, const Vec2 &t2);
        RUNEMATH_API static Vec3 GetBitangent(const Vec3 &v0, const Vec3 &v1, const Vec3 &v2, const Vec2 &t0, const Vec2 &t1, const Vec2 &t2);
        RUNEMATH_API static T GetArea(const Vec3 &a, const Vec3 &b, const Vec3 &c);
    };

    template<typename T>
    std::ostream& operator <<(std::ostream &stream, const Vec3<T> &v)
    {
        return stream << ("(" + tostrf(v.x) + ", " + tostrf(v.y) + ", " + tostrf(v.z) + ")");
    }

    template<typename T>
    std::istream& operator >>(std::istream &stream, Vec3<T> &v)
    {
        char tmp;
        return stream >> tmp >> v.x >> tmp >> v.y >> tmp >> v.z >> tmp;
    }

    typedef Vec3<float> Vec3f;
    typedef Vec3<double> Vec3d;
    typedef Vec3<int32> Vec3i;
    typedef Vec3<bool> Vec3b;
}
