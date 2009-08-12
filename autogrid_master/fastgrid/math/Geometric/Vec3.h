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
        Trirozmerny vektor
    **************************************************************************************************/
    template<typename T>
    class Vec3
    {
    public:
        T x, y, z;

        template<typename U>
        explicit Vec3(const Vec3<U> &v) { x = T(v.x); y = T(v.y); z = T(v.z); }

        typedef Rune::Vec2<T> Vec2;

        Vec3() {}
        Vec3(T f)                                       { x = f; y = f; z = f; }
        Vec3(const Vec2 &v, T Z)                        { x = v.x; y = v.y; z = Z; }
        Vec3(T X, const Vec2 &v)                        { x = X; y = v.x; z = v.y; }
        Vec3(T X, T Y, T Z)                             { x = X; y = Y; z = Z; }

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
        void operator /=(T f)                           { x /= f; y /= f; z /= f; }

        Vec3 operator +(const Vec3 &v) const            { return Vec3(x+v.x, y+v.y, z+v.z); }
        Vec3 operator +(T f) const                      { return Vec3(x+f, y+f, z+f); }
        Vec3 operator -(const Vec3 &v) const            { return Vec3(x-v.x, y-v.y, z-v.z); }
        Vec3 operator -(T f) const                      { return Vec3(x-f, y-f, z-f); }
        Vec3 operator *(const Vec3 &v) const            { return Vec3(x*v.x, y*v.y, z*v.z); }
        Vec3 operator *(T f) const                      { return Vec3(x*f, y*f, z*f); }
        Vec3 operator /(const Vec3 &v) const            { return Vec3(x/v.x, y/v.y, z/v.z); }
        Vec3 operator /(T f) const                      { return Vec3(x/f, y/f, z/f); }
        Vec3 operator -() const                         { return Vec3(-x, -y, -z); }
        T& operator [](int i)                           { return (&x)[i]; }
        T operator [](int i) const                      { return (&x)[i]; }

        T Cube() const                                  { return x*y*z; }
        T GetMax() const                                { return Math<T>::Max(Math<T>::Max(x, y), z); }
        Vec3 GetAbs() const                             { return Vec3(Math<T>::Abs(x), Math<T>::Abs(y), Math<T>::Abs(z)); }
        Vec3 GetNormalized() const                      { T f = MagnitudeInv(); return Vec3(x*f, y*f, z*f); }
        T Magnitude() const                             { return Math<T>::Sqrt(MagnitudeSqr()); }
        T MagnitudeInv() const                          { return Math<T>::Rsqrt(MagnitudeSqr()); }
        T MagnitudeSqr() const                          { return x*x + y*y + z*z; }
        Vec3 Reflect(const Vec3 &normal) const          { return *this - 2*Dot(normal, *this)*normal; }
        Vec3 Project(const Vec3 &v) const               { return *this * Dot(v, *this)/MagnitudeSqr(); }
        Vec3 Orthogonalize(const Vec3 &v) const         { return v - Project(v); }
        void Normalize()                                { operator *=(MagnitudeInv()); }
        void PackTo01()                                 { operator *=(T(0.5)); operator +=(T(0.5)); }
        void PackBack()                                 { operator *=(2); operator -=(1); }
        void Set(T X, T Y, T Z)                         { x = X; y = Y; z = Z; }
        void Rotate(const Vec3 &v, const T angle)       { Rotate(v.x, v.y, v.z, angle); }
        void Inverse()                                  { x = 1/x; y = 1/y; z = 1/z; }
        Vec3 Frac() const                               { return Vec3(Math<T>::Frac(x), Math<T>::Frac(y), Math<T>::Frac(z)); }

        bool SafeIsEqual(const Vec3 &v) const
        {
            return Math<T>::SafeIsEqual(x, v.x) && Math<T>::SafeIsEqual(y, v.y) && Math<T>::SafeIsEqual(z, v.z);
        }

        /**
            Rotuje vektor kolem rotacni osy (x,y,z) uhlem angle
        **************************************************************************************************/
        void Rotate(T x, T y, T z, T angle)
        {
            T s = Math<T>::SinDeg(angle), c = Math<T>::CosDeg(angle);
            T cr = 1-c, xx = x*x, xy = x*y, xz = x*z, yy = y*y, yz = y*z, zz = z*z, sx = s*x, sy = s*y, sz = s*z;
            Vec3<T> r1(xx+c*(1-xx), xy*cr-sz,    xz*cr+sy);
            Vec3<T> r2(xy*cr+sz,    yy+c*(1-yy), yz*cr-sx);
            Vec3<T> r3(xz*cr-sy,    yz*cr+sx,    zz+c*(1-zz));
            Set(Dot(*this, r1), Dot(*this, r2), Dot(*this, r3));
        }

        /**
            Rotuje vektor kolem osy X uhlem angle
        **************************************************************************************************/
        void RotateX(T angle)
        {
            T s = Math<T>::SinDeg(angle), c = Math<T>::CosDeg(angle);
            Set(x, y*c-z*s, y*s+z*c);
        }

        /**
            Rotuje vektor kolem osy Y uhlem angle
        **************************************************************************************************/
        void RotateY(T angle)
        {
            T s = Math<T>::SinDeg(angle), c = Math<T>::CosDeg(angle);
            Set(x*c+z*s, y, -x*s+z*c);
        }

        /**
            Rotuje vektor kolem osy Z uhlem angle
        **************************************************************************************************/
        void RotateZ(T angle)
        {
            T s = Math<T>::SinDeg(angle), c = Math<T>::CosDeg(angle);
            Set(x*c-y*s, x*s+y*c, z);
        }

        static T Dot(const Vec3 &v1, const Vec3 &v2)
        {
            return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
        }

        static Vec3 Cross(const Vec3 &v1, const Vec3 &v2)
        {
            return Vec3(v1.y*v2.z - v1.z*v2.y, v1.z*v2.x - v1.x*v2.z, v1.x*v2.y - v1.y*v2.x);
        }

        template<typename Operator>
        static Vec3 ScalarOperator(const Vec3 &v1, const Vec3 &v2, Operator &op)
        {
            return Vec3(op(v1.x, v2.x), op(v1.y, v2.y), op(v1.z, v2.z));
        }

        static T DistanceSqr(const Vec3 &v1, const Vec3 &v2)
        {
            return (v1 - v2).MagnitudeSqr();
        }

        static T Distance(const Vec3 &v1, const Vec3 &v2)
        {
            return Math<T>::Sqrt(DistanceSqr(v1, v2));
        }

        static Vec3 CalculateNormal(const Vec3 &v0, const Vec3 &v1, const Vec3 &v2)
        {
            return CalculateNormalUnnorm(v0, v1, v2).GetNormalized();
        }

        static Vec3 CalculateNormalUnnorm(const Vec3 &v0, const Vec3 &v1, const Vec3 &v2)
        {
            return Cross(v2-v1, v0-v2);
        }

        static Vec3 Center(const Vec3 &min, const Vec3 &max)
        {
            return (min + max) * T(0.5);
        }

        static T Angle(const Vec3 &a, const Vec3 &b)
        {
            return Math<T>::Acos(Dot(a, b)*a.MagnitudeInv()*b.MagnitudeInv());
        }

        static T AngleUnnorm(const Vec3 &a, const Vec3 &b)
        {
            return Math<T>::Acos(Dot(a, b));
        }

        /**
            Vraci tangent vektor (vektor mapovani textury pro normalmapping), vstupem jsou vertexy
            trojuhelniku a jeho texturove koordinaty
        **************************************************************************************************/
        static Vec3 CalculateTangent(const Vec3 &v0, const Vec3 &v1, const Vec3 &v2, const Vec2 &t0, const Vec2 &t1, const Vec2 &t2)
        {
            return ((v1 - v0) * (t2.y - t0.y) + (v2 - v0) * (t0.y - t1.y)).GetNormalized();
        }

        /**
            Vraci bitangent vektor (vektor mapovani textury pro normalmapping), vstupem jsou vertexy
            trojuhelniku a jeho texturove koordinaty
        **************************************************************************************************/
        static Vec3 CalculateBitangent(const Vec3 &v0, const Vec3 &v1, const Vec3 &v2, const Vec2 &t0, const Vec2 &t1, const Vec2 &t2)
        {
            return ((v1 - v0) * (t2.x - t0.x) + (v2 - v0) * (t0.x - t1.x)).GetNormalized();
        }

        /**
            Vraci obsah trojuhelniku
        **************************************************************************************************/
        static T CalculateAreaOfTriangle(const Vec3 &a, const Vec3 &b, const Vec3 &c)
        {
            return Cross(b - a, c - a).Magnitude() * T(0.5);
        }
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
    typedef Vec3<uint32> Vec3ui;
    typedef Vec3<bool> Vec3b;
}
