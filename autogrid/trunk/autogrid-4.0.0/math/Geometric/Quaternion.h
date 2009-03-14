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

#define RUNE_ADD_Q(func) Quaternion q; q.func; *this *= q;

namespace Rune
{
    /**
        Quaternion - uchovava rotaci bez zmeny meritka a umi mezi nimi interpolovat
    **************************************************************************************************/
    template<typename T>
    class Quaternion
    {
    public:
        Vec4<T> v;

        Quaternion() {}
        Quaternion(T x, T y, T z, T w): v(x, y, z, w) {}
        Quaternion(const Vec3<T> &v, T angle)           { SetRotated(v, angle); }
        Quaternion(const Vec3<T> &start, const Vec3<T> &end) { SetRotated(start, end); }
        Quaternion(T anglez, T anglex, T angley)        { SetYawPitchRoll(anglez, anglex, angley); }
        explicit Quaternion(const Vec4<T> &V): v(V) {}

        operator T* () const                            { return v; }
        operator const T* () const                      { return v; }
        void operator +=(T f)                           { v += f; }
        void operator +=(const Quaternion &q)           { v += q.v; }
        void operator *=(T f)                           { v *= f; }
        void operator *=(const Quaternion &q)           { operator =(*this * q); }
        void operator /=(T f)                           { v /= f; }
        Quaternion operator +(T f) const                { return Quaternion(v + f); }
        Quaternion operator +(const Quaternion &q) const{ return Quaternion(v + q.v); }
        Quaternion operator *(T f) const                { return Quaternion(v * f); }
        RUNEMATH_API Quaternion operator *(const Quaternion &q) const;
        Quaternion operator /(T f) const                { return Quaternion(v / f); }

        void Normalize()                                { v.Normalize(); }
        void Conjugate()                                { v.Set(-v.x, -v.y, -v.z, v.w); }
        void Inverse()                                  { Conjugate(); operator *=(1 / v.MagnitudeSqr()); }
        Quaternion GetInverse() const                   { Quaternion q = *this; q.Inverse(); return q; }
        void LoadIdentity()                             { v.Set(0, 0, 0, 1); }
        void SetRotated(const Vec3<T> &v, T angle)      { SetRotated(v.x, v.y, v.z, angle); }
        void Rotate(const Vec3<T> &v, T angle)          { Rotate(v.x, v.y, v.z, angle); }
        void Rotate(T x, T y, T z, T angle)             { RUNE_ADD_Q(SetRotated(x, y, z, angle)); }
        void RotateX(T angle)                           { RUNE_ADD_Q(SetRotatedX(angle)); }
        void RotateY(T angle)                           { RUNE_ADD_Q(SetRotatedY(angle)); }
        void RotateZ(T angle)                           { RUNE_ADD_Q(SetRotatedZ(angle)); }

        RUNEMATH_API Vec4<T> GetRotation() const;
        RUNEMATH_API void SetRotated(T x, T y, T z, T angle);
        RUNEMATH_API void SetRotatedX(T angle);
        RUNEMATH_API void SetRotatedY(T angle);
        RUNEMATH_API void SetRotatedZ(T angle);
        RUNEMATH_API void SetYawPitchRoll(T anglez, T anglex, T angley);
        RUNEMATH_API void SetRotated(const Vec3<T> &start, const Vec3<T> &end);
        RUNEMATH_API void SetSlerp(const Quaternion &q1, const Quaternion &q2, T t);
    };

    template<typename T>
    RUNEMATH_API Vec3<T> operator *(const Quaternion<T> &q, const Vec3<T> &v);

    template<typename T>
    RUNEMATH_API Vec3<T> operator *(const Vec3<T> &v, const Quaternion<T> &q);

    template<typename T>
    std::ostream& operator <<(std::ostream &stream, const Quaternion<T> &q)
    {
        return stream << q.v;
    }

    template<typename T>
    std::istream& operator >>(std::istream &stream, Quaternion<T> &q)
    {
        return stream >> q.v;
    }

    typedef Quaternion<float> Quaternionf;
    typedef Quaternion<double> Quaterniond;
    typedef Quaternion<int32> Quaternioni;
}
