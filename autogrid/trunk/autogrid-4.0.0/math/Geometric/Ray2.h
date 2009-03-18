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
        Trida primky ve 2D
    **************************************************************************************************/
    template<typename T>
    class Ray2
    {
    public:
        Vec2<T> normal;
        T c;

        Ray2() {}
        Ray2(const Vec2<T> &Normal, const Vec2<T> &point): normal(Normal), c(-Vec2<T>::Dot(Normal, point)) {}
        Ray2(const Vec2<T> &Normal, T C): normal(Normal), c(C) {}

        void operator *=(T f)                                   { normal *= f; c *= f; }
        Ray2 operator *(T f) const                              { return Ray2(normal * f, c * f); }

        void Set(const Vec2<T> &Normal, const Vec2<T> &point)   { normal = Normal; c = -Vec2<T>::Dot(Normal, point); }
        void Set(const Vec2<T> &Normal, T C)                    { normal = Normal; c = C; }
        void Normalize()                                        { *this *= normal.MagnitudeInv(); }
        Ray2 GetNormalized() const                              { return *this * normal.MagnitudeInv(); }
        T GetDistance(const Vec2<T> &point) const               { return Vec2<T>::Dot(normal, point) + c; }
        Vec2<T> GetMirroredPoint(const Vec2<T> &point) const    { return point + normal*GetDistance(point)*-2; }
        Vec2<T> GetMirroredVector(const Vec2<T> &vec) const     { return vec + normal*Vec2<T>::Dot(normal, vec)*-2; }
        Vec2<T> GetPoint() const                                { return normal * -c; }
        Vec2<T> GetNearestPoint(const Vec2<T> &point) const     { return point + normal*GetDistance(point)*-1; }

        RUNEMATH_API int GetSide(const Vec2<T> &point) const;
        RUNEMATH_API int GetSide(const Circle2<T> &circle) const;
    };

    typedef Ray2<float> Ray2f;
    typedef Ray2<double> Ray2d;
}
