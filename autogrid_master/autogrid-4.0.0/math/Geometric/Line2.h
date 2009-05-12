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
        Trida usecky ve 2D
    **************************************************************************************************/
    template<typename T>
    class Line2
    {
    public:
        Vec2<T> origin;     // Stred usecky
        Ray2<T> ray;        // Primka prochazejici obema body usecky
        T extent;           // Vzdalenost stredu a konce usecky

        Line2() {}
        Line2(const Vec2<T> &p1, const Vec2<T> &p2) { SetPoints(p1, p2); }
        Vec2<T> GetPoint1() const { return origin + ray.normal.GetNormal() * extent; }
        Vec2<T> GetPoint2() const { return origin - ray.normal.GetNormal() * extent; }

        RUNEMATH_API void SetPoints(const Vec2<T> &p1, const Vec2<T> &p2);
        RUNEMATH_API int GetSide(const Vec2<T> &_point) const;
    };

    typedef Line2<float> Line2f;
    typedef Line2<double> Line2d;
}
