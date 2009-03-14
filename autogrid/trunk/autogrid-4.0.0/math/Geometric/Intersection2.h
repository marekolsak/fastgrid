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
    template<typename T>
    RUNEMATH_API bool Intersect(const Ray2<T> &ray1, const Ray2<T> &ray2);

    template<typename T>
    RUNEMATH_API bool Intersect(const Ray2<T> &ray1, const Ray2<T> &ray2, Vec2<T> &point);

    template<typename T>
    RUNEMATH_API bool Intersect(const Line2<T> &line, const Ray2<T> &ray, Vec2<T> &point);

    template<typename T>
    RUNEMATH_API bool Intersect(const Line2<T> &line, const Ray2<T> &ray);

    template<typename T>
    RUNEMATH_API bool Intersect(const Line2<T> &line1, const Line2<T> &line2, Vec2<T> &point);

    template<typename T>
    RUNEMATH_API bool Intersect(const Line2<T> &line1, const Line2<T> &line2);

    template<typename T>
    RUNEMATH_API bool Intersect(const Circle2<T> &s1, const Circle2<T> &s2);

    template<typename T>
    RUNEMATH_API bool Intersect(const Circle2<T> &s1, const Circle2<T> &s2, Vec2<T> &point);

    template<typename T>
    RUNEMATH_API bool Intersect(const Circle2<T> &s, const Vec2<T> &point);

    template<typename T>
    RUNEMATH_API bool Intersect(const Rect2<T> &r, const Vec2<T> &point);
}
