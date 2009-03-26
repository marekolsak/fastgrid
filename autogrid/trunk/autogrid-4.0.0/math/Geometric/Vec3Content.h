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
    template<typename T, int spaceXY, int spaceYZ>
    class Vec3Content
    {
    public:
        union
        {
            struct
            {
                T x;
                T __unusedXY1[spaceXY];
                T y;
                T __unusedYZ1[spaceYZ];
                T z;
            };
            struct
            {
                Vec2<T, spaceXY> xy;
            };
            struct
            {
                Vec2<T, spaceXY, VEC2CONTENT_ORDER_YX> yx;
            };
            struct
            {
                T __unusedXY2[1+spaceXY];
                Vec2<T, spaceYZ> yz;
            };
            struct
            {
                T __unusedXY2[1+spaceXY];
                Vec2<T, spaceYZ, VEC2CONTENT_ORDER_YX> zy;
            };
            struct
            {
                Vec2<T, spaceXY+1+spaceYZ> xz;
            };
            struct
            {
                Vec2<T, spaceXY+1+spaceYZ, VEC2CONTENT_ORDER_YX> zx;
            };
        };
    };

    template<typename T, int spaceXY>
    class Vec3Content<T, spaceXY, 0>
    {
    public:
        union
        {
            struct
            {
                T x;
                T __unusedXY1[spaceXY];
                T y;
                T z;
            };
            struct
            {
                Vec2<T, spaceXY> xy;
            };
            struct
            {
                Vec2<T, spaceXY, VEC2CONTENT_ORDER_YX> yx;
            };
            struct
            {
                T __unusedXY2[1+spaceXY];
                Vec2<T> yz;
            };
            struct
            {
                T __unusedXY2[1+spaceXY];
                Vec2<T, 0, VEC2CONTENT_ORDER_YX> zy;
            };
            struct
            {
                Vec2<T, spaceXY+1> xz;
            };
            struct
            {
                Vec2<T, spaceXY+1, VEC2CONTENT_ORDER_YX> zx;
            };
        };
    };

    template<typename T, int spaceYZ>
    class Vec3Content<T, 0, spaceYZ>
    {
    public:
        union
        {
            struct
            {
                T x;
                T y;
                T __unusedYZ1[spaceYZ];
                T z;
            };
            struct
            {
                Vec2<T> xy;
            };
            struct
            {
                Vec2<T, 0, VEC2CONTENT_ORDER_YX> yx;
            };
            struct
            {
                T __unusedX;
                Vec2<T, spaceYZ> yz;
            };
            struct
            {
                T __unusedX;
                Vec2<T, spaceYZ, VEC2CONTENT_ORDER_YX> zy;
            };
            struct
            {
                Vec2<T, 1+spaceYZ> xz;
            };
            struct
            {
                Vec2<T, 1+spaceYZ, VEC2CONTENT_ORDER_YX> zx;
            };
        };
    };

    template<typename T>
    class Vec3Content<T, 0, 0>
    {
    public:
        union
        {
            struct
            {
                T x;
                T y;
                T z;
            };
            struct
            {
                Vec2<T> xy;
            };
            struct
            {
                Vec2<T, 0, VEC2CONTENT_ORDER_YX> yx;
            };
            struct
            {
                T __unusedX;
                Vec2<T> yz;
            };
            struct
            {
                T __unusedX;
                Vec2<T, 0, VEC2CONTENT_ORDER_YX> zy;
            };
            struct
            {
                Vec2<T, 1> xz;
            };
            struct
            {
                Vec2<T, 1, VEC2CONTENT_ORDER_YX> zx;
            };
        };
    };
}
