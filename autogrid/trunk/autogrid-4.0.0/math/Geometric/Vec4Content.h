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
    template<typename T>
    class Vec4Content
    {
        union
        {
            struct
            {
                T x,y,z,w;
            };
            struct
            {
                Vec2<T> xy, zw;
            };
            struct
            {
                Vec2<T, 0, VEC2CONTENT_ORDER_YX> yx, wz;
            };
            struct
            {
                Vec2<T, 1> xz;
            };
            struct
            {
                Vec2<T, 1, VEC2CONTENT_ORDER_YX> zx;
            };
            struct
            {
                Vec2<T, 2> xw;
            };
            struct
            {
                Vec2<T, 2, VEC2CONTENT_ORDER_YX> wx;
            };
            struct
            {
                T __unusedX1;
                Vec2<T, 0> yz;
            };
            struct
            {
                T __unusedX2;
                Vec2<T, 0, VEC2CONTENT_ORDER_YX> zy;
            };
            struct
            {
                T __unusedX1;
                Vec2<T, 1> yw;
            };
            struct
            {
                T __unusedX2;
                Vec2<T, 1, VEC2CONTENT_ORDER_YX> wy;
            };
            struct
            {
                Vec3<T> xyz;
            };
            struct
            {
                Vec3<T, 1> xzw;
            };
            struct
            {
                Vec3<T, 0, 1> xyw;
            };
            struct
            {
                T __unusedX3;
                Vec3<T> yzw;
            };
        };
    };
}
