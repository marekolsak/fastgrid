/*
    AutoGrid

    Copyright (C) 1989-2007, Garrett M. Morris, David S. Goodsell, Ruth Huey, Arthur J. Olson,
    All Rights Reserved.
    Copyright (C) 2008-2009, Marek Olsak (maraeo@gmail.com), All Rights Reserved.

    AutoGrid is a Trade Mark of The Scripps Research Institute.

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
// we do not want to have a redefinition of the following macros max,min
#ifdef _WIN32
    #undef min
    #undef max
#endif

#undef X
#undef Y
#undef Z
#include "math/All.h"
using namespace Rune;
#define X 0
#define Y 1
#define Z 2

// Length squared
template<typename T>
inline T lengthSquared(T x, T y, T z)
{
    return x*x + y*y + z*z;
}

// Round to 3 decimal places
template<typename T>
inline T round3dp(T x)
{
    return Math<T>::Round(x * T(1000)) / T(1000);
}

// Dot product
template<typename T>
inline T dotProduct(const T *a, const T *b)
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

// Cross product
template<typename T>
void crossProduct(T *r, const T *a, const T *b)
{
    r[0] = a[1] * b[2] - a[2] * b[1];
    r[1] = a[2] * b[0] - a[0] * b[2];
    r[2] = a[0] * b[1] - a[1] * b[0];
}

// Get angle between 2 vectors
template<typename T>
inline T angle(const T *a, const T *b)
{
    return Math<T>::Acos(dotProduct(a, b));
}
