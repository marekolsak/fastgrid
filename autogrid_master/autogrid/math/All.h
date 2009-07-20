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

// TODO: translate all czech comments into english

#pragma once
#define RUNEMATH_API

#if defined(_MSC_VER)
    // disable the warning: nameless struct/union
    #pragma warning (disable: 4201)
#endif

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include "Int.h"
#include "UnionCast.h"

// Basic
#include "Math.h"

// Vectors
#include "Geometric/Vec2.h"
#include "Geometric/Vec3.h"
#include "Geometric/Vec4.h"

// 3D
#include "Geometric/Sphere3.h"
#include "Geometric/AxisAlignedBox3.h"
#include "Geometric/Intersection3.h"
