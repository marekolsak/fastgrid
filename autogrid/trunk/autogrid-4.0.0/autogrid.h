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

// Required for a successful compilation on Visual C++
#if defined(_MSC_VER)
    // disable the warning: ' function ': was declared deprecated
    #pragma warning (disable: 4996)

    // Some functions in Visual C++ differ from those in the linux/unix environment
    #define isnan _isnan
    #define strncasecmp _strnicmp
    #define snprintf _snprintf

    #define inline __forceinline
#endif

#include "../autodock-4.0.1/autocomm.h"
#include <cmath>
#include <cfloat>

// OpenMP configuration

#define AG_OPENMP
#define AG_OPENMP_PARALLEL_FOR omp parallel for schedule(dynamic, 1)

// Macros

// Using a greater value of MAX_DIST might increase precision a little,
// but keep in mind the memory consumption increases linearly with respect to it
//#define MAX_DIST                (1<<14) // 2^14 = 16384 = 163.84 Angstroms. Maximum distance in 100ths of an Angstrom.
#define MAX_DIST                (1<<13) // 2^13 = 8192 = 81.92 Angstroms. Maximum distance in 100ths of an Angstrom.

#define NBCUTOFF                8       // non-bond cutoff = 8 Angstroms.
#define AG_MAX_ATOMS            (1<<15) // 2^15 = 32768. Maximum number of atoms in macromolecule.
#define A_DIVISOR               100     // Angstrom is divided by this in the look-up table.
#define MAX_LEN_AUTOGRID_TYPE   7
#define NUM_ALL_TYPES           32      // TODO: IS THIS REASONABLE???
#define NUM_RECEPTOR_TYPES      NUM_ALL_TYPES
#define INIT_NUM_GRID_PTS       -1

// Functions

template<typename T>
inline double angstrom(T i)
{
    return double(i) / A_DIVISOR;
}

template<typename T>
inline int lookup(T r)
{
    return int(r * A_DIVISOR);
}

#include "math.h"
