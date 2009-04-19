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

#if !defined(__CUDACC__)

    // Required for a successful compilation on Visual C++
    #if defined(_MSC_VER)
        // disable the warning: ' function ': was declared deprecated
        #pragma warning (disable: 4996)
        // disable the warning: conditional expression is constant
        #pragma warning (disable: 4127)

        // Some functions in Visual C++ differ from those in the linux/unix environment
        #define isnan _isnan
        #define strncasecmp _strnicmp
        #define snprintf _snprintf

        #define inline __forceinline
    #endif
#endif

#include "../autodock-4.0.1/autocomm.h"

#if !defined(__CUDACC__)
    #undef X
    #undef Y
    #undef Z
    #undef XYZ

    #include <cmath>
    #include <cfloat>
    #include <omp.h>
    #include "math.h"

    // Options
    // Do not uncomment these! Specify them in the command-line arguments of your compiler (G++, CMake) or in the project settings (VC++)

    // Enables the OpenMP support
    //#define AG_OPENMP

    // Enables the NVIDIA CUDA support
    //#define AG_CUDA

    // OpenMP configuration

    #define AG_OPENMP_PARALLEL_FOR omp parallel for schedule(dynamic, 1)

#endif

// Macros

// Using a greater value of MAX_DIST might increase precision a little,
// but keep in mind that a memory consumption increases linearly
#define MAX_DIST                (1<<13) // 2^13 = 8192 = 81.92 Angstroms. Maximum distance in 100ths of an Angstrom.
#define NBCUTOFF                8       // non-bond cutoff = 8 Angstroms.
#define AG_MAX_ATOMS            (1<<15) // 2^15 = 32768. Maximum number of atoms in macromolecule. // TODO: fix the code so that the number of atoms can be arbitrary
#define A_DIVISOR               100     // Angstrom is divided by this in the look-up table.
#define MAX_LEN_AUTOGRID_TYPE   7
#define NUM_ALL_TYPES           32      // TODO: IS THIS REASONABLE???
#define NUM_RECEPTOR_TYPES      NUM_ALL_TYPES
#define INIT_NUM_GRID_PTS       UINT_MAX

#if !defined(AG_CALLCONV)
    #define AG_CALLCONV
#endif

// Functions

#define DDD_FACTOR 332 // Used to convert between calories and SI units

// Distance-dependent dielectric ewds: Mehler and Solmajer, Prot Eng 4, 903-910.
template<typename T>
inline AG_CALLCONV T calculateDistDepDielInv(T distance)
{
    // Optimized formula is given by:
    T E = exp(T(-0.3153767175) * distance);
    return (T(DDD_FACTOR) + T(DDD_FACTOR * 7.7839) * E) / (T(78.4) + T(-66.57180475) * E);
}

template<typename Float, typename Int>
inline AG_CALLCONV Float indexToAngstrom(Int i)
{
    return Float(i) / Float(A_DIVISOR);
}

template<typename Int, typename Float>
inline AG_CALLCONV Int angstromToIndex(Float r)
{
    Int index = Int(r * Float(A_DIVISOR));
    return min(index, Int(MAX_DIST-1)); // make sure lookup index is in the table
}

#if !defined(__CUDACC__)

    inline double roundOutput(double a)
    {
        return fabs(a) < 0.0005 ? 0 : Mathd::Round(a * 1000) * 0.001; // round to 3 decimal places
    }

    inline int align(int value, int size)
    {
        return ((value - 1) / size + 1) * size;
    }

    // Useful macros - loop over all grid points

    #define FOR_EACH_GRID_POINT(gridPos, outputIndex) \
        /* Z axis */ \
        for (int z = 0; z < input->numGridPoints.z; z++) \
        { \
            /* gridPos contains the current grid point. */ \
            Vec3d gridPos; \
            gridPos.z = (z - input->numGridPointsDiv2.z) * input->gridSpacing; \
            int outputIndexZBase = z * input->numGridPoints.xy.Square(); \
    \
            /* Y axis */ \
            for (int y = 0; y < input->numGridPoints.y; y++) \
            { \
                gridPos.y = (y - input->numGridPointsDiv2.y) * input->gridSpacing; \
                int outputIndexZYBase = outputIndexZBase + y * input->numGridPoints.x; \
    \
                /* X axis */ \
                for (int x = 0; x < input->numGridPoints.x; x++) \
                { \
                    gridPos.x = (x - input->numGridPointsDiv2.x) * input->gridSpacing; \
                    int outputIndex = outputIndexZYBase + x;

    #define END_FOR() } } }

#endif
