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
    #define snscanf _snscanf
#endif

#include "../autodock-4.0.1/autocomm.h"
#include <cmath>
#include <cfloat>

/******************************************************************************/
/*      Name: autogrid.h                                                      */
/*  Function: Header file for Autogrid.                                       */
/* Copyright: (C) 1995, TSRI                                                  */
/*----------------------------------------------------------------------------*/
/*   Authors: Garrett Matthew Morris, David S. Goodsell                       */
/*                                                                            */
/*            The Scripps Research Institute                                  */
/*            Department of Molecular Biology, MB5                            */
/*            10666 North Torrey Pines Road                                   */
/*            La Jolla, CA 92037.                                             */
/*                                                                            */
/*            e-mail: garrett@scripps.edu                                     */
/*                    goodsell@scripps.edu                                    */
/*                                                                            */
/*      Date: 02/06/95  6-FEB-1995                                            */
/*----------------------------------------------------------------------------*/
/*    Inputs: None.                                                           */
/*   Returns: Parameters, Macro substitutions, Prototyped functions.          */
/*   Globals: (see 'autoglobal.h')                                            */
/*----------------------------------------------------------------------------*/
/* Modification Record                                                        */
/* Date     Inits   Comments                                                  */
/* 04/01/93 GMM     Created for use in makefile.                              */
/******************************************************************************/

#define MAX_DIST     16384   /* Maximum distance in 100ths of an Angstrom.    */
                             /*  = 163.84 Angstroms                           */
#define AG_MAX_ATOMS    32768   /* Maximum number of atoms in macromolecule.     */
/*    32768 = 2^15    */
/*    int 16-bit two's complement ranges 0-32767, 0 to (2^15 - 1)    */

#define ORDERED     0
#define CYLINDRICAL 1
#define SPHERICAL   2

#define A_DIVISOR   100    /* Angstrom is divided by this in look-up table. */

#define NBCUTOFF    8      /* non-bond cutoff = 8 Angstroms.                */

#define PRECISION   0.0001 /* fabs(Energies) less than this will be written as '0' */

/*----------------------------------------------------------------------------*/
/* Macros,                                                                    */
/*----------------------------------------------------------------------------*/

#define sq(a)               ((a) * (a))
#define sq_hyp(x,y,z)       ((x)*(x) + (y)*(y) + (z)*(z))
#define hypotenuse(x,y,z)   (sqrt(double(sq(x) + sq(y) + sq(z))))
#define equal(a,b,n)        (strncmp(a, b, size_t(n)) == 0)

// round() is a C99 function and not universally available
// Required to round %.3f consistently on different platforms
#if defined(HAVE_ROUND)
    #define round3dp(x) ((round((x)*1000.0L))/1000.0L)
#else
    #define round3dp(x) ((floor((x)*1000.0 + 0.5)) / 1000.0)
#endif

// we do not want to have a redefinition of the following macro max,min
#ifdef _WIN32
#undef min
#undef max
#endif

#define max(x,y)            (((x) > (y)) ? (x) : (y))
#define min(x,y)            (((x) < (y)) ? (x) : (y))
#define angstrom(i)         ((double(i)) / A_DIVISOR)
#define lookup(r)           (int((r) * A_DIVISOR))

#define NUM_ALL_TYPES 20   /*??? IS THIS REASONABLE???*/
#define MAX_LEN_AUTOGRID_TYPE 7

// from main()
#define NUM_RECEPTOR_TYPES  NUM_ALL_TYPES
#define INIT_NUM_GRID_PTS -1
