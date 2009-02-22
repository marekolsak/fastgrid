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
#include "parameterlibrary.h"

struct InputData
{
    // if the first char is equal to '\0', the filename is not specified
    char fldFilenameAVS[MAX_CHARS];
    char floatingGridFilename[MAX_CHARS];
    char receptorFilename[MAX_CHARS];
    char xyzFilename[MAX_CHARS];
    char parameterLibraryFilename[MAX_CHARS];   // the AD4 parameters .dat file name

    // variables for RECEPTOR:
    // each type is now at most two characters, eg 'NA\0'
    // NB: these are sparse arrays, some entries are not set
    char receptorTypes[NUM_RECEPTOR_TYPES][3];
    int numReceptorTypes; // number of different receptor atom types actually found in receptor PDBQT

    int numGridPointsPerMap;
    int numReceptorAtoms;

    double charge[AG_MAX_ATOMS];
    double vol[AG_MAX_ATOMS];
    double solpar[AG_MAX_ATOMS];
    int atomType[AG_MAX_ATOMS];
    HBondType hbond[AG_MAX_ATOMS];
    double coord[AG_MAX_ATOMS][XYZ];

    double cgridmin[XYZ];
    double center[XYZ];
    double covalentPoint[XYZ];         // Cartesian-coordinate of covalent affinity well.
    int ne[XYZ];
    int numGridPoints[XYZ];
    int nelements[XYZ];

    double epsilon[MAX_DIST];

    // for NEW3 desolvation terms
    double solparQ;   // unweighted value restored 3:9:05
    double invDielCal;
    double rSmooth;
    double spacing;     // One quarter of a C-C bond length.
    double covHalfWidth;
    double covBarrier;

    bool distDepDiel, disorderH;
};
