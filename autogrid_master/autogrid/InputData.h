/*
    FastGrid (formerly AutoGrid)

    Copyright (C) 2009 The Scripps Research Institute. All rights reserved.
    Copyright (C) 2009 Masaryk University. All rights reserved.

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
#include "ParameterLibrary.h"

struct InputData
{
    // Filenames
    // if the first char is equal to '\0', the filename is not specified
    char fldFilenameAVS[MAX_CHARS];
    char floatingGridFilename[MAX_CHARS];
    char receptorFilename[MAX_CHARS];
    char xyzFilename[MAX_CHARS];
    char parameterLibraryFilename[MAX_CHARS];   // the AD4 parameters .dat file name

    // Grid
    int numGridPointsPerMap;    // for the entire grid
    Vec3i numGridPoints;        // in one axis
    Vec3i numGridPointsDiv2;    // in one axis
    Vec3d gridCornerMin;        // corner of the grid (minimal coordinates)
    Vec3d gridCenter;           // center of mass where the grid is centered on
    double gridSpacing;         // One quarter of a C-C bond length.

    // variables for RECEPTOR:
    // each type is now at most two characters, eg 'NA\0'
    // NB: these are sparse arrays, some entries are not set
    char receptorTypes[NUM_RECEPTOR_TYPES][3];
    int numReceptorTypes; // number of different receptor atom types actually found in receptor PDBQT
    int numReceptorAtoms;

    // length of these arrays is equal to InputData::numReceptorAtoms
    double *charge;
    double *vol;
    double *solpar;
    int *atomType;
    HBondType *hbond;
    Vec4d *receptorAtom; // XYZ = coord, W = charge * coeff_estat * (distDepDiel ? 1 : invDielCal)

    double epsilon[MAX_DIST];

    // for NEW3 desolvation terms
    double solparQ;   // unweighted value restored 3:9:05
    double invDielCal;
    double rSmooth;
    double covHalfWidthSquaredInv;
    double covBarrier;
    Vec3d covalentPoint;         // Cartesian-coordinate of covalent affinity well.

    bool distDepDiel, disorderH;
    
    virtual ~InputData() {}
};

