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
#include "Interface.h"

#define STD_NUM_ATOMS_PER_KERNEL 4000

void stdSetGridMapParametersAsync(const int3 *numGridPoints, const int3 *numGridPointsDiv2,
                                    const float *gridSpacing, const float *gridSpacingCoalesced,
                                    float **deviceEnergies, cudaStream_t stream);
void stdSetDistDepDielTexture(const cudaArray *ptr, const cudaChannelFormatDesc *desc);
void stdSetDistDepDielLookUpTable(float **devicePtr, cudaStream_t stream);
void stdSetGridMapSliceParametersAsync(const int *outputOffsetZBase, cudaStream_t stream);
void stdSetGridMapKernelParametersAsync(const int *numAtoms, const float4 *atoms, cudaStream_t stream);

CudaKernelProc stdGetKernelProc(bool distDepDiel, DielectricKind dddKind, bool calcSlicesSeparately, bool unrollLoop);
void stdCallKernelAsync(CudaKernelProc kernel, const dim3 &grid, const dim3 &block, cudaStream_t stream);
