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
#include "Interface.h"
// dddcm = Distance Dependent Dielectric in Constant Memory

#define DDDCM_NUM_ATOMS_PER_KERNEL 2000

void dddcmSetGridMap(const int3 *numGridPoints, const int3 *numGridPointsDiv2,
                     const float *gridSpacing, float **deviceEnergies, cudaStream_t stream);
void dddcmSetDistDepDielTexture(const cudaArray *ptr, const cudaChannelFormatDesc *desc);
void dddcmSetDistDepDielLookUpTable(float **devicePtr, cudaStream_t stream);
void dddcmSetSlice(const int *zIndex, cudaStream_t stream);
void dddcmSetAtoms(const int *numAtoms, const float4 *atoms, cudaStream_t stream);

CudaKernelProc dddcmGetKernelProc(bool distDepDiel, DielectricKind dddKind, bool calcSlicesSeparately, bool unrollLoop);
void dddcmCallKernel(CudaKernelProc kernel, const dim3 &grid, const dim3 &block, cudaStream_t stream);
