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
#include <cuda_runtime_api.h>
#include <vector_functions.h>
#include "cudautils.h"

// Define a number of atoms per kernel
#define NUM_ATOMS_PER_KERNEL                4000
#define NUM_ATOMS_PER_KERNEL_DDD_CONSTMEM   2000

typedef __global__ void (*KernelProc)();

// Functions
void ciSetGridMapParametersAsync(const int3 *numGridPoints, const int3 *numGridPointsDiv2,
                                 const float *gridSpacing, const float *gridSpacingCoalesced,
                                 float **deviceEnergies, cudaStream_t stream);
void ciSetDistDepDielTexture(const cudaArray *ptr, const cudaChannelFormatDesc *desc);
void ciSetGridMapSliceParametersAsync(const int *outputOffsetZBase, cudaStream_t stream);
void ciSetGridMapKernelParametersAsync(const int *numAtoms, const float4 *atoms, cudaStream_t stream);

KernelProc ciGetKernelProc(bool distDepDiel, DielectricKind dddKind, bool calcSlicesSeparately, bool unrollLoop);
void ciCallKernelAsync(KernelProc kernel, const dim3 &grid, const dim3 &block, cudaStream_t stream);

///////////////////////////////////////////////////////////////////////////////////////////////////

#define CUDA_SAFE_CALL(call) ciCheckError((call), __FILE__, __LINE__, __FUNCTION__, #call)
#define CUDA_SAFE_KERNEL(call) ciCheckError(((call), cudaGetLastError()), __FILE__, __LINE__, __FUNCTION__, #call)

void ciCheckError(cudaError e, const char *file, int line, const char *func, const char *code);
