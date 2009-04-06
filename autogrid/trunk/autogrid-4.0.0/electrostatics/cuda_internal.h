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

#define NUM_ATOMS_PER_KERNEL 4000

#define CUDA_SAFE_CALL(call) checkErrorCUDA(call, __FILE__, __LINE__, __FUNCTION__, #call)
#define CUDA_SAFE_KERNEL(call) checkErrorCUDA((call, cudaGetLastError()), __FILE__, __LINE__, __FUNCTION__, #call)

void setGridMapParametersCUDA(int numGridPointsX, const int2 &numGridPointsDiv2, float gridSpacing);
void setGridMapSliceParametersCUDA(int numAtoms, const float4 *atoms, int outputIndexZBase);
void callKernelCUDA(const dim3 &grid, const dim3 &block, float *outEnergies, const float *epsilon);
void checkErrorCUDA(cudaError e, const char *file, int line, const char *func, const char *code);
