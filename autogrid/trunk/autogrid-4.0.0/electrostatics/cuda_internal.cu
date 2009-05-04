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

using namespace std;
#include "cuda_internal.h"

#define AG_CALLCONV __device__
#include "../autogrid.h"

// Grid size and spacing
static __constant__ int2 numGridPointsDiv2;
static __constant__ int numGridPointsX;
static __constant__ float gridSpacing;

// Per-slice parameters
static __constant__ int outputIndexZBase, numAtoms;
static __constant__ float4 atoms[NUM_ATOMS_PER_KERNEL]; // {x, y, (z-gridPosZ)^2, charge}

static __device__ float getDistDepDiel(float invR, const float *epsilon)
{
#if DDD_PROFILE == DDD_GLOBALMEM

    int index = angstromToIndex<int>(1.0f / invR);
    return epsilon[index];

#elif DDD_PROFILE == DDD_TEXTUREMEM

    TODO? // TODO

#elif DDD_PROFILE == DDD_INPLACE

    return calculateDistDepDielInv<float>(1.0f / invR);

#else
    Error?
#endif
}

// Generic kernel
template<int DistanceDependentDielectric>
static __global__ void calcGridPoint(float *outEnergies, const float *epsilon)
{
    int x = (blockIdx.x * blockDim.x + threadIdx.x) * NUM_GRIDPOINTS_PER_KERNEL;
    int y = blockIdx.y * blockDim.y + threadIdx.y;

    float gridPosX = (x - numGridPointsDiv2.x) * gridSpacing;
    float gridPosY = (y - numGridPointsDiv2.y) * gridSpacing;

    float energy[NUM_GRIDPOINTS_PER_KERNEL] = {0};

    //  Do all Receptor (protein, DNA, etc.) atoms...
    for (int ia = 0; ia < numAtoms; ia++)
    {
        float dy = gridPosY - atoms[ia].y;
        float dydzSq = dy*dy + atoms[ia].z;

        // Calculate dx[i]
        float dx[NUM_GRIDPOINTS_PER_KERNEL];
        dx[0] = gridPosX - atoms[ia].x;

        #pragma unroll
        for (int i = 1; i < NUM_GRIDPOINTS_PER_KERNEL; i++)
            dx[i] = dx[i-1] + gridSpacing;

        // The estat forcefield coefficient/weight is premultiplied
        if (DistanceDependentDielectric)
        {
            #pragma unroll
            for (int i = 0; i < NUM_GRIDPOINTS_PER_KERNEL; i++)
            {
                float invR = rsqrt(dx[i]*dx[i] + dydzSq);
                energy[i] += atoms[ia].w * min(invR, 2.f) * getDistDepDiel(invR, epsilon);
            }
        }
        else
        {
            #pragma unroll
            for (int i = 0; i < NUM_GRIDPOINTS_PER_KERNEL; i++)
                energy[i] += atoms[ia].w * min(rsqrt(dx[i]*dx[i] + dydzSq), 2.f);
        }
    }

    int outputIndex = outputIndexZBase + y * numGridPointsX + x;

    #pragma unroll
    for (int i = 0; i < NUM_GRIDPOINTS_PER_KERNEL; i++)
        outEnergies[outputIndex++] += energy[i];
}

void setGridMapParametersAsyncCUDA(const int *numGridPointsX, const int2 *numGridPointsDiv2XY, const float *gridSpacing, cudaStream_t stream)
{
    // Set common variables
    CUDA_SAFE_CALL(cudaMemcpyToSymbolAsync("numGridPointsDiv2", numGridPointsDiv2XY, sizeof(int2),  0, cudaMemcpyHostToDevice, stream));
    CUDA_SAFE_CALL(cudaMemcpyToSymbolAsync("numGridPointsX",    numGridPointsX,      sizeof(int),   0, cudaMemcpyHostToDevice, stream));
    CUDA_SAFE_CALL(cudaMemcpyToSymbolAsync("gridSpacing",       gridSpacing,         sizeof(float), 0, cudaMemcpyHostToDevice, stream));
}

void setGridMapSliceParametersAsyncCUDA(const int *outputIndexZBase, cudaStream_t stream)
{
    CUDA_SAFE_CALL(cudaMemcpyToSymbolAsync("outputIndexZBase", outputIndexZBase, sizeof(int), 0, cudaMemcpyHostToDevice, stream));
}

void setGridMapKernelParametersAsyncCUDA(const int *numAtoms, const float4 *atoms, cudaStream_t stream)
{
    CUDA_SAFE_CALL(cudaMemcpyToSymbolAsync("atoms",    atoms,    sizeof(float4) * *numAtoms, 0, cudaMemcpyHostToDevice, stream));
    CUDA_SAFE_CALL(cudaMemcpyToSymbolAsync("numAtoms", numAtoms, sizeof(int),                0, cudaMemcpyHostToDevice, stream));
}

void callKernelAsyncCUDA(const dim3 &grid, const dim3 &block, float *outEnergies, bool distDepDiel, const float *epsilon, cudaStream_t stream)
{
    if (distDepDiel)
        CUDA_SAFE_KERNEL((calcGridPoint<1><<<grid, block, stream>>>(outEnergies, epsilon)));
    else
        CUDA_SAFE_KERNEL((calcGridPoint<0><<<grid, block, stream>>>(outEnergies, 0)));
}

void checkErrorCUDA(cudaError e, const char *file, int line, const char *func, const char *code)
{
    if (e != cudaSuccess)
        fprintf(stderr, "CUDA error: '%s'\n"
                        "        in file '%s'\n"
                        "        in line %i\n"
                        "        in function '%s'\n"
                        "        in code '%s'\n", cudaGetErrorString(e), file, line, func, code);
}
