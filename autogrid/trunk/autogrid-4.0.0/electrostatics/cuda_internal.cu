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
#include <cstdio>

// Grid size and spacing
static __constant__ float *epsilon, *outEnergies;
static __constant__ int2 numGridPointsDiv2;
static __constant__ int numGridPointsX;
static __constant__ float gridSpacing, gridSpacingCoalesced;

// Per-slice parameters
static __constant__ int outputIndexZBase, numAtoms;
static __constant__ float4 atoms[NUM_ATOMS_PER_KERNEL]; // {x, y, (z-gridPosZ)^2, charge}

static inline __device__ float getDistDepDiel(float invR)
{
#if DDD_PROFILE == DDD_GLOBALMEM

    int index = angstromToIndex<int>(1.0f / invR);
    float e = epsilon[index];

#elif DDD_PROFILE == DDD_TEXTUREMEM

    TODO? // TODO

#elif DDD_PROFILE == DDD_INPLACE

    float e = calculateDistDepDielInv<float>(1.0f / invR);

#else
    Error?
#endif

    return min(invR, 2.f) * e;
}

// Basic kernel, no loop unrolling
template<int DistanceDependentDielectric>
static __global__ void calcGridPoints1()
{
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;

    float gridPosX = (x - numGridPointsDiv2.x) * gridSpacing;
    float gridPosY = (y - numGridPointsDiv2.y) * gridSpacing;

    float energy = 0;

    //  Do all Receptor (protein, DNA, etc.) atoms...
    for (int ia = 0; ia < numAtoms; ia++)
    {
        float dy = gridPosY - atoms[ia].y;
        float dydzSq = dy*dy + atoms[ia].z;

        // Calculate dx
        float dx = gridPosX - atoms[ia].x;

        // The estat forcefield coefficient/weight is premultiplied
        if (DistanceDependentDielectric)
            energy += atoms[ia].w * getDistDepDiel(rsqrt(dx*dx + dydzSq));
        else
            energy += atoms[ia].w * min(rsqrt(dx*dx + dydzSq), 2.f);
    }

    int outputIndex = outputIndexZBase + y * numGridPointsX + x;
    outEnergies[outputIndex] += energy;
}

// Kernel where the loop over grid points is unrolled by 4
template<int DistanceDependentDielectric>
static __global__ void calcGridPoints4()
{
    int x = (blockIdx.x * blockDim.x + threadIdx.x) * 4;
    int y = blockIdx.y * blockDim.y + threadIdx.y;

    float gridPosX = (x - numGridPointsDiv2.x) * gridSpacing;
    float gridPosY = (y - numGridPointsDiv2.y) * gridSpacing;

    float energy0 = 0, energy1 = 0, energy2 = 0, energy3 = 0;

    //  Do all Receptor (protein, DNA, etc.) atoms...
    for (int ia = 0; ia < numAtoms; ia++)
    {
        float dy = gridPosY - atoms[ia].y;
        float dydzSq = dy*dy + atoms[ia].z;

        // Calculate dx[i]
        float dx0, dx1, dx2, dx3;
        dx0 = gridPosX - atoms[ia].x;
        dx1 = dx0 + gridSpacing;
        dx2 = dx1 + gridSpacing;
        dx3 = dx2 + gridSpacing;

        // The estat forcefield coefficient/weight is premultiplied
        if (DistanceDependentDielectric)
        {
            energy0 += atoms[ia].w * getDistDepDiel(rsqrt(dx0*dx0 + dydzSq));
            energy1 += atoms[ia].w * getDistDepDiel(rsqrt(dx1*dx1 + dydzSq));
            energy2 += atoms[ia].w * getDistDepDiel(rsqrt(dx2*dx2 + dydzSq));
            energy3 += atoms[ia].w * getDistDepDiel(rsqrt(dx3*dx3 + dydzSq));
        }
        else
        {
            energy0 += atoms[ia].w * min(rsqrt(dx0*dx0 + dydzSq), 2.f);
            energy1 += atoms[ia].w * min(rsqrt(dx1*dx1 + dydzSq), 2.f);
            energy2 += atoms[ia].w * min(rsqrt(dx2*dx2 + dydzSq), 2.f);
            energy3 += atoms[ia].w * min(rsqrt(dx3*dx3 + dydzSq), 2.f);
        }
    }

    int outputIndex = outputIndexZBase + y * numGridPointsX + x;

    outEnergies[outputIndex++] += energy0;
    outEnergies[outputIndex++] += energy1;
    outEnergies[outputIndex++] += energy2;
    outEnergies[outputIndex  ] += energy3;
}

void setGridMapParametersAsyncCUDA(const int *numGridPointsX, const int2 *numGridPointsDiv2XY,
                                   const float *gridSpacing, const float *gridSpacingCoalesced,
                                   float **epsilon, float **outEnergies, cudaStream_t stream)
{
    // Set common variables
    CUDA_SAFE_CALL(cudaMemcpyToSymbolAsync("epsilon",              epsilon,              sizeof(float*), 0, cudaMemcpyHostToDevice, stream));
    CUDA_SAFE_CALL(cudaMemcpyToSymbolAsync("outEnergies",          outEnergies,          sizeof(float*), 0, cudaMemcpyHostToDevice, stream));
    CUDA_SAFE_CALL(cudaMemcpyToSymbolAsync("numGridPointsDiv2",    numGridPointsDiv2XY,  sizeof(int2),   0, cudaMemcpyHostToDevice, stream));
    CUDA_SAFE_CALL(cudaMemcpyToSymbolAsync("numGridPointsX",       numGridPointsX,       sizeof(int),    0, cudaMemcpyHostToDevice, stream));
    CUDA_SAFE_CALL(cudaMemcpyToSymbolAsync("gridSpacing",          gridSpacing,          sizeof(float),  0, cudaMemcpyHostToDevice, stream));
    CUDA_SAFE_CALL(cudaMemcpyToSymbolAsync("gridSpacingCoalesced", gridSpacingCoalesced, sizeof(float),  0, cudaMemcpyHostToDevice, stream));
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

void callKernelAsyncCUDA(const dim3 &grid, const dim3 &block, bool unrollLoop, bool distDepDiel, cudaStream_t stream)
{
    if (unrollLoop)
        if (distDepDiel)
            CUDA_SAFE_KERNEL((calcGridPoints4<1><<<grid, block, stream>>>()));
        else
            CUDA_SAFE_KERNEL((calcGridPoints4<0><<<grid, block, stream>>>()));
    else
        if (distDepDiel)
            CUDA_SAFE_KERNEL((calcGridPoints1<1><<<grid, block, stream>>>()));
        else
            CUDA_SAFE_KERNEL((calcGridPoints1<0><<<grid, block, stream>>>()));
}

