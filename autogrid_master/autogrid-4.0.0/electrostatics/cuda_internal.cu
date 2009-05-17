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
#if DDD_PROFILE == DDD_GLOBALMEM
static __constant__ float *epsilon;
#elif DDD_PROFILE == DDD_CONSTMEM
static __constant__ float epsilon[MAX_DIST];
#elif DDD_PROFILE == DDD_TEXTUREMEM
static texture<float, 1, cudaReadModeElementType> epsilonTexture;
#endif

static __constant__ float *outEnergies;
static __constant__ int2 numGridPointsDiv2;
static __constant__ int numGridPointsX;
static __constant__ float gridSpacing, gridSpacingCoalesced;

// Per-slice parameters
static __constant__ int outputIndexZBase, numAtoms;
static __constant__ float4 atoms[NUM_ATOMS_PER_KERNEL]; // {x, y, (z-gridPosZ)^2, charge}

// Kinds of dielectric
#define CONSTANT_DIEL           0
#define DISTANCE_DEPENDENT_DIEL 1

// Forward declaration
template<int DielectricKind>
static inline __device__ float dielectric(float invR);

// Constant dielectric
template<>
static inline __device__ float dielectric<CONSTANT_DIEL>(float invR)
{
    return fminf(invR, 2.f);
}

// Distance-dependent dielectric
template<>
static inline __device__ float dielectric<DISTANCE_DEPENDENT_DIEL>(float invR)
{
#if (DDD_PROFILE == DDD_GLOBALMEM) || (DDD_PROFILE == DDD_CONSTMEM)
    float e = epsilon[min(int(A_DIVISOR / invR), MAX_DIST-1)];
#elif DDD_PROFILE == DDD_TEXTUREMEM
    float e = tex1D(epsilonTexture, (float(A_DIVISOR) / (MAX_DIST-1)) / invR);
#elif DDD_PROFILE == DDD_INPLACE
    float e = calculateDistDepDielInv<float>(1.0f / invR);
#else
    Error?
#endif

    return dielectric<CONSTANT_DIEL>(invR) * e;
}

// Basic kernel, no loop unrolling, memory accesses are automatically coalesced
template<int DielectricKind>
static __global__ void calcGridPoints1()
{
    float2 gridPos;
    int outputIndex;
    
    // Calculate gridPos, outputIndex
    {
        int x = blockIdx.x * blockDim.x + threadIdx.x;
        int y = blockIdx.y * blockDim.y + threadIdx.y;

        gridPos.x = (x - numGridPointsDiv2.x) * gridSpacing;
        gridPos.y = (y - numGridPointsDiv2.y) * gridSpacing;
        outputIndex = outputIndexZBase + y * numGridPointsX + x;
    }

    float energy = 0;

    //  Do all Receptor (protein, DNA, etc.) atoms...
    for (int ia = 0; ia < numAtoms; ia++)
    {
        // Calculate dx, dy
        float dx = gridPos.x - atoms[ia].x;
        float dy = gridPos.y - atoms[ia].y;

        // The estat forcefield coefficient/weight is premultiplied in .w
        energy += atoms[ia].w * dielectric<DielectricKind>(rsqrtf(dx*dx + dy*dy + atoms[ia].z));
    }

    outEnergies[outputIndex] += energy;
}

// Kernel where the loop over grid points is unrolled 4x, memory accesses are explicitly coalesced
template<int DielectricKind>
static __global__ void calcGridPoints4()
{
    float2 gridPos;
    int outputIndex;

    // Calculate gridPos, outputIndex
    {
        int x = blockIdx.x * blockDim.x + threadIdx.x;
        int y = blockIdx.y * blockDim.y + threadIdx.y;

        gridPos.x = (x - numGridPointsDiv2.x) * gridSpacing;
        gridPos.y = (y - numGridPointsDiv2.y) * gridSpacing;
        outputIndex = outputIndexZBase + y * numGridPointsX + x;
    }

    float energy0 = 0, energy1 = 0, energy2 = 0, energy3 = 0;

    //  Do all Receptor (protein, DNA, etc.) atoms...
    for (int ia = 0; ia < numAtoms; ia++)
    {
        // Calculate dx[i]
        float dx0 = gridPos.x - atoms[ia].x;
        float dx1 = dx0 + gridSpacingCoalesced;
        float dx2 = dx1 + gridSpacingCoalesced;
        float dx3 = dx2 + gridSpacingCoalesced;

        float dy = gridPos.y - atoms[ia].y;
        float dydzSq = dy*dy + atoms[ia].z;

        // The estat forcefield coefficient/weight is premultiplied in .w
        energy0 += atoms[ia].w * dielectric<DielectricKind>(rsqrt(dx0*dx0 + dydzSq));
        energy1 += atoms[ia].w * dielectric<DielectricKind>(rsqrt(dx1*dx1 + dydzSq));
        energy2 += atoms[ia].w * dielectric<DielectricKind>(rsqrt(dx2*dx2 + dydzSq));
        energy3 += atoms[ia].w * dielectric<DielectricKind>(rsqrt(dx3*dx3 + dydzSq));
    }

    outEnergies[outputIndex] += energy0; outputIndex += 16;
    outEnergies[outputIndex] += energy1; outputIndex += 16;
    outEnergies[outputIndex] += energy2; outputIndex += 16;
    outEnergies[outputIndex] += energy3;
}

void setGridMapParametersAsyncCUDA(const int *numGridPointsX, const int2 *numGridPointsDiv2XY,
                                   const float *gridSpacing, const float *gridSpacingCoalesced,
                                   float **epsilon, float **outEnergies, cudaStream_t stream)
{
    // Set common variables
#if DDD_PROFILE == DDD_GLOBALMEM
    CUDA_SAFE_CALL(cudaMemcpyToSymbolAsync("epsilon",              epsilon,              sizeof(float*), 0, cudaMemcpyHostToDevice, stream));
#elif DDD_PROFILE == DDD_CONSTMEM
    CUDA_SAFE_CALL(cudaMemcpyToSymbolAsync("epsilon",              *epsilon,             sizeof(float) * MAX_DIST, 0, cudaMemcpyHostToDevice, stream));
#endif
    CUDA_SAFE_CALL(cudaMemcpyToSymbolAsync("outEnergies",          outEnergies,          sizeof(float*), 0, cudaMemcpyHostToDevice, stream));
    CUDA_SAFE_CALL(cudaMemcpyToSymbolAsync("numGridPointsDiv2",    numGridPointsDiv2XY,  sizeof(int2),   0, cudaMemcpyHostToDevice, stream));
    CUDA_SAFE_CALL(cudaMemcpyToSymbolAsync("numGridPointsX",       numGridPointsX,       sizeof(int),    0, cudaMemcpyHostToDevice, stream));
    CUDA_SAFE_CALL(cudaMemcpyToSymbolAsync("gridSpacing",          gridSpacing,          sizeof(float),  0, cudaMemcpyHostToDevice, stream));
    CUDA_SAFE_CALL(cudaMemcpyToSymbolAsync("gridSpacingCoalesced", gridSpacingCoalesced, sizeof(float),  0, cudaMemcpyHostToDevice, stream));
}

void setEpsilonTexture(const cudaArray *ptr, const cudaChannelFormatDesc *desc)
{
#if DDD_PROFILE == DDD_TEXTUREMEM
    epsilonTexture.normalized = true;
    epsilonTexture.filterMode = cudaFilterModePoint;
    epsilonTexture.addressMode[0] = cudaAddressModeClamp;
    
    CUDA_SAFE_CALL(cudaBindTextureToArray(&epsilonTexture, ptr, desc));
#endif
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
            CUDA_SAFE_KERNEL((calcGridPoints4<DISTANCE_DEPENDENT_DIEL><<<grid, block, stream>>>()));
        else
            CUDA_SAFE_KERNEL((calcGridPoints4<CONSTANT_DIEL><<<grid, block, stream>>>()));
    else
        if (distDepDiel)
            CUDA_SAFE_KERNEL((calcGridPoints1<DISTANCE_DEPENDENT_DIEL><<<grid, block, stream>>>()));
        else
            CUDA_SAFE_KERNEL((calcGridPoints1<CONSTANT_DIEL><<<grid, block, stream>>>()));
}

