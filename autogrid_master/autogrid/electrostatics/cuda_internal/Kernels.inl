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

#include "Interface.h"

#define AG_CALLCONV __device__
#include "../../autogrid.h"
#include <cstdio>

// Static parameters
static texture<float, 1, cudaReadModeElementType> epsilonTexture;
static __constant__ int3 numGridPoints, numGridPointsDiv2;
static __constant__ float gridSpacing, gridSpacingCoalesced;
static __constant__ float *deviceEnergies;
#if defined(USE_DDD_CONSTMEM)
static __constant__ float epsilon[MAX_DIST];
#else
static __constant__ float *epsilon;
#endif

// Atoms
static __constant__ int numAtoms;
static __constant__ float4 atoms[STD_NUM_ATOMS_PER_KERNEL]; // {x, y, z, charge}

// Per-slice parameters
static __constant__ int outputOffsetZBase;

// Forward declarations
template<int DielectricKind>
static inline __device__ float dielectric(float invR);
template<int DielectricKind>
static inline __device__ float distanceDependentDiel(float invR);
template<int CalcGranularity>
static inline __device__ void initialize(float3 &gridPos, int &outputIndex);
template<int CalcGranularity>
static inline __device__ float calc_dzSq(float gridPosZ, float atomZ);

// Get dielectric using the lookup table
static inline __device__ float lookupEpsilonTable(float invR)
{
    return epsilon[min(int(A_DIVISOR / invR), MAX_DIST-1)];
}

// Distance dependent dielectric - global memory
template<>
static inline __device__ float distanceDependentDiel<DistanceDependentDiel_GlobalMem>(float invR)
{
    return lookupEpsilonTable(invR);
}

// Distance dependent dielectric - constant memory
template<>
static inline __device__ float distanceDependentDiel<DistanceDependentDiel_ConstMem>(float invR)
{
    return lookupEpsilonTable(invR);
}

// Distance dependent dielectric - texture memory
template<>
static inline __device__ float distanceDependentDiel<DistanceDependentDiel_TextureMem>(float invR)
{
    return tex1D(epsilonTexture, (float(A_DIVISOR) / (MAX_DIST-1)) / invR);
}

// Distance dependent dielectric - in-place calculation
template<>
static inline __device__ float distanceDependentDiel<DistanceDependentDiel_InPlace>(float invR)
{
    return calculateDistDepDielInv<float>(1.0f / invR);
}

// Constant dielectric
template<>
static inline __device__ float dielectric<ConstantDiel>(float invR)
{
    return fminf(invR, 2.f);
}

// Distance-dependent dielectric - main function
template<int DielectricKind>
static inline __device__ float dielectric(float invR)
{
    float ddd = distanceDependentDiel<DielectricKind>(invR);
    return dielectric<ConstantDiel>(invR) * ddd;
}

// Initializes gridPos and outputIndex based on the thread ID
// Version for calculating one slice at a time
template<>
static inline __device__ void initialize<CalcOneSlice>(float3 &gridPos, int &outputIndex)
{
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;

    gridPos.x = (x - numGridPointsDiv2.x) * gridSpacing;
    gridPos.y = (y - numGridPointsDiv2.y) * gridSpacing;
    gridPos.z = 0; // unused anyway
    outputIndex = outputOffsetZBase + y * numGridPoints.x + x;
}

// Initializes gridPos and outputIndex based on the thread ID
// Version for calculating the entire gridmap at once
template<>
static inline __device__ void initialize<CalcEntireGrid>(float3 &gridPos, int &outputIndex)
{
    int x = (blockIdx.x / numGridPoints.z) * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;
    int z = blockIdx.x % numGridPoints.z;

    gridPos.x = (x - numGridPointsDiv2.x) * gridSpacing;
    gridPos.y = (y - numGridPointsDiv2.y) * gridSpacing;
    gridPos.z = (z - numGridPointsDiv2.z) * gridSpacing;
    outputIndex = z * numGridPoints.y * numGridPoints.x + y * numGridPoints.x + x;
}

// Returns precalculated squared distance between gridPos and the atom in Z direction
template<>
static inline __device__ float calc_dzSq<CalcOneSlice>(float gridPosZ, float atomZ)
{
    return atomZ;
}

// Calculates squared distance between gridPos and the atom in Z direction
template<>
static inline __device__ float calc_dzSq<CalcEntireGrid>(float gridPosZ, float atomZ)
{
    float dz = gridPosZ - atomZ;
    return dz*dz;
}

// The basic kernel, no loop unrolling, memory accesses are automatically coalesced
template<int DielectricKind, int CalcGranularity>
static __global__ void calcSlice_1Point()
{
    float3 gridPos;
    int outputIndex;
    initialize<CalcGranularity>(gridPos, outputIndex);

    float energy = 0;

    //  Do all Receptor (protein, DNA, etc.) atoms...
    for (int ia = 0; ia < numAtoms; ia++)
    {
        // Calculate dx, dy
        float dx = gridPos.x - atoms[ia].x;
        float dy = gridPos.y - atoms[ia].y;
        float dzSq = calc_dzSq<CalcGranularity>(gridPos.z, atoms[ia].z);

        // The estat forcefield coefficient/weight is premultiplied in .w
        energy += atoms[ia].w * dielectric<DielectricKind>(rsqrtf(dx*dx + dy*dy + dzSq));
    }

    deviceEnergies[outputIndex] += energy;
}

// The kernel where the loop over grid points is unrolled 4 times, memory accesses are explicitly coalesced
template<int DielectricKind, int CalcGranularity>
static __global__ void calcSlice_4Points()
{
    float3 gridPos;
    int outputIndex;
    initialize<CalcGranularity>(gridPos, outputIndex);

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
        float dzSq = calc_dzSq<CalcGranularity>(gridPos.z, atoms[ia].z);
        float dydzSq = dy*dy + dzSq;

        // The estat forcefield coefficient/weight is premultiplied in .w
        energy0 += atoms[ia].w * dielectric<DielectricKind>(rsqrt(dx0*dx0 + dydzSq));
        energy1 += atoms[ia].w * dielectric<DielectricKind>(rsqrt(dx1*dx1 + dydzSq));
        energy2 += atoms[ia].w * dielectric<DielectricKind>(rsqrt(dx2*dx2 + dydzSq));
        energy3 += atoms[ia].w * dielectric<DielectricKind>(rsqrt(dx3*dx3 + dydzSq));
    }

    deviceEnergies[outputIndex] += energy0; outputIndex += 16;
    deviceEnergies[outputIndex] += energy1; outputIndex += 16;
    deviceEnergies[outputIndex] += energy2; outputIndex += 16;
    deviceEnergies[outputIndex] += energy3;
}

void stdSetDistDepDielTexture(const cudaArray *ptr, const cudaChannelFormatDesc *desc)
{
    epsilonTexture.normalized = true;
    epsilonTexture.filterMode = cudaFilterModePoint;
    epsilonTexture.addressMode[0] = cudaAddressModeClamp;

    CUDA_SAFE_CALL(cudaBindTextureToArray(&epsilonTexture, ptr, desc));
}

void stdSetDistDepDielLookUpTableAsync(float **devicePtr, cudaStream_t stream)
{
#if defined(USE_DDD_CONSTMEM)
    CUDA_SAFE_CALL(cudaMemcpyToSymbolAsync(::epsilon, *devicePtr, sizeof(float) * MAX_DIST, 0, cudaMemcpyDeviceToDevice, stream));
#else
    CUDA_SAFE_CALL(cudaMemcpyToSymbolAsync(::epsilon, devicePtr,  sizeof(float*),           0, cudaMemcpyHostToDevice,   stream));
#endif
}

void stdSetGridMapParametersAsync(const int3 *numGridPoints, const int3 *numGridPointsDiv2,
                                  const float *gridSpacing, const float *gridSpacingCoalesced,
                                  float **deviceEnergies, cudaStream_t stream)
{
    CUDA_SAFE_CALL(cudaMemcpyToSymbolAsync(::numGridPoints,        numGridPoints,        sizeof(int3),   0, cudaMemcpyHostToDevice, stream));
    CUDA_SAFE_CALL(cudaMemcpyToSymbolAsync(::numGridPointsDiv2,    numGridPointsDiv2,    sizeof(int3),   0, cudaMemcpyHostToDevice, stream));
    CUDA_SAFE_CALL(cudaMemcpyToSymbolAsync(::gridSpacing,          gridSpacing,          sizeof(float),  0, cudaMemcpyHostToDevice, stream));
    CUDA_SAFE_CALL(cudaMemcpyToSymbolAsync(::gridSpacingCoalesced, gridSpacingCoalesced, sizeof(float),  0, cudaMemcpyHostToDevice, stream));
    CUDA_SAFE_CALL(cudaMemcpyToSymbolAsync(::deviceEnergies,       deviceEnergies,       sizeof(float*), 0, cudaMemcpyHostToDevice, stream));
}

void stdSetGridMapSliceParametersAsync(const int *outputOffsetZBase, cudaStream_t stream)
{
    CUDA_SAFE_CALL(cudaMemcpyToSymbolAsync(::outputOffsetZBase, outputOffsetZBase, sizeof(int), 0, cudaMemcpyHostToDevice, stream));
}

void stdSetGridMapKernelParametersAsync(const int *numAtoms, const float4 *atoms, cudaStream_t stream)
{
    CUDA_SAFE_CALL(cudaMemcpyToSymbolAsync(::atoms,    atoms,    sizeof(float4) * *numAtoms, 0, cudaMemcpyHostToDevice, stream));
    CUDA_SAFE_CALL(cudaMemcpyToSymbolAsync(::numAtoms, numAtoms, sizeof(int),                0, cudaMemcpyHostToDevice, stream));
}

void stdCallKernelAsync(CudaKernelProc kernel, const dim3 &grid, const dim3 &block, cudaStream_t stream)
{
    CUDA_SAFE_KERNEL((kernel<<<grid, block, stream>>>()));
}

CudaKernelProc stdGetKernelProc(bool distDepDiel, DielectricKind dddKind, bool calcSlicesSeparately, bool unrollLoop)
{
#if defined(USE_DDD_CONSTMEM)

    if (calcSlicesSeparately)
        if (unrollLoop)
        {
            if (distDepDiel)
                switch (dddKind)
                {
                case DistanceDependentDiel_ConstMem:
                    return calcSlice_4Points<DistanceDependentDiel_ConstMem, CalcOneSlice>;
                }
        }
        else
        {
            if (distDepDiel)
                switch (dddKind)
                {
                case DistanceDependentDiel_ConstMem:
                    return calcSlice_1Point<DistanceDependentDiel_ConstMem, CalcOneSlice>;
                }
        }
    else
        if (unrollLoop)
        {
            if (distDepDiel)
                switch (dddKind)
                {
                case DistanceDependentDiel_ConstMem:
                    return calcSlice_4Points<DistanceDependentDiel_ConstMem, CalcEntireGrid>;
                }
        }
        else
        {
            if (distDepDiel)
                switch (dddKind)
                {
                case DistanceDependentDiel_ConstMem:
                    return calcSlice_1Point<DistanceDependentDiel_ConstMem, CalcEntireGrid>;
                }
        }

#else

    if (calcSlicesSeparately)
        if (unrollLoop)
        {
            if (distDepDiel)
                switch (dddKind)
                {
                case DistanceDependentDiel_GlobalMem:
                    return calcSlice_4Points<DistanceDependentDiel_GlobalMem, CalcOneSlice>;
                case DistanceDependentDiel_TextureMem:
                    return calcSlice_4Points<DistanceDependentDiel_TextureMem, CalcOneSlice>;
                case DistanceDependentDiel_InPlace:
                    return calcSlice_4Points<DistanceDependentDiel_InPlace, CalcOneSlice>;
                }
            else
                return calcSlice_4Points<ConstantDiel, CalcOneSlice>;
        }
        else
        {
            if (distDepDiel)
                switch (dddKind)
                {
                case DistanceDependentDiel_GlobalMem:
                    return calcSlice_1Point<DistanceDependentDiel_GlobalMem, CalcOneSlice>;
                case DistanceDependentDiel_TextureMem:
                    return calcSlice_1Point<DistanceDependentDiel_TextureMem, CalcOneSlice>;
                case DistanceDependentDiel_InPlace:
                    return calcSlice_1Point<DistanceDependentDiel_InPlace, CalcOneSlice>;
                }
            else
                return calcSlice_1Point<ConstantDiel, CalcOneSlice>;
        }
    else
        if (unrollLoop)
        {
            if (distDepDiel)
                switch (dddKind)
                {
                case DistanceDependentDiel_GlobalMem:
                    return calcSlice_4Points<DistanceDependentDiel_GlobalMem, CalcEntireGrid>;
                case DistanceDependentDiel_TextureMem:
                    return calcSlice_4Points<DistanceDependentDiel_TextureMem, CalcEntireGrid>;
                case DistanceDependentDiel_InPlace:
                    return calcSlice_4Points<DistanceDependentDiel_InPlace, CalcEntireGrid>;
                }
            else
                return calcSlice_4Points<ConstantDiel, CalcEntireGrid>;
        }
        else
        {
            if (distDepDiel)
                switch (dddKind)
                {
                case DistanceDependentDiel_GlobalMem:
                    return calcSlice_1Point<DistanceDependentDiel_GlobalMem, CalcEntireGrid>;
                case DistanceDependentDiel_TextureMem:
                    return calcSlice_1Point<DistanceDependentDiel_TextureMem, CalcEntireGrid>;
                case DistanceDependentDiel_InPlace:
                    return calcSlice_1Point<DistanceDependentDiel_InPlace, CalcEntireGrid>;
                }
            else
                return calcSlice_1Point<ConstantDiel, CalcEntireGrid>;
        }

#endif

    return 0;
}
