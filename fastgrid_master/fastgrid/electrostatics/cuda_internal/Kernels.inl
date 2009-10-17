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
static __constant__ float gridSpacing;
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
static __constant__ int zIndex; // index of slice in the Z direction

// Contains a reciprocal square root and its parameter
struct RsqrtDesc
{
    float x;
    float result; // = rsqrt(x)

    // notice that: sqrt(x) = x/sqrt(x) = x*result
    // that is, we get sqrt with just rsqrt+mul
};

static inline __device__ RsqrtDesc myRsqrtf(float x)
{
    RsqrtDesc i;
    i.x = x;
    i.result = rsqrt(i.x);
    return i;
}

// Forward declarations
template<int DielectricKind>
static inline __device__ float dielectric(RsqrtDesc rsqrt);
template<int DielectricKind>
static inline __device__ float distanceDependentDiel(RsqrtDesc rsqrt);

// Get dielectric using the lookup table
static inline __device__ float lookupEpsilonTable(RsqrtDesc rsqrt)
{
    return epsilon[min(int(A_DIVISOR * rsqrt.x * rsqrt.result), MAX_DIST-1)];
}

// Distance dependent dielectric - global memory
template<>
static inline __device__ float distanceDependentDiel<DistanceDependentDiel_GlobalMem>(RsqrtDesc rsqrt)
{
    return lookupEpsilonTable(rsqrt);
}

// Distance dependent dielectric - constant memory
template<>
static inline __device__ float distanceDependentDiel<DistanceDependentDiel_ConstMem>(RsqrtDesc rsqrt)
{
    return lookupEpsilonTable(rsqrt);
}

// Distance dependent dielectric - texture memory
template<>
static inline __device__ float distanceDependentDiel<DistanceDependentDiel_TextureMem>(RsqrtDesc rsqrt)
{
    return tex1D(epsilonTexture, (float(A_DIVISOR) / (MAX_DIST-1)) * rsqrt.x * rsqrt.result);
}

// Distance dependent dielectric - in-place calculation
template<>
static inline __device__ float distanceDependentDiel<DistanceDependentDiel_InPlace>(RsqrtDesc rsqrt)
{
    return calculateDistDepDielInv<float>(rsqrt.x * rsqrt.result);
}

// Constant dielectric
template<>
static inline __device__ float dielectric<ConstantDiel>(RsqrtDesc rsqrt)
{
    return fminf(rsqrt.result, 2.f);
}

// Distance-dependent dielectric - main function
template<int DielectricKind>
static inline __device__ float dielectric(RsqrtDesc rsqrt)
{
    float ddd = distanceDependentDiel<DielectricKind>(rsqrt);
    return dielectric<ConstantDiel>(rsqrt) * ddd;
}

// Initializes gridPos and outputIndex based on the thread ID
static inline __device__ void initialize(bool unrolling, int calcGranularity, float3 &gridPos, int &outputIndex)
{
    int x, y, z;
    if (calcGranularity == CalcEntireGrid)
    {
        int gridDimZ = numGridPoints.z;
        if (unrolling) gridDimZ /= 8;

        x = (blockIdx.x / gridDimZ) * blockDim.x + threadIdx.x;
        y = blockIdx.y * blockDim.y + threadIdx.y;
        z = blockIdx.x % gridDimZ;
        if (unrolling) z *= 8;
    }
    else
    {
        x = blockIdx.x * blockDim.x + threadIdx.x;
        y = blockIdx.y * blockDim.y + threadIdx.y;
        z = zIndex;
    }

    gridPos.x = (x - numGridPointsDiv2.x) * gridSpacing;
    gridPos.y = (y - numGridPointsDiv2.y) * gridSpacing;
    gridPos.z = (z - numGridPointsDiv2.z) * gridSpacing;

    outputIndex = z * numGridPoints.y * numGridPoints.x + y * numGridPoints.x + x;
}

// The basic kernel, no loop unrolling
template<int DielectricKind, int CalcGranularity>
static __global__ void calc_1Point()
{
    float3 gridPos;
    int outputIndex;
    initialize(false, CalcGranularity, gridPos, outputIndex);

    float energy = 0;

    //  Do all Receptor (protein, DNA, etc.) atoms...
    for (int ia = 0; ia < numAtoms; ia++)
    {
        float dx = gridPos.x - atoms[ia].x;
        float dy = gridPos.y - atoms[ia].y;
        float dz = gridPos.z - atoms[ia].z;

        // The estat forcefield coefficient/weight is premultiplied in .w
        energy += atoms[ia].w * dielectric<DielectricKind>(myRsqrtf(dx*dx + dy*dy + dz*dz));
    }

    deviceEnergies[outputIndex] += energy;
}

// The kernel where the loop over grid points is unrolled by 8
template<int DielectricKind, int CalcGranularity>
static __global__ void calc_8Points()
{
    float3 gridPos;
    int outputIndex;
    initialize(true, CalcGranularity, gridPos, outputIndex);

    float energy0 = 0, energy1 = 0, energy2 = 0, energy3 = 0,
          energy4 = 0, energy5 = 0, energy6 = 0, energy7 = 0;

    //  Do all Receptor (protein, DNA, etc.) atoms...
    for (int ia = 0; ia < numAtoms; ia++)
    {
        float dx = gridPos.x - atoms[ia].x;
        float dy = gridPos.y - atoms[ia].y;
        float dxdySq = dx*dx + dy*dy;

        float dz0 = gridPos.z - atoms[ia].z;
        float dz1 = dz0 + gridSpacing;
        float dz2 = dz1 + gridSpacing;
        float dz3 = dz2 + gridSpacing;
        float dz4 = dz3 + gridSpacing;
        float dz5 = dz4 + gridSpacing;
        float dz6 = dz5 + gridSpacing;
        float dz7 = dz6 + gridSpacing;

        // The estat forcefield coefficient/weight is premultiplied in .w
        energy0 += atoms[ia].w * dielectric<DielectricKind>(myRsqrtf(dxdySq + dz0*dz0));
        energy1 += atoms[ia].w * dielectric<DielectricKind>(myRsqrtf(dxdySq + dz1*dz1));
        energy2 += atoms[ia].w * dielectric<DielectricKind>(myRsqrtf(dxdySq + dz2*dz2));
        energy3 += atoms[ia].w * dielectric<DielectricKind>(myRsqrtf(dxdySq + dz3*dz3));
        energy4 += atoms[ia].w * dielectric<DielectricKind>(myRsqrtf(dxdySq + dz4*dz4));
        energy5 += atoms[ia].w * dielectric<DielectricKind>(myRsqrtf(dxdySq + dz5*dz5));
        energy6 += atoms[ia].w * dielectric<DielectricKind>(myRsqrtf(dxdySq + dz6*dz6));
        energy7 += atoms[ia].w * dielectric<DielectricKind>(myRsqrtf(dxdySq + dz7*dz7));
    }

    int numGridPointsXMulY = numGridPoints.y * numGridPoints.x;
    deviceEnergies[outputIndex] += energy0; outputIndex += numGridPointsXMulY;
    deviceEnergies[outputIndex] += energy1; outputIndex += numGridPointsXMulY;
    deviceEnergies[outputIndex] += energy2; outputIndex += numGridPointsXMulY;
    deviceEnergies[outputIndex] += energy3; outputIndex += numGridPointsXMulY;
    deviceEnergies[outputIndex] += energy4; outputIndex += numGridPointsXMulY;
    deviceEnergies[outputIndex] += energy5; outputIndex += numGridPointsXMulY;
    deviceEnergies[outputIndex] += energy6; outputIndex += numGridPointsXMulY;
    deviceEnergies[outputIndex] += energy7;
}

void stdSetDistDepDielTexture(const cudaArray *ptr, const cudaChannelFormatDesc *desc)
{
    epsilonTexture.normalized = true;
    epsilonTexture.filterMode = cudaFilterModePoint;
    epsilonTexture.addressMode[0] = cudaAddressModeClamp;

    CUDA_SAFE_CALL(cudaBindTextureToArray(&epsilonTexture, ptr, desc));
}

void stdSetDistDepDielLookUpTable(float **devicePtr, cudaStream_t stream)
{
#if defined(USE_DDD_CONSTMEM)
    CUDA_SAFE_CALL(cudaMemcpyToSymbolAsync(::epsilon, *devicePtr, sizeof(float) * MAX_DIST, 0, cudaMemcpyDeviceToDevice, stream));
#else
    CUDA_SAFE_CALL(cudaMemcpyToSymbolAsync(::epsilon, devicePtr,  sizeof(float*),           0, cudaMemcpyHostToDevice,   stream));
#endif
}

void stdSetGridMap(const int3 *numGridPoints, const int3 *numGridPointsDiv2,
                   const float *gridSpacing, float **deviceEnergies, cudaStream_t stream)
{
    CUDA_SAFE_CALL(cudaMemcpyToSymbolAsync(::numGridPoints,        numGridPoints,        sizeof(int3),   0, cudaMemcpyHostToDevice, stream));
    CUDA_SAFE_CALL(cudaMemcpyToSymbolAsync(::numGridPointsDiv2,    numGridPointsDiv2,    sizeof(int3),   0, cudaMemcpyHostToDevice, stream));
    CUDA_SAFE_CALL(cudaMemcpyToSymbolAsync(::gridSpacing,          gridSpacing,          sizeof(float),  0, cudaMemcpyHostToDevice, stream));
    CUDA_SAFE_CALL(cudaMemcpyToSymbolAsync(::deviceEnergies,       deviceEnergies,       sizeof(float*), 0, cudaMemcpyHostToDevice, stream));
}

void stdSetSlice(const int *zIndex, cudaStream_t stream)
{
    CUDA_SAFE_CALL(cudaMemcpyToSymbolAsync(::zIndex, zIndex, sizeof(int), 0, cudaMemcpyHostToDevice, stream));
}

void stdSetAtoms(const int *numAtoms, const float4 *atoms, cudaStream_t stream)
{
    CUDA_SAFE_CALL(cudaMemcpyToSymbolAsync(::atoms,    atoms,    sizeof(float4) * *numAtoms, 0, cudaMemcpyHostToDevice, stream));
    CUDA_SAFE_CALL(cudaMemcpyToSymbolAsync(::numAtoms, numAtoms, sizeof(int),                0, cudaMemcpyHostToDevice, stream));
}

void stdCallKernel(CudaKernelProc kernel, const dim3 &grid, const dim3 &block, cudaStream_t stream)
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
                    return calc_8Points<DistanceDependentDiel_ConstMem, CalcOneSlice>;
                }
        }
        else
        {
            if (distDepDiel)
                switch (dddKind)
                {
                case DistanceDependentDiel_ConstMem:
                    return calc_1Point<DistanceDependentDiel_ConstMem, CalcOneSlice>;
                }
        }
    else
        if (unrollLoop)
        {
            if (distDepDiel)
                switch (dddKind)
                {
                case DistanceDependentDiel_ConstMem:
                    return calc_8Points<DistanceDependentDiel_ConstMem, CalcEntireGrid>;
                }
        }
        else
        {
            if (distDepDiel)
                switch (dddKind)
                {
                case DistanceDependentDiel_ConstMem:
                    return calc_1Point<DistanceDependentDiel_ConstMem, CalcEntireGrid>;
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
                    return calc_8Points<DistanceDependentDiel_GlobalMem, CalcOneSlice>;
                case DistanceDependentDiel_TextureMem:
                    return calc_8Points<DistanceDependentDiel_TextureMem, CalcOneSlice>;
                case DistanceDependentDiel_InPlace:
                    return calc_8Points<DistanceDependentDiel_InPlace, CalcOneSlice>;
                }
            else
                return calc_8Points<ConstantDiel, CalcOneSlice>;
        }
        else
        {
            if (distDepDiel)
                switch (dddKind)
                {
                case DistanceDependentDiel_GlobalMem:
                    return calc_1Point<DistanceDependentDiel_GlobalMem, CalcOneSlice>;
                case DistanceDependentDiel_TextureMem:
                    return calc_1Point<DistanceDependentDiel_TextureMem, CalcOneSlice>;
                case DistanceDependentDiel_InPlace:
                    return calc_1Point<DistanceDependentDiel_InPlace, CalcOneSlice>;
                }
            else
                return calc_1Point<ConstantDiel, CalcOneSlice>;
        }
    else
        if (unrollLoop)
        {
            if (distDepDiel)
                switch (dddKind)
                {
                case DistanceDependentDiel_GlobalMem:
                    return calc_8Points<DistanceDependentDiel_GlobalMem, CalcEntireGrid>;
                case DistanceDependentDiel_TextureMem:
                    return calc_8Points<DistanceDependentDiel_TextureMem, CalcEntireGrid>;
                case DistanceDependentDiel_InPlace:
                    return calc_8Points<DistanceDependentDiel_InPlace, CalcEntireGrid>;
                }
            else
                return calc_8Points<ConstantDiel, CalcEntireGrid>;
        }
        else
        {
            if (distDepDiel)
                switch (dddKind)
                {
                case DistanceDependentDiel_GlobalMem:
                    return calc_1Point<DistanceDependentDiel_GlobalMem, CalcEntireGrid>;
                case DistanceDependentDiel_TextureMem:
                    return calc_1Point<DistanceDependentDiel_TextureMem, CalcEntireGrid>;
                case DistanceDependentDiel_InPlace:
                    return calc_1Point<DistanceDependentDiel_InPlace, CalcEntireGrid>;
                }
            else
                return calc_1Point<ConstantDiel, CalcEntireGrid>;
        }

#endif

    return 0;
}
