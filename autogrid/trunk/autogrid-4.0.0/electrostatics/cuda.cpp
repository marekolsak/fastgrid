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

#if defined(AG_CUDA)
#include "electrostatics.h"
#include "cuda_internal.h"
#include "../exceptions.h"
#include <algorithm>
#include <windows.h>

template<typename TDst, typename TSrc>
static TDst typecast(TSrc a)
{
    return TDst(a);
}

// Same as cudaMemcpy with the exception that:
// - if cudaMemcpyHostToDevice is specified, the TSrc type is converted to the TDst type before the actual copy
// - if cudaMemcpyDeviceToHost is specified, the TSrc type is converted to the TDst type after the actual copy
// In other words, it adds an automatic conversion between types. Useful for uploading an array of doubles as floats etc.
template<typename TDst, typename TSrc>
static void myCudaMemcpy(TDst *dst, const TSrc *src, size_t numElements, cudaMemcpyKind kind)
{
    switch (kind)
    {
    case cudaMemcpyHostToDevice:
        {
            TDst *inter = new TDst[numElements];
            std::transform(src, src + numElements, inter, typecast<TDst, TSrc>);
            cudaMemcpy(dst, inter, numElements * sizeof(TDst), kind);
            delete [] inter;
        }
        break;

    case cudaMemcpyDeviceToHost:
        {
            TSrc *inter = new TSrc[numElements];
            cudaMemcpy(inter, src, numElements * sizeof(TSrc), kind);
            std::transform(inter, inter + numElements, dst, typecast<TDst, TSrc>);
            delete [] inter;
        }
        break;

    default:
        if (typeid(TDst) != typeid(TSrc))
            throw ExitProgram(0xbad);

        cudaMemcpy(dst, src, numElements * sizeof(TDst), kind);
    }
}

void calculateElectrostaticMapCUDA(const InputData *input, GridMap &elecMap)
{
    // Convert the array of Vec4ds to Vec4fs
    Vec4f *receptorAtom = new Vec4f[input->numReceptorAtoms];
    std::transform(&input->receptorAtom[0], &input->receptorAtom[0] + input->numReceptorAtoms, receptorAtom, typecast<Vec4f, Vec4d>);

    // Allocate the device memory for energies and initialize it
    float *outEnergies = 0;
    cudaMalloc((void**)&outEnergies, sizeof(float) * input->numGridPointsPerMap);
    myCudaMemcpy<float, double>(outEnergies, elecMap.energies, input->numGridPointsPerMap, cudaMemcpyHostToDevice);

    // Allocate and copy epsilon[] for a distance-dependent dielectric
    float *epsilon = 0;
    if (input->distDepDiel)
    {
        cudaMalloc((void**)&epsilon, sizeof(float) * MAX_DIST);
        myCudaMemcpy<float, double>(epsilon, input->epsilon, MAX_DIST, cudaMemcpyHostToDevice);
    }

    // Set common variables
    float gridSpacing = float(input->gridSpacing);
    setGridMapParametersCUDA(input->numGridPoints.x, make_int2(input->numGridPointsDiv2.x, input->numGridPointsDiv2.y), gridSpacing);
    OutputDebugString(cudaGetErrorString(cudaGetLastError()));
    OutputDebugString("\n");

    // This list of atoms will be moved to a constant memory
    Vec4f atomConstMem[NUM_ATOMS_PER_KERNEL];

    // Execution configuration of a kernel
    dim3 dimBlock(input->numGridPoints.x, input->numGridPoints.y);
    dim3 dimGrid(1, 1);

    // For each Z
    for (int z = 0; z < input->numGridPoints.z; z++)
    {
        // Set the base of output index
        unsigned int outputIndexZBase = z * input->numGridPoints.x * input->numGridPoints.y;

        // Calculate the Z coord of this grid point
        float gridPosZ = (z - input->numGridPointsDiv2.z) * gridSpacing;

        // For each subset NUM_ATOMS_PER_KERNEL long
        for (int iaStart = 0; iaStart < input->numReceptorAtoms; iaStart += NUM_ATOMS_PER_KERNEL)
        {
            int iaCount = Mathi::Min(input->numReceptorAtoms - iaStart, NUM_ATOMS_PER_KERNEL);

            // For each atom in the subset
            for (int ia = 0; ia < iaCount; ia++)
            {
                Vec4f &atom = receptorAtom[iaStart + ia];
                // Copy X, Y and the charge
                atomConstMem[ia].x = atom.x;
                atomConstMem[ia].y = atom.y;
                atomConstMem[ia].w = atom.w;

                // Precalculate (Z - gridPosZ)^2
                float dz = atom.z - gridPosZ;
                atomConstMem[ia].z = dz*dz;
            }

            // Move atoms to the constant memory
            setGridMapSliceParametersCUDA(iaCount, (float4*)atomConstMem, outputIndexZBase);
            OutputDebugString(cudaGetErrorString(cudaGetLastError()));
            OutputDebugString("\n");

            // Calculate the slice of the grid for the given subset of atoms
            callKernelCUDA(dimGrid, dimBlock, outEnergies, epsilon);
            OutputDebugString(cudaGetErrorString(cudaGetLastError()));
            OutputDebugString("\n");
        }
    }

    // Free the epsilon array on device
    if (epsilon)
        cudaFree(epsilon);

    // Copy output energies from device to host
    myCudaMemcpy<double, float>(elecMap.energies, outEnergies, input->numGridPointsPerMap, cudaMemcpyDeviceToHost);

    // Free device memory
    cudaFree(outEnergies);

    delete [] receptorAtom;
}

void calculateElectrostaticMap(const InputData *input, GridMap &elecMap)
{
    // TODO: is there any cuda device?
    if (1)
    {
        // TODO: divide the map to the num of devices
        for (int i = 0; i < /* num of devices */ 1; i++)
        {
            // TODO: start a new thread and call this on the subgrid:
            calculateElectrostaticMapCUDA(input, elecMap);
        }
    }
    else
        // CPU fallback
        calculateElectrostaticMapCPU(input, elecMap);
}

#endif
