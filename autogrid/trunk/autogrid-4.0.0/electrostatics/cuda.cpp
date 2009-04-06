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

template<typename TDst, typename TSrc>
static TDst typecast(TSrc a)
{
    return TDst(a);
}

// Same as cudaMemcpy with the exception that:
// - If cudaMemcpyHostToDevice or cudaMemcpyHostToHost is specified, the TSrc type is converted to the TDst type before the actual copy
// - If cudaMemcpyDeviceToHost is specified, the TSrc type is converted to the TDst type after the actual copy
// - If cudaMemcpyDeviceToDevice is specified, the ExitProgram exception is thrown
// In other words, it adds an automatic conversion between types. Useful for uploading an array of doubles as floats etc.
template<typename TDst, typename TSrc>
static cudaError_t myCudaMemcpy(TDst *dst, const TSrc *src, size_t numElements, cudaMemcpyKind kind)
{
    switch (kind)
    {
    case cudaMemcpyHostToDevice:
    case cudaMemcpyHostToHost:
        {
            TDst *inter = new TDst[numElements];
            std::transform(src, src + numElements, inter, typecast<TDst, TSrc>);
            cudaError_t e = cudaMemcpy(dst, inter, numElements * sizeof(TDst), kind);
            delete [] inter;
            return e;
        }

    case cudaMemcpyDeviceToHost:
        {
            TSrc *inter = new TSrc[numElements];
            cudaError_t e = cudaMemcpy(inter, src, numElements * sizeof(TSrc), kind);
            std::transform(inter, inter + numElements, dst, typecast<TDst, TSrc>);
            delete [] inter;
            return e;
        }

    default:
        throw ExitProgram(0xbad);
    }
}

template<typename TDst, typename TSrc>
void myCudaCopyGridMapPadded(TDst *dst, const Vec3i &numGridPointsDst, const TSrc *src, const Vec3i &numGridPointsSrc, cudaMemcpyKind kind)
{
    Vec3i numGridPointsMin = Vec3i::ScalarOperator(numGridPointsDst, numGridPointsSrc, Mathi::Min);

    for (int z = 0; z < numGridPointsMin.z; z++)
    {
        // Set the base of output indices from z
        unsigned int outputIndexZBaseDst = z * numGridPointsDst.xy.Square();
        unsigned int outputIndexZBaseSrc = z * numGridPointsSrc.xy.Square();

        for (int y = 0; y < numGridPointsMin.y; y++)
        {
            // Set the base of output indices from (z,y)
            unsigned int outputIndexZYBaseDst = outputIndexZBaseDst + y * numGridPointsDst.x;
            unsigned int outputIndexZYBaseSrc = outputIndexZBaseSrc + y * numGridPointsSrc.x;

            // Copy one row in the X axis
            CUDA_SAFE_CALL((myCudaMemcpy<TDst, TSrc>(dst + outputIndexZYBaseDst, src + outputIndexZYBaseSrc,
                                                        numGridPointsMin.x, kind)));
        }
    }
}

void calculateElectrostaticMapCUDA(const InputData *input, GridMap &elecMap)
{
    // Pad/align the grid to a size of the grid block
    dim3 dimBlock(16, 16);
    Vec3i numGridPointsPadded = input->numGridPoints;
    numGridPointsPadded.x = ((numGridPointsPadded.x - 1) / dimBlock.x + 1) * dimBlock.x;
    numGridPointsPadded.y = ((numGridPointsPadded.y - 1) / dimBlock.y + 1) * dimBlock.y;
    int numGridPointsPerMapPadded = numGridPointsPadded.Cube();
    dim3 dimGrid(numGridPointsPadded.x / dimBlock.x,
                 numGridPointsPadded.y / dimBlock.y);

    // Allocate the padded grid in a device memory
    float *outEnergies = 0;
    CUDA_SAFE_CALL(cudaMalloc((void**)&outEnergies, sizeof(float) * numGridPointsPerMapPadded));

    // Copy the initial energies from the original grid to the padded one in the device memory
    // Elements in the area of padding will be uninitialized
    myCudaCopyGridMapPadded<float, double>(outEnergies, numGridPointsPadded, elecMap.energies, input->numGridPoints, cudaMemcpyHostToDevice);

    // Allocate and copy the lookup table of epsilons for a distance-dependent dielectric
    float *epsilon = 0;
    if (input->distDepDiel)
    {
        CUDA_SAFE_CALL(cudaMalloc((void**)&epsilon, sizeof(float) * MAX_DIST));
        CUDA_SAFE_CALL((myCudaMemcpy<float, double>(epsilon, input->epsilon, MAX_DIST, cudaMemcpyHostToDevice)));
    }

    // Set common variables
    float gridSpacing = float(input->gridSpacing);
    setGridMapParametersCUDA(numGridPointsPadded.x, make_int2(input->numGridPointsDiv2.x, input->numGridPointsDiv2.y), gridSpacing);

    // Convert the array of Vec4d's to Vec4f's
    Vec4f *receptorAtom = new Vec4f[input->numReceptorAtoms];
    std::transform(&input->receptorAtom[0], &input->receptorAtom[0] + input->numReceptorAtoms, receptorAtom, typecast<Vec4f, Vec4d>);

    // Subset of atoms which will be moved to the constant memory
    Vec4f atomConstMem[NUM_ATOMS_PER_KERNEL];

    // For each Z
    for (int z = 0; z < input->numGridPoints.z; z++)
    {
        // Set the base of output index
        unsigned int outputIndexZBase = z * numGridPointsPadded.xy.Square();

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

            // Calculate the slice of the grid for the given subset of atoms
            callKernelCUDA(dimGrid, dimBlock, outEnergies, epsilon);
        }
    }

    delete [] receptorAtom;

    // Free the epsilon array on the device
    if (epsilon)
        cudaFree(epsilon);

    // Copy output energies from the device to the host
    myCudaCopyGridMapPadded<double, float>(elecMap.energies, input->numGridPoints, outEnergies, numGridPointsPadded, cudaMemcpyDeviceToHost);

    // Free device memory
    CUDA_SAFE_CALL(cudaFree(outEnergies));
}

void calculateElectrostaticMap(const InputData *input, GridMap &elecMap)
{
    // Get a device count
    int deviceCount;
    CUDA_SAFE_CALL(cudaGetDeviceCount(&deviceCount));

    if (deviceCount)
    {
        // Print a list of devices
        cudaDeviceProp prop;
        fprintf(stderr, "Found %i CUDA devices:\n", deviceCount);
        for (int i = 0; i < deviceCount; i++)
        {
            cudaGetDeviceProperties(&prop, i);
            fprintf(stderr, "%i. Name: %s, Capability: %i.%i, Total global mem: %i, Total const mem: %i\n", i+1,
                            prop.name, prop.major, prop.minor, prop.totalGlobalMem, prop.totalConstMem);
        }

        // TODO: divide the map into the num of devices
        for (int i = 0; i < /* num of devices */ 1; i++)
        {
            // TODO: start a new thread and call this on the subgrid:
            calculateElectrostaticMapCUDA(input, elecMap);
        }
    }
    else
    {
        fprintf(stderr, "No CUDA device found. Switching to the CPU version.\n");

        // CPU fallback
        calculateElectrostaticMapCPU(input, elecMap);
    }
}

#endif
