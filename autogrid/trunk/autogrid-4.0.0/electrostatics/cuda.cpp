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

// Same as cudaMemcpyAsync with the exception that:
// - If cudaMemcpyHostToDevice or cudaMemcpyHostToHost is specified, the TSrc type is converted to the TDst type before the actual copy
// - If cudaMemcpyDeviceToHost is specified, the TSrc type is converted to the TDst type after the actual copy
// - If cudaMemcpyDeviceToDevice is specified, the ExitProgram exception is thrown
// In other words, it adds an automatic conversion between types. Useful for uploading an array of doubles as floats etc.
template<typename TDst, typename TSrc>
static cudaError_t myCudaMemcpyAsync(TDst *dst, const TSrc *src, size_t numElements, cudaMemcpyKind kind, cudaStream_t stream)
{
    switch (kind)
    {
    case cudaMemcpyHostToDevice:
    case cudaMemcpyHostToHost:
        {
            TDst *inter = new TDst[numElements];
            std::transform(src, src + numElements, inter, typecast<TDst, TSrc>);
            cudaError_t e = cudaMemcpyAsync(dst, inter, numElements * sizeof(TDst), kind, stream);
            delete [] inter;
            return e;
        }

    case cudaMemcpyDeviceToHost:
        {
            TSrc *inter = new TSrc[numElements];
            cudaError_t e = cudaMemcpyAsync(inter, src, numElements * sizeof(TSrc), kind, stream);
            std::transform(inter, inter + numElements, dst, typecast<TDst, TSrc>);
            delete [] inter;
            return e;
        }

    default:
        throw ExitProgram(0xbad);
    }
}

template<typename TDst, typename TSrc>
static void myCudaCopyGridMapPaddedAsync(TDst *dst, const Vec3i &numGridPointsDst, const TSrc *src, const Vec3i &numGridPointsSrc, cudaMemcpyKind kind, cudaStream_t stream)
{
    //Vec3i numGridPointsMin = Vec3i::ScalarOperator(numGridPointsDst, numGridPointsSrc, Mathi::Min);
    Vec3i numGridPointsMin = Vec3i(Mathi::Min(numGridPointsDst.x, numGridPointsSrc.x),
                                   Mathi::Min(numGridPointsDst.y, numGridPointsSrc.y),
                                   Mathi::Min(numGridPointsDst.z, numGridPointsSrc.z));
    int numGridPointsDstXMulY = numGridPointsDst.x*numGridPointsDst.y;
    int numGridPointsSrcXMulY = numGridPointsSrc.x*numGridPointsSrc.y;

    for (int z = 0; z < numGridPointsMin.z; z++)
    {
        // Set the base of output indices from z
        int outputIndexZBaseDst = z * numGridPointsDstXMulY;
        int outputIndexZBaseSrc = z * numGridPointsSrcXMulY;

        for (int y = 0; y < numGridPointsMin.y; y++)
        {
            // Set the base of output indices from (z,y)
            int outputIndexZYBaseDst = outputIndexZBaseDst + y * numGridPointsDst.x;
            int outputIndexZYBaseSrc = outputIndexZBaseSrc + y * numGridPointsSrc.x;

            // Copy one row in the X axis
            CUDA_SAFE_CALL((myCudaMemcpyAsync<TDst, TSrc>(dst + outputIndexZYBaseDst, src + outputIndexZYBaseSrc, numGridPointsMin.x, kind, stream)));
        }
    }
}

static void calculateElectrostaticMapCUDA(const InputData *input, const ProgramParameters &programParams, GridMap &elecMap)
{
    // Create a CUDA stream
    cudaStream_t stream;
    CUDA_SAFE_CALL(cudaStreamCreate(&stream));

    // Create CUDA events
    cudaEvent_t init = 0, calc = 0, finalize = 0, end = 0;
    if (programParams.benchmarkEnabled())
    {
        CUDA_SAFE_CALL(cudaEventCreate(&init));
        CUDA_SAFE_CALL(cudaEventCreate(&calc));
        CUDA_SAFE_CALL(cudaEventCreate(&finalize));
        CUDA_SAFE_CALL(cudaEventCreate(&end));
    }

    // Pad/align the grid to a size of the grid block
    dim3 dimBlock(2*NUM_GRIDPOINTS_PER_KERNEL, 14);
    Vec3i numGridPointsPadded(align(input->numGridPoints.x, dimBlock.x),
                              align(input->numGridPoints.y, dimBlock.y),
                              input->numGridPoints.z);
    int numGridPointsPerMapPadded = numGridPointsPadded.Cube();
    dim3 dimGrid(numGridPointsPadded.x / dimBlock.x,
                 numGridPointsPadded.y / dimBlock.y);
    dimBlock.x /= NUM_GRIDPOINTS_PER_KERNEL;

    // Record the init event
    if (programParams.benchmarkEnabled())
        CUDA_SAFE_CALL(cudaEventRecord(init, stream));

    // Allocate the padded grid in global memory
    float *outEnergies = 0;
    CUDA_SAFE_CALL(cudaMalloc((void**)&outEnergies, sizeof(float) * numGridPointsPerMapPadded));

    // Allocate the lookup table for a distance-dependent dielectric
    float *epsilon = 0;
#if DDD_PROFILE == DDD_GLOBALMEM
    if (input->distDepDiel)
        CUDA_SAFE_CALL(cudaMalloc((void**)&epsilon, sizeof(float) * MAX_DIST));
#endif

    // Copy the initial energies from the original grid to the padded one in global memory
    // Elements in the area of padding will stay uninitialized
    myCudaCopyGridMapPaddedAsync<float, double>(outEnergies, numGridPointsPadded, elecMap.energies, input->numGridPoints, cudaMemcpyHostToDevice, stream);

    // Copy the lookup table for a distance-dependent dielectric to global memory
    if (input->distDepDiel && epsilon)
        CUDA_SAFE_CALL((myCudaMemcpyAsync<float, double>(epsilon, input->epsilon, MAX_DIST, cudaMemcpyHostToDevice, stream)));

    // Set some variables in constant memory
    float gridSpacing = float(input->gridSpacing);
    int2 numGridPointsDiv2XY = make_int2(input->numGridPointsDiv2.x, input->numGridPointsDiv2.y);
    setGridMapParametersAsyncCUDA(&numGridPointsPadded.x, &numGridPointsDiv2XY, &gridSpacing, stream);

    // The number of subsets we divide atoms into. Each kernel invocation contains only a subset of atoms,
    // which is limited by the size of constant memory.
    int numAtomSubsets = (input->numReceptorAtoms - 1) / NUM_ATOMS_PER_KERNEL + 1;

    // Reserve memory for each output index of a slice. We can't pass stack pointers to functions since the calls are asynchronous.
    int *outputIndexZBase = new int[input->numGridPoints.z];

    // Reserve memory for each kernel invocation. The same reason as above.
    struct AtomConstMem
    {
        int numAtoms;
        float4 atoms[NUM_ATOMS_PER_KERNEL];

    } *atomConstMem = new AtomConstMem[input->numGridPoints.z * numAtomSubsets];

    // Precalculate X*Y
    int numGridPointsPaddedXMulY = numGridPointsPadded.x*numGridPointsPadded.y;

    // Record the calc event
    if (programParams.benchmarkEnabled())
        CUDA_SAFE_CALL(cudaEventRecord(calc, stream));

    // For each Z
    for (int z = 0; z < input->numGridPoints.z; z++)
    {
        // Set the base of output index
        int *outputIndexZBasePtr = outputIndexZBase + z;
        *outputIndexZBasePtr = z * numGridPointsPaddedXMulY;
        setGridMapSliceParametersAsyncCUDA(outputIndexZBasePtr, stream);

        // Calculate the Z coord of this grid point
        double gridPosZ = (z - input->numGridPointsDiv2.z) * gridSpacing;

        // Set the pointer to memory of this slice
        AtomConstMem *atomConstMemZBase = atomConstMem + z * numAtomSubsets;

        // For each subset NUM_ATOMS_PER_KERNEL long
        for (int iaStart = 0, i = 0; iaStart < input->numReceptorAtoms; iaStart += NUM_ATOMS_PER_KERNEL, i++)
        {
            int numAtoms = Mathi::Min(input->numReceptorAtoms - iaStart, NUM_ATOMS_PER_KERNEL);

            AtomConstMem &thisAtomConstMem = atomConstMemZBase[i];
            thisAtomConstMem.numAtoms = numAtoms;

            // For each atom in the subset
            for (int ia = 0; ia < numAtoms; ia++)
            {
                const Vec4d &atom = input->receptorAtom[iaStart + ia];

                // Copy X, Y and the charge, precalculate distanceZ^2
                thisAtomConstMem.atoms[ia].x = float(atom.x);
                thisAtomConstMem.atoms[ia].y = float(atom.y);
                thisAtomConstMem.atoms[ia].z = float(Mathd::Sqr(atom.z - gridPosZ)); // Precalculate (Z - gridPosZ)^2
                thisAtomConstMem.atoms[ia].w = float(atom.w);
            }

            // Move atoms to the constant memory
            setGridMapKernelParametersAsyncCUDA(&thisAtomConstMem.numAtoms, thisAtomConstMem.atoms, stream);

            // Calculate the slice of the grid for the given subset of atoms
            callKernelAsyncCUDA(dimGrid, dimBlock, outEnergies, input->distDepDiel, epsilon, stream);
        }
    }

    // Record the finalize event
    if (programParams.benchmarkEnabled())
        CUDA_SAFE_CALL(cudaEventRecord(finalize, stream));

    // Copy output energies from the device to the host
    myCudaCopyGridMapPaddedAsync<double, float>(elecMap.energies, input->numGridPoints, outEnergies, numGridPointsPadded, cudaMemcpyDeviceToHost, stream);

    // Record the end event, synchronize
    if (programParams.benchmarkEnabled())
        CUDA_SAFE_CALL(cudaEventRecord(end, stream));
    CUDA_SAFE_CALL(cudaStreamSynchronize(stream));

    // Calculate and print times
    if (programParams.benchmarkEnabled())
    {
        float initTime, calcTime, finalizeTime;
        CUDA_SAFE_CALL(cudaEventElapsedTime(&initTime,     init,     calc));
        CUDA_SAFE_CALL(cudaEventElapsedTime(&calcTime,     calc,     finalize));
        CUDA_SAFE_CALL(cudaEventElapsedTime(&finalizeTime, finalize, end));
        double seconds = double(calcTime) / 1000;
        double atoms = input->numReceptorAtoms * double(input->numGridPointsPerMap);
        double atomsPerSec = atoms / seconds;
        fprintf(stderr, "CUDA Initialization: %i ms\n"
                        "CUDA Kernels: %i ms\n"
                        "CUDA Finalization: %i ms\n",
                        int(initTime), int(calcTime), int(finalizeTime));
        fprintf(stderr, "Electrostatics performance: %i million atoms/s\n", int(atomsPerSec / 1000000));

        // Destroy the events
        CUDA_SAFE_CALL(cudaEventDestroy(init));
        CUDA_SAFE_CALL(cudaEventDestroy(calc));
        CUDA_SAFE_CALL(cudaEventDestroy(finalize));
        CUDA_SAFE_CALL(cudaEventDestroy(end));
    }

    // Free the epsilon array on the device
    if (epsilon)
        CUDA_SAFE_CALL(cudaFree(epsilon));

    // Free device memory
    CUDA_SAFE_CALL(cudaFree(outEnergies));

    // Destroy the stream
    CUDA_SAFE_CALL(cudaStreamDestroy(stream));

    delete [] atomConstMem;
    delete [] outputIndexZBase;
}

void calculateElectrostaticMapCPU(const InputData *input, const ProgramParameters &programParams, GridMap &elecMap);

void calculateElectrostaticMap(const InputData *input, const ProgramParameters &programParams, GridMap &elecMap)
{
    if (programParams.useCUDA())
        calculateElectrostaticMapCUDA(input, programParams, elecMap);
    else
        calculateElectrostaticMapCPU(input, programParams, elecMap);
}

#endif
