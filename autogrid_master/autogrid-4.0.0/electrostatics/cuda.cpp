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
#include "../openthreads/Thread"
#include <algorithm>
#include <cstring>

template<typename TDst, typename TSrc>
static TDst typecast(TSrc a)
{
    return TDst(a);
}

void myCudaCopyGridMapPaddedAsync(float *dst, const Vec3i &numGridPointsDst, const float *src, const Vec3i &numGridPointsSrc, cudaMemcpyKind kind, cudaStream_t stream)
{
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

            // Copy one row in axis X
            CUDA_SAFE_CALL(cudaMemcpyAsync(dst + outputIndexZYBaseDst, src + outputIndexZYBaseSrc, sizeof(float) * numGridPointsMin.x, kind, stream));
        }
    }
}

static void calculateElectrostaticMapCUDA(const InputData *input, const ProgramParameters &programParams, GridMap &elecMap)
{
    // Set a device
    CUDA_SAFE_CALL(cudaSetDevice(programParams.getDeviceID()));

    // Create a CUDA stream
    cudaStream_t stream = 0;
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

    // Determine the block size
    int numGridPointsPerKernel;
    dim3 dimBlock;
    if (programParams.unrollLoopCUDA())
    {
        numGridPointsPerKernel = 4;

        if (input->distDepDiel)
        {
#if DDD_PROFILE == DDD_GLOBALMEM
            dimBlock = dim3(16, 16);
#elif DDD_PROFILE == DDD_TEXTUREMEM
            dimBlock = dim3(16, 16);
#elif DDD_PROFILE == DDD_INPLACE
            dimBlock = dim3(16, 16);
#elif DDD_PROFILE == DDD_CONSTMEM
            dimBlock = dim3(16, 16);
#else
            Error?
#endif
        }
        else
            dimBlock = dim3(16, 4);
    }
    else
    {
        numGridPointsPerKernel = 1;

        if (input->distDepDiel)
        {
#if DDD_PROFILE == DDD_GLOBALMEM
            dimBlock = dim3(16, 24);
#elif DDD_PROFILE == DDD_TEXTUREMEM
            dimBlock = dim3(16, 24);
#elif DDD_PROFILE == DDD_INPLACE
            dimBlock = dim3(16, 24);
#elif DDD_PROFILE == DDD_CONSTMEM
            dimBlock = dim3(16, 24);
#else
            Error?
#endif
        }
        else
            dimBlock = dim3(16, 24);
    }

    // Pad/align the grid to a size of the grid block
    Vec3i numGridPointsPadded(align(input->numGridPoints.x, dimBlock.x * numGridPointsPerKernel),
                              align(input->numGridPoints.y, dimBlock.y),
                              input->numGridPoints.z);
    int numGridPointsPerMapPadded = numGridPointsPadded.Cube();
    dim3 dimGrid(numGridPointsPadded.x / (dimBlock.x * numGridPointsPerKernel),
                 numGridPointsPadded.y / dimBlock.y);

    // Record the init event
    if (programParams.benchmarkEnabled())
        CUDA_SAFE_CALL(cudaEventRecord(init, stream));

    // Allocate the padded grid in global memory
    float *energiesDevice = 0;
    CUDA_SAFE_CALL(cudaMalloc((void**)&energiesDevice, sizeof(float) * numGridPointsPerMapPadded));

    // Convert doubles to floats and save them in page-locked memory
    float *energiesHost = 0;
    CUDA_SAFE_CALL(cudaMallocHost((void**)&energiesHost, sizeof(float) *input->numGridPointsPerMap));
    std::transform(elecMap.energies, elecMap.energies + input->numGridPointsPerMap, energiesHost, typecast<float, double>);

    // Allocate the lookup table for distance-dependent dielectric
    float *epsilonDevice = 0;
#if (DDD_PROFILE == DDD_GLOBALMEM) || (DDD_PROFILE == DDD_CONSTMEM) || (DDD_PROFILE == DDD_TEXTUREMEM)
    float *epsilonHost = 0;
    cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindFloat);
    cudaArray *epsilonArray;
    if (input->distDepDiel)
    {
    #if DDD_PROFILE == DDD_GLOBALMEM
        CUDA_SAFE_CALL(cudaMalloc((void**)&epsilonDevice, sizeof(float) * MAX_DIST));
    #elif DDD_PROFILE == DDD_TEXTUREMEM
        // Allocate the texture for distance-dependent dielectric
        CUDA_SAFE_CALL(cudaMallocArray(&epsilonArray, &channelDesc, MAX_DIST, 1));
    #endif
        CUDA_SAFE_CALL(cudaMallocHost((void**)&epsilonHost, sizeof(float) * MAX_DIST));
    }
#endif

    // Copy the initial energies from the original grid to the padded one in global memory
    // Elements in the area of padding will stay uninitialized
    myCudaCopyGridMapPaddedAsync(energiesDevice, numGridPointsPadded, energiesHost, input->numGridPoints, cudaMemcpyHostToDevice, stream);

#if (DDD_PROFILE == DDD_GLOBALMEM) || (DDD_PROFILE == DDD_CONSTMEM) || (DDD_PROFILE == DDD_TEXTUREMEM)
    if (input->distDepDiel)
    {
        // Convert doubles to floats
        std::transform(input->epsilon, input->epsilon + MAX_DIST, epsilonHost, typecast<float, double>);

    #if DDD_PROFILE == DDD_GLOBALMEM
        // Copy the epsilon table to global memory
        CUDA_SAFE_CALL(cudaMemcpyAsync(epsilonDevice, epsilonHost, sizeof(float) * MAX_DIST, cudaMemcpyHostToDevice, stream));
    #elif DDD_PROFILE == DDD_TEXTUREMEM
        CUDA_SAFE_CALL(cudaMemcpyToArrayAsync(epsilonArray, 0, 0, epsilonHost, sizeof(float) * MAX_DIST, cudaMemcpyHostToDevice, stream));
        setEpsilonTexture(epsilonArray, &channelDesc);
    #endif
    }
#endif

    // Set some variables in constant memory
    struct Params
    {
        float gridSpacing, gridSpacingCoalesced;
        float *epsilonParam, *energiesDevice;
        int2 numGridPointsDiv2XY;
        int numGridPointsPaddedX;
    } *params;
    CUDA_SAFE_CALL(cudaMallocHost((void**)&params, sizeof(Params)));
    params->numGridPointsPaddedX = numGridPointsPadded.x;
    params->gridSpacing = float(input->gridSpacing);
    params->gridSpacingCoalesced = params->gridSpacing * 16;
    params->epsilonParam = epsilonDevice;
#if DDD_PROFILE == DDD_CONSTMEM
    params->epsilonParam = epsilonHost;
#endif
    params->energiesDevice = energiesDevice;
    params->numGridPointsDiv2XY = make_int2(input->numGridPointsDiv2.x, input->numGridPointsDiv2.y);
    setGridMapParametersAsyncCUDA(&params->numGridPointsPaddedX, &params->numGridPointsDiv2XY, &params->gridSpacing, &params->gridSpacingCoalesced,
                                  &params->epsilonParam, &params->energiesDevice, stream);

    // The number of subsets we divide atoms into. Each kernel invocation contains only a subset of atoms,
    // which is limited by the size of constant memory.
    int numAtomSubsets = (input->numReceptorAtoms - 1) / NUM_ATOMS_PER_KERNEL + 1;

    // Reserve memory for each output index of a slice. We can't pass stack pointers to functions since the calls are asynchronous.
    int *outputIndexZBase;
    CUDA_SAFE_CALL(cudaMallocHost((void**)&outputIndexZBase, sizeof(int) * input->numGridPoints.z));

    // Reserve memory for each kernel invocation. The same reason as above.
    struct AtomConstMem
    {
        int numAtoms;
        float4 atoms[NUM_ATOMS_PER_KERNEL];

    } *atomConstMem;
    CUDA_SAFE_CALL(cudaMallocHost((void**)&atomConstMem, sizeof(AtomConstMem) * input->numGridPoints.z * numAtomSubsets));

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
        double gridPosZ = (z - input->numGridPointsDiv2.z) * params->gridSpacing;

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
            callKernelAsyncCUDA(dimGrid, dimBlock, programParams.unrollLoopCUDA(), input->distDepDiel, stream);
        }
    }

    // Record the finalize event
    if (programParams.benchmarkEnabled())
        CUDA_SAFE_CALL(cudaEventRecord(finalize, stream));

    // Copy output energies from the device to the host
    myCudaCopyGridMapPaddedAsync(energiesHost, input->numGridPoints, energiesDevice, numGridPointsPadded, cudaMemcpyDeviceToHost, stream);

    // Record the end event, synchronize
    if (programParams.benchmarkEnabled())
        CUDA_SAFE_CALL(cudaEventRecord(end, stream));
    CUDA_SAFE_CALL(cudaStreamSynchronize(stream));

    // Convert floats to doubles and save them in the elecMap object
    std::transform(energiesHost, energiesHost + input->numGridPointsPerMap, elecMap.energies, typecast<double, float>);

    // Calculate and print times
    if (programParams.benchmarkEnabled())
    {
        float initTime, calcTime, finalizeTime;
        CUDA_SAFE_CALL(cudaEventElapsedTime(&initTime,     init,     calc));
        CUDA_SAFE_CALL(cudaEventElapsedTime(&calcTime,     calc,     finalize));
        CUDA_SAFE_CALL(cudaEventElapsedTime(&finalizeTime, finalize, end));
        double seconds = double(calcTime) / 1000;
        double atoms = input->numReceptorAtoms * double(input->numGridPointsPerMap);
        double atomsPadded = input->numReceptorAtoms * double(numGridPointsPerMapPadded);
        double atomsPerSec = atoms / seconds;
        double atomsPerSecPadded = atomsPadded / seconds;
        fprintf(stderr, "CUDA Initialization: %i ms\n"
                        "CUDA Kernels: %i ms\n"
                        "CUDA Finalization: %i ms\n",
                        int(initTime), int(calcTime), int(finalizeTime));
        fprintf(stderr, "Electrostatics performance: %i million atoms/s\n", int(atomsPerSec / 1000000));
        fprintf(stderr, "Electrostatics performance: %i million atoms/s (including grid points added by padding)\n", int(atomsPerSecPadded / 1000000));

        // Destroy the events
        CUDA_SAFE_CALL(cudaEventDestroy(init));
        CUDA_SAFE_CALL(cudaEventDestroy(calc));
        CUDA_SAFE_CALL(cudaEventDestroy(finalize));
        CUDA_SAFE_CALL(cudaEventDestroy(end));
    }

#if (DDD_PROFILE == DDD_GLOBALMEM) || (DDD_PROFILE == DDD_CONSTMEM) || (DDD_PROFILE == DDD_TEXTUREMEM)
    if (input->distDepDiel)
    {
        // Free the epsilon tables
    #if DDD_PROFILE == DDD_GLOBALMEM
        CUDA_SAFE_CALL(cudaFree(epsilonDevice));
    #elif DDD_PROFILE == DDD_TEXTUREMEM
        CUDA_SAFE_CALL(cudaFreeArray(epsilonArray));
    #endif
        CUDA_SAFE_CALL(cudaFreeHost(epsilonHost));
    }
#endif

    CUDA_SAFE_CALL(cudaFreeHost(params));
    CUDA_SAFE_CALL(cudaFreeHost(atomConstMem));
    CUDA_SAFE_CALL(cudaFreeHost(outputIndexZBase));

    // Free energies
    CUDA_SAFE_CALL(cudaFree(energiesDevice));
    CUDA_SAFE_CALL(cudaFreeHost(energiesHost));

    // Destroy the stream
    CUDA_SAFE_CALL(cudaStreamDestroy(stream));
}

void checkErrorCUDA(cudaError e, const char *file, int line, const char *func, const char *code)
{
    if (e != cudaSuccess)
    {
        const char *str = cudaGetErrorString(e);
        fprintf(stderr, "CUDA error: '%s'\n"
                        "        in file '%s'\n"
                        "        in line %i\n"
                        "        in function '%s'\n"
                        "        in code '%s'\n", str, file, line, func, code);
        if (strstr(str, "launch fail"))
            throw ExitProgram(0xbad);
    }
}
///////////////////////////////////////////////////////////////////////////////////////////////////

void calculateElectrostaticMapCPU(const InputData *input, const ProgramParameters &programParams, GridMap &elecMap);

class CudaThread : public OpenThreads::Thread
{
public:
    CudaThread(const InputData *input, const ProgramParameters &programParams, GridMap &elecMap)
        : input(input), programParams(&programParams), elecMap(&elecMap) {}

    virtual void run()
    {
        calculateElectrostaticMapCUDA(input, *programParams, *elecMap);
    }

private:
    const InputData *input;
    const ProgramParameters *programParams;
    GridMap *elecMap;
};

void *calculateElectrostaticMapAsync(const InputData *input, const ProgramParameters &programParams, GridMap &elecMap)
{
    if (programParams.useCUDA())
    {
        // Create and run the thread
        CudaThread *thread = new CudaThread(input, programParams, elecMap);
        thread->start();
        return thread;
    }
    else
    {
        calculateElectrostaticMapCPU(input, programParams, elecMap);
        return 0;
    }
}

void synchronizeCalculation(void *handle)
{
    if (handle)
    {
        CudaThread *thread = (CudaThread*)handle;

        // Wait until the thread terminates
        thread->join();
        delete thread;
    }
}

#endif
