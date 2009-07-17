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
#include "cudaevents.h"
#include "cudagridmap.h"
#include "cudafloattexture1d.h"
#include "cudaconstantmemory.h"
#include "../exceptions.h"
#include "../openthreads/Thread"
#include <algorithm>

///////////////////////////////////////////////////////////////////////////////////////////////////

static void setupSliceThreadBlocks(int numThreads, int numGridPointsPerThread, const Vec3i &numGridPoints,
                                   Vec3i *numGridPointsPadded, dim3 *dimGrid, dim3 *dimBlock)
{
    if (numThreads % 16 != 0)
        throw ExitProgram(0xbad);

    // One slice at a time
    dimBlock->x = 16;
    dimBlock->y = numThreads / dimBlock->x;

    // Pad/align the grid to a size of the grid block
    numGridPointsPadded->x = align(numGridPoints.x, dimBlock->x * numGridPointsPerThread);
    numGridPointsPadded->y = align(numGridPoints.y, dimBlock->y);
    numGridPointsPadded->z = numGridPoints.z;
    dimGrid->x = numGridPointsPadded->x / (dimBlock->x * numGridPointsPerThread);
    dimGrid->y = numGridPointsPadded->y / dimBlock->y;
}

static void setupGridThreadBlocks(int numThreads, int numGridPointsPerThread, const Vec3i &numGridPoints,
                                  Vec3i *numGridPointsPadded, dim3 *dimGrid, dim3 *dimBlock)
{
    setupSliceThreadBlocks(numThreads, numGridPointsPerThread, numGridPoints, numGridPointsPadded, dimGrid, dimBlock);
    dimGrid->x *= numGridPointsPadded->z;
}

static void callKernel(CudaConstantMemory &constMem, int atomSubsetIndex, KernelProc kernelProc,
                       const dim3 &dimGrid, const dim3 &dimBlock, cudaStream_t stream)
{
    constMem.setAtomConstMem(atomSubsetIndex); // Move atoms to the constant memory

    // Calculate the entire grid for the given subset of atoms
    ciCallKernelAsync(kernelProc, dimGrid, dimBlock, stream);
}

static void calculateElectrostaticMapCUDA(const InputData *input, const ProgramParameters *programParams, GridMap *elecMap)
{
    cudaStream_t stream;

    CUDA_SAFE_CALL(cudaSetDevice(programParams->getDeviceIDCUDA()));
    CUDA_SAFE_CALL(cudaStreamCreate(&stream));
    CudaEvents events(programParams->benchmarkEnabled(), stream);

    // Determine a number of threads per block and a number of grid points calculated in each thread
    int numThreads;
    int numGridPointsPerThread;
    if (programParams->unrollLoopCUDA())
    {
        numThreads = input->distDepDiel ? 16*16 : 16*4;
        numGridPointsPerThread = 4;
    }
    else
    {
        numThreads = 16*24;
        numGridPointsPerThread = 1;
    }

    // Calculate the size of the padded grid and determine dimensions of thread blocks and the overall thread grid
    Vec3i numGridPointsPadded;
    dim3 dimGrid, dimBlock;
    if (programParams->calcSlicesSeparatelyCUDA())
        setupSliceThreadBlocks(numThreads, numGridPointsPerThread, input->numGridPoints,
                               &numGridPointsPadded, &dimGrid, &dimBlock);
    else
        setupGridThreadBlocks(numThreads, numGridPointsPerThread, input->numGridPoints,
                              &numGridPointsPadded, &dimGrid, &dimBlock);

    events.recordInitialization();

    // Create a padded gridmap on the GPU
    CudaGridMap grid(input->numGridPoints, numGridPointsPadded, elecMap->energies, stream);

    // Create a texture for distance-dependent dielectric
    CudaFloatTexture1D *texture = 0;
    if (input->distDepDiel && programParams->getDDDKindCUDA() == DistanceDependentDiel_TextureMem)
        texture = new CudaFloatTexture1D(MAX_DIST, input->epsilon, BindToKernel, stream);

    // This makes use of constant memory easier
    CudaConstantMemory constMem(stream);
    constMem.setGridMapParameters(input->numGridPointsDiv2, input->gridSpacing,
                                  numGridPointsPadded, grid.getEnergiesDevicePtr());

    KernelProc kernelProc = ciGetKernelProc(input->distDepDiel, programParams->getDDDKindCUDA(),
                                            programParams->calcSlicesSeparatelyCUDA(),
                                            programParams->unrollLoopCUDA());

    // The number of subsets we divide atoms into. Each kernel execution contains only a subset of atoms,
    // which is limited by the size of constant memory.
    int numAtomsPerKernel = programParams->getDDDKindCUDA() == DistanceDependentDiel_ConstMem ?
                            NUM_ATOMS_PER_KERNEL_DDD_CONSTMEM : NUM_ATOMS_PER_KERNEL;
    int numAtomSubsets = (input->numReceptorAtoms - 1) / numAtomsPerKernel + 1;

    // Initialize storage for atoms in page-locked system memory
    constMem.initAtoms(numAtomsPerKernel, numAtomSubsets, programParams->calcSlicesSeparatelyCUDA(),
                       input->receptorAtom, input->numReceptorAtoms);

    events.recordStartCalculation();

    // Do the gridmap calculation...
    if (programParams->calcSlicesSeparatelyCUDA())
        // For each Z
        for (int z = 0; z < input->numGridPoints.z; z++)
        {
            constMem.setZSlice(z);
            for (int i = 0; i < numAtomSubsets; i++)
                callKernel(constMem, i, kernelProc, dimGrid, dimBlock, stream);
        }
    else // !calcSlicesSeparately
        for (int i = 0; i < numAtomSubsets; i++)
            callKernel(constMem, i, kernelProc, dimGrid, dimBlock, stream);

    // Finalize
    events.recordEndCalculation();
    grid.copyFromDeviceToHost();
    events.recordFinalization();
    CUDA_SAFE_CALL(cudaStreamSynchronize(stream));
    events.printTimes(input, numGridPointsPadded.Cube());
    grid.readFromHost(elecMap->energies);

    delete texture;

    CUDA_SAFE_CALL(cudaStreamDestroy(stream));
}

void ciCheckError(cudaError e, const char *file, int line, const char *func, const char *code)
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
    CudaThread(const InputData *input, const ProgramParameters *programParams, GridMap *elecMap)
        : input(input), programParams(programParams), elecMap(elecMap) {}

    virtual void run()
    {
        calculateElectrostaticMapCUDA(input, programParams, elecMap);
    }

private:
    const InputData *input;
    const ProgramParameters *programParams;
    GridMap *elecMap;
};

void *calculateElectrostaticMapAsync(const InputData *input, const ProgramParameters &programParams, GridMap &elecMap)
{
    if (programParams.useCUDA())
        if (programParams.useCUDAThread())
        {
            // Create and start the thread
            CudaThread *thread = new CudaThread(input, &programParams, &elecMap);
            thread->start();
            return thread;
        }
        else
            calculateElectrostaticMapCUDA(input, &programParams, &elecMap);
    else
        calculateElectrostaticMapCPU(input, programParams, elecMap);
    return 0;
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
