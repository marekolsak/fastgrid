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

#if defined(AG_CUDA)
#include <cstdio>
#include "CudaEvents.h"

CudaEvents::CudaEvents(bool enabled, cudaStream_t stream): enabled(enabled), stream(stream)
{
    if (!enabled) return;
    CUDA_SAFE_CALL(cudaEventCreate(&init));
    CUDA_SAFE_CALL(cudaEventCreate(&calc));
    CUDA_SAFE_CALL(cudaEventCreate(&finalize));
    CUDA_SAFE_CALL(cudaEventCreate(&end));
}

CudaEvents::~CudaEvents()
{
    if (!enabled) return;
    CUDA_SAFE_CALL(cudaEventDestroy(init));
    CUDA_SAFE_CALL(cudaEventDestroy(calc));
    CUDA_SAFE_CALL(cudaEventDestroy(finalize));
    CUDA_SAFE_CALL(cudaEventDestroy(end));
}

void CudaEvents::recordInitialization()
{
    if (!enabled) return;
    CUDA_SAFE_CALL(cudaEventRecord(init, stream));
}

void CudaEvents::recordStartCalculation()
{
    if (!enabled) return;
    CUDA_SAFE_CALL(cudaEventRecord(calc, stream));
}

void CudaEvents::recordEndCalculation()
{
    if (!enabled) return;
    CUDA_SAFE_CALL(cudaEventRecord(finalize, stream));
}

void CudaEvents::recordFinalization()
{
    if (!enabled) return;
    CUDA_SAFE_CALL(cudaEventRecord(end, stream));
}

void CudaEvents::printTimes(const InputData *input, int numGridPointsPerMapPadded)
{
    if (!enabled) return;

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
}

#endif
