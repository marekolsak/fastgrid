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

#include <algorithm>
#include "CudaConstantMemory.h"

struct CudaConstantMemory::AtomsConstMem
{
    float4 atoms[4096];
    int numAtoms;
};

CudaConstantMemory::CudaConstantMemory(cudaStream_t stream, CudaInternalAPI *api): api(api), atomsHost(0), stream(stream),
    zIndexArray(0), numAtomSubsets(1), currentZSlice(0), epsilonHost(0)
{
    CUDA_SAFE_CALL(cudaMallocHost((void**)&paramsHost, sizeof(Params)));
}

CudaConstantMemory::~CudaConstantMemory()
{
    if (epsilonHost)
    {
        CUDA_SAFE_CALL(cudaFreeHost(epsilonHost));
        CUDA_SAFE_CALL(cudaFree(params.epsilonDevice));
    }

    CUDA_SAFE_CALL(cudaFreeHost(paramsHost));
    if (zIndexArray)
        CUDA_SAFE_CALL(cudaFreeHost(zIndexArray));
    if (atomsHost)
        CUDA_SAFE_CALL(cudaFreeHost(atomsHost));
}

void CudaConstantMemory::initDistDepDielLookUpTable(const double *epsilon)
{
    int size = sizeof(float) * MAX_DIST;
    CUDA_SAFE_CALL(cudaMallocHost((void**)&epsilonHost, size));
    CUDA_SAFE_CALL(cudaMalloc((void**)&params.epsilonDevice, size));
    paramsHost->epsilonDevice = params.epsilonDevice;

    std::transform(epsilon, epsilon + MAX_DIST, epsilonHost, typecast<float, double>);
    CUDA_SAFE_CALL(cudaMemcpyAsync(params.epsilonDevice, epsilonHost, size, cudaMemcpyHostToDevice, stream));
    api->setDistDepDielLookUpTable(&paramsHost->epsilonDevice, stream);
}

void CudaConstantMemory::setGridMapParameters(const Vec3i &numGridPointsDiv2, double gridSpacing,
                                              const Vec3i &numGridPointsPadded, float *energiesDevice)
{
    params.numGridPointsPadded = make_int3(numGridPointsPadded.x, numGridPointsPadded.y, numGridPointsPadded.z);
    params.gridSpacing = float(gridSpacing);
    params.energiesDevice = energiesDevice;
    params.numGridPointsDiv2 = make_int3(numGridPointsDiv2.x, numGridPointsDiv2.y, numGridPointsDiv2.z);
    *paramsHost = params;
    api->setGridMap(&paramsHost->numGridPointsPadded, &paramsHost->numGridPointsDiv2,
                                   &paramsHost->gridSpacing, &paramsHost->energiesDevice, stream);
}

void CudaConstantMemory::initAtoms(const Vec4d *atoms, int numAtoms, bool calculateSlicesSeparately)
{
    numAtomSubsets = (numAtoms - 1) / api->numAtomsPerKernel + 1;
    int kernelCalls = calculateSlicesSeparately ? params.numGridPointsPadded.z : 1;
    CUDA_SAFE_CALL(cudaMallocHost((void**)&atomsHost, sizeof(AtomsConstMem) * kernelCalls * numAtomSubsets));

    if (calculateSlicesSeparately) // Initialize atoms for later calculating slices separately
    {
        initZIndexArray();

        for (int z = 0; z < params.numGridPointsPadded.z; z++)
        {
            // Set the pointer to memory of this slice
            AtomsConstMem *atomConstMemZBase = atomsHost + z * numAtomSubsets;

            // For each subset numAtomsPerKernel long
            for (int iaStart = 0, i = 0; iaStart < numAtoms; iaStart += api->numAtomsPerKernel, i++)
            {
                int numAtomsInSubset = Mathi::Min(numAtoms - iaStart, api->numAtomsPerKernel);

                AtomsConstMem &thisAtomConstMem = atomConstMemZBase[i];
                thisAtomConstMem.numAtoms = numAtomsInSubset;

                // For each atom in the subset
                for (int ia = 0; ia < numAtomsInSubset; ia++)
                {
                    const Vec4d &atom = atoms[iaStart + ia];

                    // Copy X, Y, Z, and the charge
                    thisAtomConstMem.atoms[ia].x = float(atom.x);
                    thisAtomConstMem.atoms[ia].y = float(atom.y);
                    thisAtomConstMem.atoms[ia].z = float(atom.z);
                    thisAtomConstMem.atoms[ia].w = float(atom.w);
                }
            }
        }
    }
    else // Initialize atoms for later calculating the entire gridmap in one kernel call
    {
        // For each subset numAtomsPerKernel long
        for (int iaStart = 0, i = 0; iaStart < numAtoms; iaStart += api->numAtomsPerKernel, i++)
        {
            int numAtomsInSubset = Mathi::Min(numAtoms - iaStart, api->numAtomsPerKernel);

            AtomsConstMem &thisAtomConstMem = atomsHost[i];
            thisAtomConstMem.numAtoms = numAtomsInSubset;

            // For each atom in the subset
            for (int ia = 0; ia < numAtomsInSubset; ia++)
            {
                const Vec4d &atom = atoms[iaStart + ia];

                // Copy X, Y, Z, and the charge
                thisAtomConstMem.atoms[ia].x = float(atom.x);
                thisAtomConstMem.atoms[ia].y = float(atom.y);
                thisAtomConstMem.atoms[ia].z = float(atom.z);
                thisAtomConstMem.atoms[ia].w = float(atom.w);
            }
        }
    }
}

void CudaConstantMemory::initZIndexArray()
{
    int numSlices = params.numGridPointsPadded.z;

    // Reserve memory for each offseted output index of a slice and precalculate
    CUDA_SAFE_CALL(cudaMallocHost((void**)&zIndexArray, sizeof(int) * numSlices));
    for (int z = 0; z < numSlices; z++)
        zIndexArray[z] = z;
}

void CudaConstantMemory::setZSlice(int z)
{
    currentZSlice = z;

    // Set an offset of output index to constant memory
    api->setSlice(zIndexArray + z, stream);
}

void CudaConstantMemory::setAtomConstMem(int atomSubsetIndex)
{
    AtomsConstMem *thisAtomConstMem = atomsHost + currentZSlice * numAtomSubsets + atomSubsetIndex;
    api->setAtoms(&thisAtomConstMem->numAtoms, thisAtomConstMem->atoms, stream);
}
