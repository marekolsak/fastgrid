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

#include "cudaconstantmemory.h"

struct CudaConstantMemory::AtomsConstMem
{
    int numAtoms;
    float4 atoms[4096];
};

CudaConstantMemory::CudaConstantMemory(cudaStream_t stream): atomsHost(0), stream(stream),
    zOffsetArray(0), numAtomSubsets(1), currentZSlice(0)
{
    CUDA_SAFE_CALL(cudaMallocHost((void**)&paramsHost, sizeof(Params)));
}

CudaConstantMemory::~CudaConstantMemory()
{
    CUDA_SAFE_CALL(cudaFreeHost(paramsHost));
    if (zOffsetArray)
        CUDA_SAFE_CALL(cudaFreeHost(zOffsetArray));
    if (atomsHost)
        CUDA_SAFE_CALL(cudaFreeHost(atomsHost));
}

void CudaConstantMemory::setGridMapParameters(const Vec3i &numGridPointsDiv2, double gridSpacing,
                                              const Vec3i &numGridPointsPadded, float *energiesDevice)
{
    params.numGridPointsPadded = make_int3(numGridPointsPadded.x, numGridPointsPadded.y, numGridPointsPadded.z);
    params.gridSpacing = float(gridSpacing);
    params.gridSpacingCoalesced = float(gridSpacing * 16);
    params.energiesDevice = energiesDevice;
    params.numGridPointsDiv2 = make_int3(numGridPointsDiv2.x, numGridPointsDiv2.y, numGridPointsDiv2.z);
    *paramsHost = params;
    ciSetGridMapParametersAsync(&paramsHost->numGridPointsPadded, &paramsHost->numGridPointsDiv2,
                                &paramsHost->gridSpacing, &paramsHost->gridSpacingCoalesced,
                                &paramsHost->energiesDevice, stream);
}

void CudaConstantMemory::initAtoms(int numAtomsPerKernel, int numAtomSubsets, bool calculateSlicesSeparately,
                                   const Vec4d *atoms, int numAtoms)
{
    this->numAtomSubsets = numAtomSubsets;
    int kernelCalls = calculateSlicesSeparately ? params.numGridPointsPadded.z : 1;
    CUDA_SAFE_CALL(cudaMallocHost((void**)&atomsHost, sizeof(AtomsConstMem) * kernelCalls * numAtomSubsets));

    if (calculateSlicesSeparately) // Initialize atoms for later calculating slices separately
    {
        initZOffsetArray();

        for (int z = 0; z < params.numGridPointsPadded.z; z++)
        {
            // Set the pointer to memory of this slice
            AtomsConstMem *atomConstMemZBase = atomsHost + z * numAtomSubsets;

            // Calculate the Z coord of this grid point
            double gridPosZ = (z - params.numGridPointsDiv2.z) * params.gridSpacing;

            // For each subset numAtomsPerKernel long
            for (int iaStart = 0, i = 0; iaStart < numAtoms; iaStart += numAtomsPerKernel, i++)
            {
                int numAtomsInSubset = Mathi::Min(numAtoms - iaStart, numAtomsPerKernel);

                AtomsConstMem &thisAtomConstMem = atomConstMemZBase[i];
                thisAtomConstMem.numAtoms = numAtomsInSubset;

                // For each atom in the subset
                for (int ia = 0; ia < numAtomsInSubset; ia++)
                {
                    const Vec4d &atom = atoms[iaStart + ia];

                    // Copy X, Y and the charge, precalculate distanceZ^2
                    thisAtomConstMem.atoms[ia].x = float(atom.x);
                    thisAtomConstMem.atoms[ia].y = float(atom.y);
                    thisAtomConstMem.atoms[ia].z = float(Mathd::Sqr(atom.z - gridPosZ)); // Precalculate (Z - gridPosZ)^2
                    thisAtomConstMem.atoms[ia].w = float(atom.w);
                }
            }
        }
    }
    else // Initialize atoms for later calculating the entire gridmap in one kernel call
    {
        // For each subset numAtomsPerKernel long
        for (int iaStart = 0, i = 0; iaStart < numAtoms; iaStart += numAtomsPerKernel, i++)
        {
            int numAtomsInSubset = Mathi::Min(numAtoms - iaStart, numAtomsPerKernel);

            AtomsConstMem &thisAtomConstMem = atomsHost[i];
            thisAtomConstMem.numAtoms = numAtomsInSubset;

            // For each atom in the subset
            for (int ia = 0; ia < numAtomsInSubset; ia++)
            {
                const Vec4d &atom = atoms[iaStart + ia];

                // Copy X, Y, Z and the charge
                thisAtomConstMem.atoms[ia].x = float(atom.x);
                thisAtomConstMem.atoms[ia].y = float(atom.y);
                thisAtomConstMem.atoms[ia].z = float(atom.z);
                thisAtomConstMem.atoms[ia].w = float(atom.w);
            }
        }
    }
}

void CudaConstantMemory::initZOffsetArray()
{
    int numSlices = params.numGridPointsPadded.z;
    int numGridPointsPaddedXMulY = params.numGridPointsPadded.x * params.numGridPointsPadded.y;

    // Reserve memory for each offseted output index of a slice and precalculate
    CUDA_SAFE_CALL(cudaMallocHost((void**)&zOffsetArray, sizeof(int) * numSlices));
    for (int z = 0; z < numSlices; z++)
        zOffsetArray[z] = z * numGridPointsPaddedXMulY;
}

void CudaConstantMemory::setZSlice(int z)
{
    currentZSlice = z;

    // Set an offset of output index to constant memory
    ciSetGridMapSliceParametersAsync(zOffsetArray + z, stream);
}

void CudaConstantMemory::setAtomConstMem(int atomSubsetIndex)
{
    AtomsConstMem *thisAtomConstMem = atomsHost + currentZSlice * numAtomSubsets + atomSubsetIndex;
    ciSetGridMapKernelParametersAsync(&thisAtomConstMem->numAtoms, thisAtomConstMem->atoms, stream);
}
