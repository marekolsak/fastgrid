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

#pragma once
#include "Cuda_internal.h"
#include "../autogrid.h"

// This class takes care of uploading data asynchronously to constant memory on the GPU,
// which requires storing the data in page-locked system memory first.
class CudaConstantMemory
{
public:
    CudaConstantMemory(cudaStream_t stream);
    ~CudaConstantMemory();
    void setGridMapParameters(const Vec3i &numGridPointsDiv2, double gridSpacing,
                              const Vec3i &numGridPointsPadded, float *energiesDevice);
    void initAtoms(int numAtomsPerKernel, int numAtomSubsets, bool calculateSlicesSeparately,
                   const Vec4d *atoms, int numAtoms); // must be called only once and after setGridMapParameters
    void setZSlice(int z); // should be called before setAtomConstMem if used at all
    void setAtomConstMem(int atomSubsetIndex);

private:
    struct Params
    {
        int3 numGridPointsDiv2;
        int3 numGridPointsPadded;
        float gridSpacing, gridSpacingCoalesced;
        float *energiesDevice;
    };

    struct AtomsConstMem;

    Params *paramsHost, params;
    AtomsConstMem *atomsHost;
    cudaStream_t stream;
    int *zOffsetArray;
    int numAtomSubsets;
    int currentZSlice;

    void initZOffsetArray(); 
};
