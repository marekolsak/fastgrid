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
#include "calculateelecmap.h"
#include "calculateelecmap_cuda.h"
#include "../exceptions.h"
#include <algorithm>

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
    dim3 numGridPoints(input->numGridPoints.x, input->numGridPoints.y, input->numGridPoints.z);
    dim3 numGridPointsDiv2(input->numGridPointsDiv2.x, input->numGridPointsDiv2.y, input->numGridPointsDiv2.z);

    float *receptorAtomCoord = new float[input->numReceptorAtoms * 4];
    std::transform(&input->receptorAtomCoord[0].x, &input->receptorAtomCoord[0].x + input->numReceptorAtoms * 4, receptorAtomCoord, typecast<float, double>);

    // Allocate device memory for energies and initialize
    float *outEnergies = 0;
    cudaMalloc((void**)&outEnergies, sizeof(float) * input->numGridPointsPerMap);
    myCudaMemcpy<float, double>(outEnergies, elecMap.energies, input->numGridPointsPerMap, cudaMemcpyHostToDevice);

    // Allocate and copy epsilon[] for distance-dependent dielectric
    float *epsilon = 0;
    if (input->distDepDiel)
    {
        cudaMalloc((void**)&epsilon, sizeof(float) * MAX_DIST);
        myCudaMemcpy<float, double>(epsilon, input->epsilon, MAX_DIST, cudaMemcpyHostToDevice);
    }

    // Calculate the electrostatic map (distance-dependent dielectric) on the device
    cuCalculateElectrostaticMap(outEnergies, numGridPoints, numGridPointsDiv2, float(input->gridSpacing),
                                input->numReceptorAtoms, receptorAtomCoord, epsilon);

    // Free the epsilon array on device
    if (epsilon)
        cudaFree(epsilon);

    // Copy output energies from device to host
    myCudaMemcpy<double, float>(elecMap.energies, outEnergies, input->numGridPointsPerMap, cudaMemcpyDeviceToHost);

    // Free device memory
    cudaFree(outEnergies);

    delete [] receptorAtomCoord;
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

void waitForCUDA()
{
    // TODO: wait for all CUDA threads
}

#endif
