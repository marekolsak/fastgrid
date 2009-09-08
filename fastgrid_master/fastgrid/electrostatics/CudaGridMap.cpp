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
#include "CudaGridMap.h"
#include "../Exceptions.h"

CudaGridMap::CudaGridMap(const Vec3i &numGridPoints, const Vec3i &numGridPointsPadded, const double *inputEnergies, cudaStream_t stream)
    : stream(stream), numGridPoints(numGridPoints), numGridPointsPadded(numGridPointsPadded)
{
    // Allocate the padded grid in global memory
    CUDA_SAFE_CALL(cudaMalloc((void**)&energiesDevice, sizeof(float) * numGridPointsPadded.Cube()));

    // Convert doubles to floats and save them in page-locked memory
    int numGridPointsPerMap = numGridPoints.Cube();
    CUDA_SAFE_CALL(cudaMallocHost((void**)&energiesHost, sizeof(float) * numGridPointsPerMap));
    std::transform(inputEnergies, inputEnergies + numGridPointsPerMap, energiesHost, typecast<float, double>);

    // Copy the initial energies from the original grid to the padded one in global memory
    // Elements in the area of padding will stay uninitialized
    copyGridMapPadded(energiesDevice, numGridPointsPadded, energiesHost, numGridPoints, cudaMemcpyHostToDevice);
}

CudaGridMap::~CudaGridMap()
{
    CUDA_SAFE_CALL(cudaFree(energiesDevice));
    CUDA_SAFE_CALL(cudaFreeHost(energiesHost));
}

void CudaGridMap::copyFromDeviceToHost()
{
    copyGridMapPadded(energiesHost, numGridPoints, energiesDevice, numGridPointsPadded, cudaMemcpyDeviceToHost);
}

void CudaGridMap::readFromHost(double *outputEnergies)
{
    std::transform(energiesHost, energiesHost + numGridPoints.Cube(), outputEnergies, typecast<double, float>);
}

void CudaGridMap::copyGridMapPadded(float *dst,       const Vec3i &numGridPointsDst,
                                    const float *src, const Vec3i &numGridPointsSrc,
                                    cudaMemcpyKind kind)
{
    Vec3i numGridPointsMin = Vec3i(Mathi::Min(numGridPointsDst.x, numGridPointsSrc.x),
                                   Mathi::Min(numGridPointsDst.y, numGridPointsSrc.y),
                                   Mathi::Min(numGridPointsDst.z, numGridPointsSrc.z));
    int numGridPointsDstXMulY = numGridPointsDst.x * numGridPointsDst.y;
    int numGridPointsSrcXMulY = numGridPointsSrc.x * numGridPointsSrc.y;

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
