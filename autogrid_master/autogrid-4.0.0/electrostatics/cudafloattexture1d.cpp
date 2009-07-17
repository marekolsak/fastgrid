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

#include "cudafloattexture1d.h"
#include <algorithm>

CudaFloatTexture1D::CudaFloatTexture1D(int width, const double *data, CudaAction action, cudaStream_t stream)
{
    channelDesc = cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindFloat);

    // Allocate the texture on the GPU...
    CUDA_SAFE_CALL(cudaMallocArray(&deviceArray, &channelDesc, width, 1));
    // ... and in page-locked system memory
    CUDA_SAFE_CALL(cudaMallocHost((void**)&hostMem, sizeof(float) * width));

    // Convert doubles to floats and save them to page-locked system memory
    std::transform(data, data + width, hostMem, typecast<float, double>);

    // Copy floats from the page-locked memory to the GPU
    CUDA_SAFE_CALL(cudaMemcpyToArrayAsync(deviceArray, 0, 0, hostMem, sizeof(float) * width, cudaMemcpyHostToDevice, stream));

    if (action == BindToKernel)
        ciSetDistDepDielTexture(deviceArray, &channelDesc);
}

CudaFloatTexture1D::~CudaFloatTexture1D()
{
    CUDA_SAFE_CALL(cudaFreeArray(deviceArray));
    CUDA_SAFE_CALL(cudaFreeHost(hostMem));
}
