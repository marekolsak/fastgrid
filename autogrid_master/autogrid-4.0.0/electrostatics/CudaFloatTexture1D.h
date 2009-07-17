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
#include "cuda_internal/Interface.h"

// One-channel 1D float texture
class CudaFloatTexture1D
{
public:
    CudaFloatTexture1D(int width, const double *data, CudaAction action, cudaStream_t stream, CudaInternalAPI *api);
    ~CudaFloatTexture1D();
    cudaArray *getDeviceArray() { return deviceArray; }

private:
    cudaChannelFormatDesc channelDesc;
    float *hostMem;
    cudaArray *deviceArray;
};
