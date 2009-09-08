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

#pragma once
#include "cuda_internal/Interface.h"
#include "../autogrid.h"

// This class implements a gridmap representation on the GPU and takes care of padding.
class CudaGridMap
{
public:
    // The constructor creates the gridmap, and copies it to the GPU (asynchronous)
    CudaGridMap(const Vec3i &numGridPoints, const Vec3i &numGridPointsPadded, const double *inputEnergies, cudaStream_t stream);
    ~CudaGridMap();
    void copyFromDeviceToHost(); // Copies the gridmap from the GPU to page-locked system memory (asynchronous)
    void readFromHost(double *outputEnergies); // Saves the gridmap into outputEnergies
    float *getEnergiesDevicePtr() { return energiesDevice; }

private:
    cudaStream_t stream;
    Vec3i numGridPoints, numGridPointsPadded;
    float *energiesDevice, *energiesHost;

    void copyGridMapPadded(float *dst,       const Vec3i &numGridPointsDst,
                           const float *src, const Vec3i &numGridPointsSrc,
                           cudaMemcpyKind kind);
};
