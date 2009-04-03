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

#include "../autogrid.h"
#include "calculateelecmap_cuda.h"

#define NUM_ATOMS_PER_KERNEL 32

// Grid size and spacing
static __constant__ uint2 numGridPointsDiv2;
static __constant__ unsigned int numGridPointsX;
static __constant__ float gridSpacing;

static __constant__ unsigned int outputIndexZBase, numAtoms;
static __constant__ float4 atoms[NUM_ATOMS_PER_KERNEL]; // {x, y, (z-gridPosZ)^2, charge}

#define DDD

static __global__ void calculateGridPointDDD(float *outEnergies, const float *epsilon)
#include "cuda_kernel.inl"

#undef DDD

static __global__ void calculateGridPointCD(float *outEnergies)
#include "cuda_kernel.inl"


extern "C" __host__ void cuCalculateElectrostaticMap( // Distance-dependent dielectric
                         /* Output          */       float *outEnergies,
                         /* Grid parameters */       dim3 numGridPoints, dim3 numGridPointsDiv2, float gridSpacing,
                         /* Receptor atoms  */       int numReceptorAtoms, const float *receptorAtomCoord,
                         /* Dist-dep dielec.*/       const float *epsilon)
{
    float4 atoms[NUM_ATOMS_PER_KERNEL];

    dim3 dimBlock(numGridPoints.x, numGridPoints.y);
    dim3 dimGrid(1, 1);

    // Set common variables
    cudaMemcpyToSymbol("numGridPointsDiv2", &numGridPointsDiv2.x, sizeof(numGridPointsDiv2.x) * 2);
    cudaMemcpyToSymbol("numGridPointsX", &numGridPoints.x, sizeof(numGridPoints.x));
    cudaMemcpyToSymbol("gridSpacing", &gridSpacing, sizeof(gridSpacing));

    // For each Z
    for (unsigned int z = 0; z < numGridPoints.z; z++)
    {
        float gridPosZ = (z - numGridPointsDiv2.z) * gridSpacing;

        // Set the base of output index
        unsigned int outputIndexZBase = z * numGridPoints.x * numGridPoints.y;
        cudaMemcpyToSymbol("outputIndexZBase", &outputIndexZBase, sizeof(outputIndexZBase));

        // For each subset NUM_ATOMS_PER_KERNEL long
        for (int iaStart = 0; iaStart < numReceptorAtoms; iaStart += NUM_ATOMS_PER_KERNEL)
        {
            unsigned int iaCount = min(numReceptorAtoms - iaStart, NUM_ATOMS_PER_KERNEL);

            // For each atom in the subset
            for (unsigned int ia = 0; ia < iaCount; ia++)
            {
                const float *atom = receptorAtomCoord + (iaStart + ia);
                atoms[ia].x = atom[0];
                atoms[ia].y = atom[1];

                float dz = atom[2] - gridPosZ;
                atoms[ia].z = dz*dz;
                atoms[ia].w = atom[3];
            }

            // Move atoms to the constant memory
            cudaMemcpyToSymbol("atoms", &atoms[0], sizeof(float4) * iaCount);
            cudaMemcpyToSymbol("numAtoms", &iaCount, sizeof(iaCount));

            // Calculate the slice of the grid for the given subset of atoms
            if (epsilon)
                calculateGridPointDDD<<<dimGrid, dimBlock>>>(outEnergies, epsilon);
            else
                calculateGridPointCD<<<dimGrid, dimBlock>>>(outEnergies);
        }
    }
}
