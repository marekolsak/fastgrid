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

#include "calculateelecmap.h"

// Constant dielectric...
static void calculateElectrostaticMapCD(const InputData *input, GridMap &elecMap)
{
    #if defined(AG_OPENMP)
        #pragma AG_OPENMP_PARALLEL_FOR
    #endif
    FOR_EACH_GRID_POINT(gridPos, outputIndex)
    {
        double energy = elecMap.energies[outputIndex];

        //  Do all Receptor (protein, DNA, etc.) atoms...
        for (int ia = 0; ia < input->numReceptorAtoms; ia++)
        {
            // Get reciprocal of the distance from current grid point to this receptor atom (1 / |receptorAtomCoord - gridPos|)
            double invR = (input->receptorAtomCoord[ia].xyz - gridPos).MagnitudeInv();

            // Both the constant dielectric and the estat forcefield coefficient/weight are premultiplied
            energy += input->receptorAtomCoord[ia].w * Mathd::Min(invR, 2);
        }

        // Round to 3 decimal places
        elecMap.energies[outputIndex] = roundOutput(energy);
    }
    END_FOR();
}

// Distance-dependent dielectric...
static void calculateElectrostaticMapDDD(const InputData *input, GridMap &elecMap)
{
    #if defined(AG_OPENMP)
        #pragma AG_OPENMP_PARALLEL_FOR
    #endif
    FOR_EACH_GRID_POINT(gridPos, outputIndex)
    {
        double energy = elecMap.energies[outputIndex];

        //  Do all Receptor (protein, DNA, etc.) atoms...
        for (int ia = 0; ia < input->numReceptorAtoms; ia++)
        {
            // Get the distance from current grid point to this receptor atom (|receptorAtomCoord - gridPos|)
            double r = (input->receptorAtomCoord[ia].xyz - gridPos).Magnitude();
            double invR = 1 / r;

            // The estat forcefield coefficient/weight is premultiplied
            energy += input->receptorAtomCoord[ia].w * Mathd::Min(invR, 2) * input->epsilon[lookup(r)];
        }

        // Round to 3 decimal places
        elecMap.energies[outputIndex] = roundOutput(energy);
    }
    END_FOR();
}

void calculateElectrostaticMapCPU(const InputData *input, GridMap &elecMap)
{
    if (input->distDepDiel)
        calculateElectrostaticMapDDD(input, elecMap);
    else
        calculateElectrostaticMapCD(input, elecMap);
}

#if !defined(AG_CUDA)

void calculateElectrostaticMap(const InputData *input, GridMap &elecMap)
{
    calculateElectrostaticMapCPU(input, elecMap);
}

#endif
