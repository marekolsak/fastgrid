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

#include "electrostatics.h"
#include "../utils.h"

template<int DistanceDependentDielectric>
static void calculateElectrostaticMapGeneric(const InputData *input, GridMap &elecMap)
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
            // Get reciprocal of the distance from current grid point to this receptor atom (1 / |receptorAtom - gridPos|)
            double r = (input->receptorAtom[ia].xyz - gridPos).Magnitude();
            double invR = 1 / r;

            if (DistanceDependentDielectric)
                // The estat forcefield coefficient/weight is premultiplied
                energy += input->receptorAtom[ia].w * Mathd::Min(invR, 2) * input->epsilon[lookup(r)];
            else
                // Both the constant dielectric and the estat forcefield coefficient/weight are premultiplied
                energy += input->receptorAtom[ia].w * Mathd::Min(invR, 2);
        }

        // Round to 3 decimal places
        elecMap.energies[outputIndex] = roundOutput(energy);
    }
    END_FOR();
}

void calculateElectrostaticMapCPU(const InputData *input, const ProgramParameters &programParams, GridMap &elecMap)
{
    Timer *t2 = 0;
    if (programParams.benchmarkEnabled())
        t2 = Timer::startNew("Electrostatic map      ");

    if (input->distDepDiel)
        calculateElectrostaticMapGeneric<1>(input, elecMap);
    else
        calculateElectrostaticMapGeneric<0>(input, elecMap);

    if (programParams.benchmarkEnabled())
    {
        t2->stopAndLog(stderr, false);

        double seconds = t2->getReal() / double(getClocksPerSec());
        double atoms = input->numReceptorAtoms * double(input->numGridPointsPerMap);
        double atomsPerSec = atoms / seconds;
        fprintf(stderr, "Electrostatics performance: %i million atoms/s\n", int(atomsPerSec / 1000000));
    }
}

#if !defined(AG_CUDA)

void calculateElectrostaticMap(const InputData *input, const ProgramParameters &programParams, GridMap &elecMap)
{
    calculateElectrostaticMapCPU(input, programParams, elecMap);
}

#endif
