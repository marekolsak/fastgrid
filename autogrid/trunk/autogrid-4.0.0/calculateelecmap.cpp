#include "calculateelecmap.h"

// Constant dielectric...
static void calcWithConstDielectric(const InputData *input, const GridMap &elecMap)
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
            //  Get reciprocal of the distance from current grid point to this receptor atom (1 / |gridPos - receptorAtomCoord|)
            double invR = (input->receptorAtomCoord[ia] - gridPos).RMagnitude();

            // Both the constant dielectric and the estat forcefield coefficient/weight are premultiplied
            energy += input->charge_mul_coeffEstat_mulIfContDiel_invDielCal[ia] * Mathd::Min(invR, 2);
        }

        // Round to 3 decimal places
        elecMap.energies[outputIndex] = roundOutput(energy);
    }
    END_FOR();
}

// Distance-dependent dielectric...
static void calcWithDistDepDielectric(const InputData *input, const GridMap &elecMap)
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
            //  Get the distance from current grid point to this receptor atom (|gridPos - receptorAtomCoord|)
            double r = (input->receptorAtomCoord[ia] - gridPos).Magnitude();

            // Distance-dependent dielectric...
            // The estat forcefield coefficient/weight is premultiplied
            energy += input->charge_mul_coeffEstat_mulIfContDiel_invDielCal[ia] * Mathd::Min(1 / r, 2) * input->epsilon[lookup(r)];
        }

        // Round to 3 decimal places
        elecMap.energies[outputIndex] = roundOutput(energy);
    }
    END_FOR();
}

void calculateElectrostaticMap(const InputData *input, const GridMap &elecMap)
{
    if (input->distDepDiel)
    {
        calcWithDistDepDielectric(input, elecMap);
    }
    else
    {
        calcWithConstDielectric(input, elecMap);
    }
}
