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

#include "programparameters.h"
#include "exceptions.h"
#include "inputdataloader.h"
#include "utils.h"
#include "calculategridmaps.h"
#include "electrostatics/electrostatics.h"
#include <new>

void initCovalentMaps(const InputData *input, const GridMapList &gridmaps)
{
    // TODO: once covalent maps are supported, rewrite this function using NVIDIA CUDA

    for (int m = 0; m < gridmaps.getNumMaps(); m++)
        if (gridmaps[m].isCovalent)
#if defined(AG_OPENMP)
    #pragma omp parallel for schedule(dynamic, 1)
#endif
            FOR_EACH_GRID_POINT(gridPos, outputIndex)
            {
                // Distance squared from current grid point to the covalent attachment point
                double rcovSq = (input->covalentPoint - gridPos).MagnitudeSqr() * input->covHalfWidthSquaredInv;
                if (rcovSq < Mathd::Sqr(APPROX_ZERO))
                    rcovSq = Mathd::Sqr(APPROX_ZERO);

                double energy = input->covBarrier * (1 - exp(-0.69314718055994529 * rcovSq)); // -0.69314718055994529 = log(0.5)
                gridmaps[m].energies[outputIndex] = energy;
            }
            END_FOR();
}

void calculateFloatingGrid(const InputData *input, const GridMapList &gridmaps)
{
    // TODO: rewrite this function using NVIDIA CUDA

    // Calculate the so-called "Floating Grid"...
#if defined(AG_OPENMP)
    #pragma omp parallel for schedule(dynamic, 1)
#endif
    FOR_EACH_GRID_POINT(gridPos, outputIndex)
    {
        double minDistSq = Mathd::Sqr(BIG);

        //  Do all Receptor (protein, DNA, etc.) atoms...
        for (int ia = 0; ia < input->numReceptorAtoms; ia++)
        {
            //  Get distance^2 from current grid point to this receptor atom
            Vec3d distance = Vec3d(input->receptorAtom[ia]) - gridPos;
            double distSq = distance.MagnitudeSqr();
            if (distSq < minDistSq)
                minDistSq = distSq;
        }

        gridmaps.getFloatingGridMins()[outputIndex] = float(Mathd::Sqrt(minDistSq));
    }
    END_FOR();
}

// Function: Calculation of interaction energy grids for Autodock.
// Directional H_bonds from Goodford:
// Distance dependent dielectric after Mehler and Solmajer.
// Charge-based desolvation
// Copyright: (C) 2004, TSRI
//
// Authors: Garrett Matthew Morris, Ruth Huey, David S. Goodsell
//
// The Scripps Research Institute
// Department of Molecular Biology, MB5
// 10550 North Torrey Pines Road
// La Jolla, CA 92037-1000.
//
// e-mail: garrett@scripps.edu
// rhuey@scripps.edu
// goodsell@scripps.edu
//
// Helpful suggestions and advice:
// Arthur J. Olson
// Bruce Duncan, Yng Chen, Michael Pique, Victoria Roberts
// Lindy Lindstrom
//
// Inputs: Control file, receptor PDBQT file, parameter file
// Returns: Atomic affinity, desolvation and electrostatic grid maps.
void autogridMain(int argc, char **argv)
{
    // Get the time at the start of the run...
    tms tmsJobStart;
    Clock jobStart = times(&tmsJobStart);

    double versionNumber = 4.00;

    // Initialize the ProgramParameters object, which parses the command-line arguments
    ProgramParameters programParams(argc, argv);

    // Initialize the log file
    LogFile logFile(versionNumber, programParams.getProgramName(), programParams.getLogFilename());
#if defined(AG_OPENMP)
    logFile.printFormatted("Using OpenMP.\n"
                           "- Total number of cores: %i\n"
                           "- Threads available: %i\n\n", omp_get_max_threads());
#endif
#if defined(AG_CUDA)
    logFile.print("Using CUDA.\n");
#endif

    // Declaration of gridmaps, InputDataLoader::load takes care of their initialization
    GridMapList gridmaps(&logFile);

    // Initialization of free energy coefficients and atom parameters
    ParameterLibrary parameterLibrary(&logFile, programParams.getDebugLevel());

    // Reading in the grid parameter file
    InputDataLoader *inputDataLoader = new InputDataLoader(&logFile);
    inputDataLoader->load(programParams.getGridParameterFilename(), gridmaps, parameterLibrary);
    // TODO: shouldn't we put these out of the load function? :
    // - gridmaps initialization code
    // - initialization of atom parameters recIndex/mapIndex (in parameterLibrary)

    // Now we want to make the input data read-only
    const InputData *input = inputDataLoader;

    if (input->floatingGridFilename[0])
        gridmaps.enableFloatingGrid();

    // Inititializing arrays of output energies
    gridmaps.prepareGridmaps(input->numGridPointsPerMap);

    // TODO: add a smarter mechanism of checking for the available disk space, we need to know it as soon as possible.
    // the formerly implemented checks in the middle of calculations were done too late

    // Loading the parameter library from the file
    if (input->parameterLibraryFilename[0])
        parameterLibrary.load(input->parameterLibraryFilename);

    // Writing to AVS-readable gridmaps file (fld)
    saveAVSGridmapsFile(gridmaps, input, programParams, logFile);

#if defined(BOINCCOMPOUND)
    boinc_fraction_done(0.1);
#endif

    Timer *t0 = 0;
    if (programParams.benchmarkEnabled())
        t0 = Timer::startNew("Precalculations        ");

    // Calculating the lookup table of the pairwise interaction energies
    PairwiseInteractionEnergies energyLookup;
    energyLookup.calculate(gridmaps, logFile, input->numReceptorTypes, input->receptorTypes, input->rSmooth);

    // Precalculating the exponential function for receptor and ligand desolvation
    DesolvExpFunc desolvExpFunc(parameterLibrary.coeff_desolv);

    // Calculating bond vectors for directional H-bonds
    BondVectors *bondVectors = new BondVectors(input->numReceptorAtoms, &logFile);
    bondVectors->calculate(input, parameterLibrary);

    if (programParams.benchmarkEnabled())
        t0->stopAndLog(stderr);

    logFile.printFormatted("Beginning grid calculations.\n"
                           "\nCalculating %d grids over %d elements, around %d receptor atoms.\n\n"
                           /*"                    Percent   Estimated Time  Time/this plane\n"
                           "XY-plane  Z-coord   Done      Remaining       Real, User, System\n"
                           "            /Ang              /sec            /sec\n"
                           "________  ________  ________  ______________  __________________________\n\n"*/,
                           gridmaps.getNumMapsInclFloatingGrid(), input->numGridPointsPerMap, input->numReceptorAtoms);

    // TODO: rewrite writing out progress in percents
    /* Former code:

        for (all z)
        {
            tms timesGridStart;
            Clock gridStartTime = times(&timesGridStart);

            for (all y) ...

            tms timerGridEnd;
            Clock gridEndTime = times(&timerGridEnd);
            logFile.printFormatted(" %6d   %8.3lf   %5.1lf%%   ", gridCoordZ, input->gridCornerMin.z + gridPosZ, ((z+1) * 100.0) / input->numGridPoints.z);
            logFile.printTimeInHMS((gridEndTime - gridStartTime) * (input->numGridPoints.z - z));
            logFile.print("  ");
            logFile.printExecutionTimes(gridStartTime, gridEndTime, &timesGridStart, &timerGridEnd);
        }
    */

    // Covalent Atom Types are not yet supported with the new AG4/AD4 atom typing mechanism...
    Timer *t1 = 0;
    if (programParams.benchmarkEnabled())
        t1 = Timer::startNew("Covalent maps          ");
    initCovalentMaps(input, gridmaps);
    if (programParams.benchmarkEnabled())
        t1->stopAndLog(stderr);

    // Calculate the electrostatic map
    calculateElectrostaticMap(input, programParams, gridmaps.getElectrostaticMap());

    // Calculate the atom maps and the desolvation map
    calculateGridmaps(input, programParams, gridmaps, parameterLibrary, energyLookup, desolvExpFunc, bondVectors);

    // Calculate the so-called "floating grid"
    if (gridmaps.containsFloatingGrid())
    {
        Timer *t3 = 0;
        if (programParams.benchmarkEnabled())
            t3 = Timer::startNew("Floating grid          ");
        calculateFloatingGrid(input, gridmaps);
        if (programParams.benchmarkEnabled())
            t3->stopAndLog(stderr);
    }

    delete bondVectors;

#if defined(BOINCCOMPOUND)
    boinc_fraction_done(0.9);
#endif

    // Save all gridmaps
    gridmaps.saveToFiles(input, programParams.getGridParameterFilename());

    delete inputDataLoader;

    // Writing out summary
    gridmaps.logSummary();

    // Get the time at the end of the run and print the difference
    tms tmsJobEnd;
    Clock jobEnd = times(&tmsJobEnd);
    logFile.printExecutionTimesInHMS(jobStart, jobEnd, &tmsJobStart, &tmsJobEnd);
}

int main(int argc, char **argv)
{
    try
    {
        // Initialize BOINC if needed
        boincInit();

        // AutoGrid's main function
        autogridMain(argc, argv);

        // This should not return if used
        boincDone();
        return 0;
    }
    catch (ExitProgram &e)  // the ExitProgram exception is a replacement for C's exit function (we need destructors)
    {
        return e.getExitCode();
    }
    catch (std::bad_alloc &)
    {
        fprintf(stderr, "\n%s: FATAL ERROR: Not enough memory!\n", *argv);
        return 0xBADA110C; // BADALLOC
    }
}
