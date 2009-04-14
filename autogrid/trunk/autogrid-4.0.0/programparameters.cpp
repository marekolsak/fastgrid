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

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "programparameters.h"
#include "utils.h"
#include "exceptions.h"

#if defined(AG_OPENMP)
    #include <omp.h>
    #define OMP_IGNORED
#else
    #define OMP_IGNORED " (ignored)"
#endif

#if defined(AG_CUDA)
    #include "electrostatics/cuda_internal.h"
    #define CUDA_IGNORED
#else
    #define CUDA_IGNORED " (ignored)"
#endif

ProgramParameters::ProgramParameters(int argc, char **argv): debug(0), benchmark(false), nns(true), cutoffGrid(true), cuda(true)
{
    programName[0] = 0;
    gridParameterFilename[0] = 0;
    logFilename[0] = 0;

    parse(argc, argv);
}

void ProgramParameters::parse(int argc, char **argv)
{
    strncpy(programName, argv[0], MAX_CHARS);

    // Loop over arguments
    while(argc > 1)
    {
        if      (strncmp(argv[1], "-u", 2) == 0 ||
                 strncmp(argv[1], "--help", 6) == 0)
        {
            fprintf(stderr, "Usage: %s [OPTIONS]\n"
                            "\n"
                            "Stream control:\n"
                            "  -p FILE           read grid parameters from FILE\n"
                            "  -l FILE           log to FILE\n"
                            "\n"
                            "Miscellaneous:\n"
                            "  -u, --help        display this help and exit\n"
                            "  -d, --debug       print additional debug information\n"
                            "      --benchmark   print execution times to standard error output\n"
                            "      --no-nns      disable the nearest-neighbor-search optimization\n"
                            "      --no-cogrid   disable the cutoff-grid optimization\n"
                            "      --omp N       set OpenMP to use N threads at most" OMP_IGNORED "\n"
                            "      --no-cuda     disable CUDA" CUDA_IGNORED "\n"
                            "      --cuda-enum   enumerate all CUDA devices and exit" CUDA_IGNORED "\n"
                            "      --cuda-dev N  use a CUDA device number N (default: 0)" CUDA_IGNORED "\n"
                            "\n"
#if defined(AG_OPENMP)
                            "Compiled with OpenMP enabled.\n"
#else
                            "Compiled with OpenMP disabled (related options ignored).\n"
#endif
#if defined(AG_CUDA)
                            "Compiled with CUDA enabled.\n"
#else
                            "Compiled with CUDA disabled (related options ignored).\n"
#endif
                            "\n", programName);
            throw ExitProgram(0);
        }
        else if (strncmp(argv[1], "-d", 2) == 0 ||
                 strncmp(argv[1], "--debug", 7) == 0)
            ++debug;
        else if (strncmp(argv[1], "-l", 2) == 0)
        {
            strncpy(logFilename, argv[2], MAX_CHARS);
            ++argv;
            --argc;
        }
        else if (strncmp(argv[1], "-p", 2) == 0)
        {
            strncpy(gridParameterFilename, argv[2], MAX_CHARS);
            ++argv;
            --argc;
        }
        else if (strncmp(argv[1], "--benchmark", 11) == 0)
        {
            benchmark = true;
        }
        else if (strncmp(argv[1], "--no-nns", 8) == 0)
        {
            nns = false;
        }
        else if (strncmp(argv[1], "--no-cogrid", 11) == 0)
        {
            cutoffGrid = false;
        }
        else if (strncmp(argv[1], "--omp", 5) == 0)
        {
#if defined(AG_OPENMP)
            int n;
            if (sscanf(argv[2], "%i", &n) == 1)
                omp_set_num_threads(n);
            else
            {
                fprintf(stderr,"%s: '%s' is not a number\n", programName, argv[2]);
                throw ExitProgram(1);
            }
#endif
            ++argv;
            --argc;
        }
        else if (strncmp(argv[1], "--no-cuda", 9) == 0)
        {
            cuda = false;
        }
        else if (strncmp(argv[1], "--cuda-enum", 11) == 0)
        {
#if defined(AG_CUDA)
            // Get a device count
            int deviceCount;
            CUDA_SAFE_CALL(cudaGetDeviceCount(&deviceCount));

            // Print a list of devices
            cudaDeviceProp prop;
            fprintf(stderr, "Found %i CUDA devices:\n", deviceCount);
            for (int i = 0; i < deviceCount; i++)
            {
                memset(&prop, 0, sizeof(prop));
                CUDA_SAFE_CALL(cudaGetDeviceProperties(&prop, i));
                fprintf(stderr, "%i. %s, Capability: %i.%i, Global memory: %i, Constant memory: %i\n", i,
                                prop.name, prop.major, prop.minor, prop.totalGlobalMem, prop.totalConstMem);
            }
#endif
            throw ExitProgram(0);
        }
        else if (strncmp(argv[1], "--cuda-dev", 10) == 0)
        {
#if defined(AG_CUDA)
            int n;
            if (sscanf(argv[2], "%i", &n) == 1)
                CUDA_SAFE_CALL(cudaSetDevice(n));
            else
            {
                fprintf(stderr, "%s: '%s' is not a number\n", programName, argv[2]);
                throw ExitProgram(1);
            }
#endif
            ++argv;
            --argc;
        }
        else
        {
            fprintf(stderr, "%s: unknown switch %s\n", programName, argv[1]);
            throw ExitProgram(1);
        }

        --argc;
        ++argv;
    }
}
