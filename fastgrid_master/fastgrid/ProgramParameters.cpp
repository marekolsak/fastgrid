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

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "ProgramParameters.h"
#include "Utils.h"
#include "Exceptions.h"
#if !defined(_WIN32)
    #include "UnixSignals.h"
    #include <sys/wait.h>
#endif

#define ST_ENABLED "enabled"
#define ST_DISABLED "disabled (related options ignored)"

#if defined(AG_OPENMP)
    #define OMP_STATUS ST_ENABLED
#else
    #define OMP_STATUS ST_DISABLED
    #define omp_set_num_threads(n) (n)
#endif

#if defined(AG_CUDA)
    #include "electrostatics/cuda_internal/Interface.h"
    #define CUDA_STATUS ST_ENABLED
#else
    #define CUDA_STATUS ST_DISABLED
#endif

ProgramParameters::ProgramParameters(int argc, char **argv): debug(0), deviceID(0), cutoffGridMem(512), benchmark(false), nns(true),
    cutoffGrid(true), cuda(true), cudaThread(true), calcSlicesSeparately(false), v4(false), cudaUnroll(Unassigned),
    cudaDDDKind(Diel_Unassigned)
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
    while (argc > 1)
    {
        if (0);

        // Stream control:
        else if (cmp2(argv[1], "-l", "--log"))    readParamString(&argc, &argv, logFilename);
        else if (cmp2(argv[1], "-p", "--gpf"))    readParamString(&argc, &argv, gridParameterFilename);

        // Miscellaneous:
        else if (cmp(argv[1], "--benchmark"))     benchmark = true;
        else if (cmp(argv[1], "--cogrid-mem"))    cutoffGridMem = readParamInt(&argc, &argv);
        else if (cmp(argv[1], "--cuda-enum"))     cudaEnumDevicesAndExit();
        else if (cmp(argv[1], "--cuda-dev"))      deviceID = readParamInt(&argc, &argv);
        else if (cmp2(argv[1], "-d", "--debug"))  ++debug;
        else if (cmp2(argv[1], "-u", "--help"))   printHelpAndExit();
        else if (cmp(argv[1], "--no-cuda"))       cuda = false;
        else if (cmp(argv[1], "--omp"))           omp_set_num_threads(readParamInt(&argc, &argv));
        else if (cmp(argv[1], "--timeout"))       setTimeout(readParamInt(&argc, &argv));
        else if (cmp(argv[1], "--v4"))            v4 = true;

        // Advanced:
        else if (cmp(argv[1], "--cuda-ddd=g"))    cudaDDDKind = DistanceDependentDiel_GlobalMem;
        else if (cmp(argv[1], "--cuda-ddd=c"))    cudaDDDKind = DistanceDependentDiel_ConstMem;
        else if (cmp(argv[1], "--cuda-ddd=t"))    cudaDDDKind = DistanceDependentDiel_TextureMem;
        else if (cmp(argv[1], "--cuda-ddd=i"))    cudaDDDKind = DistanceDependentDiel_InPlace;
        else if (cmp(argv[1], "--cuda-unroll=y")) cudaUnroll = True;
        else if (cmp(argv[1], "--cuda-unroll=n")) cudaUnroll = False;
        else if (cmp(argv[1], "--cuda-thread=y")) cudaThread = true;
        else if (cmp(argv[1], "--cuda-thread=n")) cudaThread = false;
        else if (cmp(argv[1], "--cuda-slices=y")) calcSlicesSeparately = true;
        else if (cmp(argv[1], "--cuda-slices=n")) calcSlicesSeparately = false;
        else if (cmp(argv[1], "--no-cogrid"))     cutoffGrid = false;
        else if (cmp(argv[1], "--no-nns"))        nns = false;

        // Error:
        else
        {
            fprintf(stderr, "%s: unknown switch '%s'\n", programName, argv[1]);
            throw ExitProgram(1);
        }

        --argc;
        ++argv;
    }
}

void ProgramParameters::cudaEnumDevicesAndExit()
{
#if defined(AG_CUDA)
    // Get a device count
    int deviceCount;
    CUDA_SAFE_CALL(cudaGetDeviceCount(&deviceCount));

    // Print a list of devices
    cudaDeviceProp prop;
    fprintf(stderr, "Found %i CUDA devices:\n\n"
                    " # | Name                           | MPs | Cap | GlobalMem | ConstMem\n"
                    "-----------------------------------------------------------------------\n", deviceCount);
    for (int i = 0; i < deviceCount; i++)
    {
        memset(&prop, 0, sizeof(prop));
        CUDA_SAFE_CALL(cudaGetDeviceProperties(&prop, i));
        fprintf(stderr, "%2i | %-31s|%4i |%2i.%-2i|%6lu MiB |%5lu KiB\n", i,
                        prop.name, prop.multiProcessorCount, prop.major, prop.minor,
                        prop.totalGlobalMem >> 20, prop.totalConstMem >> 10);
    }
#endif
    throw ExitProgram(0);
}

void ProgramParameters::readParamString(int *argc, char ***argv, char *out)
{
    strncpy(out, (*argv)[2], MAX_CHARS);
    ++*argv;
    --*argc;
}

int ProgramParameters::readParamInt(int *argc, char ***argv)
{
    int n;
    if (sscanf((*argv)[2], "%i", &n) != 1)
    {
        fprintf(stderr, "%s: '%s' is not a number\n", programName, (*argv)[2]);
        throw ExitProgram(1);
    }
    ++*argv;
    --*argc;
    return n;
}

void ProgramParameters::printHelpAndExit()
{
    fprintf(stderr, "Usage: %s [OPTIONS]\n"
                    "\n"
                    "Stream control:\n"
                    "  -l, --log FILE    log to FILE, default: standard output\n"
                    "  -p, --gpf FILE    read grid parameters from FILE, default: standard input\n"
                    "\n"
                    "Miscellaneous:\n"
                    "      --benchmark   print execution times to standard error output\n"
                    "      --cogrid-mem N  reserve at most N megabytes of memory for the cutoff grid,\n"
                    "                      default: 512\n"
                    "      --cuda-enum   enumerate all CUDA devices and exit\n"
                    "      --cuda-dev N  use a CUDA device number N, default: 0\n"
                    "  -d, --debug       increment debug level\n"
                    "  -u, --help        display this help and exit\n"
                    "      --no-cuda     disable CUDA, use the CPU codepath instead\n"
                    "      --omp N       set OpenMP to use N threads at most\n"
                    "      --timeout N   terminate if calculations does not finish in N seconds\n"
                    "                    (POSIX only), this must be the first parameter if used\n"
                    "      --v4          set the AutoGrid 4.0 default parameter library\n"
                    "\n"
                    "Advanced:\n"
                    "      --cuda-ddd=g|c|t|i  set a way of calculating distance-dependent dielectric\n"
                    "                          to one of the following: g=global memory, c=constant\n"
                    "                          memory, t=texture memory, i=in-place\n"
                    "                          default: either t or i based on grid dimensions\n"
                    "      --cuda-slices=y|n  calculate gridmap slices separately i.e. only one\n"
                    "                         gridmap slice per CUDA kernel call, default: n\n"
                    "      --cuda-thread=y|n  use a separate thread for CUDA, default: y\n"
                    "      --cuda-unroll=y|n  increase performance for large gridmaps,\n"
                    "                         default: based on grid dimensions\n"
                    "      --no-cogrid   disable the cutoff-grid optimization\n"
                    "      --no-nns      disable the nearest-neighbor-search optimization\n"
                    "\n"
                    "Compiled with OpenMP " OMP_STATUS ".\n"
                    "Compiled with CUDA " CUDA_STATUS ".\n"
                    "\n", programName);
    throw ExitProgram(0);
}

void ProgramParameters::setTimeout(int seconds)
{
#if !defined(_WIN32)
    pid_t pid = fork();
    if (pid)
    {
        int status;
        alarm(seconds);
        UNIX_SIGNAL_TRY(SIGALRM)
        {
            waitpid(pid, &status, 0);
        }
        UNIX_SIGNAL_CATCH()
        {
            kill(pid, SIGKILL);
            fprintf(stderr, "Timeout! Process terminated.\n");
        }
        UNIX_SIGNAL_END()
        throw ExitProgram(0);
    }
#endif
}

bool ProgramParameters::cmp(const char *s1, const char *s2)
{
    return strcmp(s1, s2) == 0;
}

bool ProgramParameters::cmp2(const char *s1, const char *s2, const char *s3)
{
    return cmp(s1, s2) || cmp(s1, s3);
}
