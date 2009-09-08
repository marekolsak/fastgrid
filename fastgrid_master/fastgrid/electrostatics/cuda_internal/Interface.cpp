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
#include <cstring>
#include "StandardKernels.h"
#include "DDDConstMemKernels.h"
#include "../../Exceptions.h"

void getCudaInternalAPI(DielectricKind dddKind, CudaInternalAPI &api)
{
    if (dddKind != DistanceDependentDiel_ConstMem)
    {
        api.numAtomsPerKernel           = STD_NUM_ATOMS_PER_KERNEL;
        api.setDistDepDielTexture       = stdSetDistDepDielTexture;
        api.setDistDepDielLookUpTable   = stdSetDistDepDielLookUpTable;
        api.setGridMap                  = stdSetGridMap;
        api.setSlice                    = stdSetSlice;
        api.setAtoms                    = stdSetAtoms;
        api.getKernelProc               = stdGetKernelProc;
        api.callKernel                  = stdCallKernel;
    }
    else
    {
        api.numAtomsPerKernel           = DDDCM_NUM_ATOMS_PER_KERNEL;
        api.setDistDepDielTexture       = dddcmSetDistDepDielTexture;
        api.setDistDepDielLookUpTable   = dddcmSetDistDepDielLookUpTable;
        api.setGridMap                  = dddcmSetGridMap;
        api.setSlice                    = dddcmSetSlice;
        api.setAtoms                    = dddcmSetAtoms;
        api.getKernelProc               = dddcmGetKernelProc;
        api.callKernel                  = dddcmCallKernel;
    }
}

void checkCudaError(cudaError e, const char *file, int line, const char *func, const char *code)
{
    if (e != cudaSuccess)
    {
        const char *str = cudaGetErrorString(e);
        fprintf(stderr, "CUDA error: '%s'\n"
                        "        in file '%s'\n"
                        "        in line %i\n"
                        "        in function '%s'\n"
                        "        in code '%s'\n", str, file, line, func, code);
        if (strstr(str, "launch fail"))
            throw ExitProgram(0xbad);
    }
}
