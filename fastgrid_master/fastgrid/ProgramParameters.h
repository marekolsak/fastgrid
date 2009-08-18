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
#include "autogrid.h"
#include "electrostatics/CudaUtils.h"

enum ExtendedBool
{
    False = 0,
    True = 1,
    Unassigned = -1
};

class ProgramParameters
{
public:
    ProgramParameters(int argc, char **argv);
    const char *getProgramName() const              { return programName; }
    const char *getGridParameterFilename() const    { return gridParameterFilename; }
    const char *getLogFilename() const              { return logFilename; }
    int getDebugLevel() const                       { return debug; }
    bool benchmarkEnabled() const                   { return benchmark; }
    bool useNNS() const                             { return nns; }
    bool useCutoffGrid() const                      { return cutoffGrid; }
    bool useCUDA() const                            { return cuda; }
    bool useCUDAThread() const                      { return cudaThread; }
    bool calcSlicesSeparatelyCUDA() const           { return calcSlicesSeparately; }
    ExtendedBool unrollLoopCUDA() const             { return cudaUnroll; }
    int getDeviceIDCUDA() const                     { return deviceID; }
    DielectricKind getDDDKindCUDA() const           { return cudaDDDKind; }
    bool useVersion4() const                        { return v4; }
    int getCutoffGridMemoryLimit() const            { return cutoffGridMem; }

private:
    char programName[MAX_CHARS];
    char gridParameterFilename[MAX_CHARS];
    char logFilename[MAX_CHARS];
    int debug, deviceID, cutoffGridMem;
    bool benchmark, nns, cutoffGrid, cuda, cudaThread, calcSlicesSeparately, v4;
    ExtendedBool cudaUnroll;
    DielectricKind cudaDDDKind;

    void parse(int argc, char **argv);
    void cudaEnumDevicesAndExit();
    void readParamString(int *argc, char ***argv, char *out);
    int readParamInt(int *argc, char ***argv);
    void printHelpAndExit();
    void setTimeout(int seconds);
    bool cmp(const char *s1, const char *s2);
    bool cmp2(const char *s1, const char *s2, const char *s3);
};
