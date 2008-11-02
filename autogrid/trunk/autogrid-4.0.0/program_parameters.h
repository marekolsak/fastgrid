#pragma once
#include "autogrid.h"

struct ProgramParameters
{
    char programName[MAX_CHARS];
    char gridParameterFilename[MAX_CHARS];
    char logFilename[MAX_CHARS];
    int debug;

    ProgramParameters();
};

int processProgramParameters(int argc, char **argv, ProgramParameters &out);
