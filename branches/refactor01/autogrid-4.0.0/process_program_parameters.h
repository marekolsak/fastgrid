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

int process_program_parameters(int argc, char **argv, ProgramParameters &out);
