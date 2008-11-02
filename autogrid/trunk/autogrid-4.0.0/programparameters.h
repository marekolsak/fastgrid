#pragma once
#include "autogrid.h"

class ProgramParameters
{
public:
    ProgramParameters(int argc, char **argv);
    const char *getProgramName() const;
    const char *getGridParameterFilename() const;
    const char *getLogFilename() const;
    int getDebugLevel() const;

private:
    char programName[MAX_CHARS];
    char gridParameterFilename[MAX_CHARS];
    char logFilename[MAX_CHARS];
    int debug;

    int Parse(int argc, char **argv);
};
