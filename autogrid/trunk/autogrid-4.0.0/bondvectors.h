#pragma once
#include "inputdata.h"

// This class is dedicated to calculating bond vectors for directional H-bonds
class BondVectors
{
public:
    bool disorder[AG_MAX_ATOMS];
    int rexp[AG_MAX_ATOMS];
    double rvector[AG_MAX_ATOMS][XYZ];
    double rvector2[AG_MAX_ATOMS][XYZ];

    BondVectors(LogFile *logFile);
    void calculate(const InputData *input, ParameterLibrary &parameterLibrary);

private:
    LogFile *logFile;
};
