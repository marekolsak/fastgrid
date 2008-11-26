#pragma once
#include "logfile.h"

void readParameterLibrary(char FN_parameter_library[MAX_CHARS], int outputLevel, const char *programname, int debug, FILE *logFile,
                          LinearFreeEnergyModel &model, AtomParameterManager &apm);
void setupParameterLibrary(int outputLevel, const char *programname, int debug, FILE *logFile,
                           LinearFreeEnergyModel &model, AtomParameterManager &apm);
