#pragma once
#include "logfile.h"

void readParameterLibrary(char FN_parameter_library[MAX_CHARS], int outlev, const char *programname, int debug, FILE *logFile, Linear_FE_Model &AD4);
void setupParameterLibrary(int outlev, const char *programname, int debug, FILE *logFile, Linear_FE_Model &AD4);
