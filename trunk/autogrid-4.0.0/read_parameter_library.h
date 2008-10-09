#pragma once
#include <cstdio>

void read_parameter_library(char FN_parameter_library[MAX_CHARS], int outlev, const char *programname, int debug, FILE *logFile, Linear_FE_Model &AD4);
void setup_parameter_library(int outlev, const char *programname, int debug, FILE *logFile, Linear_FE_Model &AD4);
