#pragma once
#include <cstdio>

int process_program_parameters(int argc, char **argv, FILE *&GPF, FILE *&logFile, char *&programname, char (&grid_param_fn)[MAX_CHARS], int &debug, int &oldpdbq);
