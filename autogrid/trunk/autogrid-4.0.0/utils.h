#pragma once
#include <cstdio>

// BOINC
void boincInit();
void boincDone();
FILE *openFile(const char *path, const char *mode);

// Timer
void beginTimer(const char *description);
void endTimer();
