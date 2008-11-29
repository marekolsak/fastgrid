#pragma once
#include <cstdio>

// round() is a C99 function and not universally available
// Required to round %.3f consistently on different platforms
#if defined(HAVE_ROUND)
    #define round3dp(x) ((round((x)*1000.0L))/1000.0L)
#else
    #define round3dp(x) ((floor((x)*1000.0 + 0.5)) / 1000.0)
#endif

// BOINC
void boincInit();
void boincDone();
FILE *openFile(const char *path, const char *mode);

// Timer
void beginTimer(const char *description);
void endTimer();
