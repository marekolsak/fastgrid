#pragma once

// round() is a C99 function and not universally available
// Required to round %.3f consistently on different platforms
#if defined(HAVE_ROUND)
    #define round3dp(x) ((round((x)*1000.0L))/1000.0L)
#else
    #define round3dp(x) ((floor((x)*1000.0 + 0.5)) / 1000.0)
#endif

// print_error() is used with error_level where:
// error_level = one of the following:
#define FATAL_ERROR -2
#define ERROR -1
#define WARNING  0
#define INFORMATION 1
#define SUGGESTION 2

#if !defined(_WIN32)
    #if defined(CLOCKS_PER_SEC)
        #undef CLOCKS_PER_SEC
    #endif
    #define CLOCKS_PER_SEC (sysconf(_SC_CLK_TCK))
#endif

void print_error(const char *programname, FILE * fileptr, int error_level, char message[LINE_LEN]);
FILE *ag_fopen(const char *path, const char *mode);
int get_rec_index(const char key[]);
