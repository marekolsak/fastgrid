#pragma once
#include <cstdio>
#if !defined(_WIN32)
    #include <unistd.h>
#endif
#include "autogrid.h"
#include "logfile.h"

// round() is a C99 function and not universally available
// Required to round %.3f consistently on different platforms
#if defined(HAVE_ROUND)
    #define round3dp(x) ((round((x)*1000.0L))/1000.0L)
#else
    #define round3dp(x) ((floor((x)*1000.0 + 0.5)) / 1000.0)
#endif

// GPF tokens
// TODO: move this to the InputData class
enum GPFTokens
{
    GPF_NULL = 0,
    GPF_COMMENT,
    GPF_RECEPTOR,
    GPF_GRIDFLD,
    GPF_NPTS,
    GPF_SPACING,
    GPF_GRIDCENTER,
    GPF_LIGAND_TYPES,
    GPF_MAP,
    GPF_NBP_COEFFS,
    GPF_NBP_R_EPS,
    GPF_ELECMAP,
    GPF_DIEL,
    GPF_FMAP,
    GPF_DISORDER,
    GPF_SMOOTH,
    GPF_SOL_PAR,
    GPF_CONSTANT,
    GPF_COVALENTMAP,
    GPF_RECEPTOR_TYPES,
    GPF_PARAM_FILE,
    GPF_DSOLVMAP,
    GPF_QASP
};

// functions
void initBoinc();
FILE *openFile(const char *path, const char *mode);
double calculateDDDMehlerSolmajer(double distance, double aprrox_zero);
int checkSize(int nelements, char axischar, LogFile &logFile);
int parseGPFLine(char line[LINE_LEN]);
int parseTypes(char * line, char *words[], int maxwords);
int strIndex(char s[], char t[]);

void beginTimer(const char *description);
void endTimer();
