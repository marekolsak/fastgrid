#pragma once
#include <ctime>
#include <cstdio>
#include <sys/types.h>
#if !defined(_WIN32)
    #include <sys/times.h>
    #include <unistd.h>
#endif
#include "autogrid.h"
#include "structs.h"

// round() is a C99 function and not universally available
// Required to round %.3f consistently on different platforms
#if defined(HAVE_ROUND)
    #define round3dp(x) ((round((x)*1000.0L))/1000.0L)
#else
    #define round3dp(x) ((floor((x)*1000.0 + 0.5)) / 1000.0)
#endif

// WinAPI defines ERROR, we have to #undef it since we use the same name
#if defined(ERROR)
    #undef ERROR
#endif

// printError() is used with error_level where:
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
    #define CLOCKS_PER_SEC (getClocksPerSec())
#endif

// Define tokens for parsing AutoDock atomic parameter files
#define	PAR_		-1
#define	PAR_NULL	 0
#define	PAR_VDW 	 1
#define	PAR_HBOND	 2
#define	PAR_ESTAT	 3
#define	PAR_DESOLV	 4
#define	PAR_TORS	 5
#define	PAR_ATOM_PAR	 6
#define	PAR_COMMENT	 7

// GPF tokens
#define	GPF_NULL	0
#define	GPF_COMMENT	1
#define GPF_RECEPTOR	2
#define GPF_GRIDFLD	3
#define GPF_NPTS	4
#define GPF_SPACING	5
#define GPF_GRIDCENTER	7
#define GPF_LIGAND_TYPES	8
#define GPF_MAP		9
#define	GPF_NBP_COEFFS	10
#define GPF_NBP_R_EPS	11
#define GPF_ELECMAP	12
#define GPF_DIEL	13
#define GPF_FMAP	14
#define GPF_DISORDER	15
#define GPF_SMOOTH	16
#define GPF_SOL_PAR	17
#define GPF_CONSTANT	18
#define GPF_COVALENTMAP	19
#define GPF_RECEPTOR_TYPES	20
#define GPF_PARAM_FILE	21
#define GPF_DSOLVMAP	22
#define GPF_QASP	23

// tms structure
#if defined(_WIN32)
struct tms
{
    clock_t tms_utime;	// CPU time used in executing the instructions of the calling process.
    clock_t tms_stime;	// CPU time used by the system on behalf of the calling process.
    clock_t tms_cutime;	// sum of the tms_utime values and the tms_cutime values of all terminated child processes of the calling process, whose status has been reported to the parent process by wait or waitpid; see section Process Completion. In other words, it represents the total CPU time used in executing the instructions of all the terminated child processes of the calling process, excluding child processes which have not yet been reported by wait or waitpid.
    clock_t tms_cstime;	// similar to tms_cutime, but represents the total CPU time used by the system on behalf of all the terminated child processes of the calling process.
    // All of the times are given in clock ticks. These are absolute values; in a newly created process, they are all zero. See section Creating a Process.
};
#endif

// functions
void initBoinc();
FILE *openFile(const char *path, const char *mode);
char *getHostname(char *buffer, int size);
ParameterEntry *atomParameterManager_find(const char key[]);
void atomParameterManager_enter(const char key[], ParameterEntry value);
void printBanner(FILE *logFile, double versionNumber);
double calculateDDDMehlerSolmajer(double distance, double aprrox_zero);
int checkSize(int nelements, char axischar, const char *programname, FILE *logFile);
int getRecIndex(const char key[]);
int parseGPFLine(char line[LINE_LEN]);
int parseTypes(char * line, char *words[], int maxwords);
void printTimeFixed(float t, FILE *logFile);
char *getDate(int flag, char *buffer, int size);
void printError(const char *programname, FILE *fileptr, int error_level, char message[LINE_LEN]);
void printErrorFormatted(const char *programname, FILE *fileptr, int error_level, const char *format, ...);
int strIndex(char s[], char t[]);
void printClockTime(Clock duration, struct tms *start, struct tms *end, Real idct, FILE *logFile);
void printClockTimeInHMS(Clock duration, struct tms *start, struct tms *end, Real idct, FILE *logFile);

#if defined(_WIN32)
clock_t times(struct tms *buffer);
#else
int getClocksPerSec();
#endif

void beginTimer(const char *description);
void endTimer();
