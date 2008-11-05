#pragma once
#include <cstdio>
#if !defined(_WIN32)
    #include <unistd.h>
#endif
#include "autogrid.h"
#include "structs.h"
#include "logfile.h"

// round() is a C99 function and not universally available
// Required to round %.3f consistently on different platforms
#if defined(HAVE_ROUND)
    #define round3dp(x) ((round((x)*1000.0L))/1000.0L)
#else
    #define round3dp(x) ((floor((x)*1000.0 + 0.5)) / 1000.0)
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

// functions
void initBoinc();
FILE *openFile(const char *path, const char *mode);
ParameterEntry *atomParameterManager_find(const char key[]);
void atomParameterManager_enter(const char key[], ParameterEntry value);
double calculateDDDMehlerSolmajer(double distance, double aprrox_zero);
int checkSize(int nelements, char axischar, LogFile &logFile);
int getRecIndex(const char key[]);
int parseGPFLine(char line[LINE_LEN]);
int parseTypes(char * line, char *words[], int maxwords);
int strIndex(char s[], char t[]);

void beginTimer(const char *description);
void endTimer();
