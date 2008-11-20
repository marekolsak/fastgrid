#ifdef HAVE_CONFIG_H
    #include <config.h>
#endif
#include "utils.h"
#include "constants.h"
#include "exceptions.h"
#include <cstdlib>
#include <cmath>
#include <cerrno>
#include <cctype>
#include <cstring>

// initializes BOINC
void initBoinc()
{
#if defined(BOINC)
    boinc_init_diagnostics(BOINC_DIAG_DUMPCALLSTACKENABLED | BOINC_DIAG_HEAPCHECKENABLED | BOINC_DIAG_REDIRECTSTDERR | BOINC_DIAG_REDIRECTSTDOUT);

#if defined(BOINCCOMPOUND)
    BOINC_OPTIONS options;
    options.main_program = false;
    options.check_heartbeat = false;    // monitor does check heartbeat
    options.handle_trickle_ups = false;
    options.handle_trickle_downs = false;
    options.handle_process_control = false;
    options.send_status_msgs = true;    // only the worker programs (i.e. model) sends status msgs
    options.direct_process_action = true;   // monitor handles suspend/quit, but app/model doesn't

    // Initialization of Boinc
    int rc = boinc_init_options(options);   // return 0 for success
    if (rc)
    {
        fprintf(stderr, "BOINC_ERROR: boinc_init_options() failed \n");
        throw ExitProgram(rc);
    }

#else
    // All BOINC applications must initialize the BOINC interface:
    rc = boinc_init();
    if (rc)
    {
        fprintf(stderr, "BOINC_ERROR: boinc_init() failed.\n");
        throw ExitProgram(rc);
    }
#endif
#endif
}

// fopen rewrite to either use BOINC api or normal system call
FILE *openFile(const char *path, const char *mode)
{
    FILE *filep;

#if defined(BOINC)
    int rc;
    char resolved_name[512];

    rc = boinc_resolve_filename(path, resolved_name, sizeof(resolved_name));
    if (rc)
    {
        fprintf(stderr, "BOINC_ERROR: cannot open filename.%s\n", path);
        boinc_finish(rc);       // back to BOINC core
    }
    // Then open the file with boinc_fopen() not just fopen()
    filep = boinc_fopen(resolved_name, mode);
#else
    filep = fopen(path, mode);
#endif
    return filep;
}

// Atom Parameter Manager (hash, atomParameterManager_enter, atomParameterManager_find)

#define MAXKEY (256*256)

typedef ParameterEntry PE;
static PE *dictionary[MAXKEY];

static unsigned int hash(const char key[]) {
    switch (strlen(key)) {
        case 0: return 0;
        case 1: return (unsigned int)key[0];
        default: return (unsigned int)key[0] + 256*(unsigned int)key[1];
    }
}

void atomParameterManager_enter(const char key[], PE value) {
    if (dictionary[hash(key)] == 0) {
        dictionary[hash(key)] = (PE *) calloc(1, sizeof(PE));
    }
    *(dictionary[hash(key)]) = value;  // this replaces, as well as inserts
    return;
}

PE * atomParameterManager_find(const char key[]) {
    return dictionary[hash(key)];
}

double calculateDDDMehlerSolmajer(double distance, double approx_zero) {
    /*____________________________________________________________________________
     * Distance-dependent dielectric ewds: Mehler and Solmajer, Prot Eng 4, 903-910.
     *____________________________________________________________________________*/
    double epsilon = 1.0L;
    double lambda = 0.003627L;
    double epsilon0 = 78.4L;
    double A = -8.5525L;
    double B;
    B = epsilon0 - A;
    double rk= 7.7839L;
    double lambda_B;
    lambda_B = -lambda * B;

    epsilon = A + B / (1 + rk*exp(lambda_B * distance));

    if (epsilon < approx_zero) {
        epsilon = 1.0L;
    }
    return epsilon;
}

/******************************************************************************/
/*      Name: checkSize                                                       */
/*  Function: Checks that number of grid elements is valid.                   */
/* Copyright: (C) 1994, TSRI                                                  */
/*----------------------------------------------------------------------------*/
/*    Author: Garrett Morris, The Scripps Research Institute                  */
/*      Date: 13/07/92                                                        */
/*----------------------------------------------------------------------------*/
/*    Inputs: nelements, axischar                                             */
/*   Returns: nelements                                                       */
/*   Globals: MAX_GRID_PTS                                                    */
/*----------------------------------------------------------------------------*/
/* Modification Record                                                        */
/* Date     Inits   Comments                                                  */
/* 04/01/93 GMM     Created for use in makefile.                              */
/******************************************************************************/
int checkSize(int nelements, char axischar, LogFile &logFile)
{
    // nelements mustn't be negative, shouldn't be zero or larger than MAX_GRID_PTS and should be even
    if (nelements < 0)
        logFile.printErrorFormatted(FATAL_ERROR, "Negative number of %c-grid elements!  Aborting.\n\n", axischar);
    else if (nelements == 0)
        logFile.printErrorFormatted(WARNING, "0 %c-grid elements!\n\n", axischar);
    else if (nelements > MAX_GRID_PTS)
    {
        logFile.printErrorFormatted(WARNING, "Maximum number of %c-grid elements allowed is %d. Using this value.\n", axischar, MAX_GRID_PTS);
        nelements = MAX_GRID_PTS;
    }
    else if (nelements % 2 == 1)
    {
        logFile.printTitledFormatted("Number of grid elements must be even; %c-elements changed to: %d\n", axischar, nelements);
        nelements -= 1;
    }

    return nelements;
}

int getRecIndex(const char key[])
{
    ParameterEntry *foundParam;

    foundParam = atomParameterManager_find(key);
    if (foundParam != 0)
        return foundParam->recIndex;
    return -1;
}

/******************************************************************************/
/*      Name: parseGPFLine                                                       */
/*  Function: Parse the AutoGrid parameter file line                          */
/* Copyright: (C) 1995, TSRI                                                  */
/*----------------------------------------------------------------------------*/
/*    Author: Garrett Morris, The Scripps Research Institute                  */
/*      Date: 02/01/95 (1-feb-1995)                                           */
/*----------------------------------------------------------------------------*/
/*    Inputs: line                                                            */
/*   Returns: integer token describing the keyword found.                     */
/*   Globals: none.                                                           */
/*----------------------------------------------------------------------------*/
/* Modification Record                                                        */
/* Date     Inits   Comments                                                  */
/* 02/01/95 GMM     Entered code.                                             */
/******************************************************************************/

int parseGPFLine(char line[LINE_LEN])
{
    int l, i, token = -1 ;	       /* return -1 if nothing is recognized. */
    char c[LINE_LEN];

    l = (int)strIndex(line, " ");
    if (l == -1) {
        l = (int)strIndex(line, "\t");
        if (l == -1) {
            l = (int)strlen(line);
	}
    }
    for(i=0; i<l; i++) {
        c[i] = (char)tolower((int)line[i]);
    }

    if ((c[0]=='\n')||(c[0]=='\0')) {
        token = GPF_NULL;

    } else if (c[0]=='#') {
        token = GPF_COMMENT;

    } else if (equal(c,"receptor_types",14)) {
        token = GPF_RECEPTOR_TYPES;

    } else if (equal(c,"receptor",8)) {
        token = GPF_RECEPTOR;

    } else if (equal(c,"gridfld",7)) {
        token = GPF_GRIDFLD;

    } else if (equal(c,"npts",4)) {
        token = GPF_NPTS;

    } else if (equal(c,"spacing",7)) {
        token = GPF_SPACING;

    } else if (equal(c,"gridcenter",10)) {
        token = GPF_GRIDCENTER;

    } else if (equal(c,"types",5)) {
        token = GPF_LIGAND_TYPES;

    } else if (equal(c,"ligand_types",12)) {
        token = GPF_LIGAND_TYPES;


    } else if (equal(c,"map",3)) {
        token = GPF_MAP;

    } else if (equal(c,"elecmap",7)) {
        token = GPF_ELECMAP;

    } else if (equal(c,"dsolvmap",8)) {
        token = GPF_DSOLVMAP;

    } else if (equal(c,"covalentmap",11)) {
        token = GPF_COVALENTMAP;

    } else if (equal(c,"nbp_coeffs",10)) {
        token = GPF_NBP_COEFFS;

    } else if (equal(c,"nbp_r_eps",9)) {
        token = GPF_NBP_R_EPS;

    } else if (equal(c,"dielectric",10)) {
        token = GPF_DIEL;

    } else if (equal(c,"qasp",4)) {
        token = GPF_QASP;

    } else if (equal(c,"fmap",4)) {
        token = GPF_FMAP;

    } else if (equal(c,"disorder_h",10)) {
        token = GPF_DISORDER;

    } else if (equal(c,"smooth",6)) {
        token = GPF_SMOOTH;

    } else if (equal(c,"sol_par",7)) {
        token = GPF_SOL_PAR;

    } else if (equal(c,"constant",8)) {
        token = GPF_CONSTANT;

    } else if (equal(c,"parameter_file",14)) {
        token = GPF_PARAM_FILE;

    }
    return(token);
}

/******************************************************************************/
/*      Name: parseTypes                                                      */
/*  Function: Parse the AutoGrid types lines                                  */
/* Copyright: (C) 1995, TSRI                                                  */
/*----------------------------------------------------------------------------*/
/*    Author: Garrett Morris, The Scripps Research Institute                  */
/*      Date: 02/01/95 (1-feb-1995)                                           */
/*----------------------------------------------------------------------------*/
/*    Inputs: line, array of pointers, cut-off number of words                */
/*   Returns: integer, number of types found.                                 */
/*   Globals: none.                                                           */
/*----------------------------------------------------------------------------*/
/* Modification Record                                                        */
/* Date     Inits   Comments                                                  */
/* 06/02/03 RH      Entered code.                                             */
/******************************************************************************/
int parseTypes(char * line, char *words[], int maxwords)
/*utility func for parsing types*/
{
    char *char_ptr = line;
    int num_types = 0;
    /*flag for first word which is always a keyword*/
    int found_keyword = 0;
    int index = 0;

    for(;;) {
        /*skip spaces*/
        while(isspace(*char_ptr)){
            char_ptr++;
            index++;
        };
        /*done parsing when get eol 'null' character*/
        /* could get null after a space*/
        if (*char_ptr == '\0'){
            /*return number of 'types' found*/
            return num_types;
        };
        /* the first word is the keyword not a type*/
        if(found_keyword==0){
            found_keyword++;
        } else {
            /*words is a list of indicies of beginning of 1 or 2 char types*/
            words[num_types++] = char_ptr;
        };
        /*once in a type, skip possible 2nd characters up to a space or null
         * character*/
        while(!isspace(*char_ptr) && *char_ptr!='\0'){
            char_ptr++;
            index++;
        };
        /*done parsing when get eol 'null' character*/
        /* could get null after a character*/
        if(*char_ptr=='\0'){
            return num_types;
        };
        /*make each 'type' a null terminated string*/
        *char_ptr++ = '\0';
        index++;
        /*if there are too many types, return*/
        if(num_types >=maxwords){
            return num_types;
        };
    }
}

// returns index of t in s, -1 if none.
int strIndex(char s[], char t[])
{
    char *r = strstr(s, t);
    return r? int(r-s) : -1;
}

#if defined(_WIN32)


#endif

static clock_t timers[256];
static unsigned int indexOfNesting = 0;

void beginTimer(const char *description)
{
    if (indexOfNesting > 255)
    {
        fprintf(stderr, "ERROR: Cannot initiate a timer\n");
        throw ExitProgram(1);
    }

    for (unsigned int i = 0; i < indexOfNesting; i++)
        fprintf(stderr, "  ");
    fprintf(stderr, "\"%s\" {\n", description);
    tms _t;
    timers[indexOfNesting] = times(&_t);
    ++indexOfNesting;
}

void endTimer()
{
    if (indexOfNesting <= 0)
    {
        fprintf(stderr, "ERROR: Cannot terminate a timer\n");
        throw ExitProgram(1);
    }

    --indexOfNesting;
    tms _t;
    clock_t time = times(&_t) - timers[indexOfNesting];
    for (unsigned int i = 0; i < indexOfNesting; i++)
        fprintf(stderr, "  ");
    fprintf(stderr, "} took %i ms.\n", (time*1000/getClocksPerSec()));
}

// Dummy graphics API entry points.  This app does not do graphics, but it still must provide these callbacks.
#if defined(BOINC)
void app_graphics_render(int xs, int ys, double time_of_day){}
void app_graphics_reread_prefs(){}
void boinc_app_mouse_move(int x, int y, bool left, bool middle, bool right){}
void boinc_app_mouse_button(int x, int y, int which, bool is_down){}
void boinc_app_key_press(int wParam, int lParam){}
void boinc_app_key_release(int wParam, int lParam){}
#endif
