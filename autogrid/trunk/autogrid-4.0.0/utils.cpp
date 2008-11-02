#ifdef HAVE_CONFIG_H
    #include <config.h>
#endif
#if defined(_WIN32)
    #define WIN32_LEAN_AND_MEAN
    #include <windows.h>
    #include <Winsock2.h>
#endif
#include "utils.h"
#include <cstdlib>
#include <cmath>
#include <cerrno>
#include "autogrid.h"
#include "constants.h"
#include "exceptions.h"

// initializes BOINC
void ag_boinc_init()
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
        exit(rc);
    }

#else
    // All BOINC applications must initialize the BOINC interface:
    rc = boinc_init();
    if (rc)
    {
        fprintf(stderr, "BOINC_ERROR: boinc_init() failed.\n");
        exit(rc);
    }
#endif
#endif
}

// fopen rewrite to either use BOINC api or normal system call
FILE *ag_fopen(const char *path, const char *mode)
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

char *ag_gethostname(char *buffer, int size)
{
#if defined(_WIN32)
    WSADATA wsaData;
    WSAStartup(MAKEWORD(2, 2), &wsaData);
#endif

    if (gethostname(buffer, size) != 0)
        strncpy(buffer, "(gethostname returned error)", size);

#if defined(_WIN32)
    WSACleanup();
#endif
    return buffer;
}

// Atom Parameter Manager (hash, apm_enter, apm_find)

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

void apm_enter(const char key[], PE value) {
    if (dictionary[hash(key)] == 0) {
        dictionary[hash(key)] = (PE *) calloc(1, sizeof(PE));
    }
    *(dictionary[hash(key)]) = value;  // this replaces, as well as inserts
    return;
}

PE * apm_find(const char key[]) {
    return dictionary[hash(key)];
}

// Output banner...
void fprint_banner(FILE *logFile, double versionNumber)
{
    fprintf(logFile,"\n       _______________________________________________________\n");
    fprintf(logFile,"\n");
    fprintf(logFile,"__________//____________________________/////_________________/________\n");
    fprintf(logFile,"_________/__/____________/_____________/______________/_______/________\n");
    fprintf(logFile,"________/____/___________/_____________/______________________/________\n");
    fprintf(logFile,"________/____/__/_____/_/////___/////__/__////_/_///__/__////_/________\n");
    fprintf(logFile,"_______/______/_/_____/__/_____/_____/_/_____/_//___/_/_/____//________\n");
    fprintf(logFile,"_______////////_/_____/__/_____/_____/_/_____/_/______/_/_____/________\n");
    fprintf(logFile,"_______/______/_/____//__/___/_/_____/_/_____/_/______/_/____//________\n");
    fprintf(logFile,"_______/______/__////_/___///___/////___/////__/______/__////_/________\n");
    fprintf(logFile,"\n");
    fprintf(logFile,"       _______________________________________________________\n");
    fprintf(logFile,"\n");
    fprintf(logFile,"                                ______\n");
    fprintf(logFile,"                               /      \\\n");
    fprintf(logFile,"                              /        \\\n");
    fprintf(logFile,"                             /          \\\n");
    fprintf(logFile,"                             \\    /\\    /\n");
    fprintf(logFile,"                              \\  /  \\  /\n");
    fprintf(logFile,"                               \\/ /\\ \\/\n");
    fprintf(logFile,"                                 /  \\\n");
    fprintf(logFile,"                                /____\\\n");
    fprintf(logFile,"\n");
    fprintf(logFile,"\n");
    fprintf(logFile,"                ______________________________________ \n");
    fprintf(logFile,"               |                                      |\n");
    fprintf(logFile,"               |            AutoGrid %3.2lf             |\n",versionNumber);
    fprintf(logFile,"               |                                      |\n");
    fprintf(logFile,"               |        Garrett M. Morris, TSRI       |\n");
    fprintf(logFile,"               |            Ruth Huey, TSRI           |\n");
    fprintf(logFile,"               |        David S. Goodsell, TSRI       |\n");
    fprintf(logFile,"               |         Arthur J. Olson, TSRI        |\n");
    fprintf(logFile,"               |                                      |\n");
    fprintf(logFile,"               |        (c) 1989-2005, TSRI           |\n");
    fprintf(logFile,"               |   The Scripps Research Institute     |\n");
    fprintf(logFile,"               |______________________________________|\n");
    fprintf(logFile,"\n");
    fprintf(logFile,"                ______________________________________ \n");
    fprintf(logFile,"               |                                      |\n");
    fprintf(logFile,"               | Calculation of van der Waals, H-Bond,|\n");
    fprintf(logFile,"               |   Electrostatic Potential Energy, &  |\n");
    fprintf(logFile,"               |   Desolvation Free Energy Grid Maps  |\n");
    fprintf(logFile,"               |             for AutoDock             |\n");
    fprintf(logFile,"               |______________________________________|\n");
    fprintf(logFile,"\n");
    fprintf(logFile,"\n");
    fprintf(logFile,"\n");
    fprintf(logFile,"\n");
}

double calc_ddd_Mehler_Solmajer(double distance, double approx_zero) {
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
/*      Name: check_size                                                      */
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
int check_size(int nelements, char axischar, const char *programname, FILE *logFile)
{
    int oldnelements;

    if (nelements < 0) {
        fprintf(stderr, "\n%s: Error! Negative number of %c-grid elements!  Aborting.\n\n", programname, axischar);
        exit(-2);
    }
    if (nelements == 0) {
        fprintf(stderr, "\n%s: Warning! 0 %c-grid elements!\n\n", programname, axischar);
    }
    if (nelements>MAX_GRID_PTS) {
        fprintf(logFile, "%s: Warning! Maximum number of %c-grid elements allowed is %d. Using this value.\n", programname, axischar, MAX_GRID_PTS);
        nelements = MAX_GRID_PTS;
    }
    oldnelements = nelements;
    nelements = (int) ((nelements/2) * 2); // N.B.: integer divide truncates remainder.
    if (oldnelements != nelements)
        fprintf(logFile, "%s: Number of grid elements must be even; %c-elements changed to: %d\n", programname, axischar, nelements);

    return nelements;
}

#if !defined(WIN32)
int get_clocks_per_sec()
{
    long clocks = sysconf(_SC_CLK_TCK);
    if (clocks < 0)
    {
        fprintf(stderr, "\"sysconf(_SC_CLK_TCK)\" failed in get_clocks_per_sec()\n");
        exit(-1);
    }
    return clocks;
}
#endif

int get_rec_index(const char key[])
{
    ParameterEntry *found_parm;

    found_parm = apm_find(key);
    if (found_parm != 0)
        return found_parm->rec_index;
    return -1;
}

/******************************************************************************/
/*      Name: gpfparser                                                       */
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

int gpfparser(char line[LINE_LEN])
{
    int l, i, token = -1 ;	       /* return -1 if nothing is recognized. */
    char c[LINE_LEN];

    l = (int)strindex(line, " ");
    if (l == -1) {
        l = (int)strindex(line, "\t");
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
/*      Name: parsetypes                                                      */
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
int parsetypes(char * line, char *words[], int maxwords)
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

void prHMSfixed(float t, FILE *logFile)
{
    int   h, m;
    float T, s;
    float hrs = 3600., min = 60.;

    h = (int) (t/hrs);
    T = t - h*hrs;
    m = (int) (T/min);
    s = T - m*min;

    if (h == 0) {
        if (m == 0)
            fprintf(logFile,    "        %5.2fs",        s);
        else
            fprintf(logFile,   "    %2dm %05.2fs",    m, s);
    } else {
            fprintf(logFile, "%2dh %02dm %05.2fs", h, m, s);
    }
}

char *getdate(int flag, char *buffer, int size)
{
    time_t tn; /* tn = "time_now" */
    char *StringTimeDate;
    struct tm *ts;

    tn = time(&tn);
    ts = localtime(&tn);

    if (flag==1) {
        _snprintf(buffer, size, "%d:%02d %02d\" %s, %02d/%02d/%4d\n",
        ((ts->tm_hour >  12) ? (ts->tm_hour-12) : ts->tm_hour), ts->tm_min, ts->tm_sec,
        ((ts->tm_hour >= 12) ? "p.m." : "a.m."),
        (ts->tm_mon + 1), ts->tm_mday, 1900+ts->tm_year);
    } else if (flag==2) {
          StringTimeDate = ctime(&tn);
          _snprintf(buffer, size, "%s", StringTimeDate);
    } else {
        _snprintf(buffer, size, "%d:%02d %02d\" %s\n",
        ((ts->tm_hour >  12) ? (ts->tm_hour-12) : ts->tm_hour), ts->tm_min, ts->tm_sec,
        ((ts->tm_hour >= 12) ? "pm" : "am"));
    }
    return buffer;
}

void printhms(Real t, FILE *logFile)
{
    int   h,
          m;
    Real T, s;
    Real min = 60.,
	  hrs = 3600.;

    h = (int)(t/hrs);
    T = t - ((Real)h)*hrs;
    m = (int)(T/min);
    s = T - ((Real)m)*min;

    if (h == 0) {
        if (m == 0)
            fprintf(logFile,       "%.2fs",       s);
        else
            fprintf(logFile,    "%dm %05.2fs",    m, s);
    } else {
            fprintf(logFile, "%dh %02dm %05.2fs", h, m, s);
    }
}

// print an error or informational message to a file-pointer or
// standard error
void print_error(const char *programname, FILE * fileptr, int error_level, char message[LINE_LEN])
{
    char output_message[LINE_LEN];
    char tag[LINE_LEN];

    switch (error_level)
    {
    case ERROR:
    case FATAL_ERROR:
        strncpy(tag, "ERROR", LINE_LEN);
        break;
    case WARNING:
        strncpy(tag, "WARNING", LINE_LEN);
        break;
    case INFORMATION:
        strncpy(tag, "INFORMATION", LINE_LEN);
        break;
    case SUGGESTION:
        strncpy(tag, "SUGGESTION", LINE_LEN);
        break;
    }

    _snprintf(output_message, LINE_LEN, "\n%s: %s:  %s\n", programname, tag, message);

    // Records all messages in the logFile.
    fprintf(fileptr, "%s\n", output_message);

    // Only send errors, fatal errors and warnings to standard error, stderr.
    switch (error_level)
    {
    case ERROR:
    case FATAL_ERROR:
    case WARNING:
        fprintf(stderr, "%s\n", output_message);
        break;
    }

    // If this is a fatal error, exit now.
    if (error_level == FATAL_ERROR)
        exit(error_level);
}

void print_errorf(const char *programname, FILE *fileptr, int error_level, const char *format, ...)
{
    char message[LINE_LEN];
    va_list ap;
    va_start(ap, format);
    vsprintf(message, format, ap);
    va_end(ap);
    print_error(programname, fileptr, error_level, message);
}

// returns index of t in s, -1 if none.
int strindex(char s[], char t[])
{
    char *r = strstr(s, t);
    return r? int(r-s) : -1;
}

void timesys(Clock duration, struct tms *start, struct tms *end, Real idct, FILE *logFile)
{
	fprintf(logFile, "Real= %.2f,  CPU= %.2f,  System= %.2f\n",     (Real)duration * idct,
                         (Real)(end->tms_utime  - start->tms_utime) * idct,
                         (Real)(end->tms_stime  - start->tms_stime) * idct);
}

void timesyshms(Clock duration, struct tms  *start, struct tms  *end, Real idct, FILE *logFile)
{
    int h, m;
    Real t, T, s;
#ifndef USE_DOUBLE
    const Real min = 60., hrs = 3600.;
#else
    const Real min = 60.L, hrs = 3600.L;
#endif

    fprintf(logFile, "Real= ");
    t = (Real)duration * idct;
    h = (int)(t/hrs);
    T = t - ((Real)h)*hrs;
    m = (int)(T/min);
    s = T - ((Real)m)*min;
    if (h == 0) {
        if (m == 0)
#ifndef USE_DOUBLE
            fprintf(logFile,       "%.2fs",       s);
#else
            fprintf(logFile,       "%.2lfs",       s);
#endif
        else
#ifndef USE_DOUBLE
            fprintf(logFile,    "%dm %05.2fs",    m, s);
#else
            fprintf(logFile,    "%dm %05.2lfs",    m, s);
#endif
    } else {
#ifndef USE_DOUBLE
            fprintf(logFile, "%dh %02dm %05.2fs", h, m, s);
#else
            fprintf(logFile, "%dh %02dm %05.2lfs", h, m, s);
#endif
    }

    fprintf(logFile, ",  CPU= ");
    t =      (Real)((end->tms_utime  - start->tms_utime) * idct);
    h = (int)(t/hrs);
    T = t - ((Real)h)*hrs;
    m = (int)(T/min);
    s = T - ((Real)m)*min;
    if (h == 0) {
        if (m == 0)
#ifndef USE_DOUBLE
            fprintf(logFile,       "%.2fs",       s);
#else
            fprintf(logFile,       "%.2lfs",       s);
#endif
        else
#ifndef USE_DOUBLE
            fprintf(logFile,    "%dm %05.2fs",    m, s);
#else
            fprintf(logFile,    "%dm %05.2lfs",    m, s);
#endif
    } else {
#ifndef USE_DOUBLE
            fprintf(logFile, "%dh %02dm %05.2fs", h, m, s);
#else
            fprintf(logFile, "%dh %02dm %05.2lfs", h, m, s);
#endif
    }

    fprintf(logFile, ",  System= ");
    t = (Real)((end->tms_stime  - start->tms_stime) * idct);
    h = (int)(t/hrs);
    T = t - ((Real)h)*hrs;
    m = (int)(T/min);
    s = T - ((Real)m)*min;
    if (h == 0) {
        if (m == 0)
#ifndef USE_DOUBLE
            fprintf(logFile,       "%.2fs",       s);
#else
            fprintf(logFile,       "%.2lfs",       s);
#endif
        else
#ifndef USE_DOUBLE
            fprintf(logFile,    "%dm %05.2fs",    m, s);
#else
            fprintf(logFile,    "%dm %05.2lfs",    m, s);
#endif
    } else {
#ifndef USE_DOUBLE
            fprintf(logFile, "%dh %02dm %05.2fs", h, m, s);
#else
            fprintf(logFile, "%dh %02dm %05.2lfs", h, m, s);
#endif
    }

    fprintf(logFile, "\n");
}

#if defined(_WIN32)

// converts WinAPI's FILETIME to clock_t
static clock_t FileTimeToClockTime(unsigned long long fileTime)
{
    // fileTime contains the time in 100s of nanoseconds
    return clock_t((fileTime * CLOCKS_PER_SEC) / 10000000ull);
}

// there is no times(..) function on Windows so we have to implement it on our
// own
clock_t times(struct tms *buffer)
{
    if (!buffer)
    {
        _set_errno(EFAULT);
        return clock_t(-1);
    }

    unsigned long long creationTime, exitTime, kernelTime, userTime;
    GetProcessTimes(GetCurrentProcess(),
                    reinterpret_cast<FILETIME*>(&creationTime),
                    reinterpret_cast<FILETIME*>(&exitTime),
                    reinterpret_cast<FILETIME*>(&kernelTime),
                    reinterpret_cast<FILETIME*>(&userTime));

    // Fill in the tms structure
    buffer->tms_cstime = 0; // We do not use these two anyway
    buffer->tms_cutime = 0;
    buffer->tms_stime = FileTimeToClockTime(kernelTime);
    buffer->tms_utime = FileTimeToClockTime(userTime);

    // Use the high-resolution performance counter.
    // The drawback is that we cannot let this thread switch between
    // individual processors because that would give us incorrect values.
    // This can be solved by calling SetThreadAffinityMask at the beginning
    // of main(..) function in case times(..) is invoked from the main thread
    // only.
    unsigned long long freq, time;
    QueryPerformanceFrequency(reinterpret_cast<LARGE_INTEGER*>(&freq));
    QueryPerformanceCounter(reinterpret_cast<LARGE_INTEGER*>(&time));
    clock_t ret = clock_t((time * CLOCKS_PER_SEC) / (freq? freq : 1));
    return ret;
}

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
    fprintf(stderr, "} took %i ms.\n", (time*1000/CLOCKS_PER_SEC));
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
