#include <cstdlib>
#include <cstdio>
#include <cstring>
#include "msvc_compatibility.h"
#include "autocomm.h"
//#include "autoglobal.h"
#include "atom_parameter_manager.h"
#include "main_utils.h"


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
        strcpy(tag, "ERROR");
        break;
    case WARNING:
        strcpy(tag, "WARNING");
        break;
    case INFORMATION:
        strcpy(tag, "INFORMATION");
        break;
    case SUGGESTION:
        strcpy(tag, "SUGGESTION");
        break;
    }

    sprintf(output_message, "\n%s: %s:  %s\n", programname, tag, message);

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

int get_rec_index(const char key[])
{
    ParameterEntry *found_parm;

    found_parm = apm_find(key);
    if (found_parm != NULL)
        return found_parm->rec_index;
    return -1;
}

#if defined(BOINC)
// Dummy graphics API entry points.  This app does not do graphics, but it still must provide these callbacks. 

void app_graphics_render(int xs, int ys, double time_of_day){}
void app_graphics_reread_prefs(){}
void boinc_app_mouse_move(int x, int y, bool left, bool middle, bool right){}
void boinc_app_mouse_button(int x, int y, int which, bool is_down){}
void boinc_app_key_press(int wParam, int lParam){}
void boinc_app_key_release(int wParam, int lParam){}
#endif
