/*
    AutoGrid

    Copyright (C) 1989-2007, Garrett M. Morris, David S. Goodsell, Ruth Huey, Arthur J. Olson,
    All Rights Reserved.
    Copyright (C) 2008-2009, Marek Olsak (maraeo@gmail.com), All Rights Reserved.

    AutoGrid is a Trade Mark of The Scripps Research Institute.

    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation; either version 2
    of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
*/

#if defined(_MSC_VER)
    // disable the warning: ' function ': was declared deprecated
    #pragma warning (disable: 4996)
#endif

#include "utils.h"
#include "exceptions.h"
#include "times.h"

#if defined(beginTimer)
    #undef beginTimer
#endif
#if defined(endTimer)
    #undef endTimer
#endif

// the BOINC API header files
#if defined(BOINC)
    #include "diagnostics.h"
    #include "boinc_api.h"
    #include "filesys.h"    // boinc_fopen(), etc...
#endif

// initializes BOINC
void boincInit()
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

void boincDone()
{
#if defined(BOINCCOMPOUND)
    boinc_fraction_done(1);
#endif
#if defined(BOINC)
    boinc_finish(0); // should not return
#endif
}

// fopen rewrite to either use BOINC api or normal system call
FILE *boincOpenFile(const char *path, const char *mode)
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

class Timer
{
public:
    Timer(): real(0), user(0), sys(0) {}

    void start()
    {
        tms t;
        realStart = times(&t);
        userStart = t.tms_utime;
        sysStart = t.tms_stime;
    }

    void stop()
    {
        tms t;
        real += times(&t) - realStart;
        user += t.tms_utime - userStart;
        sys += t.tms_stime - sysStart;
    }

    void log() const
    {
        int cps = getClocksPerSec() ;
        fprintf(stderr, "Real: %i ms,\tCPU: %i ms,\tSys: %i ms\n", real * 1000 / cps, user * 1000 / cps, sys * 1000 / cps);
    }

    bool started() const
    {
        return real != 0;
    }

private:
    clock_t real, user, sys;
    clock_t realStart, userStart, sysStart;
};

static Timer timers[256];

void beginTimer(int id)
{
    timers[id].start();
}

void endTimer(int id)
{
    timers[id].stop();
}

void logTimers()
{
    for (int i = 0; i < 256; i++)
        if (timers[i].started())
        {
            fprintf(stderr, "Timer %i:\t", i);
            timers[i].log();
        }
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
