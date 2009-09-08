/*
    FastGrid (formerly AutoGrid)

    Copyright (C) 2009 The Scripps Research Institute. All rights reserved.
    Copyright (C) 2009 Masaryk University. All rights reserved.

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

// the BOINC API header files
#if defined(BOINC)
    #include "diagnostics.h"
    #include "boinc_api.h"
    #include "filesys.h"    // boinc_fopen(), etc...
#endif

#include "Utils.h"
#include "Exceptions.h"
#include "times.h"

void saveAVSGridmapsFile(const GridMapList &gridmaps, const InputData *input, const ProgramParameters &programParams, LogFile &logFile)
{
    FILE *fldFileAVS;
    if ((fldFileAVS = boincOpenFile(input->fldFilenameAVS, "w")) == 0)
    {
        logFile.printErrorFormatted(ERROR, "can't create grid dimensions data file %s\n", input->fldFilenameAVS);
        logFile.printError(FATAL_ERROR, "Unsuccessful completion.\n\n");
    }
    else
        logFile.printFormatted("\nCreating (AVS-readable) grid maps file : %s\n", input->fldFilenameAVS);

    int numMaps = gridmaps.getNumMapsInclFloatingGrid();
    fprintf(fldFileAVS, "# AVS field file\n#\n");
    fprintf(fldFileAVS, "# AutoDock Atomic Affinity and Electrostatic Grids\n#\n");
    fprintf(fldFileAVS, "# Created by %s.\n#\n", programParams.getProgramName());
    fprintf(fldFileAVS, "#SPACING %.3f\n", float(input->gridSpacing));
    fprintf(fldFileAVS, "#NELEMENTS %d %d %d\n", input->numGridPoints.x-1, input->numGridPoints.y-1, input->numGridPoints.z-1);
    fprintf(fldFileAVS, "#CENTER %.3lf %.3lf %.3lf\n", input->gridCenter.x, input->gridCenter.y, input->gridCenter.z);
    fprintf(fldFileAVS, "#MACROMOLECULE %s\n", input->receptorFilename);
    fprintf(fldFileAVS, "#GRID_PARAMETER_FILE %s\n#\n", programParams.getGridParameterFilename());
    fprintf(fldFileAVS, "ndim=3\t\t\t# number of dimensions in the field\n");
    fprintf(fldFileAVS, "dim1=%d\t\t\t# number of x-elements\n", input->numGridPoints.x);
    fprintf(fldFileAVS, "dim2=%d\t\t\t# number of y-elements\n", input->numGridPoints.y);
    fprintf(fldFileAVS, "dim3=%d\t\t\t# number of z-elements\n", input->numGridPoints.z);
    fprintf(fldFileAVS, "nspace=3\t\t# number of physical coordinates per point\n");
    fprintf(fldFileAVS, "veclen=%d\t\t# number of affinity values at each point\n", numMaps);
    fprintf(fldFileAVS, "data=float\t\t# data type (byte, integer, float, double)\n");
    fprintf(fldFileAVS, "field=uniform\t\t# field type (uniform, rectilinear, irregular)\n");
    for (int i = 0; i < 3; i++)
        fprintf(fldFileAVS, "coord %d file=%s filetype=ascii offset=%d\n", (i + 1), input->xyzFilename, (i * 2));
    for (int i = 0; i < gridmaps.getNumAtomMaps(); i++)
        fprintf(fldFileAVS, "label=%s-affinity\t# component label for variable %d\n", gridmaps[i].type, (i + 1));
    fprintf(fldFileAVS, "label=Electrostatics\t# component label for variable %d\n", numMaps - 2);
    fprintf(fldFileAVS, "label=Desolvation\t# component label for variable %d\n", numMaps - 1);
    if (gridmaps.containsFloatingGrid())
        fprintf(fldFileAVS, "label=Floating_Grid\t# component label for variable %d\n", numMaps);
    fprintf(fldFileAVS, "#\n# location of affinity grid files and how to read them\n#\n");

    for (int i = 0; i < gridmaps.getNumMaps(); i++)
        fprintf(fldFileAVS, "variable %d file=%s filetype=ascii skip=6\n", (i + 1), gridmaps[i].filename);

    if (gridmaps.containsFloatingGrid())
        fprintf(fldFileAVS, "variable %d file=%s filetype=ascii skip=6\n", numMaps, input->floatingGridFilename);
    fclose(fldFileAVS);
}

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

///////////////////////////////////////////////////////////////////////////////////////////////////

static Timer timers[256];

struct Timer::Private
{
private:
    friend class Timer;

    const char *name;
    clock_t real, user, sys;
    clock_t realStart, userStart, sysStart;
};

Timer *Timer::startNew(const char *name, bool start)
{
    static int i = 0;
    Timer *t = &timers[i++];
    t->p->name = name;
    if (start)
        t->start();
    return t;
}

Timer::Timer()
{
    p = new Private();
}

Timer::~Timer()
{
    delete p;
}

void Timer::start()
{
    tms t;
    p->realStart = times(&t);
    p->userStart = t.tms_utime;
    p->sysStart = t.tms_stime;
}

void Timer::stop()
{
    tms t;
    p->real += times(&t) - p->realStart;
    p->user += t.tms_utime - p->userStart;
    p->sys += t.tms_stime - p->sysStart;
}

void Timer::reset()
{
    p->real = 0;
    p->user = 0;
    p->sys = 0;
}

void Timer::log(FILE *file, bool reset)
{
    int cps = getClocksPerSec();
    fprintf(file, "%s  Real: %li ms,\tUser: %li ms,\tSys: %li ms\n", p->name, p->real * 1000 / cps, p->user * 1000 / cps, p->sys * 1000 / cps);
    fflush(file);
    if (reset)
        this->reset();
}

void Timer::stopAndLog(FILE *file, bool reset)
{
    stop();
    log(file, reset);
}

void Timer::logAll(FILE *file)
{
    for (int i = 0; i < 256; i++)
        if (timers[i].p->real != 0)
            timers[i].log(file);
}

clock_t Timer::getReal() const
{
    return p->real;
}

clock_t Timer::getUser() const
{
    return p->user;
}

clock_t Timer::getSys() const
{
    return p->sys;
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
