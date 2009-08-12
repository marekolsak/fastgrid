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

#pragma once
#include <cstdio>
#include "GridMap.h"
#include "ProgramParameters.h"

void saveAVSGridmapsFile(const GridMapList &gridmaps, const InputData *input, const ProgramParameters &programParams, LogFile &logFile);

// BOINC
void boincInit();
void boincDone();
FILE *boincOpenFile(const char *path, const char *mode);

// Timer
class Timer
{
public:
    // returns an unused timer, it's released automatically at the end of the program
    static Timer *startNew(const char *name, bool start = true);

    Timer();
    ~Timer();
    void start();
    void stop();
    void reset();
    void log(FILE *file, bool reset = true);
    void stopAndLog(FILE *file, bool reset = true);
    clock_t getReal() const;
    clock_t getUser() const;
    clock_t getSys() const;

    static void logAll(FILE *file);

private:
    struct Private;
    Private *p;
};
