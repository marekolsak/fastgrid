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
#include "autogrid.h"
#include "times.h"

// WinAPI defines ERROR, we have to #undef it since we use the same name
#if defined(ERROR)
    #undef ERROR
#endif

// LogFile::printError{Formatted}() can be used with one of the following:
enum ErrorLevel
{
    FATAL_ERROR = -2,
    ERROR = -1,
    WARNING =  0,
    INFORMATION = 1,
    SUGGESTION = 2
};

class LogFile
{
public:
    LogFile(const char *versionNumber, const char *programName, const char *filename);
    ~LogFile();

    // there precalculated values might be useful
    const char *getProgramName() const;

    // print a message to the log file
    void print(const char *msg);
    void printFormatted(const char *format, ...);

    // print program's name at the beginning of every message
    void printTitled(const char *msg);
    void printTitledFormatted(const char *format, ...);

    // print an error or informational message to the log file and (if it's an error) stderr
    // a fatal error immediately terminates the program
    void printError(ErrorLevel errorLevel, const char *msg);
    void printErrorFormatted(ErrorLevel errorLevel, const char *format, ...);

    // print a time in seconds
    void printExecutionTimes(Clock startTime, Clock endTime, tms *start, tms *end);

    // print a time in the hours-minutes-seconds format
    void printTimeInHMS(Clock time, bool fixedOutputLength = true);
    void printTimeInHMS(float time, bool fixedOutputLength = true);
    void printExecutionTimesInHMS(Clock startTime, Clock endTime, tms *start, tms *end);

private:
    FILE *file;
    char programName[MAX_CHARS];
    float invClocksPerSec;

    void printBanner(const char *versionNumber);
    void printCurrentDate(int flag);
    void printHostname();
};
