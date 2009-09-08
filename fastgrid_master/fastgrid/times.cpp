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

#if !defined(_WIN32)
    #include <unistd.h>
    #include <cstdio>
    #include "Exceptions.h"
#endif
#include "times.h"

int getClocksPerSec()
{
#if defined(_WIN32)
    return CLOCKS_PER_SEC;
#else
    int clocks = sysconf(_SC_CLK_TCK);
    if (clocks < 0)
    {
        fprintf(stderr, "\"sysconf(_SC_CLK_TCK)\" failed in getClocksPerSec()\n");
        throw ExitProgram(-1);
    }
    return clocks;
#endif
}

#if defined(_WIN32)
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <cerrno>

// converts WinAPI's FILETIME to clock_t
static clock_t FileTimeToClockTime(unsigned long long fileTime)
{
    // fileTime contains the time in 100s of nanoseconds
    return clock_t((fileTime * CLOCKS_PER_SEC) / 10000000ull);
}

// there is no times(..) function on Windows so we have to implement it ourselves
clock_t times(tms *buffer)
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
    // This can be solved by calling SetThreadAffinityMask at the start
    // of the program in case times(..) is invoked from the main thread
    // only. The InitWinSetup class takes care of that.
    unsigned long long freq, time;
    QueryPerformanceFrequency(reinterpret_cast<LARGE_INTEGER*>(&freq));
    QueryPerformanceCounter(reinterpret_cast<LARGE_INTEGER*>(&time));
    clock_t ret = clock_t((time * CLOCKS_PER_SEC) / (freq? freq : 1));
    return ret;
}

class InitWinSetup
{
public:
    InitWinSetup()
    {
        // Bound this thread to the first CPU
        SetThreadAffinityMask(GetCurrentThread(), 1);
        // Raise the process priority
        //SetPriorityClass(GetCurrentProcess(), ABOVE_NORMAL_PRIORITY_CLASS);
        SetPriorityClass(GetCurrentProcess(), HIGH_PRIORITY_CLASS);
    }
};

static InitWinSetup winSetup;

#endif
