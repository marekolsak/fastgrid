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
#include <ctime>
#if defined(WIN32)

// Our own implementation of the times function on Windows

struct tms
{
    clock_t tms_utime;    // CPU time used in executing the instructions of the calling process.
    clock_t tms_stime;    // CPU time used by the system on behalf of the calling process.

    // Not available on Windows
    clock_t tms_cutime;    // sum of the tms_utime values and the tms_cutime values of all terminated child processes of the calling process, whose status has been reported to the parent process by wait or waitpid; see section Process Completion. In other words, it represents the total CPU time used in executing the instructions of all the terminated child processes of the calling process, excluding child processes which have not yet been reported by wait or waitpid.
    clock_t tms_cstime;    // similar to tms_cutime, but represents the total CPU time used by the system on behalf of all the terminated child processes of the calling process.
    // All of the times are given in clock ticks. These are absolute values; in a newly created process, they are all zero. See section Creating a Process.
};

clock_t times(tms *buffer);

#else
#include <sys/times.h>
#endif

int getClocksPerSec();
