#if defined(_WIN32)
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <ctime>
#include <errno.h>
#include "times.h"

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

    // Use a high-resolution performance counter.
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
