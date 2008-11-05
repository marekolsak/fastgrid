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
    LogFile(double versionNumber, const char *programName, const char *filename);
    ~LogFile();

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

    // temporary type fallback
    operator FILE*() { return file; }

private:
    FILE *file;
    char programName[MAX_CHARS];
    float invClocksPerSec;

    void printBanner(double versionNumber);
    void printCurrentDate(int flag);
    void printHostname();
};
