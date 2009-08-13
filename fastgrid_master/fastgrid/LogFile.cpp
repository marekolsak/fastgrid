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

#if !defined(PACKAGE_BUGREPORT)
#define PACKAGE_BUGREPORT "maraeo@gmail.com"
#endif

#if defined(_WIN32)
    #define WIN32_LEAN_AND_MEAN
    #include <Winsock2.h>
#else
    #include <unistd.h>
#endif
#include "LogFile.h"
#include "Exceptions.h"
#include "Utils.h"
#include <cstring>
#include <cstdarg>

#define FORMATTED_MSG_MAX_SIZE (1<<14)

// Formats a message
// The reason we use a macro instead of a function is that we need to access the variable number of arguments
#define FORMAT_MESSAGE(message, messageLength, format) \
        char message[FORMATTED_MSG_MAX_SIZE]; \
        int messageLength; \
        message[sizeof(message)-1] = 0; \
        va_list ap; \
        va_start(ap, format); \
        messageLength = vsnprintf(message, sizeof(message), format, ap); \
        /* if the message buffer is not long enough or hasn't been terminated by zero */ \
        if (messageLength == -1 || message[sizeof(message)-1] != 0) \
        { \
            printError(WARNING, "The following formatted string will be truncated."); \
            message[sizeof(message)-1] = 0; \
            messageLength = sizeof(message)-1; \
        } \
        va_end(ap);

LogFile::LogFile(const char *versionNumber, const char *programName, const char *filename): file(0), invClocksPerSec(1 / float(getClocksPerSec()))
{
    strncpy(this->programName, programName, MAX_CHARS);

    // Initialize the log file
    if (!filename || !filename[0])
        file = stdout;
    else
    {
        file = boincOpenFile(filename, "w");
        if (!file)
        {
            fprintf(stderr, "\n%s: Sorry, I can't create the log file \"%s\"\n", programName, filename);
            fprintf(stderr, "\n%s: Unsuccessful Completion.\n\n", programName);
            throw ExitProgram(911);
        }
    }

    // Output basic information
    printBanner(versionNumber);
    fprintf(file, "                           $Revision: 1.71 $\n\n\n"
                  "\nMaximum number of maps that can be computed = %d (defined by MAX_MAPS in \"autocomm.h\").\n\n\n"
                  "This file was created at:\t\t\t", MAX_MAPS);
    printCurrentDate(1);
    printHostname();
    print("\n\n");
}

LogFile::~LogFile()
{
    if (file && file != stdout)
        fclose(file);
}

const char *LogFile::getProgramName() const
{
    return programName;
}

void LogFile::print(const char *msg)
{
    size_t len = strlen(msg);
    if (!fwrite(msg, len, 1, file))
        printError(FATAL_ERROR, "Not enough disk space.");
}

void LogFile::printFormatted(const char *format, ...)
{
    FORMAT_MESSAGE(message, messageLength, format);
    if (!fwrite(message, messageLength, 1, file))
        printError(FATAL_ERROR, "Not enough disk space.");
}

void LogFile::printTitled(const char *msg)
{
    fprintf(file, "%s: %s", programName, msg);
}

void LogFile::printTitledFormatted(const char *format, ...)
{
    FORMAT_MESSAGE(message, messageLength, format);
    printTitled(message);
}

void LogFile::printError(ErrorLevel errorLevel, const char *msg)
{
    const char *tags[5] = {
        "ERROR",
        "ERROR",
        "WARNING",
        "INFORMATION",
        "SUGGESTION"
    };

    char outputMessage[LINE_LEN];
    snprintf(outputMessage, LINE_LEN, "\n%s: %s:  %s\n", programName, tags[errorLevel+2], msg);
    size_t len = strlen(outputMessage);
    if (!fwrite(outputMessage, len, 1, file) && errorLevel != FATAL_ERROR)
        // if it's a fatal error, no need to care because the program will be terminated anyway
        printError(FATAL_ERROR, "Not enough disk space.");

    // Only send errors, fatal errors and warnings to standard error.
    if (errorLevel <= WARNING)
        if (!fwrite(outputMessage, len, 1, stderr) && errorLevel != FATAL_ERROR)
            // the log file *might* still be writable...
            printError(FATAL_ERROR, "Cannot write to standard error output.");

    // If this is a fatal error, exit now.
    if (errorLevel == FATAL_ERROR)
        throw ExitProgram(errorLevel);
}

void LogFile::printErrorFormatted(ErrorLevel errorLevel, const char *format, ...)
{
    FORMAT_MESSAGE(message, messageLength, format);
    printError(errorLevel, message);
}

void LogFile::printExecutionTimes(Clock startTime, Clock endTime, tms *start, tms *end)
{
    fprintf(file, "Real= %.2f,  CPU= %.2f,  System= %.2f\n",
        invClocksPerSec * (endTime - startTime),
        invClocksPerSec * (end->tms_utime - start->tms_utime),
        invClocksPerSec * (end->tms_stime - start->tms_stime));
}

void LogFile::printTimeInHMS(Clock time, bool fixedOutputLength)
{
    printTimeInHMS(invClocksPerSec * time, fixedOutputLength);
}

void LogFile::printTimeInHMS(float time, bool fixedOutputLength)
{
    float hrs = 3600, min = 60;
    int h = int(time / hrs);
    float T = time - h*hrs;
    int m = int(T / min);
    float s = T - m*min;

    if (h == 0)
        if (m == 0)
            fprintf(file, fixedOutputLength ? "        %5.2fs" : "%.2fs", s);
        else
            fprintf(file, fixedOutputLength ? "    %2dm %05.2fs" : "%dm %05.2fs", m, s);
    else
        fprintf(file, fixedOutputLength ? "%2dh %02dm %05.2fs" : "%dh %02dm %05.2fs", h, m, s);
}

void LogFile::printExecutionTimesInHMS(Clock startTime, Clock endTime, tms *start, tms *end)
{
    print("Real= ");
    printTimeInHMS(invClocksPerSec * (endTime - startTime), false);
    print(",  CPU= ");
    printTimeInHMS(invClocksPerSec * (end->tms_utime - start->tms_utime), false);
    print(",  System= ");
    printTimeInHMS(invClocksPerSec * (end->tms_stime  - start->tms_stime), false);
    print("\n");
}

// Output banner...
void LogFile::printBanner(const char *versionNumber)
{
    fprintf(file,
        "\n       _______________________________________________________\n"
        "\n"
#if defined(FASTGRID)
        "________///////_________________________/////_________________/________\n"
        "_______/_________________________/_____/______________/_______/________\n"
        "_______/_________________________/_____/______________________/________\n"
        "_______///////___/////___//////_/////__/__////_/_///__/__////_/________\n"
        "_______/______________/_/________/_____/_____/_//___/_/_/____//________\n"
        "_______/________///////__/////___/_____/_____/_/______/_/_____/________\n"
        "_______/_______/_____//_______/__/___/_/_____/_/______/_/____//________\n"
        "_______/________/////_/_//////____///___/////__/______/__////_/________\n"
#else
        "__________//____________________________/////_________________/________\n"
        "_________/__/____________/_____________/______________/_______/________\n"
        "________/____/___________/_____________/______________________/________\n"
        "________/____/__/_____/_/////___/////__/__////_/_///__/__////_/________\n"
        "_______/______/_/_____/__/_____/_____/_/_____/_//___/_/_/____//________\n"
        "_______////////_/_____/__/_____/_____/_/_____/_/______/_/_____/________\n"
        "_______/______/_/____//__/___/_/_____/_/_____/_/______/_/____//________\n"
        "_______/______/__////_/___///___/////___/////__/______/__////_/________\n"
#endif
        "\n"
        "       _______________________________________________________\n"
        "\n"
        "                                ______\n"
        "                               /      \\\n"
        "                              /        \\\n"
        "                             /          \\\n"
        "                             \\    /\\    /\n"
        "                              \\  /  \\  /\n"
        "                               \\/ /\\ \\/\n"
        "                                 /  \\\n"
        "                                /____\\\n"
        "\n"
        "\n"
        "                ______________________________________ \n"
        "               |                                      |\n"
        "               |            " APPNAME " %-8s         |\n"
#if defined(FASTGRID)
        "               |                                      |\n"
        "               |          The fork of AutoGrid        |\n"
        "               |      optimized to fully leverage     |\n"
        "               |  the power of multi-core processors  |\n"
        "               |     and NVIDIA graphics hardware.    |\n"
#endif
        "               |                                      |\n"
        "               |        Garrett M. Morris, TSRI       |\n"
        "               |            Ruth Huey, TSRI           |\n"
        "               |        David S. Goodsell, TSRI       |\n"
        "               |         Arthur J. Olson, TSRI        |\n"
        "               |       Jiri Filipovic, NCBR MUNI      |\n"
        "               |         Marek Olsak, FI MUNI         |\n"
        "               |                                      |\n"
        "               |        (C) 1989-2009, TSRI           |\n"
        "               |   The Scripps Research Institute     |\n"
        "               |                                      |\n"
        "               |           (C) 2009, MUNI             |\n"
        "               |          Masaryk University          |\n"
        "               |______________________________________|\n"
        "\n"
        "                ______________________________________ \n"
        "               |                                      |\n"
        "               | Calculation of van der Waals, H-Bond,|\n"
        "               |   Electrostatic Potential Energy, &  |\n"
        "               |   Desolvation Free Energy Grid Maps  |\n"
        "               |             for AutoDock             |\n"
        "               | For help, email %-19s |\n"
        "               |______________________________________|\n"
        "\n"
        "\n"
        "\n"
        "\n", versionNumber, PACKAGE_BUGREPORT);
}

void LogFile::printCurrentDate(int flag)
{
    time_t timeNow;
    timeNow = time(&timeNow);
    tm *ts = localtime(&timeNow);

    if (flag == 1)
        fprintf(file, "%d:%02d %02d\" %s, %02d/%02d/%4d\n",
            ((ts->tm_hour >  12) ? (ts->tm_hour-12) : ts->tm_hour), ts->tm_min, ts->tm_sec,
            ((ts->tm_hour >= 12) ? "p.m." : "a.m."),
            (ts->tm_mon + 1), ts->tm_mday, 1900+ts->tm_year);
    else if (flag == 2)
        fprintf(file, "%s", ctime(&timeNow));
    else
        fprintf(file, "%d:%02d %02d\" %s\n",
            ((ts->tm_hour >  12) ? (ts->tm_hour-12) : ts->tm_hour), ts->tm_min, ts->tm_sec,
            ((ts->tm_hour >= 12) ? "pm" : "am"));
}

void LogFile::printHostname()
{
#if defined(_WIN32)
    WSADATA wsaData;
    WSAStartup(MAKEWORD(2, 2), &wsaData);
#endif

    char buffer[MAX_CHARS];
    if (gethostname(buffer, MAX_CHARS) != 0)
        strncpy(buffer, "(gethostname returned an error)", MAX_CHARS);

#if defined(_WIN32)
    WSACleanup();
#endif

    fprintf(file, "                   using:\t\t\t\"%s\"\n", buffer);
}
