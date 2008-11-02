#include "logfile.h"
#include "utils.h"
#include "exceptions.h"

LogFile::LogFile(double versionNumber, const char *programName, const char *filename): file(0)
{
    // Initialize the log file
    if (!filename || !filename[0])
        file = stdout;
    else
    {
        file = openFile(filename, "w");
        if (!file)
        {
            fprintf(stderr, "\n%s: Sorry, I can't create the log file \"%s\"\n", programName, filename);
            fprintf(stderr, "\n%s: Unsuccessful Completion.\n\n", programName);
            throw ExitProgram(911);
        }
    }

    // Output basic information
    printBanner(file, versionNumber);
    fprintf(file, "                           $Revision: 1.58 $\n\n\n");
    fprintf(file, "\nMaximum number of maps that can be computed = %d (defined by MAX_MAPS in \"autocomm.h\").\n\n\n", MAX_MAPS);
    fprintf(file, "This file was created at:\t\t\t");
    {
        char strtmp[MAX_CHARS];
        fprintf(file, getDate(1, strtmp, MAX_CHARS));
        fprintf(file, "                   using:\t\t\t\"%s\"\n", getHostname(strtmp, MAX_CHARS));
    }
    fprintf(file, "\n\n");
}

LogFile::~LogFile()
{
    if (file && file != stdout)
        fclose(file);
}
