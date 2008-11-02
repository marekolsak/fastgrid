#pragma once
#include <cstdio>

class LogFile
{
public:
    LogFile(double versionNumber, const char *programName, const char *filename);
    ~LogFile();

    // temporary type fallback
    operator FILE*() { return file; }

private:
    FILE *file;

    void writeBanner();
};
