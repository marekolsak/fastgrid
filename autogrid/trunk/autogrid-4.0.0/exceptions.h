#pragma once
#include <exception>

class ExitProgram : public std::exception
{
public:
    ExitProgram(int exitCode): code(exitCode) {}
    int getExitCode() const { return code; }

private:
    int code;
};
