#pragma once

// The replacement for C's exit function.
// This must be used in place of exit to make sure destructors are properly called
class ExitProgram
{
public:
    ExitProgram(int exitCode): code(exitCode) {}
    int getExitCode() const { return code; }

private:
    int code;
};
