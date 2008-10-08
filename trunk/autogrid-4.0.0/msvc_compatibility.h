#pragma once
#if defined(_MSC_VER)

// disable the warning: ' function ': was declared deprecated
#pragma warning (disable: 4996)

#include <cfloat>

// Some functions in Visual C++ differ from those in the linux/unix environment
#define isnan _isnan
#define strncasecmp _strnicmp

// WinAPI defines ERROR, we have to #undef it since we use the same name for
// something else
#if defined(ERROR)
#undef ERROR
#endif

#endif
