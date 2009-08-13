#pragma once
#include <signal.h>
#include <setjmp.h>

static jmp_buf __jbuf;

static void __on_signal(int s)
{
    signal(s, SIG_DFL);
    longjmp(__jbuf, 1);
}

#define UNIX_SIGNAL_TRY(s)      \
    {                           \
        signal(s, __on_signal); \
        switch (setjmp(__jbuf)) \
        {                       \
        case 0:

#define UNIX_SIGNAL_CATCH()     \
            break;              \
        default:

#define UNIX_SIGNAL_END()       \
        }                       \
    }
