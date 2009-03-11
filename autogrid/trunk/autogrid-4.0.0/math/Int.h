#pragma once

// Definice fixnich typu
#if defined(_MSC_VER)
    #define RUNE_TYPEDEF_FIXED_INT(bits) typedef signed __int##bits int##bits; typedef unsigned __int##bits uint##bits
#else
    #define RUNE_TYPEDEF_FIXED_INT(bits) typedef int##bits##_t int##bits; typedef uint##bits##_t uint##bits
#endif

namespace Rune
{
    RUNE_TYPEDEF_FIXED_INT(8);
    RUNE_TYPEDEF_FIXED_INT(16);
    RUNE_TYPEDEF_FIXED_INT(32);
    RUNE_TYPEDEF_FIXED_INT(64);
}
