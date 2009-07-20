/*
    Linear Algebra / Math library

    Copyright (C) 2003-2009, Marek Olsak (maraeo@gmail.com), All Rights Reserved.
    Copyright (C) 2003-2005, Tomas Pastorek (tomas@tomaspastorek.cz), All Rights Reserved.

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
