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

namespace Rune
{
    // Sablona pro pretypovani pres union, vhodne pro bezpecnou reinterpretaci vsech datovych typu
    // krome pointeru a referenci (pro tyto se pouziva reinterpret_cast).
    template<typename T1, typename T2>
    inline T1 union_cast(const T2 &a)
    {
        union
        {
            T2 a;
            T1 b;
        } myUnion;
        myUnion.a = a;
        return myUnion.b;
    }
}
