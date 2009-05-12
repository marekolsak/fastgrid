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
    /**
        Trida kruhu ve 2D
    **************************************************************************************************/
    template<typename T>
    class Circle2
    {
    public:
        Vec2<T> pos;
        T radius;

        Circle2() {}
        Circle2(const Vec2<T> &position, const T circleradius): pos(position), radius(circleradius) {}
    };

    typedef Circle2<float> Circle2f;
    typedef Circle2<double> Circle2d;
}
