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

#include "../All.h"

namespace Rune
{
    /**
        Pro bod - vraci -1, kdyz je za primkou, 1 kdyz pred primkou a 0, kdyz lezi na primce
    **************************************************************************************************/
    template<typename T>
    int Ray2<T>::GetSide(const Vec2<T> &point) const
    {
        T f = GetDistance(point);
        if (f > Math<T>::Epsilon()) return 1;
        if (f < -Math<T>::Epsilon()) return -1;
        return 0;
    }

    template int Ray2<float>::GetSide(const Vec2<float> &point) const;
    template int Ray2<double>::GetSide(const Vec2<double> &point) const;

    /**
        Pro kruh - vraci -1, kdyz je za primkou, 1 kdyz pred primkou a 0, kdyz se dotyka primky
    **************************************************************************************************/
    template<typename T>
    int Ray2<T>::GetSide(const Circle2<T> &circle) const
    {
        T f = GetDistance(circle.pos);
        if (f < circle.radius) return 0;

        if (f > 0) f -= circle.radius;
        else f += circle.radius;

        if (f > Math<T>::Epsilon()) return 1;
        if (f < -Math<T>::Epsilon()) return -1;
        return 0;
    }

    template int Ray2<float>::GetSide(const Circle2<float> &circle) const;
    template int Ray2<double>::GetSide(const Circle2<double> &circle) const;
}
