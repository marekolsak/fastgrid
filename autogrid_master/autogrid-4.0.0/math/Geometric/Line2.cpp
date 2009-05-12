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
        Nastavi body usecky
    **************************************************************************************************/
    template<typename T>
    void Line2<T>::SetPoints(const Vec2<T> &p1, const Vec2<T> &p2)
    {
        origin = Vec2<T>::Center(p1, p2);
        ray.Set((p2-p1).GetNormalized().GetNormal(), origin);
        extent = Vec2<T>::Distance(origin, p1);
    }

    template void Line2<float>::SetPoints(const Vec2<float> &p1, const Vec2<float> &p2);
    template void Line2<double>::SetPoints(const Vec2<double> &p1, const Vec2<double> &p2);

    /**
        Pro bod - vraci -1, kdyz je za primkou, 1 kdyz pred primkou a 0, kdyz lezi na primce
    **************************************************************************************************/
    template<typename T>
    int Line2<T>::GetSide(const Vec2<T> &_point) const
    {
        T f = ray.GetDistance(_point);
        if (f > Math<T>::Epsilon()) return 1;
        if (f < -Math<T>::Epsilon()) return -1;
        return 0;
    }

    template int Line2<float>::GetSide(const Vec2<float> &_point) const;
    template int Line2<double>::GetSide(const Vec2<double> &_point) const;
}
