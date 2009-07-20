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
    template<typename T>
    RUNEMATH_API bool Vec4<T>::SafeIsEqual(const Vec4 &v) const
    {
        return Math<T>::SafeIsEqual(x, v.x) && Math<T>::SafeIsEqual(y, v.y) &&
               Math<T>::SafeIsEqual(z, v.z) && Math<T>::SafeIsEqual(w, v.w);
    }

    template RUNEMATH_API bool Vec4<float>::SafeIsEqual(const Vec4 &v) const;
    template RUNEMATH_API bool Vec4<double>::SafeIsEqual(const Vec4 &v) const;

    template<typename T>
    RUNEMATH_API T Vec4<T>::Dot(const Vec4<T> &v1, const Vec4<T> &v2)
    {
        return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z + v1.w*v2.w;
    }

    template RUNEMATH_API float Vec4<float>::Dot(const Vec4<float> &v1, const Vec4<float> &v2);
    template RUNEMATH_API double Vec4<double>::Dot(const Vec4<double> &v1, const Vec4<double> &v2);

    template<typename T>
    RUNEMATH_API Vec4<T> Vec4<T>::Cross(const Vec4<T> &a, const Vec4<T> &b, const Vec4<T> &c)
    {
        // BTW nejsou matice jednodussi? :)
        return Vec4<T>(a.w*b.z*c.y + a.z*b.w*c.y + a.w*b.y*c.z - a.y*b.w*c.z - a.z*b.y*c.w + a.y*b.z*c.w,
                       a.w*b.z*c.x - a.z*b.w*c.x - a.w*b.x*c.z + a.x*b.w*c.z + a.z*b.x*c.w - a.x*b.z*c.w,
                      -a.w*b.y*c.x + a.y*b.w*c.x + a.w*b.x*c.y - a.x*b.w*c.y - a.y*b.x*c.w + a.x*b.y*c.w,
                       a.z*b.y*c.x - a.y*b.z*c.x - a.z*b.x*c.y + a.x*b.z*c.y + a.y*b.x*c.z - a.x*b.y*c.z);
    }

    template RUNEMATH_API Vec4<float> Vec4<float>::Cross(const Vec4<float> &a, const Vec4<float> &b, const Vec4<float> &c);
    template RUNEMATH_API Vec4<double> Vec4<double>::Cross(const Vec4<double> &a, const Vec4<double> &b, const Vec4<double> &c);
}
