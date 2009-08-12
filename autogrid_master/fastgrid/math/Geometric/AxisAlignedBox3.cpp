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
        Vraci vsechny rohy AABB
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API void AxisAlignedBox3<T>::GetVertices(Vec3<T> *vertices) const
    {
        vertices[0] = min;
        vertices[1].Set(max.x, min.y, min.z);
        vertices[2].Set(min.x, max.y, min.z);
        vertices[3].Set(min.x, min.y, max.z);
        vertices[4].Set(max.x, max.y, min.z);
        vertices[5].Set(min.x, max.y, max.z);
        vertices[6].Set(max.x, min.y, max.z);
        vertices[7] = max;
    }

    template RUNEMATH_API void AxisAlignedBox3<float>::GetVertices(Vec3<float> *vertices) const;
    template RUNEMATH_API void AxisAlignedBox3<double>::GetVertices(Vec3<double> *vertices) const;

    /**
        Spocita AABB ze stredu a vektoru od stredu
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API void AxisAlignedBox3<T>::SetCenterAndExtents(const Vec3<T> &center, const Vec3<T> &extents)
    {
        min = center - extents;
        max = center + extents;
    }

    template RUNEMATH_API void AxisAlignedBox3<float>::SetCenterAndExtents(const Vec3<float> &center, const Vec3<float> &extents);
    template RUNEMATH_API void AxisAlignedBox3<double>::SetCenterAndExtents(const Vec3<double> &center, const Vec3<double> &extents);

    /**
        Spocita AABB z vertexu
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API void AxisAlignedBox3<T>::Approximate(const Vec3<T> *vertices, int count)
    {
        assert(count);
        min = *vertices;
        max = *vertices;
        const Vec3<T> *it, *last = vertices+count;
        for (it = vertices+1; it != last; ++it)
        {
            if (it->x < min.x)      min.x = it->x;
            else if (it->x > max.x) max.x = it->x;
            if (it->y < min.y)      min.y = it->y;
            else if (it->y > max.y) max.y = it->y;
            if (it->z < min.z)      min.z = it->z;
            else if (it->z > max.z) max.z = it->z;
        }
    }

    template RUNEMATH_API void AxisAlignedBox3<float>::Approximate(const Vec3<float> *vertices, int count);
    template RUNEMATH_API void AxisAlignedBox3<double>::Approximate(const Vec3<double> *vertices, int count);

    /**
        Spocita AABB z nekolika AABB
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API void AxisAlignedBox3<T>::Approximate(const AxisAlignedBox3<T> *boxes, int count)
    {
        min = boxes->min;
        max = boxes->max;
        const AxisAlignedBox3<T> *it, *last = boxes+count;
        for (it = boxes+1; it != last; ++it)
        {
            if (it->min.x < min.x)      min.x = it->min.x;
            else if (it->max.x > max.x) max.x = it->max.x;
            if (it->min.y < min.y)      min.y = it->min.y;
            else if (it->max.y > max.y) max.y = it->max.y;
            if (it->min.z < min.z)      min.z = it->min.z;
            else if (it->max.z > max.z) max.z = it->max.z;
        }
    }

    template RUNEMATH_API void AxisAlignedBox3<float>::Approximate(const AxisAlignedBox3<float> *boxes, int count);
    template RUNEMATH_API void AxisAlignedBox3<double>::Approximate(const AxisAlignedBox3<double> *boxes, int count);
}
