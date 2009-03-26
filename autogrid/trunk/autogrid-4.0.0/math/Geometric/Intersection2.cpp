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
        Zjisti, zda 2 primky koliduji, kvuli nepresnosti floatu je mozne, ze budou kolidovat vzdy
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API bool Intersect(const Ray2<T> &ray1, const Ray2<T> &ray2)
    {
        if (!Math<T>::SafeIsZero(ray1.normal.x - ray2.normal.x) ||
            !Math<T>::SafeIsZero(ray1.normal.y - ray2.normal.y)) return true;
        return Math<T>::SafeIsZero(ray1.c - ray2.c);
    }

    template RUNEMATH_API bool Intersect(const Ray2<float> &ray1, const Ray2<float> &ray2);
    template RUNEMATH_API bool Intersect(const Ray2<double> &ray1, const Ray2<double> &ray2);

    /**
        Zjisti, zda 2 primky koliduji, a vraci prusecik
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API bool Intersect(const Ray2<T> &ray1, const Ray2<T> &ray2, Vec2<T> &point)
    {
        if (Intersect(ray1, ray2))
        {
            // Spocita prusecik
            T f = 1 / (ray1.normal.x * ray2.normal.y - ray2.normal.x * ray1.normal.y);
            point.x =  (ray1.normal.y * ray2.c - ray2.normal.y * ray1.c) * f;
            point.y = -(ray1.normal.x * ray2.c - ray2.normal.x * ray1.c) * f;
            return true;
        }
        return false;
    }

    template RUNEMATH_API bool Intersect(const Ray2<float> &ray1, const Ray2<float> &ray2, Vec2<float> &point);
    template RUNEMATH_API bool Intersect(const Ray2<double> &ray1, const Ray2<double> &ray2, Vec2<double> &point);

    /**
        Spocita kolizi usecky a primky a vraci prusecik
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API bool Intersect(const Line2<T> &line, const Ray2<T> &ray, Vec2<T> &point)
    {
        if (Intersect(line.ray, ray, point))
            if (Vec2<T>::DistanceSqr(line.origin, point) <= line.extent*line.extent) return true;
        return false;
    }

    template RUNEMATH_API bool Intersect(const Line2<float> &line, const Ray2<float> &ray, Vec2<float> &point);
    template RUNEMATH_API bool Intersect(const Line2<double> &line, const Ray2<double> &ray, Vec2<double> &point);

    /**
        Spocita kolizi usecky a primky
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API bool Intersect(const Line2<T> &line, const Ray2<T> &ray)
    {
        Vec2<T> tmp;
        return Intersect<T>(line, ray, tmp);
    }

    template RUNEMATH_API bool Intersect(const Line2<float> &line, const Ray2<float> &ray);
    template RUNEMATH_API bool Intersect(const Line2<double> &line, const Ray2<double> &ray);

    /**
        Spocita kolizi dvou usecek a vraci prusecik
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API bool Intersect(const Line2<T> &line1, const Line2<T> &line2, Vec2<T> &point)
    {
        if (Intersect(line1.ray, line2.ray, point))
            if (Vec2<T>::DistanceSqr(line1.origin, point) <= line1.extent*line1.extent &&
                Vec2<T>::DistanceSqr(line2.origin, point) <= line2.extent*line2.extent) return true;
        return false;
    }

    template RUNEMATH_API bool Intersect(const Line2<float> &line1, const Line2<float> &line2, Vec2<float> &point);
    template RUNEMATH_API bool Intersect(const Line2<double> &line1, const Line2<double> &line2, Vec2<double> &point);

    /**
        Spocita kolizi dvou usecek
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API bool Intersect(const Line2<T> &line1, const Line2<T> &line2)
    {
        Vec2<T> tmp;
        return Intersect(line1, line2, tmp);
    }

    template RUNEMATH_API bool Intersect(const Line2<float> &line1, const Line2<float> &line2);
    template RUNEMATH_API bool Intersect(const Line2<double> &line1, const Line2<double> &line2);

    /**
        Spocita kolizi dvou kruhu
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API bool Intersect(const Circle2<T> &s1, const Circle2<T> &s2)
    {
        T x = s1.pos.x-s2.pos.x, y = s1.pos.y-s2.pos.y, r = s1.radius+s2.radius;
        return x*x + y*y <= r*r;
    }

    template RUNEMATH_API bool Intersect(const Circle2<float> &s1, const Circle2<float> &s2);
    template RUNEMATH_API bool Intersect(const Circle2<double> &s1, const Circle2<double> &s2);

    /**
        Spocita kolizi dvou kruhu, vraci navic bod dotyku kruhu
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API bool Intersect(const Circle2<T> &s1, const Circle2<T> &s2, Vec2<T> &point)
    {
        if (Intersect(s1, s2))
        {
            Vec2<T> dir = s1.pos - s2.pos;
            dir.Normalize();
            point = s2.pos + dir * s2.radius;
            return true;
        }
        return false;
    }

    template RUNEMATH_API bool Intersect(const Circle2<float> &s1, const Circle2<float> &s2, Vec2<float> &point);
    template RUNEMATH_API bool Intersect(const Circle2<double> &s1, const Circle2<double> &s2, Vec2<double> &point);

    /**
        Spocita kolizi kruhu a bodu
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API bool Intersect(const Circle2<T> &s, const Vec2<T> &point)
    {
        return Vec2<T>::DistanceSqr(s.pos, point) < s.radius*s.radius;
    }

    template RUNEMATH_API bool Intersect(const Circle2<float> &s, const Vec2<float> &point);
    template RUNEMATH_API bool Intersect(const Circle2<double> &s, const Vec2<double> &point);

    /**
        Spocita kolizi obdelnika a bodu
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API bool Intersect(const Rect2<T> &r, const Vec2<T> &point)
    {
        return r.leftTop <= point && r.rightBottom >= point;
    }

    template RUNEMATH_API bool Intersect(const Rect2<float> &r, const Vec2<float> &point);
    template RUNEMATH_API bool Intersect(const Rect2<double> &r, const Vec2<double> &point);
    template RUNEMATH_API bool Intersect(const Rect2<int32> &r, const Vec2<int32> &point);
}
