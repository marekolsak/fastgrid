/*
    Auxiliary Math library

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
// AABB = axis-aligned bounding box

namespace Rune
{
    /**
        Spocita kolizi roviny a primky
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API bool Intersect(const Plane3<T> &p, const Ray3<T> &ray, Vec3<T> &point, T &tparam)
    {
        T dp = Vec3<T>::Dot(ray.direction, p.normal);
        if (Math<T>::SafeIsZero(dp)) return false;
        T t = -p.GetDistance(ray.origin) / dp;
        if (t < 0) return false;
        point = (ray.origin + ray.direction * t);
        tparam = t;
        return true;
    }

    template RUNEMATH_API bool Intersect(const Plane3<float> &p, const Ray3<float> &ray, Vec3<float> &point, float &tparam);
    template RUNEMATH_API bool Intersect(const Plane3<double> &p, const Ray3<double> &ray, Vec3<double> &point, double &tparam);

    /**
        Spocita kolizi AABB a bodu
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API bool Intersect(const AxisAlignedBox3<T> &b, const Vec3<T> &p)
    {
        return p >= b.min && p <= b.max;
    }

    template RUNEMATH_API bool Intersect(const AxisAlignedBox3<float> &b, const Vec3<float> &p);
    template RUNEMATH_API bool Intersect(const AxisAlignedBox3<double> &b, const Vec3<double> &p);

    /**
        Spocita kolizi dvou AABB
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API bool Intersect(const AxisAlignedBox3<T> &b1, const AxisAlignedBox3<T> &b2)
    {
        if (b1.max.x < b2.min.x || b1.min.x > b2.max.x) return false;
        if (b1.max.y < b2.min.y || b1.min.y > b2.max.y) return false;
        if (b1.max.z < b2.min.z || b1.min.z > b2.max.z) return false;
        return true;
    }

    template RUNEMATH_API bool Intersect(const AxisAlignedBox3<float> &b1, const AxisAlignedBox3<float> &b2);
    template RUNEMATH_API bool Intersect(const AxisAlignedBox3<double> &b1, const AxisAlignedBox3<double> &b2);

    /**
        Spocita kolizi AABB a koule
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API bool Intersect(const AxisAlignedBox3<T> &b, const Sphere3<T> &s)
    {
        T dist = 0;

        if (s.pos.x < b.min.x)      dist += Math<T>::Sqr(s.pos.x - b.min.x);
        else if (s.pos.x > b.max.x) dist += Math<T>::Sqr(s.pos.x - b.max.x);
        if (s.pos.y < b.min.y)      dist += Math<T>::Sqr(s.pos.y - b.min.y);
        else if (s.pos.y > b.max.y) dist += Math<T>::Sqr(s.pos.y - b.max.y);
        if (s.pos.z < b.min.z)      dist += Math<T>::Sqr(s.pos.z - b.min.z);
        else if (s.pos.z > b.max.z) dist += Math<T>::Sqr(s.pos.z - b.max.z);

        return dist <= s.radius*s.radius;
    }

    template RUNEMATH_API bool Intersect(const AxisAlignedBox3<float> &b, const Sphere3<float> &s);
    template RUNEMATH_API bool Intersect(const AxisAlignedBox3<double> &b, const Sphere3<double> &s);

    /**
        Spocita kolizi AABB a koule
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API bool Intersect(const AxisAlignedBox3<T> &b, Interior bi, const Sphere3<T> &s, Interior si)
    {
        switch (si)
        {
        case SOLID:
            switch (bi)
            {
            case SOLID:
                return Intersect(b, s);

            case HOLLOW:
                {
                    T dist = 0;
                    bool face = false;

                    if (s.pos.x < b.min.x)      { dist += Math<T>::Sqr(s.pos.x - b.min.x); face = true; }
                    else if (s.pos.x > b.max.x) { dist += Math<T>::Sqr(s.pos.x - b.max.x); face = true; }
                    else if (s.pos.x - b.min.x <= s.radius) face = true;
                    else if (b.max.x - s.pos.x <= s.radius) face = true;

                    if (s.pos.y < b.min.y)      { dist += Math<T>::Sqr(s.pos.y - b.min.y); face = true; }
                    else if (s.pos.y > b.max.y) { dist += Math<T>::Sqr(s.pos.y - b.max.y); face = true; }
                    else if (s.pos.y - b.min.y <= s.radius) face = true;
                    else if (b.max.y - s.pos.y <= s.radius) face = true;

                    if (s.pos.z < b.min.z)      { dist += Math<T>::Sqr(s.pos.z - b.min.z); face = true; }
                    else if (s.pos.z > b.max.z) { dist += Math<T>::Sqr(s.pos.z - b.max.z); face = true; }
                    else if (s.pos.z - b.min.z <= s.radius) face = true;
                    else if (b.max.z - s.pos.z <= s.radius) face = true;

                    return face && dist <= s.radius*s.radius;
                }
            }

        case HOLLOW:
            {
                T minX = Math<T>::Sqr(s.pos.x - b.min.x);
                T maxX = Math<T>::Sqr(s.pos.x - b.max.x);
                T minY = Math<T>::Sqr(s.pos.y - b.min.y);
                T maxY = Math<T>::Sqr(s.pos.y - b.max.y);
                T minZ = Math<T>::Sqr(s.pos.z - b.min.z);
                T maxZ = Math<T>::Sqr(s.pos.z - b.max.z);

                bool outsideX = s.pos.x < b.min.x || s.pos.x > b.max.x;
                bool outsideY = s.pos.y < b.min.y || s.pos.y > b.max.y;
                bool outsideZ = s.pos.z < b.min.z || s.pos.z > b.max.z;

                T x = Math<T>::Min(minX, maxX);
                T y = Math<T>::Min(minY, maxY);
                T z = Math<T>::Min(minZ, maxZ);

                T distMin = 0;
                if (outsideX) distMin += x;
                if (outsideY) distMin += y;
                if (outsideZ) distMin += z;

                T rSq = s.radius*s.radius;
                if (distMin > rSq)
                    return false;

                T distMax = Math<T>::Max(minX, maxX) +
                            Math<T>::Max(minY, maxY) +
                            Math<T>::Max(minZ, maxZ);

                switch (bi)
                {
                case SOLID:
                    return rSq <= distMax;

                case HOLLOW:
                    return rSq <= distMax && (outsideX || outsideY || outsideZ || x <= rSq || y <= rSq || z <= rSq);
                }
            }
        }
        return false;
    }

    template RUNEMATH_API bool Intersect(const AxisAlignedBox3<float> &b, Interior bi, const Sphere3<float> &s, Interior si);
    template RUNEMATH_API bool Intersect(const AxisAlignedBox3<double> &b, Interior bi, const Sphere3<double> &s, Interior si);

    /**
        Spocita kolizi AABB a primky
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API bool Intersect(const AxisAlignedBox3<T> &b, const Ray3<T> &ray)
    {
        Vec3<T> sides = b.GetExtents();
        Vec3<T> diff = ray.origin - b.GetCenter();

        if ((Math<T>::Abs(diff.x) > sides.x && diff.x * ray.direction.x >= 0) ||
            (Math<T>::Abs(diff.y) > sides.y && diff.y * ray.direction.y >= 0) ||
            (Math<T>::Abs(diff.z) > sides.z && diff.z * ray.direction.z >= 0)) return false;

        Vec3<T> rayDirAbs = ray.direction.GetAbs();
        Vec3<T> crossAbs = Vec3<T>::Cross(ray.direction, diff).GetAbs();

        return !(crossAbs.x > sides.y * rayDirAbs.z + sides.z * rayDirAbs.y ||
                 crossAbs.y > sides.x * rayDirAbs.z + sides.z * rayDirAbs.x ||
                 crossAbs.z > sides.x * rayDirAbs.y + sides.y * rayDirAbs.x);
    }

    template RUNEMATH_API bool Intersect(const AxisAlignedBox3<float> &b, const Ray3<float> &ray);
    template RUNEMATH_API bool Intersect(const AxisAlignedBox3<double> &b, const Ray3<double> &ray);

    /**
        Spocita kolizi dvou kouli
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API bool Intersect(const Sphere3<T> s1, const Sphere3<T> &s2)
    {
        T d = Math<T>::Sqr(s1.pos.x-s2.pos.x) + Math<T>::Sqr(s1.pos.y-s2.pos.y) + Math<T>::Sqr(s1.pos.z-s2.pos.z);
        return d <= Math<T>::Sqr(s1.radius+s2.radius);
    }

    template RUNEMATH_API bool Intersect(const Sphere3<float> s1, const Sphere3<float> &s2);
    template RUNEMATH_API bool Intersect(const Sphere3<double> s1, const Sphere3<double> &s2);

    /**
        Spocita kolizi dvou kouli, vraci navic bod dotyku kouli
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API bool Intersect(const Sphere3<T> &s1, const Sphere3<T> &s2, Vec3<T> &point)
    {
        if (Intersect(s1, s2))
        {
            Vec3<T> dir = s1.pos - s2.pos;
            dir.Normalize();
            point = s2.pos + dir * s2.radius;
            return true;
        }
        return false;
    }

    template RUNEMATH_API bool Intersect(const Sphere3<float> &s1, const Sphere3<float> &s2, Vec3<float> &point);
    template RUNEMATH_API bool Intersect(const Sphere3<double> &s1, const Sphere3<double> &s2, Vec3<double> &point);

    /**
        Spocita kolizi koule a bodu
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API bool Intersect(const Sphere3<T> &s, const Vec3<T> &point)
    {
        return Vec3<T>::GetDistanceSqr(s.pos, point) < s.radius * s.radius;
    }

    template RUNEMATH_API bool Intersect(const Sphere3<float> &s, const Vec3<float> &point);
    template RUNEMATH_API bool Intersect(const Sphere3<double> &s, const Vec3<double> &point);

    /**
        Spocita kolizi koule a primky
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API bool Intersect(const Sphere3<T> &s, const Ray3<T> &ray)
    {
        Vec3<T> rd = s.pos - ray.origin;
        T a = Vec3<T>::Dot(ray.direction, ray.direction);
        T b = Vec3<T>::Dot(rd, ray.direction) * 2;
        T c = Vec3<T>::Dot(rd, rd) - s.radius*s.radius;
        T d = Math<T>::Discriminant(a, b, c);
        return d >= 0 && b * Math<T>::Abs(b) + d > 0;
    }

    template RUNEMATH_API bool Intersect(const Sphere3<float> &s, const Ray3<float> &ray);
    template RUNEMATH_API bool Intersect(const Sphere3<double> &s, const Ray3<double> &ray);

    /**
        Spocita kolizi koule a primky, vraci navic bod dotyku a vzdalenost
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API bool Intersect(const Sphere3<T> &s, const Ray3<T> &ray, Vec3<T> &point, T &distance)
    {
        Vec3<T> rd = s.pos - ray.origin;
        T a = Vec3<T>::Dot(ray.direction, ray.direction);
        T b = Vec3<T>::Dot(rd, ray.direction) * 2;
        T c = Vec3<T>::Dot(rd, rd) - s.radius*s.radius;
        T d = Math<T>::Discriminant(a, b, c);

        if (d >= 0 && b * Math<T>::Abs(b) + d > 0)
        {
            T t = -Math<T>::QuadraticEquationRoot(a, b, d, 1);
            point = (ray.origin + ray.direction * t);
            distance = t;
            return true;
        }
        return false;
    }

    template RUNEMATH_API bool Intersect(const Sphere3<float> &s, const Ray3<float> &ray, Vec3<float> &point, float &distance);
    template RUNEMATH_API bool Intersect(const Sphere3<double> &s, const Ray3<double> &ray, Vec3<double> &point, double &distance);

    /**
        Spocita kolizi komoleho jehlanu a bodu
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API bool Intersect(const Frustum3<T> &f, const Vec3<T> &point)
    {
        for (int i = 0; i < 6; i++)
            if (f.planes[i].GetSide(point) != 1) return false;
        return true;
    }

    template RUNEMATH_API bool Intersect(const Frustum3<float> &f, const Vec3<float> &point);
    template RUNEMATH_API bool Intersect(const Frustum3<double> &f, const Vec3<double> &point);

    /**
        Spocita kolizi komoleho jehlanu a koule
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API bool Intersect(const Frustum3<T> &f, const Sphere3<T> &s)
    {
        for (int i = 0; i < 6; i++)
            if (f.planes[i].GetDistance(s.pos) < -s.radius) return false;
        return true;
    }

    template RUNEMATH_API bool Intersect(const Frustum3<float> &f, const Sphere3<float> &s);
    template RUNEMATH_API bool Intersect(const Frustum3<double> &f, const Sphere3<double> &s);

    /**
        Spocita kolizi komoleho jehlanu a AABB
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API bool Intersect(const Frustum3<T> &f, const AxisAlignedBox3<T> &box)
    {
        Vec3<T> corners[8];
        box.GetVertices(corners);
        for (int i = 0; i < 6; i++)
        {
            bool res = false;
            for (int j = 0; j < 8; j++)
                if (f.planes[i].GetSide(corners[j]) == 1)
                {
                    res = true;
                    break;
                }
            if (!res)
                return false;
        }
        return true;
    }

    template RUNEMATH_API bool Intersect(const Frustum3<float> &f, const AxisAlignedBox3<float> &box);
    template RUNEMATH_API bool Intersect(const Frustum3<double> &f, const AxisAlignedBox3<double> &box);
}
