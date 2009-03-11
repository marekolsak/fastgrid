#pragma once

namespace Rune
{
    enum Interior
    {
        SOLID = 0,
        HOLLOW
    };

    template<typename T>
    RUNEMATH_API bool Intersect(const Plane3<T> &p, const Ray3<T> &ray, Vec3<T> &point, T &tparam);
    template<typename T>
    RUNEMATH_API bool Intersect(const AxisAlignedBox3<T> &b, const Vec3<T> &p);
    template<typename T>
    RUNEMATH_API bool Intersect(const AxisAlignedBox3<T> &b1, const AxisAlignedBox3<T> &b2);
    template<typename T>
    RUNEMATH_API bool Intersect(const AxisAlignedBox3<T> &b, const Sphere3<T> &s);
    template<typename T>
    RUNEMATH_API bool Intersect(const AxisAlignedBox3<T> &b, Interior bi, const Sphere3<T> &s, Interior si);
    template<typename T>
    RUNEMATH_API bool Intersect(const AxisAlignedBox3<T> &b, const Ray3<T> &ray);
    template<typename T>
    RUNEMATH_API bool Intersect(const Sphere3<T> s1, const Sphere3<T> &s2);
    template<typename T>
    RUNEMATH_API bool Intersect(const Sphere3<T> &s1, const Sphere3<T> &s2, Vec3<T> &point);
    template<typename T>
    RUNEMATH_API bool Intersect(const Sphere3<T> &s, const Vec3<T> &point);
    template<typename T>
    RUNEMATH_API bool Intersect(const Sphere3<T> &s, const Ray3<T> &ray);
    template<typename T>
    RUNEMATH_API bool Intersect(const Sphere3<T> &s, const Ray3<T> &ray, Vec3<T> &point, T &distance);
    template<typename T>
    RUNEMATH_API bool Intersect(const Frustum3<T> &f, const Vec3<T> &point);
    template<typename T>
    RUNEMATH_API bool Intersect(const Frustum3<T> &f, const Sphere3<T> &s);
    template<typename T>
    RUNEMATH_API bool Intersect(const Frustum3<T> &f, const AxisAlignedBox3<T> &box);
}
