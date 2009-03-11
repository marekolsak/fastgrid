#pragma once

namespace Rune
{
    template<typename T>
    RUNEMATH_API bool Intersect(const Ray2<T> &ray1, const Ray2<T> &ray2);

    template<typename T>
    RUNEMATH_API bool Intersect(const Ray2<T> &ray1, const Ray2<T> &ray2, Vec2<T> &point);

    template<typename T>
    RUNEMATH_API bool Intersect(const Line2<T> &line, const Ray2<T> &ray, Vec2<T> &point);

    template<typename T>
    RUNEMATH_API bool Intersect(const Line2<T> &line, const Ray2<T> &ray);

    template<typename T>
    RUNEMATH_API bool Intersect(const Line2<T> &line1, const Line2<T> &line2, Vec2<T> &point);

    template<typename T>
    RUNEMATH_API bool Intersect(const Line2<T> &line1, const Line2<T> &line2);

    template<typename T>
    RUNEMATH_API bool Intersect(const Circle2<T> &s1, const Circle2<T> &s2);

    template<typename T>
    RUNEMATH_API bool Intersect(const Circle2<T> &s1, const Circle2<T> &s2, Vec2<T> &point);

    template<typename T>
    RUNEMATH_API bool Intersect(const Circle2<T> &s, const Vec2<T> &point);

    template<typename T>
    RUNEMATH_API bool Intersect(const Rect2<T> &r, const Vec2<T> &point);
}
