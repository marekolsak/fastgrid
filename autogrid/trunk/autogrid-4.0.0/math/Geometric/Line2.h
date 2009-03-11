#pragma once

namespace Rune
{
    /**
        Trida usecky ve 2D
    **************************************************************************************************/
    template<typename T>
    class Line2
    {
    public:
        Vec2<T> origin;     // Stred usecky
        Ray2<T> ray;        // Primka prochazejici obema body usecky
        T extent;           // Vzdalenost stredu a konce usecky

        Line2() {}
        Line2(const Vec2<T> &p1, const Vec2<T> &p2) { SetPoints(p1, p2); }
        Vec2<T> GetPoint1() const { return origin + ray.normal.GetNormal() * extent; }
        Vec2<T> GetPoint2() const { return origin - ray.normal.GetNormal() * extent; }

        RUNEMATH_API void SetPoints(const Vec2<T> &p1, const Vec2<T> &p2);
        RUNEMATH_API int GetSide(const Vec2<T> &_point) const;
    };

    typedef Line2<float> Line2f;
    typedef Line2<double> Line2d;
}
