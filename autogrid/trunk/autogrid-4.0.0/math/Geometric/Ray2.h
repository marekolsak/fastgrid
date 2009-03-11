#pragma once

namespace Rune
{
    /**
        Trida primky ve 2D
    **************************************************************************************************/
    template<typename T>
    class Ray2
    {
    public:
        Vec2<T> normal;
        T c;

        Ray2() {}
        Ray2(const Vec2<T> &Normal, const Vec2<T> &point): normal(Normal), c(-Vec2<T>::Dot(Normal, point)) {}
        Ray2(const Vec2<T> &Normal, T C): normal(Normal), c(C) {}

        void operator *=(T f)                                   { normal *= f; c *= f; }
        Ray2 operator *(T f) const                              { return Ray2(normal * f, c * f); }

        void Set(const Vec2<T> &Normal, const Vec2<T> &point)   { normal = Normal; c = -Vec2<T>::Dot(Normal, point); }
        void Set(const Vec2<T> &Normal, T C)                    { normal = Normal; c = C; }
        void Normalize()                                        { *this *= normal.RMagnitude(); }
        Ray2 GetNormalized() const                              { return *this * normal.RMagnitude(); }
        T GetDistance(const Vec2<T> &point) const               { return Vec2<T>::Dot(normal, point) + c; }
        Vec2<T> GetMirroredPoint(const Vec2<T> &point) const    { return point + normal*GetDistance(point)*-2; }
        Vec2<T> GetMirroredVector(const Vec2<T> &vec) const     { return vec + normal*Vec2<T>::Dot(normal, vec)*-2; }
        Vec2<T> GetPoint() const                                { return normal * -c; }
        Vec2<T> GetNearestPoint(const Vec2<T> &point) const     { return point + normal*GetDistance(point)*-1; }

        RUNEMATH_API int GetSide(const Vec2<T> &point) const;
        RUNEMATH_API int GetSide(const Circle2<T> &circle) const;
    };

    typedef Ray2<float> Ray2f;
    typedef Ray2<double> Ray2d;
}
