#pragma once

namespace Rune
{
    /**
        Trida primky ve 3D
    **************************************************************************************************/
    template<typename T>
    class Ray3
    {
    public:
        Vec3<T> origin;
        Vec3<T> direction;

        Ray3() {}
        Ray3(const Vec3<T> &Origin, const Vec3<T> &Direction): origin(Origin), direction(Direction) {}

        RUNEMATH_API void SetSecondPoint(Vec3<T> &p);
    };

    template<typename T>
    RUNEMATH_API Ray3<T> operator *(const Matrix3<T> &m, const Ray3<T> &r);
    template<typename T>
    RUNEMATH_API Ray3<T> operator *(const Ray3<T> &r, const Matrix3<T> &m);
    template<typename T>
    RUNEMATH_API Ray3<T> operator *(const Matrix4<T> &m, const Ray3<T> &r);
    template<typename T>
    RUNEMATH_API Ray3<T> operator *(const Ray3<T> &r, const Matrix4<T> &m);

    typedef Ray3<float> Ray3f;
    typedef Ray3<double> Ray3d;
}
