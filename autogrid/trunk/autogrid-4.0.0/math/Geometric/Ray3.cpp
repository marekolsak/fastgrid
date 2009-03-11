#include "../All.h"

namespace Rune
{
    template<typename T>
    RUNEMATH_API void Ray3<T>::SetSecondPoint(Vec3<T> &p)
    {
        direction = (p - origin).GetNormalized();
    }

    template RUNEMATH_API void Ray3<float>::SetSecondPoint(Vec3<float> &p);
    template RUNEMATH_API void Ray3<double>::SetSecondPoint(Vec3<double> &p);

    /**
        Transformuje primku v prostoru
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API Ray3<T> operator *(const Matrix3<T> &m, const Ray3<T> &r)
    {
        return Ray3<T>(m * r.origin, m * r.direction);
    }

    template RUNEMATH_API Ray3<float> operator *(const Matrix3<float> &m, const Ray3<float> &r);
    template RUNEMATH_API Ray3<double> operator *(const Matrix3<double> &m, const Ray3<double> &r);

    /**
        Transformuje primku v prostoru transponovanou matici
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API Ray3<T> operator *(const Ray3<T> &r, const Matrix3<T> &m)
    {
        return Ray3<T>(r.origin * m, r.direction * m);
    }

    template RUNEMATH_API Ray3<float> operator *(const Ray3<float> &r, const Matrix3<float> &m);
    template RUNEMATH_API Ray3<double> operator *(const Ray3<double> &r, const Matrix3<double> &m);

    /**
        Transformuje primku v prostoru
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API Ray3<T> operator *(const Matrix4<T> &m, const Ray3<T> &r)
    {
        return Ray3<T>(m * Vertex3<T>(r.origin), m * r.direction);
    }

    template RUNEMATH_API Ray3<float> operator *(const Matrix4<float> &m, const Ray3<float> &r);
    template RUNEMATH_API Ray3<double> operator *(const Matrix4<double> &m, const Ray3<double> &r);

    /**
        Transformuje primku v prostoru transponovanou matici
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API Ray3<T> operator *(const Ray3<T> &r, const Matrix4<T> &m)
    {
        return Ray3<T>(Vertex3<T>(r.origin) * m, r.direction * m);
    }

    template RUNEMATH_API Ray3<float> operator *(const Ray3<float> &r, const Matrix4<float> &m);
    template RUNEMATH_API Ray3<double> operator *(const Ray3<double> &r, const Matrix4<double> &m);
}
