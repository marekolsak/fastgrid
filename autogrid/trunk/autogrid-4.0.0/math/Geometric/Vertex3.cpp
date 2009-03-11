#include "../All.h"

namespace Rune
{
    /**
        Transformuje vertex
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API Vertex3<T> operator *(const Matrix3<T> &m, const Vertex3<T> &v)
    {
        return m * Vec3<T>(v);
    }

    template RUNEMATH_API Vertex3<float> operator *(const Matrix3<float> &m, const Vertex3<float> &v);
    template RUNEMATH_API Vertex3<double> operator *(const Matrix3<double> &m, const Vertex3<double> &v);

    /**
        Transformuje vertex transponovanou matici
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API Vertex3<T> operator *(const Vertex3<T> &v, const Matrix3<T> &m)
    {
        return Vec3<T>(v) * m;
    }

    template RUNEMATH_API Vertex3<float> operator *(const Vertex3<float> &v, const Matrix3<float> &m);
    template RUNEMATH_API Vertex3<double> operator *(const Vertex3<double> &v, const Matrix3<double> &m);

    /**
        Transformuje vertex (je ekvivalentni transformaci vektoru (x, y, z, 1))
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API Vertex3<T> operator *(const Matrix4<T> &m, const Vertex3<T> &v)
    {
        return m * Vec3<T>(v) + Vec3<T>(m[12], m[13], m[14]);
    }

    template RUNEMATH_API Vertex3<float> operator *(const Matrix4<float> &m, const Vertex3<float> &v);
    template RUNEMATH_API Vertex3<double> operator *(const Matrix4<double> &m, const Vertex3<double> &v);

    /**
        Transformuje vertex transponovanou matici (je ekvivalentni transformaci vektoru (x, y, z, 1))
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API Vertex3<T> operator *(const Vertex3<T> &v, const Matrix4<T> &m)
    {
        return Vec3<T>(v) * m + Vec3<T>(m[3], m[7], m[11]);
    }

    template RUNEMATH_API Vertex3<float> operator *(const Vertex3<float> &v, const Matrix4<float> &m);
    template RUNEMATH_API Vertex3<double> operator *(const Vertex3<double> &v, const Matrix4<double> &m);
}
