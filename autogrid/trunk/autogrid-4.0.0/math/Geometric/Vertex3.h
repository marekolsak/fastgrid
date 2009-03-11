#pragma once

namespace Rune
{
    /**
        Bod v prostoru, obalujici trida pro typove odliseni nasobeni vektoru a vertexu matici
    **************************************************************************************************/
    template<typename T>
    class Vertex3 : public Vec3<T>
    {
    public:
        typedef Rune::Vec3<T> Vec3;

        Vertex3() {}
        Vertex3(const Vec3 &V): Vec3(V) {}
    };

    template<typename T>
    RUNEMATH_API Vertex3<T> operator *(const Matrix3<T> &m, const Vertex3<T> &v);
    template<typename T>
    RUNEMATH_API Vertex3<T> operator *(const Vertex3<T> &v, const Matrix3<T> &m);
    template<typename T>
    RUNEMATH_API Vertex3<T> operator *(const Matrix4<T> &m, const Vertex3<T> &v);
    template<typename T>
    RUNEMATH_API Vertex3<T> operator *(const Vertex3<T> &v, const Matrix4<T> &m);

    typedef Vertex3<float> Vertex3f;
    typedef Vertex3<double> Vertex3d;
    typedef Vertex3<int32> Vertex3i;
}
