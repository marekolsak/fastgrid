#pragma once

namespace Rune
{
    /**
        Trida komoleho jehlanu
    **************************************************************************************************/
    template<typename T>
    class Frustum3
    {
    public:
        Plane3<T> planes[6];

        /**
            Nastavi 6 ploch definujici komoly jehlan pohledu kamery za predpokladu, ze "m"
            obsahuje projekcni matici, pripadne muze taky obsahovat soucin projekcni a modelovaci matice
        **************************************************************************************************/
        RUNEMATH_API void Calculate(const Matrix4<T> &m);
    };

    template<typename T>
    RUNEMATH_API Frustum3<T> operator *(const Matrix3<T> &m, const Frustum3<T> &f);
    template<typename T>
    RUNEMATH_API Frustum3<T> operator *(const Frustum3<T> &f, const Matrix3<T> &m);
    template<typename T>
    RUNEMATH_API Frustum3<T> operator *(const Matrix4<T> &m, const Frustum3<T> &f);
    template<typename T>
    RUNEMATH_API Frustum3<T> operator *(const Frustum3<T> &f, const Matrix4<T> &m);

    typedef Frustum3<float> Frustum3f;
    typedef Frustum3<double> Frustum3d;
}
