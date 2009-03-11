#pragma once

namespace Rune
{
    /**
        Trida koule ve 3D
    **************************************************************************************************/
    template<typename T>
    class RUNEMATH_API Sphere3
    {
    public:
        Vec3<T> pos;
        T radius;

        Sphere3() {}
        Sphere3(const Vec3<T> &position, T sphereradius): pos(position), radius(sphereradius) {}

        void ApproximateBestOf3(const Vec3<T> *vertices, int count);
        void ApproximateFromAABB(const Vec3<T> *vertices, int count);
        void ApproximateFromAABB(const Sphere3<T> *spheres, int count);
        void ApproximateFromAverage(const Vec3<T> *vertices, int count);
        void ApproximateFromMostDistantVerts(const Vec3<T> *vertices, int count);

    private:
        static Vec3<T> FindAABBCenter(const Vec3<T> *vertices, int count);
        static Vec3<T> FindAverageCenter(const Vec3<T> *vertices, int count);
        static Vec3<T> FindCenterFromMostDistantVerts(const Vec3<T> *vertices, int count);
        static T CalculateRadius(const Vec3<T> &center, const Vec3<T> *vertices, int count);
    };

    template<typename T>
    RUNEMATH_API Sphere3<T> operator *(const Matrix3<T> &m, const Sphere3<T> &s);
    template<typename T>
    RUNEMATH_API Sphere3<T> operator *(const Sphere3<T> &s, const Matrix3<T> &m);
    template<typename T>
    RUNEMATH_API Sphere3<T> operator *(const Matrix4<T> &m, const Sphere3<T> &s);
    template<typename T>
    RUNEMATH_API Sphere3<T> operator *(const Sphere3<T> &s, const Matrix4<T> &m);

    typedef Sphere3<float> Sphere3f;
    typedef Sphere3<double> Sphere3d;
}
