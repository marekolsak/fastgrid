#include "../All.h"

namespace Rune
{
    /**
        Spocita bounding sphere, pouziva 3 ruzne zpusoby, vypocet neni uplne presny
    **************************************************************************************************/
    template<typename T>
    void Sphere3<T>::ApproximateBestOf3(const Vec3<T> *vertices, int count)
    {
        Sphere3<T> aabb, av, md;
        aabb.ApproximateFromAABB(vertices, count);
        av.ApproximateFromAverage(vertices, count);
        md.ApproximateFromMostDistantVerts(vertices, count);

        // Nastavi tu nejmensi
        if (aabb.radius < av.radius && aabb.radius < md.radius)
            *this = aabb;
        else if (av.radius < md.radius)
            *this = av;
        else
            *this = md;
    }

    template<typename T>
    void Sphere3<T>::ApproximateFromAABB(const Vec3<T> *vertices, int count)
    {
        pos = FindAABBCenter(vertices, count);
        radius = CalculateRadius(pos, vertices, count);
    }

    template<typename T>
    void Sphere3<T>::ApproximateFromAABB(const Sphere3<T> *spheres, int count)
    {
        // TODO: tato aproximace je hodne nepresna, pro presnejsi aproximaci by bylo lepsi vyuzit knihovnu CGAL (je pod LGPL)
        Vec3<T> *vertices = new Vec3<T>[count];
        T maxR = 0;
        for (int i = 0; i < count; i++)
        {
            vertices[i] = spheres[i].pos;
            if (maxR < spheres[i].radius)
                maxR = spheres[i].radius;
        }

        ApproximateFromAABB(vertices, count);
        delete [] vertices;
        radius += maxR;
    }

    template<typename T>
    void Sphere3<T>::ApproximateFromAverage(const Vec3<T> *vertices, int count)
    {
        pos = FindAverageCenter(vertices, count);
        radius = CalculateRadius(pos, vertices, count);
    }

    template<typename T>
    void Sphere3<T>::ApproximateFromMostDistantVerts(const Vec3<T> *vertices, int count)
    {
        pos = FindCenterFromMostDistantVerts(vertices, count);
        radius = CalculateRadius(pos, vertices, count);
    }

    template<typename T>
    Vec3<T> Sphere3<T>::FindAABBCenter(const Vec3<T> *vertices, int count)
    {
        AxisAlignedBox3<T> box;
        box.Approximate(vertices, count);
        return Vec3<T>::GetCenter(box.min, box.max);
    }

    template<typename T>
    Vec3<T> Sphere3<T>::FindAverageCenter(const Vec3<T> *vertices, int count)
    {
        Vec3<T> result = *vertices;
        const Vec3<T> *it, *last = vertices+count;
        for (it = vertices+1; it != last; ++it)
            result += *it;
        result *= 1 / T(count);
        return result;
    }

    template<typename T>
    Vec3<T> Sphere3<T>::FindCenterFromMostDistantVerts(const Vec3<T> *vertices, int count)
    {
        T maxDist = 0;
        Vec3<T> result;
        const Vec3<T> *it, *it2, *last = vertices+count;
        for (it = vertices; it != last; ++it)
            for (it2 = it; it2 != last; ++it2)
            {
                T dist = Vec3<T>::GetDistanceSqr(*it, *it2);
                if (dist > maxDist)
                {
                    maxDist = dist;
                    result = Vec3<T>::GetCenter(*it, *it2);
                }
            }
        return result;
    }

    template<typename T>
    T Sphere3<T>::CalculateRadius(const Vec3<T> &center, const Vec3<T> *vertices, int count)
    {
        T maxRadius = 0;
        const Vec3<T> *it, *last = vertices+count;
        for (it = vertices; it != last; ++it)
        {
            T radius = Vec3<T>::GetDistanceSqr(*it, center);
            if (radius > maxRadius) maxRadius = radius;
        }
        return Math<T>::Sqrt(maxRadius);
    }

    template
    class RUNEMATH_API Sphere3<float>;
    template
    class RUNEMATH_API Sphere3<double>;

    /**
        Transformuje kouli
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API Sphere3<T> operator *(const Matrix3<T> &m, const Sphere3<T> &s)
    {
        return Sphere3<T>(m * s.pos, (m * Vec3<T>(s.radius, 0, 0)).Magnitude());
    }

    template RUNEMATH_API Sphere3<float> operator *(const Matrix3<float> &m, const Sphere3<float> &s);
    template RUNEMATH_API Sphere3<double> operator *(const Matrix3<double> &m, const Sphere3<double> &s);

    /**
        Transformuje kouli transponovanou matici
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API Sphere3<T> operator *(const Sphere3<T> &s, const Matrix3<T> &m)
    {
        return Sphere3<T>(s.pos * m, (Vec3<T>(s.radius, 0, 0) * m).Magnitude());
    }
    
    template RUNEMATH_API Sphere3<float> operator *(const Sphere3<float> &s, const Matrix3<float> &m);
    template RUNEMATH_API Sphere3<double> operator *(const Sphere3<double> &s, const Matrix3<double> &m);

    /**
        Transformuje kouli
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API Sphere3<T> operator *(const Matrix4<T> &m, const Sphere3<T> &s)
    {
        return Sphere3<T>(m * Vertex3<T>(s.pos), (m * Vec3<T>(s.radius, 0, 0)).Magnitude());
    }

    template RUNEMATH_API Sphere3<float> operator *(const Matrix4<float> &m, const Sphere3<float> &s);
    template RUNEMATH_API Sphere3<double> operator *(const Matrix4<double> &m, const Sphere3<double> &s);

    /**
        Transformuje kouli transponovanou matici
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API Sphere3<T> operator *(const Sphere3<T> &s, const Matrix4<T> &m)
    {
        return Sphere3<T>(Vertex3<T>(s.pos) * m, (Vec3<T>(s.radius, 0, 0) * m).Magnitude());
    }

    template RUNEMATH_API Sphere3<float> operator *(const Sphere3<float> &s, const Matrix4<float> &m);
    template RUNEMATH_API Sphere3<double> operator *(const Sphere3<double> &s, const Matrix4<double> &m);
}
