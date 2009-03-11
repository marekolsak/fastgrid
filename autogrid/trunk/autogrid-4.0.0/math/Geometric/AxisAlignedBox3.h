#pragma once

namespace Rune
{
    /**
        Trida kvadru ve 3D, ktery je rovnobezny s osami (AABB = axis-aligned bounding box)
    **************************************************************************************************/
    template<typename T>
    class AxisAlignedBox3
    {
    public:
        Vec3<T> min, max;

        AxisAlignedBox3() {}
        AxisAlignedBox3(const Vec3<T> &Min, const Vec3<T> &Max): min(Min), max(Max) {}

        Vec3<T> GetCenter() const   { return Vec3<T>::GetCenter(min, max); }
        Vec3<T> GetExtents() const  { return (max - min) * T(0.5); }

        RUNEMATH_API void GetVertices(Vec3<T> *vertices) const;
        RUNEMATH_API void SetCenterAndExtents(const Vec3<T> &center, const Vec3<T> &extents);
        RUNEMATH_API void Approximate(const Vec3<T> *vertices, int count);
        RUNEMATH_API void Approximate(const AxisAlignedBox3 *boxes, int count);
    };

    template<typename T>
    RUNEMATH_API AxisAlignedBox3<T> operator *(const Matrix3<T> &m, const AxisAlignedBox3<T> &b);
    template<typename T>
    RUNEMATH_API AxisAlignedBox3<T> operator *(const AxisAlignedBox3<T> &b, const Matrix3<T> &m);
    template<typename T>
    RUNEMATH_API AxisAlignedBox3<T> operator *(const Matrix4<T> &m, const AxisAlignedBox3<T> &b);
    template<typename T>
    RUNEMATH_API AxisAlignedBox3<T> operator *(const AxisAlignedBox3<T> &b, const Matrix4<T> &m);

    typedef AxisAlignedBox3<float> AxisAlignedBox3f;
    typedef AxisAlignedBox3<double> AxisAlignedBox3d;
}
