#pragma once

namespace Rune
{
    /**
        Trida roviny ve 3D
    **************************************************************************************************/
    template<typename T>
    class Plane3
    {
    public:
        Vec3<T> normal;
        T d;

        Plane3() {}
        Plane3(const Vec3<T> &in_normal, const Vec3<T> &point): normal(in_normal), d(-Vec3<T>::Dot(normal, point)) {}
        Plane3(const Vec3<T> &in_normal, T D): normal(in_normal), d(D) {}
        Plane3(T a, T b, T c, T D): normal(a, b, c), d(D) {}
        Plane3(const Vec3<T> &p1, const Vec3<T> &p2, const Vec3<T> &p3): normal(Vec3<T>::GetNormal(p1, p2, p3)), d(-Vec3<T>::Dot(normal, p1)) {}

        void operator *=(T f)                               { normal *= f; d *= f; }
        Plane3 operator *(T f) const                        { return Plane3(normal * f, d * f); }
        bool operator ==(const Plane3 &p) const             { return normal == p.normal && d == p.d; }
        bool operator !=(const Plane3 &p) const             { return normal != p.normal && d != p.d; }

        Vec4<T> &GetVec4()                                  { return *reinterpret_cast<Vec4<T>*>(&normal.x); }
        const Vec4<T> &GetVec4() const                      { return *reinterpret_cast<const Vec4<T>*>(&normal.x); }
        void Set(T a, T b, T c, T D)                        { normal.Set(a, b, c); d = D; }
        void Set(Vec3<T> &v, T D)                           { normal = v; d = D; }
        void Normalize()                                    { *this *= normal.RMagnitude(); }
        Plane3 GetNormalized() const                        { return *this * normal.RMagnitude(); }
        T GetDistance(const Vec3<T> &point) const           { return Vec3<T>::Dot(normal, point) + d; }
        Vec3<T> GetMirroredPoint(const Vec3<T> &point) const{ return point + normal*GetDistance(point)*-2; }
        Vec3<T> GetMirroredVector(const Vec3<T> &vec) const { return vec + normal*Vec3<T>::Dot(normal, vec)*-2; }
        Vec3<T> GetPoint() const                            { return normal * -d; }
        Vec3<T> GetNearestPoint(const Vec3<T> &point) const { return point + normal*GetDistance(point)*-1; }

        RUNEMATH_API int GetSide(const Vec3<T> &point) const;
        RUNEMATH_API int GetSide(const Sphere3<T> &sphere) const;
        RUNEMATH_API int GetSide(const AxisAlignedBox3<T> &box) const;
    };

    template<typename T>
    RUNEMATH_API Plane3<T> operator *(const Matrix3<T> &m, const Plane3<T> &p);
    template<typename T>
    RUNEMATH_API Plane3<T> operator *(const Plane3<T> &p, const Matrix3<T> &m);
    template<typename T>
    RUNEMATH_API Plane3<T> operator *(const Matrix4<T> &m, const Plane3<T> &p);
    template<typename T>
    RUNEMATH_API Plane3<T> operator *(const Plane3<T> &p, const Matrix4<T> &m);

    template<typename T>
    std::ostream& operator <<(std::ostream &stream, const Plane3<T> &p)
    {
        return operator <<(stream, p.GetVec4());
    }

    template<typename T>
    std::istream& operator >>(std::istream &stream, Plane3<T> &p)
    {
        return stream >> p.GetVec4();
    }

    typedef Plane3<float> Plane3f;
    typedef Plane3<double> Plane3d;
}
