#pragma once

namespace Rune
{
    /**
        Ctyrrozmerny vektor
    **************************************************************************************************/
    template<typename T>
    class Vec4
    {
    public:
        T x,y,z,w;

        typedef Rune::Vec2<T> Vec2;
        typedef Rune::Vec3<T> Vec3;

        Vec4() {}
        Vec4(T f): x(f), y(f), z(f), w(f) {}
        Vec4(T X, T Y, T Z, T W): x(X), y(Y), z(Z), w(W) {}
        Vec4(T X, T Y, const Vec2 &v): x(X), y(Y), z(v.x), w(v.y) {}
        Vec4(T X, const Vec2 &v, T W): x(X), y(v.x), z(v.y), w(W) {}
        Vec4(const Vec2 &v, T Z, T W): x(v.x), y(v.y), z(Z), w(W) {}
        Vec4(const Vec2 &v1, const Vec2 &v2): x(v1.x), y(v1.y), z(v2.x), w(v2.y) {}
        Vec4(const Vec3 &v, T W): x(v.x), y(v.y), z(v.z), w(W) {}
        Vec4(T X, const Vec3 &v): x(X), y(v.x), z(v.y), w(v.z) {}

        bool operator == (const Vec4 &v) const          { return x == v.x && y == v.y && z == v.z && w == v.w; }
        bool operator != (const Vec4 &v) const          { return x != v.x && y != v.y && z != v.z && w != v.w; }
        bool operator <(const Vec4 &v) const            { return x < v.x && y < v.y && z < v.z && w < v.w; }
        bool operator >(const Vec4 &v) const            { return x > v.x && y > v.y && z > v.z && w > v.w; }
        bool operator <=(const Vec4 &v) const           { return x <= v.x && y <= v.y && z <= v.z && w <= v.w; }
        bool operator >=(const Vec4 &v) const           { return x >= v.x && y >= v.y && z >= v.z && w >= v.w; }

        bool operator ==(T f) const                     { return x == f && y == f && z == f && w == f; }
        bool operator !=(T f) const                     { return x != f && y != f && z != f && w != f; }
        bool operator <(T f) const                      { return x < f && y < f && z < f && w < f; }
        bool operator >(T f) const                      { return x > f && y > f && z > f && w > f; }
        bool operator <=(T f) const                     { return x <= f && y <= f && z <= f && w <= f; }
        bool operator >=(T f) const                     { return x >= f && y >= f && z >= f && w >= f; }

        void operator = (const Vec4 &v)                 { x = v.x; y = v.y; z = v.z; w = v.w; }
        void operator +=(const Vec4 &v)                 { x += v.x; y += v.y; z += v.z; w += v.w; }
        void operator +=(T f)                           { x += f; y += f; z += f; w += f; }
        void operator -=(const Vec4 &v)                 { x -= v.x; y -= v.y; z -= v.z; w -= v.w; }
        void operator -=(T f)                           { x -= f; y -= f; z -= f; w -= f; }
        void operator *=(const Vec4 &v)                 { x *= v.x; y *= v.y; z *= v.z; w *= v.w; }
        void operator *=(T f)                           { x *= f; y *= f; z *= f; w *= f; }
        void operator /=(T f)                           { f = 1/f; x *= f; y *= f; z *= f; w *= f; }

        Vec4 operator +(const Vec4 &v) const            { return Vec4(x+v.x, y+v.y, z+v.z, w+v.w); }
        Vec4 operator +(T f) const                      { return Vec4(x+f, y+f, z+f, w+f); }
        Vec4 operator -(const Vec4 &v) const            { return Vec4(x-v.x, y-v.y, z-v.z, w-v.w); }
        Vec4 operator -(T f) const                      { return Vec4(x-f, y-f, z-f, w-f); }
        Vec4 operator *(const Vec4 &v) const            { return Vec4(x*v.x, y*v.y, z*v.z, w*v.w); }
        Vec4 operator *(T f) const                      { return Vec4(x*f, y*f, z*f, w*f); }
        Vec4 operator /(T f) const                      { f = 1/f; return Vec4(x*f, y*f, z*f, w*f); }
        Vec4 operator -() const                         { return Vec4(-x, -y, -z, -w); }
        T& operator[] (int i)                           { return (&x)[i]; }
        T operator[] (int i) const                      { return (&x)[i]; }

        operator T* ()                                  { return &x; }
        operator const T* () const                      { return &x; }

        Vec3 GetVec3() const                            { return Vec3(x, y, z); }
        Vec4 GetAbs() const                             { return Vec4(Math<T>::Abs(x), Math<T>::Abs(y), Math<T>::Abs(z), Math<T>::Abs(w)); }
        Vec4 GetNormalized() const                      { T f = RMagnitude(); return Vec4(x*f, y*f, z*f, w*f); }
        T Magnitude() const                             { return Math<T>::Sqrt(MagnitudeSqr()); }
        T RMagnitude() const                            { return Math<T>::Rsqrt(MagnitudeSqr()); }
        T MagnitudeSqr() const                          { return x*x + y*y + z*z + w*w; }
        void Normalize()                                { operator *=(RMagnitude()); }
        bool IsEmpty() const                            { return x == 0 && y == 0 && z == 0 && w == 0; }
        void PackTo01()                                 { operator *=(0.5f); operator +=(0.5f); }
        void PackBack()                                 { operator -=(0.5f); operator *=(2); }
        void Set(const T X, const T Y, const T Z, const T W) { x = X; y = Y; z = Z; w = W; }
        void Inverse()                                  { x = 1/x; y = 1/y; z = 1/z; w = 1/w; }
        Vec4 Frac() const                               { return Vec4(Math<T>::Frac(x), Math<T>::Frac(y), Math<T>::Frac(z), Math<T>::Frac(w)); }

        RUNEMATH_API bool SafeIsEqual(const Vec4 &v) const;
        RUNEMATH_API static T Dot(const Vec4 &v1, const Vec4 &v2);
        RUNEMATH_API static Vec4 Cross(const Vec4 &a, const Vec4 &b, const Vec4 &c);
    };

    template<typename T>
    std::ostream& operator <<(std::ostream &stream, const Vec4<T> &v)
    {
        return stream << ("(" + tostrf(v.x) + ", " + tostrf(v.y) + ", " + tostrf(v.z) + ", " + tostrf(v.w) + ")");
    }

    template<typename T>
    std::istream& operator >>(std::istream &stream, Vec4<T> &v)
    {
        char tmp;
        return stream >> tmp >> v.x >> tmp >> v.y >> tmp >> v.z >> tmp >> v.w >> tmp;
    }

    typedef Vec4<float> Vec4f;
    typedef Vec4<double> Vec4d;
    typedef Vec4<int32> Vec4i;
    typedef Vec4<bool> Vec4b;
}
