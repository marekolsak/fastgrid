#pragma once

namespace Rune
{
    /**
        Dvourozmerny vektor
    **************************************************************************************************/
    template<typename T>
    class Vec2
    {
    public:
        T x,y;

        Vec2() {}
        Vec2(T f): x(f), y(f) {}
        Vec2(T X, T Y): x(X), y(Y) {}

        bool operator ==(const Vec2 &v) const           { return x == v.x && y == v.y; }
        bool operator !=(const Vec2 &v) const           { return x != v.x && y != v.y; }
        bool operator <(const Vec2 &v) const            { return x < v.x && y < v.y; }
        bool operator >(const Vec2 &v) const            { return x > v.x && y > v.y; }
        bool operator <=(const Vec2 &v) const           { return x <= v.x && y <= v.y; }
        bool operator >=(const Vec2 &v) const           { return x >= v.x && y >= v.y; }

        bool operator ==(T f) const                     { return x == f && y == f; }
        bool operator !=(T f) const                     { return x != f && y != f; }
        bool operator <(T f) const                      { return x < f && y < f; }
        bool operator >(T f) const                      { return x > f && y > f; }
        bool operator <=(T f) const                     { return x <= f && y <= f; }
        bool operator >=(T f) const                     { return x >= f && y >= f; }

        void operator =(const Vec2 &v)                  { x = v.x; y = v.y; }
        void operator +=(const Vec2 &v)                 { x += v.x; y += v.y; }
        void operator +=(T f)                           { x += f; y += f; }
        void operator -=(const Vec2 &v)                 { x -= v.x; y -= v.y; }
        void operator -=(T f)                           { x -= f; y -= f; }
        void operator *=(T f)                           { x *= f; y *= f; }
        void operator /=(T f)                           { x /= f; y /= f; }

        Vec2 operator +(const Vec2 &v) const            { return Vec2(x+v.x, y+v.y); }
        Vec2 operator +(T f) const                      { return Vec2(x+f, y+f); }
        Vec2 operator -(const Vec2 &v) const            { return Vec2(x-v.x, y-v.y); }
        Vec2 operator -(T f) const                      { return Vec2(x-f, y-f); }
        Vec2 operator *(T f) const                      { return Vec2(x*f, y*f); }
        Vec2 operator *(const Vec2 &v) const            { return Vec2(x*v.x, y*v.y); }
        Vec2 operator /(T f) const                      { return Vec2(x/f, y/f); }
        Vec2 operator -() const                         { return Vec2(-x, -y); }
        T& operator [](int i)                           { return (&x)[i]; }
        T operator [](int i) const                      { return (&x)[i]; }

        bool SafeIsEqual(const Vec2 &v) const           { return Math<T>::SafeIsEqual(x, v.x) && Math<T>::SafeIsEqual(y, v.y); }
        Vec2 GetAbs() const                             { return Vec2(Math<T>::Abs(x), Math<T>::Abs(y)); }
        Vec2 GetNormalized() const                      { T f = RMagnitude(); return Vec2(x*f, y*f); }
        Vec2 GetNormal() const                          { return Vec2(-y, x); }
        T Magnitude() const                             { return Math<T>::Sqrt(MagnitudeSqr()); }
        T RMagnitude() const                            { return Math<T>::Rsqrt(MagnitudeSqr()); }
        T MagnitudeSqr() const                          { return x*x + y*y; }
        void Normalize()                                { operator *=(RMagnitude()); }
        bool IsEmpty() const                            { return x == 0 && y == 0; }
        void PackTo01()                                 { operator *=(T(0.5)); operator +=(T(0.5)); }
        void PackBack()                                 { operator -=(T(0.5)); operator *=(2); }
        void Set(T X, T Y)                              { x = X; y = Y; }
        void Inverse()                                  { x = 1/x; y = 1/y; }
        Vec2 Frac() const                               { return Vec2(Math<T>::Frac(x), Math<T>::Frac(y)); }

        operator T* ()                                  { return &x; }
        operator const T* () const                      { return &x; }

        static T Dot(const Vec2 &v1, const Vec2 &v2)                { return v1.x*v2.x + v1.y*v2.y; }
        static T GetDistanceSqr(const Vec2 &v1, const Vec2 &v2)     { return (v1 - v2).MagnitudeSqr(); }
        static T GetDistance(const Vec2 &v1, const Vec2 &v2)        { return Math<T>::Sqrt(GetDistanceSqr(v1, v2)); }
        static T GetAngle(const Vec2 &a, const Vec2 &b)             { return Math<T>::Acos(Dot(a, b)*a.RMagnitude()*b.RMagnitude()); }
        static Vec2 GetCenter(const Vec2 &min, const Vec2 &max)     { return (min + max) / 2; }
        static Vec2 Reflect(const Vec2 &impact, const Vec2 &normal) { return impact - 2*Dot(normal, impact)*normal; }
    };

    template<typename T>
    std::ostream& operator <<(std::ostream &stream, const Vec2<T> &v)
    {
        return stream << ("(" + tostrf(v.x) + ", " + tostrf(v.y) + ")");
    }

    template<typename T>
    std::istream& operator >>(std::istream &stream, Vec2<T> &v)
    {
        char tmp;
        return stream >> tmp >> v.x >> tmp >> v.y >> tmp;
    }

    typedef Vec2<float> Vec2f;
    typedef Vec2<double> Vec2d;
    typedef Vec2<int32> Vec2i;
    typedef Vec2<bool> Vec2b;
}
