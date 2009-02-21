#pragma once

// Define fixed-size integer types
#if defined(_MSC_VER)
    #define TYPEDEF_FIXED_INT(bits) typedef signed __int##bits int##bits; typedef unsigned __int##bits uint##bits
#else
    #define TYPEDEF_FIXED_INT(bits) typedef int##bits##_t int##bits; typedef uint##bits##_t uint##bits
#endif

TYPEDEF_FIXED_INT(8);   // int8, uint8
TYPEDEF_FIXED_INT(16);  // int16, uint16
TYPEDEF_FIXED_INT(32);  // int32, uint32
TYPEDEF_FIXED_INT(64);  // int64, uint64

// Helper function for data type reinterpretation
template<typename T1, typename T2>
inline T1 union_cast(T2 a)
{
    union
    {
        T2 a;
        T1 b;
    } myUnion;
    myUnion.a = a;
    return myUnion.b;
}

inline float rsqrtApprox(float a)
{
    return union_cast<float>(0x5f3759df - (union_cast<int32>(a) >> 1));
}

// Commented out because it's too slow and doesn't give us any advantage over standard 1/sqrt
// Isn't the 64bit integer arithmetic little sluggish on 32bit architecture?
/*inline double rsqrtApprox(double a)
{
    return union_cast<double>(0x5fe6ec85e7de30daLL - (union_cast<int64>(a) >> 1));
}*/

// Fast approximation of 1/SQRT(a)
template<typename T>
inline T rsqrt(T a)
{
    //return 1 / sqrt(a);

    T g = rsqrtApprox(float(a));

    // Improve the accuracy using Newton's method
    T halfa = 0.5 * a;
    g *= 1.5 - halfa*g*g;
    g *= 1.5 - halfa*g*g;

    return g;
}

// Square
template<typename T>
inline T sq(T a)
{
    return a*a;
}

// Cube
template<typename T>
inline T cube(T a)
{
    return a*a*a;
}

// Hypotenuse squared
template<typename T>
inline T hypotenuseSq(T x, T y, T z)
{
    return x*x + y*y + z*z;
}

// Hypotenuse
template<typename T>
inline T hypotenuse(T x, T y, T z)
{
    return 1 / rsqrt(hypotenuseSq<T>(x, y, z));
}

// Reciprocal of hypotenuse
template<typename T>
inline T hypotenuseInv(T x, T y, T z)
{
    return rsqrt<T>(hypotenuseSq<T>(x, y, z));
}

// Normalize the vector
template<typename T>
inline void normalize(T &x, T &y, T &z)
{
    T hypInv = hypotenuseInv<T>(x, y, z);
    x *= hypInv;
    y *= hypInv;
    z *= hypInv;
}

// we do not want to have a redefinition of the following macro max,min
#ifdef _WIN32
    #undef min
    #undef max
#endif

// Max
template<typename T>
inline T max(T x, T y)
{
    return x > y ? x : y;
}

// Min
template<typename T>
inline T min(T x, T y)
{
    return x < y ? x : y;
}

// Make sure the value is between min and max
template<typename T>
inline T clamp(T x, T min, T max)
{
    return x < min ? min :
           x > max ? max :
           x;
}
