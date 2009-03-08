#pragma once

// Define fixed-sized integer types
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
#if 0
    //return 1 / sqrt(a);

#else
    T g = T(rsqrtApprox(float(a)));

    // Improve the accuracy using Newton's method
    T halfa = T(0.5) * a;
    for (int i = 0; i < 4; i++)
        g *= T(1.5) - halfa*g*g;

    return g;
#endif
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
inline T lengthSquared(T x, T y, T z)
{
    return x*x + y*y + z*z;
}

// Length of a vector
/*template<typename T>
inline T length(T x, T y, T z)
{
    return sqrt(lengthSquared<T>(x, y, z));
}*/

// Reciprocal of length
template<typename T>
inline T lengthInv(T x, T y, T z)
{
    return rsqrt<T>(lengthSquared<T>(x, y, z));
}

// Normalize the vector
template<typename T>
inline void normalize(T &x, T &y, T &z)
{
    T hypInv = lengthInv<T>(x, y, z);
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

// Acos clamped to [-1, 1]
template<typename T>
inline T acos_clamped(T x)
{
    return acos(clamp(x, T(-1), T(1)));
}

// round() is a C99 function and not universally available
// Required to round consistently on different platforms
#if !defined(HAVE_ROUND)
template<typename T>
inline T round(T x)
{
    return floor(x + T(0.5));
}
#endif

// Round to 3 decimal places
template<typename T>
inline T round3dp(T x)
{
    return round(x * T(1000)) / T(1000);
}

// Subtract two vectors
template<typename T>
inline void subtractVectors(T *r, const T *a, const T *b)
{
    r[0] = a[0] - b[0];
    r[1] = a[1] - b[1];
    r[2] = a[2] - b[2];
}

// Dot product
template<typename T>
inline T dotProduct(const T *a, const T *b)
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

// Multiply a vector by a scalar
template<typename T>
inline void scalarProduct(T *r, const T *a, T s)
{
    r[0] = a[0] * s;
    r[1] = a[1] * s;
    r[2] = a[2] * s;
}

// Get angle between 2 vectors
template<typename T>
inline T angle(const T *a, const T *b)
{
    return acos_clamped(dotProduct(a, b));
}

// Cross product
template<typename T>
void crossProduct(T *r, const T *a, const T *b)
{
    r[0] = a[1] * b[2] - a[2] * b[1];
    r[1] = a[2] * b[0] - a[0] * b[2];
    r[2] = a[0] * b[1] - a[1] * b[0];
}
