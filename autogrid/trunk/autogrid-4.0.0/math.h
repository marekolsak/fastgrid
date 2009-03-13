#pragma once
// we do not want to have a redefinition of the following macro max,min
#ifdef _WIN32
    #undef min
    #undef max
#endif

#undef X
#undef Y
#undef Z
#include "math/All.h"
using namespace Rune;
#define X 0
#define Y 1
#define Z 2

// Length squared
template<typename T>
inline T lengthSquared(T x, T y, T z)
{
    return x*x + y*y + z*z;
}

// Round to 3 decimal places
template<typename T>
inline T round3dp(T x)
{
    return Math<T>::Round(x * T(1000)) / T(1000);
}

// Dot product
template<typename T>
inline T dotProduct(const T *a, const T *b)
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

// Cross product
template<typename T>
void crossProduct(T *r, const T *a, const T *b)
{
    r[0] = a[1] * b[2] - a[2] * b[1];
    r[1] = a[2] * b[0] - a[0] * b[2];
    r[2] = a[0] * b[1] - a[1] * b[0];
}

// Get angle between 2 vectors
template<typename T>
inline T angle(const T *a, const T *b)
{
    return Math<T>::Acos(dotProduct(a, b));
}
