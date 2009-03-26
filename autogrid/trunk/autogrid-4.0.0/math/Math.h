/*
    Linear Algebra / Math library

    Copyright (C) 2003-2009, Marek Olsak (maraeo@gmail.com), All Rights Reserved.
    Copyright (C) 2003-2005, Tomas Pastorek (tomas@tomaspastorek.cz), All Rights Reserved.

    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation; either version 2
    of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
*/

#pragma once

namespace Rune
{
    /**
        Definuje matematicke konstanty a funkce, vhodne pro pouziti v sablonach. Vetsina z nich je
        designovana pro typy realnych cisel.
    **************************************************************************************************/
    template<typename T>
    class Math
    {
    public:
        // Konstanty
        static T Pi()                   { return T(3.1415926535897932384626433832795); }
        static T Pi2()                  { return T(6.283185307179586476925286766559); }
        static T PiHalf()               { return T(1.5707963267948966192313216916398); }
        static T Deg()                  { return T(57.295779513082320876798154814105); }
        static T Rad()                  { return T(0.017453292519943295769236907684886); }
        static T Epsilon()              { return T(0.0001); }
        static T Epsilon2()             { return T(0.000001); }
        static T Sqrt2()                { return T(1.4142135623730950488016887242097); }
        static T Sqrt2Inv()             { return T(0.70710678118654752440084436210485); }

        // Funkce
        static T Sin(T f)               { return T(sin(f)); }
        static T SinDeg(T f)            { return Sin(f*Rad()); }
        static T Asin(T f)              { return T(asin(f)); }
        static T Cos(T f)               { return T(cos(f)); }
        static T CosDeg(T f)            { return Cos(f*Rad()); }
        static T Acos(T f)              { return (f > 1)? 0 : (f < -1)? Pi() : T(acos(f)); }
        static T Tan(T f)               { return T(tan(f)); }
        static T TanDeg(T f)            { return Tan(f*Rad()); }
        static T Atan(T f)              { return T(atan(f)); }
        static T Cotan(T f)             { return 1 / Tan(f); }
        static T CotanDeg(T f)          { return Cotan(f*Rad()); }
        static T Atan2(T f, T g)        { return T(atan2(f, g)); }

        static T Sqr(T f)               { return f*f; }
        static T Sqrt(T f)              { return T(sqrt(f)); }
        static T Cube(T f)              { return f*f*f; }

        static T Abs(T f)               { return abs(f); }
        static T Min(T f, T g)          { return (f < g)? f : g; }
        static T Max(T f, T g)          { return (f > g)? f : g; }
        static T Clamp(T f, T min, T max) { return Min(Max(f, min), max); }
        static T Saturate(T f)          { return Clamp(f, 0, 1); }
        static T Sign(T f)              { return f < 0 ? -1 : (f > 0 ?  1 : 0); }

        static T Floor(T f)             { return T(floor(f)); }
        static T Ceil(T f)              { return T(ceil(f)); }
        static T Round(T f)             { return Floor(f + T(0.5)); }
        static T Frac(T f)              { return T(fmod(f, 1)); }

        static bool SafeIsZero(T f)     { return f < Epsilon() && f > -Epsilon(); }
        static bool SafeIsEqual(T f, T g) { return f < g+Epsilon() && f > g-Epsilon(); }

        static T Discriminant(T a, T b, T c) { return b*b - 4*a*c; }
        static T QuadraticEquationRoot(T a, T b, T d, T sign) { return (-b + sign * Sqrt(d)) / (2*a); }

        static float RsqrtApprox(float a)
        {
            return union_cast<float>(0x5f3759df - (union_cast<int32>(a) >> 1));
        }

        // Commented out because it's too slow and doesn't give us any advantage over standard 1/sqrt
        // Isn't the 64bit integer arithmetic little sluggish on 32bit architecture?
        /*inline double RsqrtApprox(double a)
        {
            return union_cast<double>(0x5fe6ec85e7de30daLL - (union_cast<int64>(a) >> 1));
        }*/

#if 1
        // 1 / Sqrt(f).
        static T Rsqrt(T a)
        {
            return 1 / Sqrt(a);
        }
#else
        // Fast approximation of 1 / Sqrt(f).
        static T Rsqrt(T a, int iterations = 4)
        {
            T g = T(RsqrtApprox(float(a)));

            // Improve the accuracy using Newton's method
            T halfa = T(0.5) * a;
            for (int i = 0; i < iterations; i++)
                g *= T(1.5) - halfa*g*g;

            return g;
        }
#endif
    };

    typedef Math<float> Mathf;
    typedef Math<double> Mathd;
    typedef Math<int32> Mathi;
}
