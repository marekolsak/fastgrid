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

#include "../All.h"

namespace Rune
{
    /**
        Vynasobi 2 quaterniony mezi sebou
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API Quaternion<T> Quaternion<T>::operator *(const Quaternion<T> &q) const
    {
        return Quaternion<T>(v.w*q.v.x + v.x*q.v.w + v.y*q.v.z - v.z*q.v.y,
                             v.w*q.v.y + v.y*q.v.w + v.z*q.v.x - v.x*q.v.z,
                             v.w*q.v.z + v.z*q.v.w + v.x*q.v.y - v.y*q.v.x,
                             v.w*q.v.w - v.x*q.v.x - v.y*q.v.y - v.z*q.v.z);
    }

    template RUNEMATH_API Quaternion<float> Quaternion<float>::operator *(const Quaternion<float> &q) const;
    template RUNEMATH_API Quaternion<double> Quaternion<double>::operator *(const Quaternion<double> &q) const;

    /**
        Vraci rotacni osu jako (x,y,z) a uhel ve stupnich jako w
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API Vec4<T> Quaternion<T>::GetRotation() const
    {
        T rm = v.MagnitudeInv(), w = v.w * rm, s = Math<T>::Rsqrt(1 - w*w) * rm;
        return Vec4<T>(v.x*s, v.y*s, v.z*s, Math<T>::Acos(w) * Math<T>::Deg() * 2);
    }

    template RUNEMATH_API Vec4<float> Quaternion<float>::GetRotation() const;
    template RUNEMATH_API Vec4<double> Quaternion<double>::GetRotation() const;

    /**
        Spocita rotaci z rotacni osy a uhlu
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API void Quaternion<T>::SetRotated(T x, T y, T z, T angle)
    {
        T halfAngle = angle * T(0.5), s = Math<T>::SinDeg(halfAngle);
        v.Set(x*s, y*s, z*s, Math<T>::CosDeg(halfAngle));
    }

    template RUNEMATH_API void Quaternion<float>::SetRotated(float x, float y, float z, float angle);
    template RUNEMATH_API void Quaternion<double>::SetRotated(double x, double y, double z, double angle);

    /**
        Spocita rotaci kolem osy X
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API void Quaternion<T>::SetRotatedX(T angle)
    {
        T halfAngle = angle * T(0.5);
        v.Set(Math<T>::SinDeg(halfAngle), 0, 0, Math<T>::CosDeg(halfAngle));
    }

    template RUNEMATH_API void Quaternion<float>::SetRotatedX(float angle);
    template RUNEMATH_API void Quaternion<double>::SetRotatedX(double angle);

    /**
        Spocita rotaci kolem osy Y
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API void Quaternion<T>::SetRotatedY(T angle)
    {
        T halfAngle = angle * T(0.5);
        v.Set(0, Math<T>::SinDeg(halfAngle), 0, Math<T>::CosDeg(halfAngle));
    }

    template RUNEMATH_API void Quaternion<float>::SetRotatedY(float angle);
    template RUNEMATH_API void Quaternion<double>::SetRotatedY(double angle);

    /**
        Spocita rotaci kolem osy Z
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API void Quaternion<T>::SetRotatedZ(T angle)
    {
        T halfAngle = angle * T(0.5);
        v.Set(0, 0, Math<T>::SinDeg(halfAngle), Math<T>::CosDeg(halfAngle));
    }

    template RUNEMATH_API void Quaternion<float>::SetRotatedZ(float angle);
    template RUNEMATH_API void Quaternion<double>::SetRotatedZ(double angle);

    /**
        Spocita rotaci z eulerovych uhlu v poradi yaw-pitch-roll (ZXY)
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API void Quaternion<T>::SetYawPitchRoll(T anglez, T anglex, T angley)
    {
        SetRotatedZ(anglez);
        RotateX(anglex);
        RotateY(angley);
    }
    
    template RUNEMATH_API void Quaternion<float>::SetYawPitchRoll(float anglez, float anglex, float angley);
    template RUNEMATH_API void Quaternion<double>::SetYawPitchRoll(double anglez, double anglex, double angley);

    /**
        Spocita takovou rotaci, ktera je potreba k orotovani vektoru start do vektoru end
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API void Quaternion<T>::SetRotated(const Vec3<T> &start, const Vec3<T> &end)
    {
        T halfAngle = Math<T>::Acos(Vec3<T>::Dot(start, end)) * T(0.5);
        Vec3<T> cross = Vec3<T>::Cross(start, end) * Math<T>::Sin(halfAngle);
        v.Set(cross.x, cross.y, cross.z, Math<T>::Cos(halfAngle));
    }

    template RUNEMATH_API void Quaternion<float>::SetRotated(const Vec3<float> &start, const Vec3<float> &end);
    template RUNEMATH_API void Quaternion<double>::SetRotated(const Vec3<double> &start, const Vec3<double> &end);

    /**
        Spocita sferickou linearni interpolaci mezi dvema quaterniony, t musi byt v intervalu <0, 1>
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API void Quaternion<T>::SetSlerp(const Quaternion<T> &q1, const Quaternion<T> &q2, T t)
    {
        T angle = Math<T>::Acos(Vec4<T>::Dot(q1.v, q2.v));
        if (Math<T>::SafeIsZero(angle)) operator =(q1);
        else operator =((q1*Math<T>::Sin((1-t)*angle) + q2*Math<T>::Sin(t*angle)) / Math<T>::Sin(angle));
    }

    template RUNEMATH_API void Quaternion<float>::SetSlerp(const Quaternion<float> &q1, const Quaternion<float> &q2, float t);
    template RUNEMATH_API void Quaternion<double>::SetSlerp(const Quaternion<double> &q1, const Quaternion<double> &q2, double t);

    /**
        Transformuje vektor
        Vzhledem k narocnosti vypoctu je mnohem rychlejsi prevest quaternion na matici a nasledne s ni
        transformovat ten vektor. Tato metoda tu je uvedena jen pro demonstraci, ze to jde i primo.
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API Vec3<T> operator *(const Quaternion<T> &q, const Vec3<T> &v)
    {
        Vec4<T> v4 = (q * Quaternion<T>(v.x, v.y, v.z, 1) * q.GetInverted()).v;
        return v4.xyz / v4.w;
    }

    template RUNEMATH_API Vec3<float> operator *(const Quaternion<float> &q, const Vec3<float> &v);
    template RUNEMATH_API Vec3<double> operator *(const Quaternion<double> &q, const Vec3<double> &v);

    /**
        Transformuje vektor, vypocet je shodny s transformovanim transponovanou matici
        Vzhledem k narocnosti vypoctu je mnohem rychlejsi prevest quaternion na matici a nasledne s ni
        transformovat ten vektor. Tato metoda tu je uvedena jen pro demonstraci, ze to jde i primo.
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API Vec3<T> operator *(const Vec3<T> &v, const Quaternion<T> &q)
    {
        Vec4<T> v4 = (q.GetInverted() * Quaternion<T>(v.x, v.y, v.z, 1) * q).v;
        return v4.xyz / v4.w;
    }

    template RUNEMATH_API Vec3<float> operator *(const Vec3<float> &v, const Quaternion<float> &q);
    template RUNEMATH_API Vec3<double> operator *(const Vec3<double> &v, const Quaternion<double> &q);
}
