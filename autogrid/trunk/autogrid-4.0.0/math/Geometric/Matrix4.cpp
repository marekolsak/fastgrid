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
        Prevede quaternion na matici
    **************************************************************************************************/
    template<typename T>
    void Matrix4<T>::operator =(const Quaternion<T> &q)
    {
        T xx = q.v.x*q.v.x, xy = q.v.x*q.v.y, xz = q.v.x*q.v.z, xw = q.v.x*q.v.w, yy = q.v.y*q.v.y;
        T yz = q.v.y*q.v.z, yw = q.v.y*q.v.w, zz = q.v.z*q.v.z, zw = q.v.z*q.v.w;
        RUNE_SET_M4(m,  1-2*(yy+zz),    2*(xy+zw),      2*(xz-yw),      0,
                        2*(xy-zw),      1-2*(xx+zz),    2*(yz+xw),      0,
                        2*(xz+yw),      2*(yz-xw),      1-2*(xx+yy),    0,
                        0,              0,              0,              1);
    }

    /**
        Prevede rotacni cast matice do quaternionu a vrati jej
    **************************************************************************************************/
    template<typename T>
    Matrix4<T>::operator Quaternion<T>() const
    {
        T t = 1 + m[0] + m[5] + m[10];
        if (t > Math<T>::Epsilon())
        {
            T s = Math<T>::Sqrt(t)*2, sr = 1/s;
            return Quaternion<T>((m[6]-m[9]) * sr, (m[8]-m[2]) * sr, (m[1]-m[4]) * sr, T(0.25) * s);
        }
        if (m[0] > m[5] && m[0] > m[10])
        {
            T s = Math<T>::Sqrt(1 + m[0] - m[5] - m[10]) * 2, sr = 1/s;
            return Quaternion<T>(T(0.25) * s, (m[1] + m[4]) * sr, (m[8] + m[2]) * sr, (m[6] - m[9]) * sr);
        }
        if (m[5] > m[10])
        {
            T s = Math<T>::Sqrt(1 + m[5] - m[0] - m[10]) * 2, sr = 1/s;
            return Quaternion<T>((m[1] + m[4]) * sr, T(0.25) * s, (m[6] + m[9]) * sr, (m[8] - m[2]) * sr);
        }
        T s = Math<T>::Sqrt(1 + m[10] - m[0] - m[5]) * 2, sr = 1/s;
        return Quaternion<T>((m[8] + m[2]) * sr, (m[6] + m[9]) * sr, T(0.25) * s, (m[1] - m[4]) * sr);
    }

    /**
        Spocita determinant
    **************************************************************************************************/
    template<typename T>
    T Matrix4<T>::Determinant() const
    {
        return (m[0]  * m[5]  - m[1]  * m[4]) * (m[10] * m[15] - m[11] * m[14]) -
               (m[0]  * m[6]  - m[2]  * m[4]) * (m[9]  * m[15] - m[11] * m[13]) +
               (m[0]  * m[7]  - m[3]  * m[4]) * (m[9]  * m[14] - m[10] * m[13]) +
               (m[1]  * m[6]  - m[2]  * m[5]) * (m[8]  * m[15] - m[11] * m[12]) -
               (m[1]  * m[7]  - m[3]  * m[5]) * (m[8]  * m[14] - m[10] * m[12]) +
               (m[2]  * m[7]  - m[3]  * m[6]) * (m[8]  * m[13] - m[9]  * m[12]);
    }

    /**
        Nastavi perspektivni matici pro projekci
    **************************************************************************************************/
    template<typename T>
    void Matrix4<T>::SetPerspective(T fov, T aspect, T near, T far)
    {
        T h = Math<T>::CotanDeg(fov * T(0.5)), w = aspect * h, nearFar = near - far;
        RUNE_SET_M4(m, w, 0, 0, 0, 0, h, 0, 0, 0, 0, far/nearFar, -1, 0, 0, near*far/nearFar, 0);
    }

    /**
        Nastavi ortho matici pro projekci
    **************************************************************************************************/
    template<typename T>
    void Matrix4<T>::SetOrtho(T left, T top, T right, T bottom, T near, T far)
    {
        T w = right - left, h = top - bottom, farNear = far + near;
        RUNE_SET_M4(m,  2 / w,                0,                      0,                          0,
                        0,                    2 / h,                  0,                          0,
                        0,                    0,                      2 / farNear,                0,
                        -(right + left) / w,  -(top + bottom) / h,    (far - near) / -farNear,    1);
    }

    /**
        Spocita soucin dvou matic
    **************************************************************************************************/
    template<typename T>
    void Matrix4<T>::SetProduct(const Matrix4<T> &m1, const Matrix4<T> &m2)
    {
        m[0]  = m1[0] * m2[0]  + m1[4] * m2[1]  + m1[8]  * m2[2]  + m1[12] * m2[3];
        m[1]  = m1[1] * m2[0]  + m1[5] * m2[1]  + m1[9]  * m2[2]  + m1[13] * m2[3];
        m[2]  = m1[2] * m2[0]  + m1[6] * m2[1]  + m1[10] * m2[2]  + m1[14] * m2[3];
        m[3]  = m1[3] * m2[0]  + m1[7] * m2[1]  + m1[11] * m2[2]  + m1[15] * m2[3];
        m[4]  = m1[0] * m2[4]  + m1[4] * m2[5]  + m1[8]  * m2[6]  + m1[12] * m2[7];
        m[5]  = m1[1] * m2[4]  + m1[5] * m2[5]  + m1[9]  * m2[6]  + m1[13] * m2[7];
        m[6]  = m1[2] * m2[4]  + m1[6] * m2[5]  + m1[10] * m2[6]  + m1[14] * m2[7];
        m[7]  = m1[3] * m2[4]  + m1[7] * m2[5]  + m1[11] * m2[6]  + m1[15] * m2[7];
        m[8]  = m1[0] * m2[8]  + m1[4] * m2[9]  + m1[8]  * m2[10] + m1[12] * m2[11];
        m[9]  = m1[1] * m2[8]  + m1[5] * m2[9]  + m1[9]  * m2[10] + m1[13] * m2[11];
        m[10] = m1[2] * m2[8]  + m1[6] * m2[9]  + m1[10] * m2[10] + m1[14] * m2[11];
        m[11] = m1[3] * m2[8]  + m1[7] * m2[9]  + m1[11] * m2[10] + m1[15] * m2[11];
        m[12] = m1[0] * m2[12] + m1[4] * m2[13] + m1[8]  * m2[14] + m1[12] * m2[15];
        m[13] = m1[1] * m2[12] + m1[5] * m2[13] + m1[9]  * m2[14] + m1[13] * m2[15];
        m[14] = m1[2] * m2[12] + m1[6] * m2[13] + m1[10] * m2[14] + m1[14] * m2[15];
        m[15] = m1[3] * m2[12] + m1[7] * m2[13] + m1[11] * m2[14] + m1[15] * m2[15];
    }

    /**
        Spocita rotacni matici podle rotacni osy a uhlu
    **************************************************************************************************/
    template<typename T>
    void Matrix4<T>::SetRotated(T x, T y, T z, T angle)
    {
        T s = Math<T>::SinDeg(angle), c = Math<T>::CosDeg(angle), t = 1-c;
        T xx = x*x, xy = x*y, xz = x*z, yy = y*y, yz = y*z, zz = z*z, xs = x*s, ys = y*s, zs = z*s;
        RUNE_SET_M4(m,  xx * t + c,     xy * t + zs,    xz * t - ys,    0,
                        xy * t - zs,    yy * t + c,     yz * t + xs,    0,
                        xz * t + ys,    yz * t - xs,    zz * t + c,     0,
                        0,              0,              0,              1);
    }

    /**
        Spocita rotacni matici podle Eulerovych uhlu v poradi pitch-roll-yaw (= XYZ)
    **************************************************************************************************/
    template<typename T>
    void Matrix4<T>::SetPitchRollYaw(T anglex, T angley, T anglez)
    {
        T cx = Math<T>::CosDeg(anglex), sx = Math<T>::SinDeg(anglex), cy = Math<T>::CosDeg(angley);
        T sy = Math<T>::SinDeg(angley), cz = Math<T>::CosDeg(anglez), sz = Math<T>::SinDeg(anglez);
        RUNE_SET_M4(m,  cy * cz,                    cy * sz,                    -sy,        0,
                        sx * sy * cz - cx * sz,     sx * sy * sz + cx * cz,     sx * cy,    0,
                        cx * sy * cz + sx * sz,     cx * sy * sz - sx * cz,     cx * cy,    0,
                        0,                          0,                          0,          1);
    }

    /**
        Spocita rotacni matici, rotacni osa je X
    **************************************************************************************************/
    template<typename T>
    void Matrix4<T>::SetRotatedX(T angle)
    {
        T s = Math<T>::SinDeg(angle), c = Math<T>::CosDeg(angle);
        RUNE_SET_M4(m, 1, 0, 0, 0, 0, c, s, 0, 0, -s, c, 0, 0, 0, 0, 1);
    }

    /**
        Spocita rotacni matici, rotacni osa je Y
    **************************************************************************************************/
    template<typename T>
    void Matrix4<T>::SetRotatedY(T angle)
    {
        T s = Math<T>::SinDeg(angle), c = Math<T>::CosDeg(angle);
        RUNE_SET_M4(m, c, 0, -s, 0, 0, 1, 0, 0, s, 0, c, 0, 0, 0, 0, 1);
    }

    /**
        Spocita rotacni matici, rotacni osa je Z
    **************************************************************************************************/
    template<typename T>
    void Matrix4<T>::SetRotatedZ(T angle)
    {
        T s = Math<T>::SinDeg(angle), c = Math<T>::CosDeg(angle);
        RUNE_SET_M4(m, c, s, 0, 0, -s, c, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1);
    }

    /**
        Spocita matici pro transformaci do souradnicoveho systemu trojuhelniku (normala shodna s osou Z)
    **************************************************************************************************/
    template<typename T>
    void Matrix4<T>::SetMatrixToFaceSpace(const Vec3<T> &v0, const Vec3<T> &v1, const Vec3<T> &v2)
    {
        Vec3<T> normal(Vec3<T>::CalculateNormal(v0, v1, v2));
        Vec3<T> axis(Vec3<T>::Cross(normal, Vec3<T>(0, 0, 1)));
        axis.Normalize();
        SetRotated(axis, Math<T>::Acos(Vec3<T>::Dot(normal, Vec3<T>(0, 0, 1))) * Math<T>::Deg());
    }

    /**
        Spocita matici pohledu
    **************************************************************************************************/
    template<typename T>
    void Matrix4<T>::SetLookAt(const Vec3<T> &eye, const Vec3<T> &center, const Vec3<T> &up)
    {
        Vec3<T> f(center-eye);
        f.Normalize();
        Vec3<T> s = Vec3<T>::Cross(f, up.GetNormalized());
        s.Normalize();
        Vec3<T> u = Vec3<T>::Cross(s, f);
        RUNE_SET_M4(m, s.x, u.x, -f.x, 0, s.y, u.y, -f.y, 0, s.z, u.z, -f.z, 0, 0, 0, 0, 1);
        Translate(-eye);
    }

    /**
        Vynasobi matici, ktera ma na diagonale vektor a vsude jinde nuly
    **************************************************************************************************/
    template<typename T>
    void Matrix4<T>::MultipleDiagonal(const Vec4<T> &v)
    {
        RUNE_SET_M4(m,  m[0]  * v.x, m[1]  * v.y, m[2]  * v.z, m[3]  * v.w,
                        m[4]  * v.x, m[5]  * v.y, m[6]  * v.z, m[7]  * v.w,
                        m[8]  * v.x, m[9]  * v.y, m[10] * v.z, m[11] * v.w,
                        m[12] * v.x, m[13] * v.y, m[14] * v.z, m[15] * v.w);
    }

    /**
        Z druhe strany vynasobi matici, ktera ma na diagonale vektor a vsude jinde nuly
    **************************************************************************************************/
    template<typename T>
    void Matrix4<T>::MultipleDiagonalInv(const Vec4<T> &v)
    {
        RUNE_SET_M4(m,  m[0]  * v.x, m[1]  * v.x, m[2]  * v.x, m[3]  * v.x,
                        m[4]  * v.y, m[5]  * v.y, m[6]  * v.y, m[7]  * v.y,
                        m[8]  * v.z, m[9]  * v.z, m[10] * v.z, m[11] * v.z,
                        m[12] * v.w, m[13] * v.w, m[14] * v.w, m[15] * v.w);
    }

    /**
        Spocita inverzni matici
    **************************************************************************************************/
    template<typename T>
    void Matrix4<T>::Inverse()
    {
        T f00 = m[0]  * m[5]  - m[1]  * m[4];
        T f01 = m[0]  * m[6]  - m[2]  * m[4];
        T f02 = m[0]  * m[7]  - m[3]  * m[4];
        T f03 = m[1]  * m[6]  - m[2]  * m[5];
        T f04 = m[1]  * m[7]  - m[3]  * m[5];
        T f05 = m[2]  * m[7]  - m[3]  * m[6];
        T f10 = m[8]  * m[13] - m[9]  * m[12];
        T f11 = m[8]  * m[14] - m[10] * m[12];
        T f12 = m[8]  * m[15] - m[11] * m[12];
        T f13 = m[9]  * m[14] - m[10] * m[13];
        T f14 = m[9]  * m[15] - m[11] * m[13];
        T f15 = m[10] * m[15] - m[11] * m[14];

        T determinant = f00*f15 - f01*f14 + f02*f13 + f03*f12 - f04*f11 + f05*f10;
        if (Math<T>::SafeIsZero(determinant)) return;

        Matrix4<T> tmp;
        tmp[0]  = + m[5]  * f15 - m[6]  * f14 + m[7]  * f13;
        tmp[1]  = - m[1]  * f15 + m[2]  * f14 - m[3]  * f13;
        tmp[2]  = + m[13] * f05 - m[14] * f04 + m[15] * f03;
        tmp[3]  = - m[9]  * f05 + m[10] * f04 - m[11] * f03;
        tmp[4]  = - m[4]  * f15 + m[6]  * f12 - m[7]  * f11;
        tmp[5]  = + m[0]  * f15 - m[2]  * f12 + m[3]  * f11;
        tmp[6]  = - m[12] * f05 + m[14] * f02 - m[15] * f01;
        tmp[7]  = + m[8]  * f05 - m[10] * f02 + m[11] * f01;
        tmp[8]  = + m[4]  * f14 - m[5]  * f12 + m[7]  * f10;
        tmp[9]  = - m[0]  * f14 + m[1]  * f12 - m[3]  * f10;
        tmp[10] = + m[12] * f04 - m[13] * f02 + m[15] * f00;
        tmp[11] = - m[8]  * f04 + m[9]  * f02 - m[11] * f00;
        tmp[12] = - m[4]  * f13 + m[5]  * f11 - m[6]  * f10;
        tmp[13] = + m[0]  * f13 - m[1]  * f11 + m[2]  * f10;
        tmp[14] = - m[12] * f03 + m[13] * f01 - m[14] * f00;
        tmp[15] = + m[8]  * f03 - m[9]  * f01 + m[10] * f00;

        *this = tmp;
        MultipleScalar(1 / determinant);
    }

    /**
        Spocita inverzni matici rychlejsim zpusobem, ale musi byt ortogonalni (pouze posun a rotace)
    **************************************************************************************************/
    template<typename T>
    void Matrix4<T>::FastInverse()
    {
        Matrix4<T> temp;
        RUNE_SET_M4(temp,   m[0], m[4], m[8],  0,
                            m[1], m[5], m[9],  0,
                            m[2], m[6], m[10], 0, // Transponovani rotace
                            -m[12]*m[0] - m[13]*m[1] - m[14]*m[2],    // Prepocitani posunu
                            -m[12]*m[4] - m[13]*m[5] - m[14]*m[6],
                            -m[12]*m[8] - m[13]*m[9] - m[14]*m[10], 1);
        *this = temp;
    }

    /**
        Spocita transponovanou matici
    **************************************************************************************************/
    template<typename T>
    void Matrix4<T>::Transpose()
    {
        T tmp;
        tmp = m[1];  m[1]  = m[4];  m[4]  = tmp;
        tmp = m[2];  m[2]  = m[8];  m[8]  = tmp;
        tmp = m[3];  m[3]  = m[12]; m[12] = tmp;
        tmp = m[6];  m[6]  = m[9];  m[9]  = tmp;
        tmp = m[7];  m[7]  = m[13]; m[13] = tmp;
        tmp = m[11]; m[11] = m[14]; m[14] = tmp;
    }

    template
    class RUNEMATH_API Matrix4<float>;
    template
    class RUNEMATH_API Matrix4<double>;

    /**
        Transformuje prostorovy vektor (je ekvivalentni transformaci vektoru (x, y, z, 0))
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API Vec3<T> operator *(const Matrix4<T> &m, const Vec3<T> &v)
    {
        return Vec3<T>(v.x * m[0] + v.y * m[4] + v.z * m[8],
                       v.x * m[1] + v.y * m[5] + v.z * m[9],
                       v.x * m[2] + v.y * m[6] + v.z * m[10]);
    }

    template RUNEMATH_API Vec3<float> operator *(const Matrix4<float> &m, const Vec3<float> &v);
    template RUNEMATH_API Vec3<double> operator *(const Matrix4<double> &m, const Vec3<double> &v);

    /**
        Transformuje prostorovy vektor transponovanou matici (je ekvivalentni transformaci vektoru (x, y, z, 0))
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API Vec3<T> operator *(const Vec3<T> &v, const Matrix4<T> &m)
    {
        return Vec3<T>(v.x * m[0] + v.y * m[1] + v.z * m[2],
                       v.x * m[4] + v.y * m[5] + v.z * m[6],
                       v.x * m[8] + v.y * m[9] + v.z * m[10]);
    }

    template RUNEMATH_API Vec3<float> operator *(const Vec3<float> &v, const Matrix4<float> &m);
    template RUNEMATH_API Vec3<double> operator *(const Vec3<double> &v, const Matrix4<double> &m);

    /**
        Transformuje prostorovy vektor v homogennich souradnicich
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API Vec4<T> operator *(const Matrix4<T> &m, const Vec4<T> &v)
    {
        return Vec4<T>(v.x * m[0] + v.y * m[4] + v.z * m[8]  + v.w * m[12],
                       v.x * m[1] + v.y * m[5] + v.z * m[9]  + v.w * m[13],
                       v.x * m[2] + v.y * m[6] + v.z * m[10] + v.w * m[14],
                       v.x * m[3] + v.y * m[7] + v.z * m[11] + v.w * m[15]);
    }

    template RUNEMATH_API Vec4<float> operator *(const Matrix4<float> &m, const Vec4<float> &v);
    template RUNEMATH_API Vec4<double> operator *(const Matrix4<double> &m, const Vec4<double> &v);

    /**
        Transformuje prostorovy vektor v homogennich souradnicich transponovanou matici
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API Vec4<T> operator *(const Vec4<T> &v, const Matrix4<T> &m)
    {
        return Vec4<T>(v.x * m[0] + v.y * m[1] + v.z * m[2]  + v.w * m[3],
                       v.x * m[4] + v.y * m[5] + v.z * m[6]  + v.w * m[7],
                       v.x * m[8] + v.y * m[9] + v.z * m[10] + v.w * m[11],
                       v.x * m[12] + v.y * m[13] + v.z * m[14] + v.w * m[15]);
    }

    template RUNEMATH_API Vec4<float> operator *(const Vec4<float> &v, const Matrix4<float> &m);
    template RUNEMATH_API Vec4<double> operator *(const Vec4<double> &v, const Matrix4<double> &m);
}
