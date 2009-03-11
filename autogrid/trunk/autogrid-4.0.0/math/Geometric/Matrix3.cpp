#include "../All.h"

namespace Rune
{
    /**
        Prevede quaternion na matici
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API void Matrix3<T>::operator =(const Quaternion<T> &q)
    {
        T xx = q.v.x*q.v.x, xy = q.v.x*q.v.y, xz = q.v.x*q.v.z, xw = q.v.x*q.v.w, yy = q.v.y*q.v.y;
        T yz = q.v.y*q.v.z, yw = q.v.y*q.v.w, zz = q.v.z*q.v.z, zw = q.v.z*q.v.w;
        RUNE_SET_M3(m,  1-2*(yy+zz),    2*(xy+zw),      2*(xz-yw),
                        2*(xy-zw),      1-2*(xx+zz),    2*(yz+xw),
                        2*(xz+yw),      2*(yz-xw),      1-2*(xx+yy));
    }

    template RUNEMATH_API void Matrix3<float>::operator =(const Quaternion<float> &q);
    template RUNEMATH_API void Matrix3<double>::operator =(const Quaternion<double> &q);

    /**
        Prevede matici do quaternionu a vrati jej
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API Matrix3<T>::operator Quaternion<T>()
    {
        T t = 1 + m[0] + m[4] + m[8];
        if (t > Math<T>::Epsilon())
        {
            T s = Math<T>::Sqrt(t)*2, sr = 1/s;
            return Quaternion<T>((m[5]-m[7]) * sr, (m[6]-m[2]) * sr, (m[1]-m[3]) * sr, T(0.25) * s);
        }
        if (m[0] > m[4] && m[0] > m[8])
        {
            T s = Math<T>::Sqrt(1 + m[0] - m[4] - m[8]) * 2, sr = 1/s;
            return Quaternion<T>(T(0.25) * s, (m[1] + m[3]) * sr, (m[6] + m[2]) * sr, (m[5] - m[7]) * sr);
        }
        if (m[4] > m[8])
        {
            T s = Math<T>::Sqrt(1 + m[4] - m[0] - m[8]) * 2, sr = 1/s;
            return Quaternion<T>((m[1] + m[3]) * sr, T(0.25) * s, (m[5] + m[7]) * sr, (m[6] - m[2]) * sr);
        }
        T s = Math<T>::Sqrt(1 + m[8] - m[0] - m[4]) * 2, sr = 1/s;
        return Quaternion<T>((m[6] + m[2]) * sr, (m[5] + m[7]) * sr, T(0.25) * s, (m[1] - m[3]) * sr);
    }

    template RUNEMATH_API Matrix3<float>::operator Quaternion<float>();
    template RUNEMATH_API Matrix3<double>::operator Quaternion<double>();

    /**
        Spocita determinant
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API T Matrix3<T>::Determinant() const
    {
        return m[0] * (m[4] * m[8] - m[5] * m[7]) +
               m[1] * (m[5] * m[6] - m[3] * m[8]) +
               m[2] * (m[3] * m[7] - m[4] * m[6]);
    }

    template RUNEMATH_API float Matrix3<float>::Determinant() const;
    template RUNEMATH_API double Matrix3<double>::Determinant() const;

    /**
        Spocita soucin dvou matic
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API void Matrix3<T>::SetProduct(const Matrix3<T> &m1, const Matrix3<T> &m2)
    {
        m[0] = m1[0] * m2[0] + m1[3] * m2[1] + m1[6] * m2[2];
        m[1] = m1[1] * m2[0] + m1[4] * m2[1] + m1[7] * m2[2];
        m[2] = m1[2] * m2[0] + m1[5] * m2[1] + m1[8] * m2[2];
        m[3] = m1[0] * m2[3] + m1[3] * m2[4] + m1[6] * m2[5];
        m[4] = m1[1] * m2[3] + m1[4] * m2[4] + m1[7] * m2[5];
        m[5] = m1[2] * m2[3] + m1[5] * m2[4] + m1[8] * m2[5];
        m[6] = m1[0] * m2[6] + m1[3] * m2[7] + m1[6] * m2[8];
        m[7] = m1[1] * m2[6] + m1[4] * m2[7] + m1[7] * m2[8];
        m[8] = m1[2] * m2[6] + m1[5] * m2[7] + m1[8] * m2[8];
    }

    template RUNEMATH_API void Matrix3<float>::SetProduct(const Matrix3<float> &m1, const Matrix3<float> &m2);
    template RUNEMATH_API void Matrix3<double>::SetProduct(const Matrix3<double> &m1, const Matrix3<double> &m2);

    /**
        Spocita soucin dvou matic, pricemz u druhe se zada akorat diagonala (byva vetsinou scale)
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API void Matrix3<T>::SetProductDiagonal(const Matrix3<T> &m1, const Vec3<T> &v)
    {
        RUNE_SET_M3(m,  m1[0] * v.x, m1[1] * v.y, m1[2] * v.z,
                        m1[3] * v.x, m1[4] * v.y, m1[5] * v.z,
                        m1[6] * v.x, m1[7] * v.y, m1[8] * v.z);
    }

    template RUNEMATH_API void Matrix3<float>::SetProductDiagonal(const Matrix3<float> &m1, const Vec3<float> &v);
    template RUNEMATH_API void Matrix3<double>::SetProductDiagonal(const Matrix3<double> &m1, const Vec3<double> &v);

    /**
        Spocita soucin dvou matic, pricemz u druhe se zada akorat diagonala, ktera se vynasobi z druhe strany (byva vetsinou scale)
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API void Matrix3<T>::SetProductDiagonalInv(const Matrix3<T> &m1, const Vec3<T> &v)
    {
        RUNE_SET_M3(m,  m1[0] * v.x, m1[1] * v.x, m1[2] * v.x,
                        m1[3] * v.y, m1[4] * v.y, m1[5] * v.y,
                        m1[6] * v.z, m1[7] * v.z, m1[8] * v.z);
    }

    template RUNEMATH_API void Matrix3<float>::SetProductDiagonalInv(const Matrix3<float> &m1, const Vec3<float> &v);
    template RUNEMATH_API void Matrix3<double>::SetProductDiagonalInv(const Matrix3<double> &m1, const Vec3<double> &v);

    /**
        Spocita soucin dvou matic floatem (byva vetsinou scale stejny pro vsechny osy)
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API void Matrix3<T>::SetProductScalar(const Matrix3<T> &m1, T f)
    {
        RUNE_SET_M3(m,  m1[0] * f, m1[1] * f, m1[2] * f,
                        m1[3] * f, m1[4] * f, m1[5] * f,
                        m1[6] * f, m1[7] * f, m1[8] * f);
    }

    template RUNEMATH_API void Matrix3<float>::SetProductScalar(const Matrix3<float> &m1, float f);
    template RUNEMATH_API void Matrix3<double>::SetProductScalar(const Matrix3<double> &m1, double f);

    /**
        Spocita rotacni matici podle rotacni osy a uhlu
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API void Matrix3<T>::SetRotated(T x, T y, T z, T angle)
    {
        T s = Math<T>::SinDeg(angle), c = Math<T>::CosDeg(angle), t = 1-c;
        T xx = x*x, xy = x*y, xz = x*z, yy = y*y, yz = y*z, zz = z*z, xs = x*s, ys = y*s, zs = z*s;
        RUNE_SET_M3(m,  xx * t + c,     xy * t + zs,    xz * t - ys,
                        xy * t - zs,    yy * t + c,     yz * t + xs,
                        xz * t + ys,    yz * t - xs,    zz * t + c);
    }

    template RUNEMATH_API void Matrix3<float>::SetRotated(float x, float y, float z, float angle);
    template RUNEMATH_API void Matrix3<double>::SetRotated(double x, double y, double z, double angle);

    /**
        Spocita rotacni matici podle Eulerovych uhlu v poradu pitch-roll-yaw (= XYZ)
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API void Matrix3<T>::SetPitchRollYaw(T anglex, T angley, T anglez)
    {
        T cx = Math<T>::CosDeg(anglex), sx = Math<T>::SinDeg(anglex), cy = Math<T>::CosDeg(angley);
        T sy = Math<T>::SinDeg(angley), cz = Math<T>::CosDeg(anglez), sz = Math<T>::SinDeg(anglez);
        RUNE_SET_M3(m,   cy * cz,                    cy * sz,                    -sy,
                    sx * sy * cz - cx * sz,     sx * sy * sz + cx * cz,     sx * cy,
                    cx * sy * cz + sx * sz,     cx * sy * sz - sx * cz,     cx * cy);
    }

    template RUNEMATH_API void Matrix3<float>::SetPitchRollYaw(float anglex, float angley, float anglez);
    template RUNEMATH_API void Matrix3<double>::SetPitchRollYaw(double anglex, double angley, double anglez);

    /**
        Spocita rotacni matici, rotacni osa je X
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API void Matrix3<T>::SetRotatedX(T angle)
    {
        T s = Math<T>::SinDeg(angle), c = Math<T>::CosDeg(angle);
        RUNE_SET_M3(m, 1, 0, 0, 0, c, s, 0, -s, c);
    }

    template RUNEMATH_API void Matrix3<float>::SetRotatedX(float angle);
    template RUNEMATH_API void Matrix3<double>::SetRotatedX(double angle);

    /**
        Spocita rotacni matici, rotacni osa je Y
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API void Matrix3<T>::SetRotatedY(T angle)
    {
        T s = Math<T>::SinDeg(angle), c = Math<T>::CosDeg(angle);
        RUNE_SET_M3(m, c, 0, -s, 0, 1, 0, s, 0, c);
    }

    template RUNEMATH_API void Matrix3<float>::SetRotatedY(float angle);
    template RUNEMATH_API void Matrix3<double>::SetRotatedY(double angle);

    /**
        Spocita rotacni matici, rotacni osa je Z
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API void Matrix3<T>::SetRotatedZ(T angle)
    {
        T s = Math<T>::SinDeg(angle), c = Math<T>::CosDeg(angle);
        RUNE_SET_M3(m, c, s, 0, -s, c, 0, 0, 0, 1);
    }

    template RUNEMATH_API void Matrix3<float>::SetRotatedZ(float angle);
    template RUNEMATH_API void Matrix3<double>::SetRotatedZ(double angle);

    /**
        Spocita matici pro transformaci do souradnicoveho systemu trojuhelniku (normala shodna s osou Z)
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API void Matrix3<T>::SetMatrixToFaceSpace(const Vec3<T> &v0, const Vec3<T> &v1, const Vec3<T> &v2)
    {
        Vec3<T> normal(Vec3<T>::GetNormal(v0, v1, v2));
        Vec3<T> axis(Vec3<T>::Cross(normal, Vec3<T>(0, 0, 1)));
        axis.Normalize();
        SetRotated(axis, Math<T>::Acos(Vec3<T>::Dot(normal, Vec3<T>(0, 0, 1)))*Math<T>::Deg());
    }

    template RUNEMATH_API void Matrix3<float>::SetMatrixToFaceSpace(const Vec3<float> &v0, const Vec3<float> &v1, const Vec3<float> &v2);
    template RUNEMATH_API void Matrix3<double>::SetMatrixToFaceSpace(const Vec3<double> &v0, const Vec3<double> &v1, const Vec3<double> &v2);

    /**
        Vynasobi matici, ktera ma na diagonale vektor a vsude jinde nuly ("v" byva vetsinou scale)
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API void Matrix3<T>::MultipleDiagonal(const Vec3<T> &v)
    {
        RUNE_SET_M3(m,  m[0] * v.x, m[1] * v.y, m[2] * v.z,
                        m[3] * v.x, m[4] * v.y, m[5] * v.z,
                        m[6] * v.x, m[7] * v.y, m[8] * v.z);
    }

    template RUNEMATH_API void Matrix3<float>::MultipleDiagonal(const Vec3<float> &v);
    template RUNEMATH_API void Matrix3<double>::MultipleDiagonal(const Vec3<double> &v);

    /**
        Z druhe strany vynasobi matici, ktera ma na diagonale vektor a vsude jinde nuly ("v" byva vetsinou scale)
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API void Matrix3<T>::MultipleDiagonalInv(const Vec3<T> &v)
    {
        RUNE_SET_M3(m,  m[0] * v.x, m[1] * v.x, m[2] * v.x,
                        m[3] * v.y, m[4] * v.y, m[5] * v.y,
                        m[6] * v.z, m[7] * v.z, m[8] * v.z);
    }

    template RUNEMATH_API void Matrix3<float>::MultipleDiagonalInv(const Vec3<float> &v);
    template RUNEMATH_API void Matrix3<double>::MultipleDiagonalInv(const Vec3<double> &v);

    /**
        Spocita transponovanou matici
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API void Matrix3<T>::Transpose()
    {
        T tmp;
        tmp = m[1]; m[1] = m[3]; m[3]  = tmp;
        tmp = m[2]; m[2] = m[6]; m[6]  = tmp;
        tmp = m[5]; m[5] = m[7]; m[7]  = tmp;
    }

    template RUNEMATH_API void Matrix3<float>::Transpose();
    template RUNEMATH_API void Matrix3<double>::Transpose();

    /**
        Spocita adjungovanou matici a nastavi ji do parametru
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API void Matrix3<T>::GetAdjoint(Matrix3<T> &adjn) const
    {
        RUNE_SET_M3(adjn,   m[4] * m[8] - m[5] * m[7],     m[2] * m[7] - m[1] * m[8],      m[1] * m[5] - m[2] * m[4],
                            m[5] * m[6] - m[3] * m[8],     m[0] * m[8] - m[2] * m[6],      m[2] * m[3] - m[0] * m[5],
                            m[3] * m[7] - m[4] * m[6],     m[1] * m[6] - m[0] * m[7],      m[0] * m[4] - m[1] * m[3]);
    }

    template RUNEMATH_API void Matrix3<float>::GetAdjoint(Matrix3<float> &adjn) const;
    template RUNEMATH_API void Matrix3<double>::GetAdjoint(Matrix3<double> &adjn) const;

    /**
        Spocita inverzni matici
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API void Matrix3<T>::Inverse()
    {
        Matrix3<T> temp;
        GetAdjoint(temp);

        T determinant = m[0] * temp[0] + m[1] * temp[3] + m[2] * temp[6];
        if (Math<T>::SafeIsZero(determinant)) return;

        *this = temp;
        MultipleScalar(1 / determinant);
    }

    template RUNEMATH_API void Matrix3<float>::Inverse();
    template RUNEMATH_API void Matrix3<double>::Inverse();

    /**
        Vraci Eulerovy uhly v poradi yaw-pitch-roll (= ZXY)
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API Vec3<T> Matrix3<T>::GetYawPitchRoll(const Matrix3<T> &m)
    {
        Vec3<T> s(m.GetScale());
        s.Inverse();
        return Vec3<T>(Math<T>::Deg() * Math<T>::Asin(-m[7] * s.y),
                       Math<T>::Deg() * Math<T>::Atan2(m[6] * s.x, m[8] * s.z),
                       Math<T>::Deg() * Math<T>::Atan2(m[1] * s.y, m[4] * s.y));
    }

    template RUNEMATH_API Vec3<float> Matrix3<float>::GetYawPitchRoll(const Matrix3<float> &m);
    template RUNEMATH_API Vec3<double> Matrix3<double>::GetYawPitchRoll(const Matrix3<double> &m);

    /**
        Transformuje prostorovy vektor
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API Vec3<T> operator *(const Matrix3<T> &m, const Vec3<T> &v)
    {
        return Vec3<T>(v.x * m[0] + v.y * m[3] + v.z * m[6],
                       v.x * m[1] + v.y * m[4] + v.z * m[7],
                       v.x * m[2] + v.y * m[5] + v.z * m[8]);
    }

    template RUNEMATH_API Vec3<float> operator *(const Matrix3<float> &m, const Vec3<float> &v);
    template RUNEMATH_API Vec3<double> operator *(const Matrix3<double> &m, const Vec3<double> &v);

    /**
        Transformuje prostorovy vektor transponovanou matici
    **************************************************************************************************/
    template<typename T>
    RUNEMATH_API Vec3<T> operator *(const Vec3<T> &v, const Matrix3<T> &m)
    {
        return Vec3<T>(v.x * m[0] + v.y * m[1] + v.z * m[2],
                       v.x * m[3] + v.y * m[4] + v.z * m[5],
                       v.x * m[6] + v.y * m[7] + v.z * m[8]);
    }

    template RUNEMATH_API Vec3<float> operator *(const Vec3<float> &v, const Matrix3<float> &m);
    template RUNEMATH_API Vec3<double> operator *(const Vec3<double> &v, const Matrix3<double> &m);
}
