#pragma once

// Usetrime si par radku psani :)
#define RUNE_SET_M4(m, m0, m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13, m14, m15) \
    m[0] = m0, m[1] = m1, m[2] = m2,   m[3] = m3,   m[4] = m4,   m[5] = m5,   m[6] = m6,   m[7] = m7,   \
    m[8] = m8, m[9] = m9, m[10] = m10, m[11] = m11, m[12] = m12, m[13] = m13, m[14] = m14, m[15] = m15;
#define RUNE_ADD_M4(func) Matrix4 m1 = *this, m2; m2.func; SetProduct(m1, m2);

namespace Rune
{
    /**
        Transformacni matice 4x4
    **************************************************************************************************/
    template<typename T>
    class RUNEMATH_API Matrix4
    {
        T m[16];

    public:
        Matrix4() {}
        Matrix4(const Quaternion<T> &q)                 { operator =(q); }

        // Operatory * a *= pro nasobeni matic nejsou. Pro (m = a * b) pouzijte m.SetProduct(a, b) a pro (m *= a) pouzijte m.Multiple(a);
        T &operator [](int i)                           { return m[i]; }
        T operator [](int i) const                      { return m[i]; }
        operator T*()                                   { return m; }
        operator const T*() const                       { return m; }

        // Jednodussi funkce
        Vec3<T> GetTranslation() const                  { return Vec3<T>(m[12], m[13], m[14]); }
        Vec4<T> GetRow(int i) const                     { return Vec4<T>(m[i], m[i+4], m[i+8], m[i+12]); }
        Vec4<T> GetColumn(int i) const                  { int i4 = i*4; return Vec4<T>(m[i4], m[i4+1], m[i4+2], m[i4+3]); }
        void GetMatrix3(Matrix3<T> &m2) const           { RUNE_SET_M3(m2,  m[0], m[1], m[2], m[4], m[5], m[6], m[8], m[9], m[10]); }
        void SetRow(const Vec4<T> &v, int i)            { m[i] = v.x; m[i+4] = v.y; m[i+8] = v.z; m[i+12] = v.w; }
        void SetColumn(const Vec4<T> &v, int i)         { int i4 = i*4; m[i4] = v.x; m[i4+1] = v.y; m[i4+2] = v.z; m[i4+3] = v.w; }
        void SetMatrix3(const Matrix3<T> &m3)           { RUNE_SET_M4(m, m3[0], m3[1], m3[2], 0, m3[3], m3[4], m3[5], 0, m3[6], m3[7], m3[8], 0, 0, 0, 0, 1); }
        void SetSum(const Matrix4 &m1, const Matrix4 &m2) { for (int i = 0; i < 16; i++) m[i] = m1[i] + m2[i]; }
        void SetIdentity()                              { RUNE_SET_M4(m, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1); }
        void SetTranslated(T x, T y, T z)               { RUNE_SET_M4(m, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, x, y, z, 1); }
        void SetTranslated(const Vec3<T> &v)            { SetTranslated(v.x, v.y, v.z); }
        void SetRotated(const Vec3<T> &v, T angle)      { SetRotated(v.x, v.y, v.z, angle); }
        void SetScaled(T x, T y, T z)                   { RUNE_SET_M4(m, x, 0, 0, 0, 0, y, 0, 0, 0, 0, z, 0, 0, 0, 0, 1); }
        void SetScaled(const Vec3<T> &v)                { SetScaled(v.x, v.y, v.z); }
        void SetScaled(T f)                             { SetScaled(f, f, f); }
        void SetOrtho2D(T width, T height)              { SetOrtho(0, 0, width, height, 0, 1); }
        void SetYawPitchRoll(T angZ, T angX, T angY)    { SetRotatedZ(angZ); RotateX(angX); RotateY(angY); }
        void Multiple(const Matrix4 &m)                 { Matrix4 tmp = *this; SetProduct(tmp, m); }
        void Multiple(const Matrix3<T> &m)              { RUNE_ADD_M4(SetMatrix3(m)); }
        void Translate(T x, T y, T z)                   { RUNE_ADD_M4(SetTranslated(x, y, z)); }
        void Translate(const Vec3<T> &v)                { Translate(v.x, v.y, v.z); }
        void Rotate(T x, T y, T z, T angle)             { RUNE_ADD_M4(SetRotated(x, y, z, angle)); }
        void Rotate(const Vec3<T> &v, T angle)          { Rotate(v.x, v.y, v.z, angle); }
        void RotateX(T angle)                           { RUNE_ADD_M4(SetRotatedX(angle)); }
        void RotateY(T angle)                           { RUNE_ADD_M4(SetRotatedY(angle)); }
        void RotateZ(T angle)                           { RUNE_ADD_M4(SetRotatedZ(angle)); }
        void Scale(T x, T y, T z)                       { MultipleDiagonal(Vec4<T>(x, y, z, 1)); }
        void Scale(const Vec3<T> &v)                    { MultipleDiagonal(Vec4<T>(v, 1)); }
        void Scale(T f)                                 { MultipleDiagonal(Vec4<T>(f, f, f, 1)); }
        void MultipleScalar(T f)                        { for (int i = 0; i < 16; i++) m[i] *= f; }

        void operator =(const Quaternion<T> &q);
        operator Quaternion<T>() const;
        T Determinant() const;
        void SetPerspective(T fov, T aspect, T near, T far);
        void SetOrtho(T left, T top, T right, T bottom, T near, T far);
        void SetProduct(const Matrix4 &m1, const Matrix4 &m2);
        void SetRotated(T x, T y, T z, T angle);
        void SetPitchRollYaw(T anglex, T angley, T anglez);
        void SetRotatedX(T angle);
        void SetRotatedY(T angle);
        void SetRotatedZ(T angle);
        void SetMatrixToFaceSpace(const Vec3<T> &v0, const Vec3<T> &v1, const Vec3<T> &v2);
        void SetLookAt(const Vec3<T> &eye, const Vec3<T> &center, const Vec3<T> &up);
        void MultipleDiagonal(const Vec4<T> &v);
        void MultipleDiagonalInv(const Vec4<T> &v);
        void Inverse();
        void FastInverse();
        void Transpose();
    };

    template<typename T>
    RUNEMATH_API Vec3<T> operator *(const Matrix4<T> &m, const Vec3<T> &v);
    template<typename T>
    RUNEMATH_API Vec3<T> operator *(const Vec3<T> &v, const Matrix4<T> &m);
    template<typename T>
    RUNEMATH_API Vec4<T> operator *(const Matrix4<T> &m, const Vec4<T> &v);
    template<typename T>
    RUNEMATH_API Vec4<T> operator *(const Vec4<T> &v, const Matrix4<T> &m);

    template<typename T>
    std::ostream& operator <<(std::ostream &stream, const Matrix4<T> &m)
    {
        stream << m.GetColumn(0) << "\n";
        stream << m.GetColumn(1) << "\n";
        stream << m.GetColumn(2) << "\n";
        stream << m.GetColumn(3) << "\n";
        return stream;
    }

    typedef Matrix4<float> Matrix4f;
    typedef Matrix4<double> Matrix4d;
}
