#pragma once

// Usetrime si par radku psani :)
#define RUNE_SET_M3(m, m0, m1, m2, m3, m4, m5, m6, m7, m8) \
    m[0] = m0, m[1] = m1, m[2] = m2, m[3] = m3, m[4] = m4, m[5] = m5, m[6] = m6, m[7] = m7, m[8] = m8;
#define RUNE_ADD_M3(func) Matrix3 m1 = *this, m2; m2.func; SetProduct(m1, m2);

namespace Rune
{
    /**
        Transformacni matice 3x3
    **************************************************************************************************/
    template<typename T>
    class Matrix3
    {
        T m[9];

    public:
        Matrix3() {}
        Matrix3(const Quaternion<T> &q)                 { operator =(q); }

        // Operatory * a *= pro nasobeni matic nejsou, protoze dochazi ke zbytecnemu kopirovani cele
        // matice. Pro (m = a * b) pouzijte m.SetProduct(a, b) a pro (m *= a) pouzijte m.Multiple(a);

        T &operator [](int i)      { return m[i]; }
        T operator [](int i) const { return m[i]; }
        operator T*()              { return m; }
        operator const T*() const  { return m; }

        // Jednodussi funkce
        Vec3<T> GetRow(int i) const                     { return Vec3<T>(m[i], m[i+3], m[i+6]); }
        Vec3<T> GetColumn(int i) const                  { int i3 = i*3; return Vec3<T>(m[i3], m[i3+1], m[i3+2]); }
        Vec3<T> GetScale() const                        { return Vec3<T>(GetRow(0).Magnitude(), GetRow(1).Magnitude(), GetRow(2).Magnitude()); }
        T GetUniformScale() const                       { return GetRow(0).Magnitude(); }
        void SetRow(const Vec3<T> &v, int i)            { m[i] = v.x; m[i+3] = v.y; m[i+6] = v.z; }
        void SetColumn(const Vec3<T> &v, int i)         { int i3 = i*3; m[i3] = v.x; m[i3+1] = v.y; m[i3+2] = v.z; }
        void SetSum(const Matrix3 &m1, const Matrix3 &m2) { for (int i = 0; i < 9; i++) m[i] = m1[i] + m2[i]; }
        void SetIdentity()                              { RUNE_SET_M3(m, 1, 0, 0, 0, 1, 0, 0, 0, 1); }
        void SetRotated(const Vec3<T> &v, T angle)      { SetRotated(v.x, v.y, v.z, angle); }
        void SetScaled(const T x, const T y, const T z) { RUNE_SET_M3(m, x, 0, 0, 0, y, 0, 0, 0, z); }
        void SetScaled(const Vec3<T> &v)                { SetScaled(v.x, v.y, v.z); }
        void SetScaled(const T f)                       { SetScaled(f, f, f); }
        void SetYawPitchRoll(T angZ, T angX, T angY)    { SetRotatedZ(angZ); RotateX(angX); RotateY(angY); }
        void Multiple(const Matrix3 &m)                 { Matrix3 tmp = *this; SetProduct(tmp, m); }
        void Rotate(T x, T y, T z, T angle)             { RUNE_ADD_M3(SetRotated(x, y, z, angle)); }
        void Rotate(const Vec3<T> &v, T angle)          { Rotate(v.x, v.y, v.z, angle); }
        void Rotate(T anglex, T angley, T anglez)       { RUNE_ADD_M3(SetRotated(anglex, angley, anglez)); }
        void RotateX(T angle)                           { RUNE_ADD_M3(SetRotatedX(angle)); }
        void RotateY(T angle)                           { RUNE_ADD_M3(SetRotatedY(angle)); }
        void RotateZ(T angle)                           { RUNE_ADD_M3(SetRotatedZ(angle)); }
        void Scale(const T x, const T y, const T z)     { MultipleDiagonal(Vec3<T>(x, y, z)); }
        void Scale(const Vec3<T> &v)                    { MultipleDiagonal(v); }
        void Scale(const T f)                           { MultipleScalar(f); }
        void MultipleScalar(T f)                        { for (int i = 0; i < 9; i++) m[i] *= f; }
        void Adjoint()                                  { Matrix3 temp; GetAdjoint(temp); *this = temp; }

        RUNEMATH_API void operator =(const Quaternion<T> &q);
        RUNEMATH_API operator Quaternion<T>();
        RUNEMATH_API T Determinant() const;
        RUNEMATH_API void SetProduct(const Matrix3 &m1, const Matrix3 &m2);
        RUNEMATH_API void SetProductDiagonal(const Matrix3 &m1, const Vec3<T> &v);
        RUNEMATH_API void SetProductDiagonalInv(const Matrix3 &m1, const Vec3<T> &v);
        RUNEMATH_API void SetProductScalar(const Matrix3 &m1, T f);
        RUNEMATH_API void SetRotated(T x, T y, T z, T angle);
        RUNEMATH_API void SetPitchRollYaw(T anglex, T angley, T anglez);
        RUNEMATH_API void SetRotatedX(T angle);
        RUNEMATH_API void SetRotatedY(T angle);
        RUNEMATH_API void SetRotatedZ(T angle);
        RUNEMATH_API void SetMatrixToFaceSpace(const Vec3<T> &v0, const Vec3<T> &v1, const Vec3<T> &v2);
        RUNEMATH_API void MultipleDiagonal(const Vec3<T> &v);
        RUNEMATH_API void MultipleDiagonalInv(const Vec3<T> &v);
        RUNEMATH_API void Transpose();
        RUNEMATH_API void GetAdjoint(Matrix3 &adjn) const;
        RUNEMATH_API void Inverse();
        RUNEMATH_API Vec3<T> GetYawPitchRoll(const Matrix3 &m);
    };

    template<typename T>
    RUNEMATH_API Vec3<T> operator *(const Matrix3<T> &m, const Vec3<T> &v);

    template<typename T>
    RUNEMATH_API Vec3<T> operator *(const Vec3<T> &v, const Matrix3<T> &m);

    template<typename T>
    std::ostream& operator <<(std::ostream &stream, const Matrix3<T> &m)
    {
        stream << m.GetColumn(0) << "\n";
        stream << m.GetColumn(1) << "\n";
        stream << m.GetColumn(2) << "\n";
        return stream;
    }

    typedef Matrix3<float> Matrix3f;
    typedef Matrix3<double> Matrix3d;
}
