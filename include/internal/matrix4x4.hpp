// Copyright 2023 Suraj Shettigar
// SPDX-License-Identifier: Apache-2.0

#ifndef CGMATH_INTERNAL_MATRIX4x4_HPP
#define CGMATH_INTERNAL_MATRIX4x4_HPP

#include "matrix2x2.hpp"

namespace cgmath::internal {

/// 4x4 Matrix using Column Major Mathematical Notation and stored in a Row
/// Major memory layout. It uses block matrix notation for storing the matrix.
/// i.e it uses 4 2x2 matrices to store a single 4x4 matrices.
/// \n
/// | a b c d |\n
/// | e f g h |\n
/// | i j k l |\n
/// | m n o p |\n
/// ==\n
/// | A B |\n
/// | C D |\n
class Matrix4x4 {
 public:
  // Constructors / Destructors

  /// Constructs an identity 4x4 matrix.
  explicit Matrix4x4()
#if defined(USE_INTRINSICS) && defined(__AVX512F__)
  // TODO: Matrix4x4 Identity initialization for AVX512
#else
      : m_value{Matrix2x2{1.0, 0.0, 0.0, 1.0}, Matrix2x2{0.0, 0.0, 0.0, 0.0},
                Matrix2x2{0.0, 0.0, 0.0, 0.0}, Matrix2x2{1.0, 0.0, 0.0, 1.0}}
#endif
      {};

  /// Constructs a 4x4 matrix with x, y, z and w (translation in the case of 3D
  /// affine transform matrix) column values.
  explicit Matrix4x4(FLOAT x_x, FLOAT x_y, FLOAT x_z, FLOAT x_w, FLOAT y_x,
                     FLOAT y_y, FLOAT y_z, FLOAT y_w, FLOAT z_x, FLOAT z_y,
                     FLOAT z_z, FLOAT z_w, FLOAT t_x, FLOAT t_y, FLOAT t_z,
                     FLOAT t_w)
#if defined(USE_INTRINSICS) && defined(__AVX512F__)
  // TODO: Matrix4x4 Identity initialization for AVX512
#else
      : m_value{Matrix2x2{x_x, x_y, y_x, y_y}, Matrix2x2{z_x, z_y, t_x, t_y},
                Matrix2x2{x_z, x_w, y_z, y_w}, Matrix2x2{z_z, z_w, t_z, t_w}}
#endif
      {};

  /// Constructs a 4x4 matrix with 3D x, y, z and w (translation in the case of
  /// 3D affine transform matrix) array values.
  explicit Matrix4x4(const std::array<FLOAT, 4> &x_axis,
                     const std::array<FLOAT, 4> &y_axis,
                     const std::array<FLOAT, 4> &z_axis,
                     const std::array<FLOAT, 4> &translation)
#if defined(USE_INTRINSICS) && defined(__AVX512F__)
  // TODO: Matrix4x4 Identity initialization for AVX512
#else
      : m_value{Matrix2x2{x_axis[0], x_axis[1], y_axis[0], y_axis[1]},
                Matrix2x2{z_axis[0], z_axis[1], translation[0], translation[1]},
                Matrix2x2{x_axis[2], x_axis[3], y_axis[2], y_axis[3]},
                Matrix2x2{z_axis[2], z_axis[3], translation[2], translation[3]}}
#endif
      {};

  /// Constructs a 4x4 matrix with 3D x, y, z and w (translation in the case of
  /// 3D affine transform matrix) vector values.
  explicit Matrix4x4(const Vector &x_axis, const Vector &y_axis,
                     const Vector &z_axis, const Vector &translation)
#if defined(USE_INTRINSICS) && defined(__AVX512F__)
      // TODO: Matrix4x4 Identity initialization for AVX512
      {};
#else
  {
#ifdef USE_INTRINSICS
#ifdef USE_DOUBLE
#if defined(__AVX2__) || defined(__AVX__)
    // [x_x, y_x, x_y, y_y]
    m_value[0].m_value.m_value =
        _mm256_shuffle_pd(x_axis.m_value, y_axis.m_value, 0b0000);
    // [z_x, t_x, z_y, t_y]
    m_value[1].m_value.m_value =
        _mm256_shuffle_pd(z_axis.m_value, translation.m_value, 0b0000);
    // [x_z, y_z, x_w, y_w]
    m_value[2].m_value.m_value =
        _mm256_shuffle_pd(x_axis.m_value, y_axis.m_value, 0b1111);
    // [z_z, t_z, z_w, t_w]
    m_value[3].m_value.m_value =
        _mm256_shuffle_pd(z_axis.m_value, translation.m_value, 0b1111);
#else
    // [x_x, y_x]
    m_value[0].m_value.m_value[0] =
        _mm_shuffle_pd(x_axis.m_value[0], y_axis.m_value[0], 0b00);
    // [x_y, y_y]
    m_value[0].m_value.m_value[1] =
        _mm_shuffle_pd(x_axis.m_value[1], y_axis.m_value[1], 0b00);
    // [z_x, t_x]
    m_value[1].m_value.m_value[0] =
        _mm256_shuffle_pd(z_axis.m_value[0], translation.m_value[0], 0b00);
    // [z_y, t_y]
    m_value[1].m_value.m_value[1] =
        _mm256_shuffle_pd(z_axis.m_value[1], translation.m_value[1], 0b00);
    // [x_z, y_z]
    m_value[2].m_value.m_value[0] =
        _mm256_shuffle_pd(x_axis.m_value[0], y_axis.m_value[0], 0b11);
    // [x_w, y_w]
    m_value[2].m_value.m_value[1] =
        _mm256_shuffle_pd(x_axis.m_value[1], y_axis.m_value[1], 0b11);
    // [z_z, t_z]
    m_value[3].m_value.m_value[0] =
        _mm256_shuffle_pd(z_axis.m_value[0], translation.m_value[0], 0b11);
    // [z_w, t_w]
    m_value[3].m_value.m_value[1] =
        _mm256_shuffle_pd(z_axis.m_value[1], translation.m_value[1], 0b11);
#endif  // AVX INTRINSICS
#else
    // [x_x, x_y, y_x, y_y]
    m_value[0].m_value.m_value =
        _mm_shuffle_ps(x_axis.m_value, y_axis.m_value, _MM_SHUFFLE(1, 0, 1, 0));
    // [x_x, y_x, x_y, y_y]
    m_value[0].m_value.m_value =
        _mm_shuffle_ps(m_value[0].m_value.m_value, m_value[0].m_value.m_value,
                       _MM_SHUFFLE(3, 1, 2, 0));

    // [z_x, z_y, t_x, t_y]
    m_value[1].m_value.m_value = _mm_shuffle_ps(
        z_axis.m_value, translation.m_value, _MM_SHUFFLE(1, 0, 1, 0));
    // [z_x, t_x, z_y, t_y]
    m_value[1].m_value.m_value =
        _mm_shuffle_ps(m_value[1].m_value.m_value, m_value[1].m_value.m_value,
                       _MM_SHUFFLE(3, 1, 2, 0));

    // [x_z, x_w, y_z, y_w]
    m_value[2].m_value.m_value =
        _mm_shuffle_ps(x_axis.m_value, y_axis.m_value, _MM_SHUFFLE(3, 2, 3, 2));
    // [x_z, y_z, x_w, y_w]
    m_value[2].m_value.m_value =
        _mm_shuffle_ps(m_value[2].m_value.m_value, m_value[2].m_value.m_value,
                       _MM_SHUFFLE(3, 1, 2, 0));

    // [z_z, z_w, t_z, t_w]
    m_value[3].m_value.m_value = _mm_shuffle_ps(
        z_axis.m_value, translation.m_value, _MM_SHUFFLE(3, 2, 3, 2));
    // [z_z, t_z, z_w, t_w]
    m_value[3].m_value.m_value =
        _mm_shuffle_ps(m_value[3].m_value.m_value, m_value[3].m_value.m_value,
                       _MM_SHUFFLE(3, 1, 2, 0));
#endif  // USE_DOUBLE
#else
    m_value[0].m_value.setX(x_axis[0]);
    m_value[0].m_value.setY(y_axis[0]);
    m_value[0].m_value.setZ(x_axis[1]);
    m_value[0].m_value.setW(y_axis[1]);

    m_value[1].m_value.setX(z_axis[0]);
    m_value[1].m_value.setY(translation[0]);
    m_value[1].m_value.setZ(z_axis[1]);
    m_value[1].m_value.setW(translation[1]);

    m_value[2].m_value.setX(x_axis[2]);
    m_value[2].m_value.setY(y_axis[2]);
    m_value[2].m_value.setZ(x_axis[3]);
    m_value[2].m_value.setW(y_axis[3]);

    m_value[3].m_value.setX(z_axis[2]);
    m_value[3].m_value.setY(translation[2]);
    m_value[3].m_value.setZ(z_axis[3]);
    m_value[3].m_value.setW(translation[3]);
#endif  // USE_INTRINSICS
  };
#endif

  ~Matrix4x4() = default;

  // Copy
  Matrix4x4(const Matrix4x4 &matrix) {
    if (this != &matrix) m_value = matrix.m_value;
  }
  Matrix4x4 &operator=(const Matrix4x4 &matrix) {
    if (this != &matrix) m_value = matrix.m_value;
    return *this;
  }

  // Getters / Setters
  /**
   * Returns the matrix row.
   * @param index Row index.
   * @return 4D vector containing given row values. (Eg: [x_x, y_x, z_x, t_x])
   */
  [[nodiscard]] inline Vector operator[](size_t index) const {
#if defined(USE_INTRINSICS) && defined(__AVX512F__)
    // TODO: Matrix4x4 [] operator for AVX512.
#else
    // 0 -> 0, 1 -> 0, 2 -> 2, 3 -> 2
    const auto i =
        static_cast<size_t>(std::floor(static_cast<FLOAT>(index) / 2.0) * 2.0);
    // 0 -> 0,1,0,1 ; 1 -> 2,3,2,3 ; 2 -> 0,1,0,1; 3 -> 2,3,2,3
    const auto j = (index * 2) % 4;
    return Vector{m_value[i].m_value[j], m_value[i].m_value[j + 1],
                  m_value[i + 1].m_value[j], m_value[i + 1].m_value[j + 1]};
#endif
  }

  /**
   * Returns the matrix row.
   * @param index Row index.
   * @return 4D vector containing given row values. (Eg: [x_x, y_x, z_x, t_x])
   */
  [[nodiscard]] inline Vector getRow(size_t index) const {
    return (*this)[index];
  }

  /**
   * Returns the matrix column.
   * @param index Column index.
   * @return 4D vector containing given column values. (Eg: [x_x, x_y, x_z,
   * x_w])
   */
  [[nodiscard]] inline Vector getColumn(size_t index) const {
#if defined(USE_INTRINSICS) && defined(__AVX512F__)
// TODO: Matrix4x4 getColumn() for AVX512.
#else
    // 0 -> 0, 1 -> 0, 2 -> 1, 3 -> 1
    const auto i =
        static_cast<size_t>(std::floor(static_cast<FLOAT>(index) / 2.0));
    // 0 -> 0,2,0,2 ; 1 -> 1,3,1,3 ; 2 -> 0,2,0,2; 3 -> 1,3,1,3
    const auto j = index % 2;
    return Vector{m_value[i].m_value[j], m_value[i].m_value[j + 2],
                  m_value[i + 2].m_value[j], m_value[i + 2].m_value[j + 2]};
#endif
  }

  /// Returns x-axis (column 0) of the matrix.
  [[nodiscard]] inline Vector getXAxis() const { return getColumn(0); }
  /// Returns y-axis (column 1) of the matrix.
  [[nodiscard]] inline Vector getYAxis() const { return getColumn(1); }
  /// Returns z-axis (column 2) of the matrix.
  [[nodiscard]] inline Vector getZAxis() const { return getColumn(2); }
  /// Returns the translation values (w-axis) (column 3) of the matrix.
  [[nodiscard]] inline Vector getTranslation() const { return getColumn(3); }

  /// Sets the x-axis column with the given values.
  inline void setXAxis(FLOAT x_x, FLOAT x_y, FLOAT x_z, FLOAT x_w) {
#if defined(USE_INTRINSICS) && defined(__AVX512F__)
// TODO: Matrix4x4 setXAxis for AVX512.
#else
    m_value[0].m_value.setX(x_x);
    m_value[0].m_value.setZ(x_y);
    m_value[2].m_value.setX(x_z);
    m_value[2].m_value.setZ(x_w);
#endif
  }
  /// Sets the x-axis column with the given 4D array values.
  inline void setXAxis(const std::array<FLOAT, 4> &x_axis) {
    setXAxis(x_axis[0], x_axis[1], x_axis[2], x_axis[3]);
  }
  /// Sets the x-axis column with the given 4D vector values.
  inline void setXAxis(const Vector &x_axis) {
    // TODO: Matrix4x4 write efficient setXAxis function for simd vectors.
    setXAxis(x_axis[0], x_axis[1], x_axis[2], x_axis[3]);
  }
  /// Sets the y-axis column with the given values.
  inline void setYAxis(FLOAT y_x, FLOAT y_y, FLOAT y_z, FLOAT y_w) {
#if defined(USE_INTRINSICS) && defined(__AVX512F__)
    // TODO: Matrix4x4 setYAxis for AVX512.
#else
    m_value[0].m_value.setY(y_x);
    m_value[0].m_value.setW(y_y);
    m_value[2].m_value.setY(y_z);
    m_value[2].m_value.setW(y_w);
#endif
  }
  /// Sets the y-axis column with the given 4D array values.
  inline void setYAxis(const std::array<FLOAT, 4> &y_axis) {
    setYAxis(y_axis[0], y_axis[1], y_axis[2], y_axis[3]);
  }
  /// Sets the y-axis column with the given 4D vector values.
  inline void setYAxis(const Vector &y_axis) {
    // TODO: Matrix4x4 write efficient setYAxis function for simd vectors.
    setYAxis(y_axis[0], y_axis[1], y_axis[2], y_axis[3]);
  }
  /// Sets the z-axis column with the given values.
  inline void setZAxis(FLOAT z_x, FLOAT z_y, FLOAT z_z, FLOAT z_w) {
#if defined(USE_INTRINSICS) && defined(__AVX512F__)
// TODO: Matrix4x4 setZAxis for AVX512.
#else
    m_value[1].m_value.setX(z_x);
    m_value[1].m_value.setZ(z_y);
    m_value[3].m_value.setX(z_z);
    m_value[3].m_value.setZ(z_w);
#endif
  }
  /// Sets the z-axis column with the given 4D array values.
  inline void setZAxis(const std::array<FLOAT, 4> &z_axis) {
    setZAxis(z_axis[0], z_axis[1], z_axis[2], z_axis[3]);
  }
  /// Sets the z-axis column with the given 4D vector values.
  inline void setZAxis(const Vector &z_axis) {
    // TODO: Matrix4x4 write efficient setZAxis function for simd vectors.
    setZAxis(z_axis[0], z_axis[1], z_axis[2], z_axis[3]);
  }
  /// Sets the translation column (w-axis) with the given values.
  inline void setTranslation(FLOAT t_x, FLOAT t_y, FLOAT t_z, FLOAT t_w) {
#if defined(USE_INTRINSICS) && defined(__AVX512F__)
// TODO: Matrix4x4 setTranslation for AVX512.
#else
    m_value[1].m_value.setY(t_x);
    m_value[1].m_value.setW(t_y);
    m_value[3].m_value.setY(t_z);
    m_value[3].m_value.setW(t_w);
#endif
  }
  /// Sets the translation column (w-axis) with the given 4D array values.
  inline void setTranslation(const std::array<FLOAT, 4> &translation) {
    setTranslation(translation[0], translation[1], translation[2],
                   translation[3]);
  }
  /// Sets the translation column (w-axis) with the given 4D vector values.
  inline void setTranslation(const Vector &translation) {
    // TODO: Matrix4x4 write efficient setTranslation function for simd vectors.
    setTranslation(translation[0], translation[1], translation[2],
                   translation[3]);
  }

  // Operators
  inline bool operator==(const Matrix4x4 &rhs) const {
#if defined(USE_INTRINSICS) && defined(__AVX512F__)
    // TODO: Matrix4x4 approx equal comparision for AVX512.
#else
    return m_value[0] == rhs.m_value[0] && m_value[1] == rhs.m_value[1] &&
           m_value[2] == rhs.m_value[2] && m_value[3] == rhs.m_value[3];
#endif
  }

  inline bool operator!=(const Matrix4x4 &rhs) const {
#if defined(USE_INTRINSICS) && defined(__AVX512F__)
// TODO: Matrix4x4 approx not equal comparision for AVX512.
#else
    return m_value[0] != rhs.m_value[0] || m_value[1] != rhs.m_value[1] ||
           m_value[2] != rhs.m_value[2] || m_value[3] != rhs.m_value[3];
#endif
  }

  /// Matrix-Matrix multiplication.
  friend inline Matrix4x4 operator*(const Matrix4x4 &lhs,
                                    const Matrix4x4 &rhs) {
#if defined(USE_INTRINSICS) && defined(__AVX512F__)
    // TODO: Matrix4x4  matrix-matrix mulitplication for AVX512.
#else
    // Matrix block multiplication method.
    // Refer:
    // https://en.wikipedia.org/wiki/Block_matrix#Block_matrix_multiplication
    return Matrix4x4{std::array<Matrix2x2, 4>{
        lhs.m_value[0] * rhs.m_value[0] + lhs.m_value[1] * rhs.m_value[2],
        lhs.m_value[0] * rhs.m_value[1] + lhs.m_value[1] * rhs.m_value[3],
        lhs.m_value[2] * rhs.m_value[0] + lhs.m_value[3] * rhs.m_value[2],
        lhs.m_value[2] * rhs.m_value[1] + lhs.m_value[3] * rhs.m_value[3]}};
#endif
  }

  friend inline Vector operator*(const Matrix4x4 &lhs, const Vector &rhs) {
#if defined(USE_INTRINSICS) && defined(__AVX512F__)
    // TODO: Matrix4x4  matrix-vector mulitplication for AVX512.
#else
    // Matrix block multiplication method.
    // Refer:
    // https://en.wikipedia.org/wiki/Block_matrix#Block_matrix_multiplication
    // | A  B |  | X |  ==> | AX + BY |
    // | C  D |  | Y |      | CX + DY |
    //  Where,
    // X = | vec_x 0 |  and Y = | vec_z 0 |
    //     | vec_y 0 |          | vec_t 0 |
#ifdef USE_INTRINSICS
#ifdef USE_DOUBLE
#if defined(__AVX2__) || defined(__AVX__)
    __m256d zero = _mm256_setzero_pd();
    __m256d mat_x =
        _mm256_shuffle_pd(rhs.m_value, zero, _MM_SHUFFLE(3, 2, 1, 0));
    // [z, t, 0, 0]
    __m256d mat_y =
        _mm256_permute2f128_pd(rhs.m_value, zero, _MM_SHUFFLE(0, 2, 0, 1));
    mat_y = _mm256_permute4x64_pd(mat_y, _MM_SHUFFLE(3, 1, 2, 0));
    Matrix2x2 X{Vector{mat_x}};
    Matrix2x2 Y{Vector{mat_y}};
#else
    __m128d zero = _mm_setzero_pd();
    __m128d mat_x0 = _mm_shuffle_pd(rhs.m_value[0], zero, _MM_SHUFFLE2(1, 0));
    __m128d mat_x1 = _mm_shuffle_pd(rhs.m_value[0], zero, _MM_SHUFFLE2(1, 1));
    __m128d mat_y0 = _mm_shuffle_pd(rhs.m_value[1], zero, _MM_SHUFFLE2(1, 0));
    __m128d mat_y1 = _mm_shuffle_pd(rhs.m_value[1], zero, _MM_SHUFFLE2(1, 1));
    Matrix2x2 X{Vector{mat_x0, mat_x1}};
    Matrix2x2 Y{Vector{mat_y0, mat_y1}};
#endif  // AVX INTRINSICS
#else
    __m128 zero = _mm_setzero_ps();
    __m128 mat_x = _mm_shuffle_ps(rhs.m_value, zero, _MM_SHUFFLE(1, 0, 1, 0));
    mat_x = _mm_shuffle_ps(mat_x, mat_x, _MM_SHUFFLE(3, 1, 2, 0));
    __m128 mat_y = _mm_shuffle_ps(rhs.m_value, zero, _MM_SHUFFLE(1, 0, 3, 2));
    mat_y = _mm_shuffle_ps(mat_y, mat_y, _MM_SHUFFLE(3, 1, 2, 0));
    Matrix2x2 X{Vector{mat_x}};
    Matrix2x2 Y{Vector{mat_y}};
#endif  // USE_DOUBLE
#else
    Matrix2x2 X{Vector{rhs[0], 0.0, rhs[1], 0.0}};
    Matrix2x2 Y{Vector{rhs[2], 0.0, rhs[3], 0.0}};
#endif  // USE_INTRINSICS
    Matrix2x2 res_x = lhs.m_value[0] * X + lhs.m_value[1] * Y;
    Matrix2x2 res_y = lhs.m_value[2] * X + lhs.m_value[3] * Y;
    return Vector{res_x.m_value[0], res_x.m_value[2], res_y.m_value[0],
                  res_y.m_value[2]};
#endif
  }

  // Functions

  static inline Matrix4x4 transpose(const Matrix4x4 &matrix) {
#if defined(USE_INTRINSICS) && defined(__AVX512F__)
    // TODO: Matrix4x4  transpose for AVX512.
#else
    // Block matrix transpose method
    // Refer:
    // https://en.wikipedia.org/wiki/Block_matrix#Block_transpose
    // | A  B |  ==> | T(A) T(C) |
    // | C  D |      | T(B) T(D) |
    return Matrix4x4{
        std::array<Matrix2x2, 4>{Matrix2x2::transpose(matrix.m_value[0]),
                                 Matrix2x2::transpose(matrix.m_value[2]),
                                 Matrix2x2::transpose(matrix.m_value[1]),
                                 Matrix2x2::transpose(matrix.m_value[3])}};
#endif
  }
  [[nodiscard]] inline FLOAT determinant() const {
#if defined(USE_INTRINSICS) && defined(__AVX512F__)
    // TODO: Matrix4x4  determinant for AVX512.
#else
    // Block matrix determinant method
    // Refer:
    // https://en.wikipedia.org/wiki/Block_matrix#Block_matrix_determinant
    // | a    b    c    d |     | A   B |
    // | e    f    g    h | ==  | C   D |
    // | i    j    k    l |
    // | m    n    o    p |
    // Determinant = Determinant(A)*Determinant(D- C*Inverse(A)*B);
    //             = Determinant(AD - BC)
    return (m_value[0] * m_value[3] - m_value[1] * m_value[2]).determinant();
#endif
  }

  static inline Matrix4x4 inverse(const Matrix4x4 &matrix) {
#if defined(USE_INTRINSICS) && defined(__AVX512F__)
    // TODO: Matrix4x4  inverse for AVX512.
#else
    // Block matrix method
    // Refer: https://en.wikipedia.org/wiki/Block_matrix#Block_matrix_inversion
    // | a    b    c    d |     | A   B |
    // | e    f    g    h | ==  | C   D |
    // | i    j    k    l |
    // | m    n    o    p |
    //
    // | A    B | == | X    Y |
    // | C    D |    | Z    W |
    // Where,
    // X = 1 / Det(matrix) * Adj(Det(D)*A - B*(Adj(D)*C))
    // Y = 1 / Det(matrix) * Adj(Det(B)*C - D*Adj(Adj(A)*B))
    // Z = 1 / Det(matrix) * Adj(Det(C)*B - A*Adj(Adj(D)*C))
    // W = 1 / Det(matrix) * Adj(Det(A)*D - C*(Adj(A)*B))
    const FLOAT r_det = static_cast<FLOAT>(1.0) / matrix.determinant();

    const Matrix2x2 X =
        r_det * Matrix2x2::adjoint(
                    matrix.m_value[3].determinant() * matrix.m_value[0] -
                    matrix.m_value[1] * Matrix2x2::adjoint(matrix.m_value[3]) *
                        matrix.m_value[2]);
    const Matrix2x2 Y =
        r_det *
        Matrix2x2::adjoint(
            matrix.m_value[1].determinant() * matrix.m_value[2] -
            matrix.m_value[3] *
                Matrix2x2::adjoint(Matrix2x2::adjoint(matrix.m_value[0]) *
                                   matrix.m_value[1]));
    const Matrix2x2 Z =
        r_det *
        Matrix2x2::adjoint(
            matrix.m_value[2].determinant() * matrix.m_value[1] -
            matrix.m_value[0] *
                Matrix2x2::adjoint(Matrix2x2::adjoint(matrix.m_value[3]) *
                                   matrix.m_value[2]));
    const Matrix2x2 W =
        r_det * Matrix2x2::adjoint(
                    matrix.m_value[0].determinant() * matrix.m_value[3] -
                    matrix.m_value[2] * Matrix2x2::adjoint(matrix.m_value[0]) *
                        matrix.m_value[1]);
    return Matrix4x4{std::array<Matrix2x2, 4>{X, Y, Z, W}};
#endif
  }

 private:
#if defined(USE_INTRINSICS) && defined(__AVX512F__)
#ifdef USE_DOUBLE
  std::array<__m512d, 2> m_value;
  explicit Matrix4x4(const std::array<__m512d, 2> &value) : m_value{value} {};
#else
  __m512 m_value;
  explicit Matrix4x4(__m512 value) : m_value{value} {};
#endif
#else
  std::array<Matrix2x2, 4> m_value{
      Matrix2x2{1.0, 0.0, 0.0, 1.0}, Matrix2x2{0.0, 0.0, 0.0, 0.0},
      Matrix2x2{0.0, 0.0, 0.0, 0.0}, Matrix2x2{1.0, 0.0, 0.0, 1.0}};
  explicit Matrix4x4(const std::array<Matrix2x2, 4> &value) : m_value{value} {};
#endif
};

}  // namespace cgmath::internal

#endif  // CGMATH_INTERNAL_MATRIX4x4_HPP
