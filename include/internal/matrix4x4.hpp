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
    // [x_x, y_x, x_z, y_z]
    __m256d tmp0 = _mm256_shuffle_pd(x_axis.m_value, y_axis.m_value, 0b0000);
    // [z_x, t_x, z_z, t_z]
    __m256d tmp1 =
        _mm256_shuffle_pd(z_axis.m_value, translation.m_value, 0b0000);
    // [x_y, y_y, x_w, y_w]
    __m256d tmp2 = _mm256_shuffle_pd(x_axis.m_value, y_axis.m_value, 0b1111);
    // [z_y, t_y, z_w, t_w]
    __m256d tmp3 =
        _mm256_shuffle_pd(z_axis.m_value, translation.m_value, 0b1111);

    // [x_x, y_x, x_y, y_y]
    m_value[0].m_value.m_value =
        _mm256_permute2f128_pd(tmp0, tmp2, _MM_SHUFFLE(0, 2, 0, 0));
    // [z_x, t_x, z_y, t_y]
    m_value[1].m_value.m_value =
        _mm256_permute2f128_pd(tmp1, tmp3, _MM_SHUFFLE(0, 2, 0, 0));
    // [x_z, y_z, x_w, y_w]
    m_value[2].m_value.m_value =
        _mm256_permute2f128_pd(tmp0, tmp2, _MM_SHUFFLE(0, 3, 0, 1));
    // [z_z, t_z, z_w, t_w]
    m_value[3].m_value.m_value =
        _mm256_permute2f128_pd(tmp1, tmp3, _MM_SHUFFLE(0, 3, 0, 1));
#else
    // [x_x, y_x]
    m_value[0].m_value.m_value[0] =
        _mm_shuffle_pd(x_axis.m_value[0], y_axis.m_value[0], 0b00);
    // [x_y, y_y]
    m_value[0].m_value.m_value[1] =
        _mm_shuffle_pd(x_axis.m_value[0], y_axis.m_value[0], 0b11);
    // [z_x, t_x]
    m_value[1].m_value.m_value[0] =
        _mm_shuffle_pd(z_axis.m_value[0], translation.m_value[0], 0b00);
    // [z_y, t_y]
    m_value[1].m_value.m_value[1] =
        _mm_shuffle_pd(z_axis.m_value[0], translation.m_value[0], 0b11);
    // [x_z, y_z]
    m_value[2].m_value.m_value[0] =
        _mm_shuffle_pd(x_axis.m_value[1], y_axis.m_value[1], 0b00);
    // [x_w, y_w]
    m_value[2].m_value.m_value[1] =
        _mm_shuffle_pd(x_axis.m_value[1], y_axis.m_value[1], 0b11);
    // [z_z, t_z]
    m_value[3].m_value.m_value[0] =
        _mm_shuffle_pd(z_axis.m_value[1], translation.m_value[1], 0b00);
    // [z_w, t_w]
    m_value[3].m_value.m_value[1] =
        _mm_shuffle_pd(z_axis.m_value[1], translation.m_value[1], 0b11);
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

  inline Vector operator*(const Vector &rhs) const {
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
    // [x, 0, z, 0]
    __m256d mat_x = _mm256_shuffle_pd(rhs.m_value, zero, 0b0000);
    // [y, 0, t, 0]
    __m256d mat_y = _mm256_shuffle_pd(rhs.m_value, zero, 0b0101);

    Matrix2x2 X{
        Vector{_mm256_permute2f128_pd(mat_x, mat_y, _MM_SHUFFLE(0, 2, 0, 0))}};
    Matrix2x2 Y{
        Vector{_mm256_permute2f128_pd(mat_x, mat_y, _MM_SHUFFLE(0, 3, 0, 1))}};
#else
    __m128d zero = _mm_setzero_pd();
    __m128d mat_x0 = _mm_shuffle_pd(rhs.m_value[0], zero, _MM_SHUFFLE2(1, 0));
    __m128d mat_x1 = _mm_shuffle_pd(rhs.m_value[0], zero, _MM_SHUFFLE2(1, 1));
    __m128d mat_y0 = _mm_shuffle_pd(rhs.m_value[1], zero, _MM_SHUFFLE2(1, 0));
    __m128d mat_y1 = _mm_shuffle_pd(rhs.m_value[1], zero, _MM_SHUFFLE2(1, 1));
    Matrix2x2 X{Vector{{mat_x0, mat_x1}}};
    Matrix2x2 Y{Vector{{mat_y0, mat_y1}}};
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
    Matrix2x2 res_x = m_value[0] * X + m_value[1] * Y;
    Matrix2x2 res_y = m_value[2] * X + m_value[3] * Y;
    return Vector{res_x.m_value[0], res_x.m_value[2], res_y.m_value[0],
                  res_y.m_value[2]};
#endif
  }

  // Type-Conversions
  /// Convert to a human-readable string value.
  explicit inline operator std::string() const {
#ifdef USE_INTRINSICS
#if defined(__AVX512F__)
    // TODO: Matrix4x4 string conversion for AVX512.
#else
    std::array<FLOAT, 4> x{}, y{}, z{}, w{};
    m_value[0].m_value.get(x.data());
    m_value[1].m_value.get(y.data());
    m_value[2].m_value.get(z.data());
    m_value[3].m_value.get(w.data());
    return "[" + std::to_string(x[0]) + ", " + std::to_string(x[1]) + ", " +
           std::to_string(y[0]) + ", " + std::to_string(y[1]) + "\n" +
           std::to_string(x[2]) + ", " + std::to_string(x[3]) + ", " +
           std::to_string(y[2]) + ", " + std::to_string(y[3]) + "\n" +
           std::to_string(z[0]) + ", " + std::to_string(z[1]) + ", " +
           std::to_string(w[0]) + ", " + std::to_string(w[1]) + "\n" +
           std::to_string(z[2]) + ", " + std::to_string(z[3]) + ", " +
           std::to_string(w[2]) + ", " + std::to_string(w[3]) + "]";
#endif
#else
    return "[" + std::to_string(m_value[0].m_value[0]) + ", " +
           std::to_string(m_value[0].m_value[1]) + ", " +
           std::to_string(m_value[1].m_value[0]) + ", " +
           std::to_string(m_value[1].m_value[1]) + "\n" +
           std::to_string(m_value[0].m_value[2]) + ", " +
           std::to_string(m_value[0].m_value[3]) + ", " +
           std::to_string(m_value[1].m_value[2]) + ", " +
           std::to_string(m_value[1].m_value[3]) + "\n" +
           std::to_string(m_value[2].m_value[0]) + ", " +
           std::to_string(m_value[2].m_value[1]) + ", " +
           std::to_string(m_value[3].m_value[0]) + ", " +
           std::to_string(m_value[3].m_value[1]) + "\n" +
           std::to_string(m_value[2].m_value[2]) + ", " +
           std::to_string(m_value[2].m_value[3]) + ", " +
           std::to_string(m_value[3].m_value[2]) + ", " +
           std::to_string(m_value[3].m_value[3]) + "]";
#endif  // USE_INTRINSICS
  }

  friend inline std::ostream &operator<<(std::ostream &out,
                                         const Matrix4x4 &mat) {
    return out << static_cast<std::string>(mat);
  }

  // Functions

  /// Returns the transpose of the given 4x4 matrix.
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

  /// Calculates the determinant of the given 4x4 matrix.
  [[nodiscard]] inline FLOAT determinant() const {
#if defined(USE_INTRINSICS) && defined(__AVX512F__)
    // TODO: Matrix4x4  determinant for AVX512.
#else
    // TODO: Matrix4x4 optimize determinant calculation.
    // Block matrix determinant method
    // Refer:
    // https://en.wikipedia.org/wiki/Block_matrix#Block_matrix_determinant
    // | a    b    c    d |     | A   B |
    // | e    f    g    h | ==  | C   D |
    // | i    j    k    l |
    // | m    n    o    p |
    // Determinant = Determinant(A)*Determinant(D- C*Inverse(A)*B)
    //             = Determinant(D)*Determinant(A- C*Inverse(D)*B)
    //             = Det(A)*Det(D) + Det(B)*Det(C) - Trace(
    //             ((Adj(A)*B)*((Adj(D)*C))
    return m_value[0].determinant() * m_value[3].determinant() +
           m_value[1].determinant() * m_value[2].determinant() -
           Matrix2x2::trace((Matrix2x2::adjoint(m_value[0]) * m_value[1]) *
                            (Matrix2x2::adjoint(m_value[3]) * m_value[2]));
#endif
  }

  /// Returns the inverse of the given 4x4 matrix.
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

  /**
   * Creates a 3D scaling matrix from the given x, y and z axis scales.
   * @param scale 3D vector containing the x, y and z axis scales.
   * @return 4x4 Affine transform matrix containing only the 3D scale component.
   */
  static inline Matrix4x4 scale(const Vector &scale) {
    // | x  0  0  0 |
    // | 0  y  0  0 |
    // | 0  0  z  0 |
    // | 0  0  0  1 |
#if defined(USE_INTRINSICS) && defined(__AVX512F__)
    // TODO: Matrix4x4 scale for AVX512.
#else

#if defined(USE_INTRINSICS)
    std::array<FLOAT, 4> val{};
    scale.get(val.data());
    return Matrix4x4{std::array<Matrix2x2, 4>{
        Matrix2x2{Vector{val[0], 0.0, 0.0, val[1]}}, Matrix2x2{Vector{0.0}},
        Matrix2x2{Vector{0.0}}, Matrix2x2{Vector{val[2], 0.0, 0.0, 1.0}}}};
#else
    return Matrix4x4{std::array<Matrix2x2, 4>{
        Matrix2x2{Vector{scale[0], 0.0, 0.0, scale[1]}}, Matrix2x2{Vector{0.0}},
        Matrix2x2{Vector{0.0}}, Matrix2x2{Vector{scale[2], 0.0, 0.0, 1.0}}}};
#endif
#endif
  }

  /**
   * Creates a 3D x-axis rotation matrix from the given angle in degree.
   * @param angle A scalar angle for x-axis given in degree. The angle is
   * assumed positive for counter-clockwise rotation and negative for clockwise
   * rotation.
   * @return 4x4 Affine transform matrix containing only the x-axis rotational
   * component.
   */
  static inline Matrix4x4 rotationX(FLOAT angle) {
    // | 1   0      0    0 |
    // | 0  cosθ  -sinθ  0 |
    // | 0  sinθ   cosθ  0 |
    // | 0   0      0    1 |
#if defined(USE_INTRINSICS) && defined(__AVX512F__)
    // TODO: Matrix4x4 scale for AVX512.
#else
#ifdef USE_INTRINSICS
#ifdef USE_DOUBLE
#if defined(__AVX2__) || defined(__AVX__)

#else

#endif  // AVX INTRINSICS
#else
    __m128 theta = _mm_set_ps(0.0f, angle, 90.0f, angle);
    // [cos, 0, cos, 1]
    __m128 cos = _mm_cosd_ps(theta);
    // [sin, 1, sin, 0]
    __m128 sin = _mm_sind_ps(theta);
    // [sin, 1, -sin, 0]
    sin = _mm_xor_ps(sin, _mm_set_ps(0.0f, -0.0f, 0.0f, 0.0f));

    __m128 mat0 = _mm_shuffle_ps(cos, cos, _MM_SHUFFLE(0, 1, 1, 3));
    __m128 mat1 = _mm_shuffle_ps(sin, sin, _MM_SHUFFLE(3, 2, 3, 3));
    __m128 mat2 = _mm_shuffle_ps(sin, sin, _MM_SHUFFLE(3, 3, 0, 3));
    __m128 mat3 = _mm_shuffle_ps(cos, cos, _MM_SHUFFLE(3, 1, 1, 0));

    return Matrix4x4{std::array<Matrix2x2, 4>{
        Matrix2x2{Vector{mat0}}, Matrix2x2{Vector{mat1}},
        Matrix2x2{Vector{mat2}}, Matrix2x2{Vector{mat3}}}};
#endif  // USE_DOUBLE
#else
    const FLOAT cos = std::cos(radian(angle));
    const FLOAT sin = std::sin(radian(angle));
    return Matrix4x4{
        std::array<Matrix2x2, 4>{Matrix2x2{Vector{1.0, 0.0, 0.0, cos}},
                                 Matrix2x2{Vector{0.0, 0.0, -sin, 0.0}},
                                 Matrix2x2{Vector{0.0, sin, 0.0, 0.0}},
                                 Matrix2x2{Vector{cos, 0.0, 0.0, 1.0}}}};
#endif  // USE_INTRINSICS
#endif
  }

  /**
   * Creates a 3D y-axis rotation matrix from the given angle in degree.
   * @param angle A scalar angle for y-axis given in degree. The angle is
   * assumed positive for counter-clockwise rotation and negative for clockwise
   * rotation.
   * @return 4x4 Affine transform matrix containing only the y-axis rotational
   * component.
   */
  static inline Matrix4x4 rotationY(FLOAT angle) {
    // |  cosθ  0  sinθ  0 |
    // |   0    1   0    0 |
    // | -sinθ  0  cosθ  0 |
    // |   0    0   0    1 |
#if defined(USE_INTRINSICS) && defined(__AVX512F__)
    // TODO: Matrix4x4 scale for AVX512.
#else
#ifdef USE_INTRINSICS
#ifdef USE_DOUBLE
#if defined(__AVX2__) || defined(__AVX__)

#else

#endif  // AVX INTRINSICS
#else
    __m128 theta = _mm_set_ps(0.0f, angle, 90.0f, angle);
    // [cos, 0, cos, 1]
    __m128 cos = _mm_cosd_ps(theta);
    // [sin, 1, sin, 0]
    __m128 sin = _mm_sind_ps(theta);

    __m128 mat0 = _mm_shuffle_ps(cos, cos, _MM_SHUFFLE(3, 1, 1, 0));
    __m128 mat1 = _mm_shuffle_ps(sin, sin, _MM_SHUFFLE(3, 3, 3, 0));
    __m128 mat2 = _mm_xor_ps(mat1, _mm_set_ps(0.0f, 0.0f, 0.0f, -0.0f));

    return Matrix4x4{std::array<Matrix2x2, 4>{
        Matrix2x2{Vector{mat0}}, Matrix2x2{Vector{mat1}},
        Matrix2x2{Vector{mat2}}, Matrix2x2{Vector{mat0}}}};
#endif  // USE_DOUBLE
#else
    const FLOAT cos = std::cos(radian(angle));
    const FLOAT sin = std::sin(radian(angle));
    return Matrix4x4{
        std::array<Matrix2x2, 4>{Matrix2x2{Vector{cos, 0.0, 0.0, 1.0}},
                                 Matrix2x2{Vector{sin, 0.0, 0.0, 0.0}},
                                 Matrix2x2{Vector{-sin, 0.0, 0.0, 0.0}},
                                 Matrix2x2{Vector{cos, 0.0, 0.0, 1.0}}}};
#endif  // USE_INTRINSICS
#endif
  }

  /**
   * Creates a 3D z-axis rotation matrix from the given angle in degree.
   * @param angle A scalar angle for z-axis given in degree. The angle is
   * assumed positive for counter-clockwise rotation and negative for clockwise
   * rotation.
   * @return 4x4 Affine transform matrix containing only the z-axis rotational
   * component.
   */
  static inline Matrix4x4 rotationZ(FLOAT angle) {
    // | cosθ  -sinθ  0  0 |
    // | sinθ   cosθ  0  0 |
    // |  0      0    1  0 |
    // |  0      0    0  1 |
#if defined(USE_INTRINSICS) && defined(__AVX512F__)
    // TODO: Matrix4x4 scale for AVX512.
#else
#ifdef USE_INTRINSICS
#ifdef USE_DOUBLE
#if defined(__AVX2__) || defined(__AVX__)

#else

#endif  // AVX INTRINSICS
#else
    __m128 theta = _mm_set1_ps(angle);
    __m128 zero = _mm_setzero_ps();
    // [cos, cos, cos, cos]
    __m128 cos = _mm_cosd_ps(theta);
    // [sin, sin, sin, sin]
    __m128 sin = _mm_sind_ps(theta);

    // [cos, cos, sin, sin]
    __m128 mat0 = _mm_shuffle_ps(cos, sin, _MM_SHUFFLE(3, 2, 1, 0));
    // [cos, sin, sin, cos]
    mat0 = _mm_shuffle_ps(mat0, mat0, _MM_SHUFFLE(1, 3, 2, 0));
    // [cos, -sin, sin, cos]
    mat0 = _mm_xor_ps(mat0, _mm_set_ps(0.0f, 0.0f, -0.0f, 0.0f));

    return Matrix4x4{std::array<Matrix2x2, 4>{
        Matrix2x2{Vector{mat0}}, Matrix2x2{Vector{zero}},
        Matrix2x2{Vector{zero}},
        Matrix2x2{Vector{_mm_set_ps(1.0f, 0.0f, 0.0f, 1.0f)}}}};
#endif  // USE_DOUBLE
#else
    const FLOAT cos = std::cos(radian(angle));
    const FLOAT sin = std::sin(radian(angle));
    return Matrix4x4{std::array<Matrix2x2, 4>{
        Matrix2x2{Vector{cos, -sin, sin, cos}}, Matrix2x2{Vector{0.0}},
        Matrix2x2{Vector{0.0}}, Matrix2x2{Vector{1.0, 0.0, 0.0, 1.0}}}};
#endif  // USE_INTRINSICS
#endif
  }

  /**
   * Creates a 3D rotation matrix from the given euler angles in degree.
   * @param angles A 3D vector containing euler angles for each axis. The angle
   * is assumed positive for counter-clockwise rotation and negative for
   * clockwise rotation.
   * @return 4x4 Affine transform matrix containing only the 3D rotational
   * component.
   */
  static inline Matrix4x4 rotation(const Vector &euler_angles) {
    // R = roll*pitch*head (yaw)
    // R = R(z)*R(x)*R(y)

    // | cosZ*cosY-sinZ*sinX*sinY  -sinZ*cosX  cosZ*sinY+sinZ*sinX*cosY  0 |
    // | sinZ*cosY+cosZ*sinX*sinY   cosZ*cosX  sinZ*sinY-cosZ*sinX*cosY  0 |
    // |        -cosX*sinY            sinX            cosX*cosY          0 |
    // |             0                 0                  0              1 |

#if defined(USE_INTRINSICS) && defined(__AVX512F__)
    // TODO: Matrix4x4 scale for AVX512.
#else
#ifdef USE_INTRINSICS
#ifdef USE_DOUBLE
#if defined(__AVX2__) || defined(__AVX__)

#else

#endif  // AVX INTRINSICS
#else
    // [cosX, cosY, cosZ, 1]
    __m128 cos = _mm_cosd_ps(euler_angles.m_value);
    // [sinX, sinY, sinZ, 0]
    __m128 sin = _mm_sind_ps(euler_angles.m_value);

    // [cosZ, cosZ, sinZ, sinZ]
    __m128 tmp = _mm_shuffle_ps(cos, sin, _MM_SHUFFLE(2, 2, 2, 2));
    // [cosZ, sinZ, sinZ, cosZ]
    tmp = _mm_shuffle_ps(tmp, tmp, _MM_SHUFFLE(1, 3, 2, 0));
    // [cosZ, -sinZ, sinZ, -cosZ]
    tmp = _mm_xor_ps(tmp, _mm_set_ps(-0.0f, 0.0f, -0.0f, 0.0f));
    // [cosZ*sinX, -sinZ*sinX, sinZ*sinX, -cosZ*sinX]
    tmp = _mm_mul_ps(tmp, _mm_shuffle_ps(sin, sin, _MM_SHUFFLE(0, 0, 0, 0)));
    // [cosZ*sinX*sinY, -sinZ*sinX*sinY, sinZ*sinX*cosY, -cosZ*sinX*cosY]
    tmp = _mm_mul_ps(tmp, _mm_shuffle_ps(sin, cos, _MM_SHUFFLE(1, 1, 1, 1)));

    // [-sinZ*sinX*sinY, cosZ*sinX*sinY, sinZ*sinX*cosY, -cosZ*sinX*cosY]
    tmp = _mm_shuffle_ps(tmp, tmp, _MM_SHUFFLE(3, 2, 0, 1));

    // [cosZ, cosZ, sinZ, sinZ]
    __m128 tmp0 = _mm_shuffle_ps(cos, sin, _MM_SHUFFLE(2, 2, 2, 2));
    // [cosZ, sinZ, cosZ, sinZ]
    tmp0 = _mm_shuffle_ps(tmp0, tmp0, _MM_SHUFFLE(3, 1, 2, 0));

    // [cosY, cosY, sinY, sinY]
    __m128 tmp1 = _mm_shuffle_ps(cos, sin, _MM_SHUFFLE(1, 1, 1, 1));

    // [cosZ*cosY, sinZ*cosY, cosZ*sinY, sinZ*sinY]
    __m128 ele0 = _mm_mul_ps(tmp0, tmp1);
    // [cosZ*cosY-sinZ*sinX*sinY, sinZ*cosY+cosZ*sinX*sinY,
    // cosZ*sinY+sinZ*sinX*cosY, sinZ*sinY-cosZ*sinX*cosY]
    ele0 = _mm_add_ps(ele0, tmp);

    // [cosX, cosX, cosX, cosX]
    tmp0 = _mm_shuffle_ps(cos, cos, _MM_SHUFFLE(0, 0, 0, 0));
    // [cosX, cosX, -cosX, -cosX]
    tmp0 = _mm_xor_ps(tmp0, _mm_set_ps(-0.0f, -0.0f, 0.0f, 0.0f));

    // [cosY, cosZ, sinY, sinZ]
    tmp1 = _mm_shuffle_ps(cos, sin, _MM_SHUFFLE(2, 1, 2, 1));

    // [cosX*cosY, cosX*cosZ, -cosX*sinY, -cosX*sinZ]
    __m128 ele1 = _mm_mul_ps(tmp0, tmp1);

    // [cosZ*cosY-sinZ*sinX*sinY, sinZ*cosY+cosZ*sinX*sinY, -cosX*sinZ,
    // cosX*cosZ]
    tmp0 = _mm_shuffle_ps(ele0, ele1, _MM_SHUFFLE(1, 3, 1, 0));
    // [cosZ*cosY-sinZ*sinX*sinY, -cosX*sinZ, sinZ*cosY+cosZ*sinX*sinY,
    // cosX*cosZ]
    tmp0 = _mm_shuffle_ps(tmp0, tmp0, _MM_SHUFFLE(3, 1, 2, 0));

    // [cosZ*sinY+sinZ*sinX*cosY, sinZ*sinY-cosZ*sinX*cosY, 0, 0]
    tmp1 = _mm_shuffle_ps(ele0, _mm_setzero_ps(), _MM_SHUFFLE(0, 0, 3, 2));
    // [cosZ*sinY+sinZ*sinX*cosY, 0, sinZ*sinY-cosZ*sinX*cosY, 0]
    tmp1 = _mm_shuffle_ps(tmp1, tmp1, _MM_SHUFFLE(3, 1, 2, 0));

    // [-cosX*sinY, -cosX*sinY, sinX, 0]
    ele0 = _mm_shuffle_ps(ele1, sin, _MM_SHUFFLE(3, 0, 2, 2));
    // [-cosX*sinY, sinX, 0, 0]
    ele0 = _mm_shuffle_ps(ele0, ele0, _MM_SHUFFLE(3, 3, 2, 0));

    // [cosX*cosY, 0, 0, 0]
    ele1 = _mm_move_ss(_mm_setzero_ps(), ele1);
    // [cosX*cosY, 0, 1, 1]
    ele1 = _mm_shuffle_ps(ele1, cos, _MM_SHUFFLE(3, 3, 1, 0));
    // [cosX*cosY, 0, 0, 1]
    ele1 = _mm_shuffle_ps(ele1, ele1, _MM_SHUFFLE(3, 1, 1, 0));

    return Matrix4x4{std::array<Matrix2x2, 4>{
        Matrix2x2{Vector{tmp0}}, Matrix2x2{Vector{tmp1}},
        Matrix2x2{Vector{ele0}}, Matrix2x2{Vector{ele1}}}};
#endif  // USE_DOUBLE
#else
    const std::array<FLOAT, 3> cos{std::cos(radian(euler_angles[0])),
                                   std::cos(radian(euler_angles[1])),
                                   std::cos(radian(euler_angles[2]))};
    const std::array<FLOAT, 3> sin{std::sin(radian(euler_angles[0])),
                                   std::sin(radian(euler_angles[1])),
                                   std::sin(radian(euler_angles[2]))};
    return Matrix4x4{std::array<Matrix2x2, 4>{
        Matrix2x2{Vector{
            cos[2] * cos[1] - sin[2] * sin[0] * sin[1], -sin[2] * cos[0],
            sin[2] * cos[1] + cos[2] * sin[0] * sin[1], cos[2] * cos[0]}},
        Matrix2x2{Vector{cos[2] * sin[1] + sin[2] * sin[0] * cos[1], 0.0,
                         sin[2] * sin[1] - cos[2] * sin[0] * cos[1], 0.0}},
        Matrix2x2{Vector{-cos[0] * sin[1], sin[0], 0.0, 0.0}},
        Matrix2x2{Vector{cos[0] * cos[1], 0.0, 0.0, 1.0}}}};
#endif  // USE_INTRINSICS
#endif
  }

  /**
   * Creates a 3D rotation matrix describing a rotation over an arbitrary axis.
   * @param angle The scalar angle by which to rotate about the given axis. The
   * angle is assumed positive for counter-clockwise rotation and negative for
   * clockwise rotation.
   * @param axis A 3D normalized vector describing the axis around which to
   * rotate.
   * @return 4x4 Affine transform matrix containing only the 3D rotational
   * component.
   */
  static inline Matrix4x4 rotationOverAxis(FLOAT angle, const Vector &axis) {
    // Refer: https://dl.acm.org/doi/10.5555/90767.90908
    // |    cos+mcos*a_x²       mcos*a_x*a_y-a_z*sin mcos*a_x*a_z+a_y*sin 0 |
    // | mcos*a_x*a_y + a_z*sin    cos+mcos*a_y²     mcos*a_y*a_z-a_x*sin 0 |
    // | mcos*a_x*a_z - a_y*sin mcos*a_y*a_z+a_x*sin    cos+mcos*a_z²     0 |
    // |           0                     0                  0             1 |
#if defined(USE_INTRINSICS) && defined(__AVX512F__)
    // TODO: Matrix4x4 scale for AVX512.
#else
#ifdef USE_INTRINSICS
#ifdef USE_DOUBLE
#if defined(__AVX2__) || defined(__AVX__)

#else

#endif  // AVX INTRINSICS
#else
    __m128 theta = _mm_set1_ps(angle);
    // [cos, cos, cos, cos]
    __m128 cos = _mm_cosd_ps(theta);
    // [sin, sin, sin, sin]
    __m128 sin = _mm_sind_ps(theta);
    // [1 - cos, 1 - cos, 1 - cos, 1 - cos]
    __m128 om_cos = _mm_sub_ps(_mm_set1_ps(1.0f), cos);

    // [a_x², a_y², a_z², 0]
    __m128 diag = _mm_mul_ps(axis.m_value, axis.m_value);
    // [(1 - cos)a_x², (1 - cos)a_y², (1 - cos)a_z², 0]
    diag = _mm_mul_ps(om_cos, diag);
    // [cos + (1 - cos)a_x², cos + (1 - cos)a_y², cos + (1 - cos)a_z², 0]
    diag = _mm_add_ps(cos, diag);

    // [(1-cos)a_x, (1-cos)a_y, (1-cos)a_z, 0]
    om_cos = _mm_mul_ps(om_cos, axis.m_value);
    // [(1-cos)a_x, (1-cos)a_y, (1-cos)a_z, 0]
    om_cos = _mm_shuffle_ps(om_cos, om_cos, _MM_SHUFFLE(3, 2, 1, 0));
    // [(1-cos)a_x*a_y, (1-cos)a_y*a_z, (1-cos)a_x*a_z, 0]
    om_cos = _mm_mul_ps(om_cos, _mm_shuffle_ps(axis.m_value, axis.m_value,
                                               _MM_SHUFFLE(3, 0, 2, 1)));

    // [a_x*sin, a_y*sin, a_z*sin, 0]
    cos = _mm_mul_ps(axis.m_value, sin);
    // [a_z*sin, a_x*sin, a_y*sin, 0]
    cos = _mm_shuffle_ps(cos, cos, _MM_SHUFFLE(3, 1, 0, 2));
    // [-a_z*sin, -a_x*sin, -a_y*sin, 0]
    sin = _mm_xor_ps(cos, _mm_set_ps(0.0f, -0.0f, -0.0f, -0.0f));

    // [(1-cos)a_x*a_y + a_z*sin, (1-cos)a_y*a_z + a_x*sin, (1-cos)a_x*a_z +
    // a_y*sin, 0]
    cos = _mm_add_ps(om_cos, cos);
    // [(1-cos)a_x*a_y - a_z*sin, (1-cos)a_y*a_z - a_x*sin, (1-cos)a_x*a_z -
    // a_y*sin, 0]
    sin = _mm_add_ps(om_cos, sin);

    // [cos + (1 - cos)a_x², cos + (1 - cos)a_y², (1-cos)a_x*a_y + a_z*sin, 0]
    theta = _mm_shuffle_ps(diag, cos, _MM_SHUFFLE(3, 0, 1, 0));
    // [0, cos + (1 - cos)a_x², cos + (1 - cos)a_y², (1-cos)a_x*a_y + a_z*sin]
    theta = _mm_shuffle_ps(theta, theta, _MM_SHUFFLE(2, 1, 0, 3));
    // [(1-cos)a_x*a_y - a_z*sin, cos + (1 - cos)a_x², cos + (1 - cos)a_y²,
    // (1-cos)a_x*a_y + a_z*sin]
    theta = _mm_move_ss(theta, sin);
    // [cos + (1 - cos)a_x², (1-cos)a_x*a_y - a_z*sin, (1-cos)a_x*a_y + a_z*sin,
    // cos + (1 - cos)a_y²]
    theta = _mm_shuffle_ps(theta, theta, _MM_SHUFFLE(2, 3, 0, 1));

    // [(1-cos)a_x*a_z + a_y*sin, 0, (1-cos)a_y*a_z - a_x*sin, 0]
    om_cos = _mm_shuffle_ps(cos, sin, _MM_SHUFFLE(3, 1, 3, 2));

    // [(1-cos)a_y*a_z + a_x*sin, 0, (1-cos)a_x*a_z - a_y*sin, 0]
    cos = _mm_shuffle_ps(cos, sin, _MM_SHUFFLE(3, 2, 3, 1));
    // [(1-cos)a_x*a_z - a_y*sin, (1-cos)a_y*a_z + a_x*sin, 0, 0]
    cos = _mm_shuffle_ps(cos, cos, _MM_SHUFFLE(3, 1, 0, 2));

    // [cos + (1 - cos)a_z², 0, 0, 0]
    diag = _mm_shuffle_ps(diag, diag, _MM_SHUFFLE(3, 3, 3, 2));
    // [cos + (1 - cos)a_z², 0, 0, 1]
    diag = _mm_shuffle_ps(diag, _mm_set_ps(1.0f, 0.0f, 0.0f, 0.0f),
                          _MM_SHUFFLE(3, 2, 1, 0));

    return Matrix4x4{std::array<Matrix2x2, 4>{
        Matrix2x2{Vector{theta}}, Matrix2x2{Vector{om_cos}},
        Matrix2x2{Vector{cos}}, Matrix2x2{Vector{diag}}}};
#endif  // USE_DOUBLE
#else
    const FLOAT cos = std::cos(angle);
    const FLOAT om_cos = static_cast<FLOAT>(1.0) - cos;
    const FLOAT sin = std::sin(angle);
    return Matrix4x4{std::array<Matrix2x2, 4>{
        Matrix2x2{Vector{cos + om_cos * axis[0] * axis[0],
                         om_cos * axis[0] * axis[1] - axis[2] * sin,
                         om_cos * axis[0] * axis[1] + axis[2] * sin,
                         cos + om_cos * axis[1] * axis[1]}},
        Matrix2x2{Vector{om_cos * axis[0] * axis[2] + axis[1] * sin, 0.0,
                         om_cos * axis[1] * axis[2] - axis[0] * sin, 0.0}},
        Matrix2x2{Vector{om_cos * axis[0] * axis[2] - axis[1] * sin,
                         om_cos * axis[1] * axis[2] + axis[0] * sin, 0.0, 0.0}},
        Matrix2x2{Vector{cos + om_cos * axis[2] * axis[2], 0.0, 0.0, 1.0}}}};
#endif  // USE_INTRINSICS
#endif
  }

  /**
   * Creates a 3D translation matrix from the given x, y and z axis values.
   * @param translate 3D vector containing the x, y and z axis translation
   * values.
   * @return 4x4 Affine transform matrix containing only the 3D translation
   * component.
   */
  static inline Matrix4x4 translation(const Vector &translate) {
    // | 1  0  0  x |
    // | 0  1  0  y |
    // | 0  0  1  z |
    // | 0  0  0  1 |
#if defined(USE_INTRINSICS) && defined(__AVX512F__)
    // TODO: Matrix4x4 scale for AVX512.
#else
#ifdef USE_INTRINSICS
    std::array<FLOAT, 4> val{};
    translate.get(val.data());
    return Matrix4x4{std::array<Matrix2x2, 4>{
        Matrix2x2{Vector{1.0, 0.0, 0.0, 1.0}},
        Matrix2x2{Vector{0.0, val[0], 0.0, val[1]}}, Matrix2x2{Vector{0.0}},
        Matrix2x2{Vector{1.0, val[2], 0.0, 1.0}}}};
#else
    return Matrix4x4{std::array<Matrix2x2, 4>{
        Matrix2x2{Vector{1.0, 0.0, 0.0, 1.0}},
        Matrix2x2{Vector{0.0, translate[0], 0.0, translate[1]}},
        Matrix2x2{Vector{0.0}},
        Matrix2x2{Vector{1.0, translate[2], 0.0, 1.0}}}};
#endif  // USE_INTRINSICS
#endif
  }

  /**
   * Returns the inverse of the given 4x4 affine transformation matrix.
   * The matrix should have only translation, rotation and scaling components.
   * @param mat The transform matrix to invert.
   * @return 4x4 inverse affine transformation matrix
   */
  static inline Matrix4x4 transformInverse(const Matrix4x4 &matrix) {
    // | x_x  y_x  z_x  t_x |
    // | x_y  y_y  z_y  t_y |
    // | x_z  y_z  z_z  t_z |
    // |  0    0    0    1  |
    // Inverse =
    // | 1/|x|²*x_x  1/|x|²*x_y  1/|x|²*x_z  dot(-t, 1/|x|²*x) |
    // | 1/|y|²*y_x  1/|y|²*y_y  1/|y|²*y_z  dot(-t, 1/|y|²*y) |
    // | 1/|z|²*z_x  1/|z|²*z_y  1/|z|²*z_z  dot(-t, 1/|z|²*z) |
    // |     0           0           0                1        |
    // =
    // | x_x'  x_y'  x_z' dot(-t, x') |
    // | y_x'  y_y'  y_z' dot(-t, y') |
    // | z_x'  z_y'  z_z' dot(-t, z') |
    // |  0     0     0        1      |
#if defined(USE_INTRINSICS) && defined(__AVX512F__)
    // TODO: Matrix4x4 scale for AVX512.
#else
#ifdef USE_INTRINSICS
#ifdef USE_DOUBLE
#if defined(__AVX2__) || defined(__AVX__)

#else

#endif  // AVX INTRINSICS
#else
    // [x_x, y_x, z_x, t_x]
    __m128 x = _mm_shuffle_ps(matrix.m_value[0].m_value.m_value,
                              matrix.m_value[1].m_value.m_value,
                              _MM_SHUFFLE(1, 0, 1, 0));
    // [x_y, y_y, z_y, t_y]
    __m128 y = _mm_shuffle_ps(matrix.m_value[0].m_value.m_value,
                              matrix.m_value[1].m_value.m_value,
                              _MM_SHUFFLE(3, 2, 3, 2));
    // [x_z, y_z, z_z, t_z]
    __m128 z = _mm_shuffle_ps(matrix.m_value[2].m_value.m_value,
                              matrix.m_value[3].m_value.m_value,
                              _MM_SHUFFLE(1, 0, 1, 0));

    // [x_x², y_x², z_x², t_x²]
    __m128 scale = _mm_mul_ps(x, x);
    // [x_x² + x_y², y_x² + y_u², z_x² + z_y², t_x² + t_y²]
    scale = _mm_add_ps(scale, _mm_mul_ps(y, y));
    // [|x|², |y|², |z|², |t|²]
    scale = _mm_add_ps(scale, _mm_mul_ps(z, z));
    scale = _mm_blendv_ps(scale, _mm_set1_ps(1.0f),
                          _mm_cmplt_ps(scale, _mm_set1_ps(EPSILON)));
    // [1/|x|², 1/|y|², 1/|z|², 0]
    scale = _mm_div_ps(_mm_set_ps(0.0f, 1.0f, 1.0f, 1.0f), scale);

    // [x_x', y_x', z_x', 0]
    x = _mm_mul_ps(scale, x);
    // [x_y', y_y', z_y', 0]
    y = _mm_mul_ps(scale, y);
    // [x_z', y_z', z_z', 0]
    z = _mm_mul_ps(scale, z);

    // [t_x, t_y, t_z, 0]
    __m128 t = _mm_shuffle_ps(matrix.m_value[1].m_value.m_value,
                              matrix.m_value[3].m_value.m_value,
                              _MM_SHUFFLE(2, 1, 3, 1));
    // [-t_x, -t_y, -t_z, 0]
    t = _mm_xor_ps(t, _mm_set1_ps(-0.0f));

    // [-t_x*x_x', -t_x*y_x', -t_x*z_x', 0]
    __m128 tmp = _mm_mul_ps(x, _mm_shuffle_ps(t, t, _MM_SHUFFLE(3, 0, 0, 0)));
    // [-t_x*x_x' + -t_y*x_y', -t_x*y_x' + -t_y*y_y', -t_x*z_x' + -t_y*z_y', 0]
    tmp = _mm_add_ps(
        tmp, _mm_mul_ps(y, _mm_shuffle_ps(t, t, _MM_SHUFFLE(3, 1, 1, 1))));
    // [dot(-t, x'), dot(-t, y'), dot(-t, z'), 0]
    t = _mm_add_ps(
        tmp, _mm_mul_ps(z, _mm_shuffle_ps(t, t, _MM_SHUFFLE(3, 2, 2, 2))));

    // [x_x', y_x', x_y', y_y']
    scale = _mm_shuffle_ps(x, y, _MM_SHUFFLE(1, 0, 1, 0));
    tmp = x;
    // [x_x', x_y', y_x', y_y']
    x = _mm_shuffle_ps(scale, scale, _MM_SHUFFLE(3, 1, 2, 0));

    // [z_x', 0, z_y', 0]
    scale = _mm_shuffle_ps(tmp, y, _MM_SHUFFLE(3, 2, 3, 2));
    tmp = z;
    // [z_x', z_y, 0, 0]
    z = _mm_shuffle_ps(scale, scale, _MM_SHUFFLE(3, 1, 2, 0));

    // [x_z', y_z', dot(-t, x'), dot(-t, y')]
    scale = _mm_shuffle_ps(tmp, t, _MM_SHUFFLE(1, 0, 1, 0));
    // [x_z', dot(-t, x'), y_z', dot(-t, y')]
    y = _mm_shuffle_ps(scale, scale, _MM_SHUFFLE(3, 1, 2, 0));

    // [dot(-t, z'), 0, z_z', 0]
    t = _mm_shuffle_ps(t, tmp, _MM_SHUFFLE(3, 2, 3, 2));
    // [z_z', dot(-t, z'), 0, 0]
    t = _mm_shuffle_ps(t, t, _MM_SHUFFLE(3, 1, 0, 2));
    // [z_z', dot(-t, z'), 0, 1]
    t = _mm_shuffle_ps(t, _mm_set_ps(1.0f, 0.0f, 0.0f, 0.0f),
                       _MM_SHUFFLE(3, 2, 1, 0));

    return Matrix4x4{
        std::array<Matrix2x2, 4>{Matrix2x2{Vector{x}}, Matrix2x2{Vector{y}},
                                 Matrix2x2{Vector{z}}, Matrix2x2{Vector{t}}}};
#endif  // USE_DOUBLE
#else
    Vector x_axis{matrix.m_value[0].m_value[0], matrix.m_value[0].m_value[2],
                  matrix.m_value[2].m_value[0], 0.0};
    Vector y_axis{matrix.m_value[0].m_value[1], matrix.m_value[0].m_value[3],
                  matrix.m_value[2].m_value[1], 0.0};
    Vector z_axis{matrix.m_value[1].m_value[0], matrix.m_value[1].m_value[2],
                  matrix.m_value[3].m_value[0], 0.0};
    Vector minus_t{-matrix.m_value[1].m_value[1], -matrix.m_value[1].m_value[3],
                   -matrix.m_value[3].m_value[1], 1.0};
    const FLOAT x_inv_scale = static_cast<FLOAT>(1.0) / x_axis.lengthSquared();
    const FLOAT y_inv_scale = static_cast<FLOAT>(1.0) / y_axis.lengthSquared();
    const FLOAT z_inv_scale = static_cast<FLOAT>(1.0) / z_axis.lengthSquared();

    x_axis *= x_inv_scale;
    y_axis *= y_inv_scale;
    z_axis *= z_inv_scale;
    minus_t = Vector{Vector::dot(minus_t, x_axis), Vector::dot(minus_t, y_axis),
                     Vector::dot(minus_t, z_axis), 1.0};

    return Matrix4x4{std::array<Matrix2x2, 4>{
        Matrix2x2{Vector{x_axis[0], x_axis[1], y_axis[0], y_axis[1]}},
        Matrix2x2{Vector{x_axis[2], minus_t[0], y_axis[2], minus_t[1]}},
        Matrix2x2{Vector{z_axis[0], z_axis[1], 0.0, 0.0}},
        Matrix2x2{Vector{z_axis[2], minus_t[2], 0.0, 1.0}}}};
#endif  // USE_INTRINSICS
#endif
  }

  /**
   * Returns the inverse of the given 4x4 affine transformation matrix with unit
   * scale. The matrix should have only translation and rotation components.
   * @param mat The transform matrix to invert.
   * @return 4x4 inverse affine transformation matrix
   */
  static inline Matrix4x4 transformInverseUnitScale(const Matrix4x4 &matrix) {
// | x_x  y_x  z_x  t_x |
// | x_y  y_y  z_y  t_y |
// | x_z  y_z  z_z  t_z |
// |  0    0    0    1  |
// Inverse =
// | x_x  x_y  x_z dot(-t, x) |
// | y_x  y_y  y_z dot(-t, y) |
// | z_x  z_y  z_z dot(-t, z) |
// |  0    0    0       1     |
#if defined(USE_INTRINSICS) && defined(__AVX512F__)
// TODO: Matrix4x4 scale for AVX512.
#else
#ifdef USE_INTRINSICS
#ifdef USE_DOUBLE
#if defined(__AVX2__) || defined(__AVX__)

#else

#endif  // AVX INTRINSICS
#else
    // [x_x, y_x, z_x, t_x]
    __m128 x = _mm_shuffle_ps(matrix.m_value[0].m_value.m_value,
                              matrix.m_value[1].m_value.m_value,
                              _MM_SHUFFLE(1, 0, 1, 0));
    // [x_y, y_y, z_y, t_y]
    __m128 y = _mm_shuffle_ps(matrix.m_value[0].m_value.m_value,
                              matrix.m_value[1].m_value.m_value,
                              _MM_SHUFFLE(3, 2, 3, 2));
    // [x_z, y_z, z_z, t_z]
    __m128 z = _mm_shuffle_ps(matrix.m_value[2].m_value.m_value,
                              matrix.m_value[3].m_value.m_value,
                              _MM_SHUFFLE(1, 0, 1, 0));
    // [t_x, t_y, t_z, 0]
    __m128 t = _mm_shuffle_ps(matrix.m_value[1].m_value.m_value,
                              matrix.m_value[3].m_value.m_value,
                              _MM_SHUFFLE(2, 1, 3, 1));
    // [-t_x, -t_y, -t_z, 0]
    t = _mm_xor_ps(t, _mm_set1_ps(-0.0f));

    // [-t_x*x_x, -t_x*y_x, -t_x*z_x, 0]
    __m128 tmp = _mm_mul_ps(x, _mm_shuffle_ps(t, t, _MM_SHUFFLE(3, 0, 0, 0)));
    // [-t_x*x_x + -t_y*x_y, -t_x*y_x + -t_y*y_y, -t_x*z_x + -t_y*z_y, 0]
    tmp = _mm_add_ps(tmp, _mm_mul_ps(y, _mm_shuffle_ps(t, t, _MM_SHUFFLE(3, 1, 1, 1))));
    // [dot(-t, x), dot(-t, y), dot(-t, z), 0]
    t = _mm_add_ps(tmp, _mm_mul_ps(z, _mm_shuffle_ps(t, t, _MM_SHUFFLE(3, 2, 2, 2))));

    tmp = x;
    // [x_x, y_x, x_y, y_y]
    x = _mm_shuffle_ps(x, y, _MM_SHUFFLE(1, 0, 1, 0));
    // [x_x, x_y, y_x, y_y]
    x = _mm_shuffle_ps(x, x, _MM_SHUFFLE(3, 1, 2, 0));

    // [z_x, t_x, z_y, t_y]
    tmp = _mm_shuffle_ps(tmp, y, _MM_SHUFFLE(3, 2, 3, 2));
    // [z_x, z_y, t_x, t_y]
    tmp = _mm_shuffle_ps(tmp, tmp, _MM_SHUFFLE(3, 1, 2, 0));
    // [z_x, z_y, 0, 0]
    tmp = _mm_movelh_ps(tmp, _mm_setzero_ps());

    // [x_z, y_z, dot(-t, x), dot(-t, y)]
    y = _mm_shuffle_ps(z, t, _MM_SHUFFLE(1, 0, 1, 0));
    // [x_z, dot(-t, x), y_z, dot(-t, y)]
    y = _mm_shuffle_ps(y, y, _MM_SHUFFLE(3, 1, 2, 0));

    // [dot(-t, z), 0, z_z, t_z]
    t = _mm_shuffle_ps(t, z, _MM_SHUFFLE(3, 2, 3, 2));
    // [z_z, dot(-t, z), 0, 0]
    t = _mm_shuffle_ps(t, t, _MM_SHUFFLE(1, 1, 0, 2));
    // [z_z, dot(-t, z), 0, 1]
    t = _mm_shuffle_ps(t, _mm_set_ps(1.0f, 0.0f, 0.0f, 0.0f),
                       _MM_SHUFFLE(3, 2, 1, 0));

    // [z_x, 0, z_y, 0]
    z = tmp;

    return Matrix4x4{
        std::array<Matrix2x2, 4>{Matrix2x2{Vector{x}}, Matrix2x2{Vector{y}},
                                 Matrix2x2{Vector{z}}, Matrix2x2{Vector{t}}}};
#endif  // USE_DOUBLE
#else
    Vector x_axis{matrix.m_value[0].m_value[0], matrix.m_value[0].m_value[2],
                  matrix.m_value[2].m_value[0], 0.0};
    Vector y_axis{matrix.m_value[0].m_value[1], matrix.m_value[0].m_value[3],
                  matrix.m_value[2].m_value[1], 0.0};
    Vector z_axis{matrix.m_value[1].m_value[0], matrix.m_value[1].m_value[2],
                  matrix.m_value[3].m_value[0], 0.0};
    Vector minus_t{-matrix.m_value[1].m_value[1], -matrix.m_value[1].m_value[3],
                   -matrix.m_value[3].m_value[1], 1.0};

    minus_t = Vector{Vector::dot(minus_t, x_axis), Vector::dot(minus_t, y_axis),
                     Vector::dot(minus_t, z_axis), 1.0};

    return Matrix4x4{std::array<Matrix2x2, 4>{
        Matrix2x2{Vector{x_axis[0], x_axis[1], y_axis[0], y_axis[1]}},
        Matrix2x2{Vector{x_axis[2], minus_t[0], y_axis[2], minus_t[1]}},
        Matrix2x2{Vector{z_axis[0], z_axis[1], 0.0, 0.0}},
        Matrix2x2{Vector{z_axis[2], minus_t[2], 0.0, 1.0}}}};
#endif  // USE_INTRINSICS
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
