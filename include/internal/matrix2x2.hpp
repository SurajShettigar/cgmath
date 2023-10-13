// Copyright 2023 Suraj Shettigar
// SPDX-License-Identifier: Apache-2.0

#ifndef CGMATH_INTERNAL_MATRIX2X2_HPP
#define CGMATH_INTERNAL_MATRIX2X2_HPP

#include "vector.hpp"

namespace cgmath::internal {

/// 2x2 Matrix using Column Major Mathematical Notation and stored in a Row
/// Major memory layout.
class Matrix2x2 {
 public:
  // Constructors / Destructors

  /// Constructs an identity 2x2 matrix.
  explicit Matrix2x2() : m_value{Vector{1.0, 0.0, 0.0, 1.0}} {};
  /// Constructs a 2x2 matrix with x-axis and y-axis column values.
  explicit Matrix2x2(FLOAT x_x, FLOAT x_y, FLOAT y_x, FLOAT y_y)
      : m_value{x_x, y_x, x_y, y_y} {};
  /// Constructs a 2x2 matrix with 2D x-axis and y-axis array values.
  explicit Matrix2x2(const std::array<FLOAT, 2> &x_axis,
                     const std::array<FLOAT, 2> &y_axis)
      : m_value{x_axis[0], y_axis[0], x_axis[1], y_axis[1]} {};
  /// Constructs a 2x2 matrix with 2D x-axis and y-axis array values.
  explicit Matrix2x2(const Vector &x_axis, const Vector &y_axis) {
#ifdef USE_INTRINSICS
#ifdef USE_DOUBLE
#if defined(__AVX512F__) || defined(__AVX2__) || defined(__AVX__)
    __m128d tmp = _mm256_castpd256_pd128(y_axis.m_value);
    m_value.m_value = _mm256_insertf128_pd(x_axis.m_value, tmp, 0b1);
    m_value.m_value = _mm256_permute4x64_pd (m_value.m_value, _MM_SHUFFLE(3, 1, 2, 0));
#else
    m_value.m_value[0] = _mm_shuffle_pd(x_axis.m_value[0], y_axis.m_value[0], _MM_SHUFFLE2(0, 0));
    m_value.m_value[1] = _mm_shuffle_pd(x_axis.m_value[0], y_axis.m_value[0], _MM_SHUFFLE2(1, 1));
#endif  // AVX INTRINSICS
#else
    m_value.m_value = _mm_shuffle_ps(x_axis.m_value, y_axis.m_value,
                                     _MM_SHUFFLE(1, 0, 1, 0));
    m_value.m_value = _mm_shuffle_ps(m_value.m_value, m_value.m_value,
                                     _MM_SHUFFLE(3, 1, 2, 0));
#endif  // USE_DOUBLE
#else
    m_value = Vector{x_axis[0], y_axis[0], x_axis[1], y_axis[1]};
#endif  // USE_INTRINSICS
  };

  ~Matrix2x2() = default;

  // Copy
  Matrix2x2(const Matrix2x2 &matrix) {
    if (this != &matrix) m_value = matrix.m_value;
  }
  Matrix2x2 &operator=(const Matrix2x2 &matrix) {
    if (this != &matrix) m_value = matrix.m_value;
    return *this;
  }

  // Getters / Setters
  /**
   * Returns the matrix row. [x_x, y_x, 0, 0] or [x_y, y_y, 0, 0]
   * @param index Row index.
   * @return 2D vector containing given row values.
   */
  [[nodiscard]] inline Vector operator[](size_t index) const {
    std::array<FLOAT, 4> val{};
    m_value.get(val.data());
    return Vector{val[index * 2], val[index * 2 + 1], 0.0, 0.0};
  }

  /**
   * Returns the matrix row. [x_x, y_x, 0, 0] or [x_y, y_y, 0, 0]
   * @param index Row index.
   * @return 2D vector containing given row values.
   */
  [[nodiscard]] inline Vector getRow(size_t index) const {
    return (*this)[index];
  }

  /**
   * Returns the matrix column. [x_x, x_y, 0, 0] or [y_x, y_y, 0, 0]
   * @param index Column index.
   * @return 2D vector containing given column values.
   */
  [[nodiscard]] inline Vector getColumn(size_t index) const {
    std::array<FLOAT, 4> val{};
    m_value.get(val.data());
    return Vector{val[index], val[index + 2], 0.0, 0.0};
  }

  /// Returns x-axis (column 0) of the matrix.
  [[nodiscard]] inline Vector getXAxis() const { return getColumn(0); }
  /// Returns y-axis (column 1) of the matrix.
  [[nodiscard]] inline Vector getYAxis() const { return getColumn(1); }

  /// Sets the x-axis column with the given values.
  inline void setXAxis(FLOAT x_x, FLOAT x_y) {
    m_value.setX(x_x);
    m_value.setZ(x_y);
  }
  /// Sets the x-axis column with the given 2D array values.
  inline void setXAxis(const std::array<FLOAT, 2> &x_axis) {
    setXAxis(x_axis[0], x_axis[1]);
  }
  /// Sets the x-axis column with the given 2D vector values.
  inline void setXAxis(const Vector &x_axis) {
#ifdef USE_INTRINSICS
#ifdef USE_DOUBLE
#if defined(__AVX512F__) || defined(__AVX2__) || defined(__AVX__)
    __m256d tmp =
        _mm256_permute4x64_pd(x_axis.m_value, _MM_SHUFFLE(3, 1, 2, 0));
    m_value.m_value = _mm256_shuffle_pd(tmp, m_value.m_value, 0b1010);
#else
    m_value.m_value[0] = _mm_shuffle_pd(x_axis.m_value[0], m_value.m_value[0],
                                        _MM_SHUFFLE2(1, 0));
    m_value.m_value[1] = _mm_shuffle_pd(x_axis.m_value[0], m_value.m_value[1],
                                        _MM_SHUFFLE2(1, 1));
#endif  // AVX INTRINSICS
#else
    m_value.m_value = _mm_shuffle_ps(m_value.m_value, m_value.m_value,
                                     _MM_SHUFFLE(2, 0, 3, 1));
    m_value.m_value = _mm_movelh_ps(x_axis.m_value, m_value.m_value);
    m_value.m_value = _mm_shuffle_ps(m_value.m_value, m_value.m_value,
                                     _MM_SHUFFLE(3, 1, 2, 0));
#endif  // USE_DOUBLE
#else
    setXAxis(x_axis[0], x_axis[1]);
#endif  // USE_INTRINSICS
  }

  /// Sets the y-axis column with the given values.
  inline void setYAxis(FLOAT y_x, FLOAT y_y) {
    m_value.setY(y_x);
    m_value.setW(y_y);
  }
  /// Sets the y-axis column with the given 2D array values.
  inline void setYAxis(const std::array<FLOAT, 2> &y_axis) {
    setYAxis(y_axis[0], y_axis[1]);
  }
  /// Sets the y-axis column with the given 2D vector values.
  inline void setYAxis(const Vector &y_axis) {
#ifdef USE_INTRINSICS
#ifdef USE_DOUBLE
#if defined(__AVX512F__) || defined(__AVX2__) || defined(__AVX__)
    __m256d tmp =
        _mm256_permute4x64_pd(y_axis.m_value, _MM_SHUFFLE(3, 1, 2, 0));
    m_value.m_value = _mm256_shuffle_pd(m_value.m_value, tmp, 0b0000);
#else
    m_value.m_value[0] = _mm_shuffle_pd(m_value.m_value[0], y_axis.m_value[0],
                                        _MM_SHUFFLE2(0, 0));
    m_value.m_value[1] = _mm_shuffle_pd(m_value.m_value[1], y_axis.m_value[0],
                                        _MM_SHUFFLE2(1, 0));
#endif  // AVX INTRINSICS
#else
    m_value.m_value = _mm_shuffle_ps(m_value.m_value, m_value.m_value,
                                     _MM_SHUFFLE(3, 1, 2, 0));
    m_value.m_value = _mm_movelh_ps(y_axis.m_value, m_value.m_value);
    m_value.m_value = _mm_shuffle_ps(m_value.m_value, m_value.m_value,
                                     _MM_SHUFFLE(1, 3, 0, 2));
#endif  // USE_DOUBLE
#else
    setYAxis(y_axis[0], y_axis[1]);
#endif  // USE_INTRINSICS
  }

  // Operators
  inline bool operator==(const Matrix2x2 &rhs) const {
    return m_value == rhs.m_value;
  }

  inline bool operator!=(const Matrix2x2 &rhs) const {
    return m_value != rhs.m_value;
  }

  /// Matrix-Matrix multiplication.
  friend inline Matrix2x2 operator*(const Matrix2x2 &lhs,
                                    const Matrix2x2 &rhs) {
    Matrix2x2 temp = Matrix2x2::transpose(rhs);
    return Matrix2x2{Vector{Vector::dot(lhs.getRow(0), temp.getRow(0)),
                            Vector::dot(lhs.getRow(0), temp.getRow(1)),
                            Vector::dot(lhs.getRow(1), temp.getRow(0)),
                            Vector::dot(lhs.getRow(1), temp.getRow(1))}};
  }

  /// Matrix-2D Column vector multiplication.
  friend inline Vector operator*(const Matrix2x2 &lhs, const Vector &rhs) {
    return Vector{Vector::dot(lhs.getRow(0), rhs),
                  Vector::dot(lhs.getRow(1), rhs), 0.0, 0.0};
  }

  // Type-Conversions
  /// Convert to a human-readable string value.
  explicit inline operator std::string() const {
#ifdef USE_INTRINSICS
    std::array<FLOAT, 4> val{};
    m_value.get(val.data());
    return "[" + std::to_string(val[0]) + ", " + std::to_string(val[1]) + "\n" +
           std::to_string(val[2]) + ", " + std::to_string(val[3]) + "]";
#else
    return "[" + std::to_string(m_value[0]) + ", " +
           std::to_string(m_value[1]) + "\n" + std::to_string(m_value[2]) +
           ", " + std::to_string(m_value[3]) + "]";
#endif  // USE_INTRINSICS
  }

  friend inline std::ostream &operator<<(std::ostream &out,
                                         const Matrix2x2 &mat) {
    return out << static_cast<std::string>(mat);
  }

  // Functions
  /// Returns the determinant of the 2x2 matrix.
  [[nodiscard]] inline FLOAT determinant() const {
#ifdef USE_INTRINSICS
#ifdef USE_DOUBLE
#if defined(__AVX512F__) || defined(__AVX2__) || defined(__AVX__)
    __m256d tmp =
        _mm256_permute4x64_pd(m_value.m_value, _MM_SHUFFLE(0, 1, 2, 3));
    tmp = _mm256_mul_pd(m_value.m_value, tmp);
    return _mm256_cvtsd_f64(_mm256_hsub_pd(tmp, tmp));
#else
    __m128d tmp = _mm_shuffle_pd(m_value.m_value[1], m_value.m_value[1],
                                 _MM_SHUFFLE2(0, 1));
    tmp = _mm_mul_pd(m_value.m_value[0], tmp);
    return _mm_cvtsd_f64(_mm_hsub_pd(tmp, tmp));
#endif  // AVX INTRINSICS
#else
    __m128 tmp = _mm_shuffle_ps(m_value.m_value, m_value.m_value,
                                _MM_SHUFFLE(0, 1, 2, 3));
    tmp = _mm_mul_ps(m_value.m_value, tmp);
    return _mm_cvtss_f32(_mm_hsub_ps(tmp, tmp));
#endif  // USE_DOUBLE
#else
    return m_value[0] * m_value[3] - m_value[1] * m_value[2];
#endif  // USE_INTRINSICS
  }

  /// Transposes the given 2x2 matrix.
  static inline Matrix2x2 transpose(const Matrix2x2 &matrix) {
#ifdef USE_INTRINSICS
#ifdef USE_DOUBLE
#if defined(__AVX512F__) || defined(__AVX2__) || defined(__AVX__)
    return Matrix2x2{Vector{_mm256_permute4x64_pd(matrix.m_value.m_value,
                                                  _MM_SHUFFLE(3, 1, 2, 0))}};
#else
    return Matrix2x2{Vector{
        {_mm_shuffle_pd(matrix.m_value.m_value[0], matrix.m_value.m_value[1],
                        _MM_SHUFFLE2(0, 0)),
         _mm_shuffle_pd(matrix.m_value.m_value[0], matrix.m_value.m_value[1],
                        _MM_SHUFFLE2(1, 1))}}};
#endif  // AVX INTRINSICS
#else
    return Matrix2x2{
        Vector{_mm_shuffle_ps(matrix.m_value.m_value, matrix.m_value.m_value,
                              _MM_SHUFFLE(3, 1, 2, 0))}};
#endif  // USE_DOUBLE
#else
    return Matrix2x2{matrix.m_value[0], matrix.m_value[1], matrix.m_value[2],
                     matrix.m_value[3]};
#endif  // USE_INTRINSICS
  }

  /// Returns the inverse of the given 2x2 matrix, if it exists.
  /// Otherwise returns the original matrix as is.
  static inline Matrix2x2 inverse(const Matrix2x2 &matrix) {
    FLOAT determinant = matrix.determinant();
    if (approxEqual(determinant, 0.0)) return matrix;
    determinant = static_cast<FLOAT>(1.0) / determinant;
    Vector adj = matrix.m_value;
#ifdef USE_INTRINSICS
#ifdef USE_DOUBLE
#if defined(__AVX512F__) || defined(__AVX2__) || defined(__AVX__)
    __m256d neg_mask = _mm256_set_pd(0.0, -0.0, -0.0, 0.0);
    adj.m_value = _mm256_permute4x64_pd(adj.m_value, _MM_SHUFFLE(0, 2, 1, 3));
    adj.m_value = _mm256_xor_pd(adj.m_value, neg_mask);
#else
    __m128d neg_mask = _mm_set_pd(0.0, -0.0);
    adj.m_value[0] =
        _mm_shuffle_pd(adj.m_value[1], adj.m_value[0], _MM_SHUFFLE2(1, 1));
    adj.m_value[1] =
        _mm_shuffle_pd(adj.m_value[0], adj.m_value[1], _MM_SHUFFLE2(0, 0));
    adj.m_value[0] = _mm_xor_pd(adj.m_value[0], neg_mask);
    adj.m_value[1] = _mm_xor_pd(adj.m_value[1], neg_mask);
    adj.m_value[1] =
        _mm_shuffle_pd(adj.m_value[1], adj.m_value[1], _MM_SHUFFLE2(0, 1));
#endif  // AVX INTRINSICS
#else
    __m128 neg_mask = _mm_set_ps(0.0f, -0.0f, -0.0f, 0.0f);
    adj.m_value =
        _mm_shuffle_ps(adj.m_value, adj.m_value, _MM_SHUFFLE(0, 2, 1, 3));
    adj.m_value = _mm_xor_ps(adj.m_value, neg_mask);
#endif  // USE_DOUBLE
#else
    FLOAT tmp = adj[0];
    adj.setX(adj[3]);
    adj.setY(-adj[1]);
    adj.setZ(-adj[2]);
    adj.setW(tmp);
#endif  // USE_INTRINSICS
    return Matrix2x2(determinant * adj);
  }

 private:
  Vector m_value;

  /// Construct a 2x2 matrix a 4D vector which has matrix values in the
  /// memory layout order. i.e [x_x, y_x, x_y, y_y].
  explicit Matrix2x2(const Vector &value)
  : m_value {value} {};
};
}  // namespace cgmath::internal

#endif  // CGMATH_MATRIX2X2_HPP
