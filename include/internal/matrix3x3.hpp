// Copyright 2023 Suraj Shettigar
// SPDX-License-Identifier: Apache-2.0

#ifndef CGMATH_INTERNAL_MATRIX3X3_HPP
#define CGMATH_INTERNAL_MATRIX3X3_HPP

#include "vector.hpp"

namespace cgmath::internal {

/// 3x3 Matrix using Column Major Mathematical Notation and stored in a Row
/// Major memory layout.
class Matrix3x3 {
 public:
  // Constructors / Destructors
  /// Constructs an identity 3x3 matrix.
  explicit Matrix3x3()
      : m_value{Vector{1.0, 0.0, 0.0}, Vector{0.0, 1.0, 0.0},
                Vector{0.0, 0.0, 1.0}} {};
  /// Constructs a 3x2 matrix with x-axis, y-axis and z-axis column values.
  explicit Matrix3x3(FLOAT x_x, FLOAT x_y, FLOAT x_z, FLOAT y_x, FLOAT y_y,
                     FLOAT y_z, FLOAT z_x, FLOAT z_y, FLOAT z_z)
      : m_value{Vector{x_x, y_x, z_x}, Vector{x_y, y_y, z_y},
                Vector{x_z, y_z, z_z}} {};
  /// Constructs a 3x3 matrix with 3D x-axis, y-axis and z-axis array values.
  explicit Matrix3x3(const std::array<FLOAT, 3> &x_axis,
                     const std::array<FLOAT, 3> &y_axis,
                     const std::array<FLOAT, 3> &z_axis)
      : m_value{Vector{x_axis[0], y_axis[0], z_axis[0]},
                Vector{x_axis[1], y_axis[1], z_axis[1]},
                Vector{x_axis[2], y_axis[2], z_axis[2]}} {};
  /// Constructs a 3x3 matrix with 3D x-axis, y-axis and z-axis vector values.
  explicit Matrix3x3(const Vector &x_axis, const Vector &y_axis,
                     const Vector &z_axis)
      : m_value{x_axis, y_axis, z_axis} {
    *this = Matrix3x3::transpose(*this);
  };

  ~Matrix3x3() = default;

  // Copy
  Matrix3x3(const Matrix3x3 &matrix) {
    if (this != &matrix) m_value = matrix.m_value;
  }
  Matrix3x3 &operator=(const Matrix3x3 &matrix) {
    if (this != &matrix) m_value = matrix.m_value;
    return *this;
  }

  // Getters / Setters
  /**
   * Returns the matrix row.
   * @param index Row index.
   * @return 2D vector containing given row values. (Eg: [x_x, y_x, z_x, 0])
   */
  [[nodiscard]] inline Vector operator[](size_t index) const {
    return m_value[index];
  }

  /**
   * Returns the matrix row.
   * @param index Row index.
   * @return 3D vector containing given row values. (Eg: [x_x, y_x, z_x, 0])
   */
  [[nodiscard]] inline Vector getRow(size_t index) const {
    return (*this)[index];
  }

  /**
   * Returns the matrix column.
   * @param index Column index.
   * @return 3D vector containing given column values. (Eg: [x_x, x_y, x_z, 0])
   */
  [[nodiscard]] inline Vector getColumn(size_t index) const {
    Matrix3x3 val = Matrix3x3::transpose(*this);
    return val.m_value[index];
  }

  /// Returns x-axis (column 0) of the matrix.
  [[nodiscard]] inline Vector getXAxis() const { return getColumn(0); }
  /// Returns y-axis (column 1) of the matrix.
  [[nodiscard]] inline Vector getYAxis() const { return getColumn(1); }
  /// Returns y-axis (column 2) of the matrix.
  [[nodiscard]] inline Vector getZAxis() const { return getColumn(2); }

  /// Sets the x-axis column with the given values.
  inline void setXAxis(FLOAT x_x, FLOAT x_y, FLOAT x_z) {
    m_value[0].setX(x_x);
    m_value[1].setX(x_y);
    m_value[2].setX(x_z);
  }
  /// Sets the x-axis column with the given 3D array values.
  inline void setXAxis(const std::array<FLOAT, 3> &x_axis) {
    setXAxis(x_axis[0], x_axis[1], x_axis[2]);
  }
  /// Sets the x-axis column with the given 3D vector values.
  inline void setXAxis(const Vector &x_axis) {
    std::array<FLOAT, 4> val{};
    x_axis.get(val.data());
    setXAxis(val[0], val[1], val[2]);
  }

  /// Sets the y-axis column with the given values.
  inline void setYAxis(FLOAT y_x, FLOAT y_y, FLOAT y_z) {
    m_value[0].setY(y_x);
    m_value[1].setY(y_y);
    m_value[2].setY(y_z);
  }
  /// Sets the y-axis column with the given 3D array values.
  inline void setYAxis(const std::array<FLOAT, 3> &y_axis) {
    setYAxis(y_axis[0], y_axis[1], y_axis[2]);
  }
  /// Sets the y-axis column with the given 3D vector values.
  inline void setYAxis(const Vector &y_axis) {
    std::array<FLOAT, 4> val{};
    y_axis.get(val.data());
    setYAxis(val[0], val[1], val[2]);
  }

  /// Sets the z-axis column with the given values.
  inline void setZAxis(FLOAT z_x, FLOAT z_y, FLOAT z_z) {
    m_value[0].setZ(z_x);
    m_value[1].setZ(z_y);
    m_value[2].setZ(z_z);
  }
  /// Sets the z-axis column with the given 3D array values.
  inline void setZAxis(const std::array<FLOAT, 3> &z_axis) {
    setZAxis(z_axis[0], z_axis[1], z_axis[2]);
  }
  /// Sets the z-axis column with the given 3D vector values.
  inline void setZAxis(const Vector &z_axis) {
    std::array<FLOAT, 4> val{};
    z_axis.get(val.data());
    setZAxis(val[0], val[1], val[2]);
  }

  // Operators
  inline bool operator==(const Matrix3x3 &rhs) const {
    return m_value[0] == rhs.m_value[0] && m_value[1] == rhs.m_value[1] &&
           m_value[2] == rhs.m_value[2];
  }

  inline bool operator!=(const Matrix3x3 &rhs) const {
    return m_value[0] != rhs.m_value[0] || m_value[1] != rhs.m_value[1] ||
           m_value[2] != rhs.m_value[2];
  }

  /// Matrix-Matrix multiplication.
  friend inline Matrix3x3 operator*(const Matrix3x3 &lhs,
                                    const Matrix3x3 &rhs) {
    Matrix3x3 temp = Matrix3x3::transpose(rhs);
    return Matrix3x3{
        Vector::dot(lhs.getRow(0), temp.getRow(0)),
        Vector::dot(lhs.getRow(1), temp.getRow(0)),
        Vector::dot(lhs.getRow(2), temp.getRow(0)),
        Vector::dot(lhs.getRow(0), temp.getRow(1)),
        Vector::dot(lhs.getRow(1), temp.getRow(1)),
        Vector::dot(lhs.getRow(2), temp.getRow(1)),
        Vector::dot(lhs.getRow(0), temp.getRow(2)),
        Vector::dot(lhs.getRow(1), temp.getRow(2)),
        Vector::dot(lhs.getRow(2), temp.getRow(2)),
    };
  }

  /// Matrix-3D Column vector multiplication.
  friend inline Vector operator*(const Matrix3x3 &lhs, const Vector &rhs) {
    return Vector{Vector::dot(lhs.getRow(0), rhs),
                  Vector::dot(lhs.getRow(1), rhs),
                  Vector::dot(lhs.getRow(2), rhs), 0.0};
  }

  // Type-Conversions
  /// Convert to a human-readable string value.
  explicit inline operator std::string() const {
#ifdef USE_INTRINSICS
    std::array<FLOAT, 4> x{}, y{}, z{};
    m_value[0].get(x.data());
    m_value[1].get(y.data());
    m_value[2].get(z.data());
    return "[" + std::to_string(x[0]) + ", " + std::to_string(x[1]) + ", " +
           std::to_string(x[2]) + "\n" + std::to_string(y[0]) + ", " +
           std::to_string(y[1]) + ", " + std::to_string(y[2]) + "\n" +
           std::to_string(z[0]) + ", " + std::to_string(z[1]) + ", " +
           std::to_string(z[2]) + "]";
#else
    return "[" + std::to_string(m_value[0][0]) + ", " +
           std::to_string(m_value[0][1]) + ", " +
           std::to_string(m_value[0][2]) + "\n" +
           std::to_string(m_value[1][0]) + ", " +
           std::to_string(m_value[1][1]) + ", " +
           std::to_string(m_value[1][2]) + "\n" +
           std::to_string(m_value[2][0]) + ", " +
           std::to_string(m_value[2][1]) + ", " +
           std::to_string(m_value[2][2]) + "]";
#endif  // USE_INTRINSICS
  }

  friend inline std::ostream &operator<<(std::ostream &out,
                                         const Matrix3x3 &mat) {
    return out << static_cast<std::string>(mat);
  }

  // Functions
  /// Returns the determinant of the 3x3 matrix.
  [[nodiscard]] inline FLOAT determinant() const {
// | a    b   c |
// | d    e   f |
// | g    h   i |
// Determinant = aei + bfg + cdh - afh - bdi - ceg
#ifdef USE_INTRINSICS
#ifdef USE_DOUBLE
#if defined(__AVX512F__) || defined(__AVX2__) || defined(__AVX__)
    // [a, b, c] * [e, f, d] * [i, g, h]
    __m256d tmp1 = m_value[0].m_value;
    __m256d tmp2 =
        _mm256_permute4x64_pd(m_value[1].m_value, _MM_SHUFFLE(3, 0, 2, 1));
    __m256d tmp3 =
        _mm256_permute4x64_pd(m_value[2].m_value, _MM_SHUFFLE(3, 1, 0, 2));
    __m256d tmp4 = _mm256_mul_pd(tmp1, tmp2);
    tmp4 = _mm256_mul_pd(tmp3, tmp4);

    // [a, b, c] * [f, d, e] * [h, i, g]
    tmp2 = _mm256_permute4x64_pd(m_value[1].m_value, _MM_SHUFFLE(3, 1, 0, 2));
    tmp3 = _mm256_permute4x64_pd(m_value[2].m_value, _MM_SHUFFLE(3, 0, 2, 1));
    __m256d tmp5 = _mm256_mul_pd(tmp1, tmp2);
    tmp5 = _mm256_mul_pd(tmp3, tmp5);

    // [aei, bfg, cdh] - [afh, bdi, ceg] => [x, y, z]
    tmp5 = _mm256_sub_pd(tmp4, tmp5);
    // [x + y, z + 0, x + y, z + 0] => [xy, z0, xy, z0]
    tmp5 = _mm256_hadd_pd(tmp5, tmp5);
    // [xy + z0, xy + z0, xy + z0, xy + z0] => [xyz0, xyz0, xyz0, xyz0]
    tmp5 = _mm256_hadd_pd(tmp5, tmp5);
    return _mm256_cvtsd_f64(tmp5);
#else
    // TODO: Matrix3x3 SSE Double determinant.
#endif  // AVX INTRINSICS
#else
    // [a, b, c] * [e, f, d] * [i, g, h]
    __m128 tmp1 = m_value[0].m_value;
    __m128 tmp2 = _mm_shuffle_ps(m_value[1].m_value, m_value[1].m_value,
                                 _MM_SHUFFLE(3, 0, 2, 1));
    __m128 tmp3 = _mm_shuffle_ps(m_value[2].m_value, m_value[2].m_value,
                                 _MM_SHUFFLE(3, 1, 0, 2));
    __m128 tmp4 = _mm_mul_ps(tmp1, tmp2);
    tmp4 = _mm_mul_ps(tmp3, tmp4);

    // [a, b, c] * [f, d, e] * [h, i, g]
    tmp2 = _mm_shuffle_ps(m_value[1].m_value, m_value[1].m_value,
                          _MM_SHUFFLE(3, 1, 0, 2));
    tmp3 = _mm_shuffle_ps(m_value[2].m_value, m_value[2].m_value,
                          _MM_SHUFFLE(3, 0, 2, 1));
    __m128 tmp5 = _mm_mul_ps(tmp1, tmp2);
    tmp5 = _mm_mul_ps(tmp3, tmp5);

    // [aei, bfg, cdh] - [afh, bdi, ceg] => [x, y, z]
    tmp5 = _mm_sub_ps(tmp4, tmp5);
    // [x + y, z + 0, x + y, z + 0] => [xy, z0, xy, z0]
    tmp5 = _mm_hadd_ps(tmp5, tmp5);
    // [xy + z0, xy + z0, xy + z0, xy + z0] => [xyz0, xyz0, xyz0, xyz0]
    tmp5 = _mm_hadd_ps(tmp5, tmp5);
    return _mm_cvtss_f32(tmp5);
#endif  // USE_DOUBLE
#else
    return m_value[0][0] * m_value[1][1] * m_value[2][2] +
           m_value[0][1] * m_value[1][2] * m_value[2][0] +
           m_value[0][2] * m_value[1][0] * m_value[2][1] -
           m_value[0][0] * m_value[1][2] * m_value[2][1] -
           m_value[0][1] * m_value[1][0] * m_value[2][2] -
           m_value[0][2] * m_value[1][1] * m_value[2][0];
#endif  // USE_INTRINSICS
  }

  /// Transposes the given 3x3 matrix.
  static inline Matrix3x3 transpose(const Matrix3x3 &matrix) {
    // TODO: Matrix3x3 transpose.
    return matrix;
  }

  /// Returns the inverse of the given 3x3 matrix, if it exists.
  /// Otherwise returns the original matrix as is.
  static inline Matrix3x3 inverse(const Matrix3x3 &matrix) {
    // TODO: Matrix3x3 inverse.
    return matrix;
  }

 private:
  std::array<Vector, 3> m_value{};
};
}  // namespace cgmath::internal

#endif  // CGMATH_INTERNAL_MATRIX3X3_HPP
