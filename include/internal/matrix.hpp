// Copyright 2023 Suraj Shettigar
// SPDX-License-Identifier: Apache-2.0

#ifndef CGMATH_INTERNAL_MATRIX_HPP
#define CGMATH_INTERNAL_MATRIX_HPP

#include "matrix2x2.hpp"

namespace cgmath::internal {

/// 4x4 Matrix using Column Major Mathematical Notation and stored in a Row
/// Major memory layout.
class Matrix {
 public:
  // Constructors / Destructors
  explicit Matrix()
      : m_value{Vector{1.0, 0.0, 0.0, 0.0}, Vector{0.0, 1.0, 0.0, 0.0},
                Vector{0.0, 0.0, 1.0, 0.0}, Vector{0.0, 0.0, 0.0, 1.0}} {};

  explicit Matrix(FLOAT x_x, FLOAT x_y, FLOAT x_z, FLOAT x_w, FLOAT y_x,
                  FLOAT y_y, FLOAT y_z, FLOAT y_w, FLOAT z_x, FLOAT z_y,
                  FLOAT z_z, FLOAT z_w, FLOAT t_x, FLOAT t_y, FLOAT t_z,
                  FLOAT t_w)
      : m_value{Vector{x_x, y_x, z_x, t_x}, Vector{x_y, y_y, z_y, t_y},
                Vector{x_z, y_z, z_z, t_z}, Vector{x_w, y_w, z_w, t_w}} {};

  explicit Matrix(const std::array<FLOAT, 4> &x_axis,
                  const std::array<FLOAT, 4> &y_axis,
                  const std::array<FLOAT, 4> &z_axis,
                  const std::array<FLOAT, 4> &translation)
      : m_value{Vector{x_axis[0], y_axis[0], z_axis[0], translation[0]},
                Vector{x_axis[1], y_axis[1], z_axis[1], translation[1]},
                Vector{x_axis[2], y_axis[2], z_axis[2], translation[2]},
                Vector{x_axis[3], y_axis[3], z_axis[3], translation[3]}} {};

  explicit Matrix(const Vector &x_axis, const Vector &y_axis,
                  const Vector &z_axis, const Vector &translation)
      : m_value{x_axis, y_axis, z_axis, translation} {};

  ~Matrix() = default;

  // Copy
  Matrix(const Matrix &matrix) {
    if (this != &matrix) m_value = matrix.m_value;
  }
  Matrix &operator=(const Matrix &matrix) {
    if (this != &matrix) m_value = matrix.m_value;
    return *this;
  }

  // Getters / Setters
  [[nodiscard]] inline Vector operator[](size_t index) const {
    return m_value[index];
  }
  [[nodiscard]] inline Vector getRow(size_t index) const {
    return m_value[index];
  }
  [[nodiscard]] inline Vector getColumn(size_t index) const {
#ifdef USE_INTRINSICS
    Matrix temp = Matrix::transpose(*this);
    return temp[index];
#else
    return Vector{m_value[0][index], m_value[1][index], m_value[2][index],
                  m_value[3][index]};
#endif
  }
  [[nodiscard]] inline Vector getXAxis() const { return getColumn(0); }
  [[nodiscard]] inline Vector getYAxis() const { return getColumn(1); }
  [[nodiscard]] inline Vector getZAxis() const { return getColumn(2); }
  [[nodiscard]] inline Vector getTranslation() const { return getColumn(3); }

  inline void setXAxis(FLOAT x_x, FLOAT x_y, FLOAT x_z, FLOAT x_w) {
    m_value[0].setX(x_x);
    m_value[1].setX(x_y);
    m_value[2].setX(x_z);
    m_value[3].setX(x_w);
  }
  inline void setXAxis(const std::array<FLOAT, 4> &x_axis) {
    m_value[0].setX(x_axis[0]);
    m_value[1].setX(x_axis[1]);
    m_value[2].setX(x_axis[2]);
    m_value[3].setX(x_axis[3]);
  }
  inline void setXAxis(const Vector &x_axis) {
    m_value[0].setX(x_axis[0]);
    m_value[1].setX(x_axis[1]);
    m_value[2].setX(x_axis[2]);
    m_value[3].setX(x_axis[3]);
  }
  inline void setYAxis(FLOAT y_x, FLOAT y_y, FLOAT y_z, FLOAT y_w) {
    m_value[0].setY(y_x);
    m_value[1].setY(y_y);
    m_value[2].setY(y_z);
    m_value[3].setY(y_w);
  }
  inline void setYAxis(const std::array<FLOAT, 4> &y_axis) {
    m_value[0].setY(y_axis[0]);
    m_value[1].setY(y_axis[1]);
    m_value[2].setY(y_axis[2]);
    m_value[3].setY(y_axis[3]);
  }
  inline void setYAxis(const Vector &y_axis) {
    m_value[0].setY(y_axis[0]);
    m_value[1].setY(y_axis[1]);
    m_value[2].setY(y_axis[2]);
    m_value[3].setY(y_axis[3]);
  }
  inline void setZAxis(FLOAT z_x, FLOAT z_y, FLOAT z_z, FLOAT z_w) {
    m_value[0].setZ(z_x);
    m_value[1].setZ(z_y);
    m_value[2].setZ(z_z);
    m_value[3].setZ(z_w);
  }
  inline void setZAxis(const std::array<FLOAT, 4> &z_axis) {
    m_value[0].setZ(z_axis[0]);
    m_value[1].setZ(z_axis[1]);
    m_value[2].setZ(z_axis[2]);
    m_value[3].setZ(z_axis[3]);
  }
  inline void setZAxis(const Vector &z_axis) {
    m_value[0].setZ(z_axis[0]);
    m_value[1].setZ(z_axis[1]);
    m_value[2].setZ(z_axis[2]);
    m_value[3].setZ(z_axis[3]);
  }
  inline void setTranslation(FLOAT t_x, FLOAT t_y, FLOAT t_z, FLOAT t_w) {
    m_value[0].setW(t_x);
    m_value[1].setW(t_y);
    m_value[2].setW(t_z);
    m_value[3].setW(t_w);
  }
  inline void setTranslation(const std::array<FLOAT, 4> &translation) {
    m_value[0].setW(translation[0]);
    m_value[1].setW(translation[1]);
    m_value[2].setW(translation[2]);
    m_value[3].setW(translation[3]);
  }
  inline void setTranslation(const Vector &translation) {
    m_value[0].setW(translation[0]);
    m_value[1].setW(translation[1]);
    m_value[2].setW(translation[2]);
    m_value[3].setW(translation[3]);
  }

  // Operators
  inline bool operator==(const Matrix &rhs) const {
    return m_value[0] == rhs.m_value[0] && m_value[1] == rhs.m_value[1] &&
           m_value[2] == rhs.m_value[2] && m_value[3] == rhs.m_value[3];
  }

  inline bool operator!=(const Matrix &rhs) const {
    return m_value[0] != rhs.m_value[0] || m_value[1] != rhs.m_value[1] ||
           m_value[2] != rhs.m_value[2] || m_value[3] != rhs.m_value[3];
  }

  friend inline Matrix operator*(const Matrix &lhs, const Matrix &rhs) {
    Matrix temp = Matrix::transpose(rhs);
    return Matrix{Vector{Vector::dot(lhs.getRow(0), temp.getRow(0)),
                         Vector::dot(lhs.getRow(0), temp.getRow(1)),
                         Vector::dot(lhs.getRow(0), temp.getRow(2)),
                         Vector::dot(lhs.getRow(0), temp.getRow(3))},
                  Vector{Vector::dot(lhs.getRow(1), temp.getRow(0)),
                         Vector::dot(lhs.getRow(1), temp.getRow(1)),
                         Vector::dot(lhs.getRow(1), temp.getRow(2)),
                         Vector::dot(lhs.getRow(1), temp.getRow(3))},
                  Vector{Vector::dot(lhs.getRow(2), temp.getRow(0)),
                         Vector::dot(lhs.getRow(2), temp.getRow(1)),
                         Vector::dot(lhs.getRow(2), temp.getRow(2)),
                         Vector::dot(lhs.getRow(2), temp.getRow(3))},
                  Vector{Vector::dot(lhs.getRow(3), temp.getRow(0)),
                         Vector::dot(lhs.getRow(3), temp.getRow(1)),
                         Vector::dot(lhs.getRow(3), temp.getRow(2)),
                         Vector::dot(lhs.getRow(3), temp.getRow(3))}};
  }

  friend inline Vector operator*(const Matrix &lhs, const Vector &rhs) {
    return Vector{
        Vector::dot(lhs.getRow(0), rhs), Vector::dot(lhs.getRow(1), rhs),
        Vector::dot(lhs.getRow(2), rhs), Vector::dot(lhs.getRow(3), rhs)};
  }

  // Functions
  inline FLOAT determinant() const {
    FLOAT determinant = 0.0;
    return determinant;
  }
  static inline Matrix transpose(const Matrix &matrix) {
#ifdef USE_INTRINSICS
#ifdef USE_DOUBLE
#if defined(__AVX512F__)
    // TODO: Transpose AVX512 matrix.
    return matrix;
#elif defined(__AVX2__) || defined(__AVX__)
    // x_x, y_x, x_y, y_y
    __m256d row1 = _mm256_permute2f128_pd(matrix.m_value[0].m_value,
                                          matrix.m_value[1].m_value,
                                          _MM_SHUFFLE(1, 0, 1, 0));
    // x_x, x_y, y_x, y_y
    row1 = _mm256_permute4x64_pd(row1, _MM_SHUFFLE(3, 1, 2, 0));

    // z_x, t_x, z_y, t_y
    __m256d row2 = _mm256_permute2f128_pd(matrix.m_value[0].m_value,
                                          matrix.m_value[1].m_value,
                                          _MM_SHUFFLE(3, 2, 3, 2));
    // z_x, z_y, t_x, t_y
    row2 = _mm256_permute4x64_pd(row2, _MM_SHUFFLE(3, 1, 2, 0));

    // x_z, y_z, x_w, y_w
    __m256d row3 = _mm256_permute2f128_pd(matrix.m_value[2].m_value,
                                          matrix.m_value[3].m_value,
                                          _MM_SHUFFLE(1, 0, 1, 0));
    // x_z, x_w, y_z, y_w
    row3 = _mm256_permute4x64_pd(row3, _MM_SHUFFLE(3, 1, 2, 0));

    // z_z, t_z, z_w, t_w
    __m256d row4 = _mm256_permute2f128_pd(matrix.m_value[2].m_value,
                                          matrix.m_value[3].m_value,
                                          _MM_SHUFFLE(3, 2, 3, 2));
    // z_z, z_w, t_z, t_w
    row4 = _mm256_permute4x64_pd(row4, _MM_SHUFFLE(3, 1, 2, 0));

    return Matrix{
        // x_x, x_y, x_z, x_w
        Vector{_mm256_permute2f128_pd(row1, row3, _MM_SHUFFLE(1, 0, 1, 0))},
        // y_x, y_y, y_z, y_w
        Vector{_mm256_permute2f128_pd(row1, row3, _MM_SHUFFLE(3, 2, 3, 2))},
        // z_x, z_y, z_z, z_w
        Vector{_mm256_permute2f128_pd(row2, row4, _MM_SHUFFLE(1, 0, 1, 0))},
        // t_x, t_y, t_z, t_w
        Vector{_mm256_permute2f128_pd(row2, row4, _MM_SHUFFLE(3, 2, 3, 2))}};
#else
    // x_x, x_y
    __m128d row11 = _mm_shuffle_pd(matrix.m_value[0].m_value[0],
                                   matrix.m_value[1].m_value[0], 0b00);
    // x_z, x_w
    __m128d row12 = _mm_shuffle_pd(matrix.m_value[2].m_value[0],
                                   matrix.m_value[3].m_value[0], 0b00);
    // y_x, y_y
    __m128d row21 = _mm_shuffle_pd(matrix.m_value[0].m_value[0],
                                   matrix.m_value[1].m_value[0], 0b11);
    // y_z, y_w
    __m128d row22 = _mm_shuffle_pd(matrix.m_value[2].m_value[0],
                                   matrix.m_value[3].m_value[0], 0b11);
    // z_x, z_y
    __m128d row31 = _mm_shuffle_pd(matrix.m_value[0].m_value[1],
                                   matrix.m_value[1].m_value[1], 0b00);
    // z_z, z_w
    __m128d row32 = _mm_shuffle_pd(matrix.m_value[2].m_value[1],
                                   matrix.m_value[3].m_value[1], 0b00);
    // t_x, t_y
    __m128d row41 = _mm_shuffle_pd(matrix.m_value[0].m_value[1],
                                   matrix.m_value[1].m_value[1], 0b11);
    // t_z, t_w
    __m128d row42 = _mm_shuffle_pd(matrix.m_value[2].m_value[1],
                                   matrix.m_value[3].m_value[1], 0b11);

    return Matrix{Vector{{row11, row12}}, Vector{{row21, row22}},
                  Vector{{row31, row32}}, Vector{{row41, row42}}};
#endif  // AVX INTRINSICS
#else
    // x_x, y_x, x_y, y_y
    __m128 row1 =
        _mm_shuffle_ps(matrix.m_value[0].m_value, matrix.m_value[1].m_value,
                       _MM_SHUFFLE(1, 0, 1, 0));
    // z_x, t_x, z_y, t_y
    __m128 row2 =
        _mm_shuffle_ps(matrix.m_value[0].m_value, matrix.m_value[1].m_value,
                       _MM_SHUFFLE(3, 2, 3, 2));
    // x_z, y_z, x_w, y_w
    __m128 row3 =
        _mm_shuffle_ps(matrix.m_value[2].m_value, matrix.m_value[3].m_value,
                       _MM_SHUFFLE(1, 0, 1, 0));
    // z_z, t_z, z_w, t_w
    __m128 row4 =
        _mm_shuffle_ps(matrix.m_value[2].m_value, matrix.m_value[3].m_value,
                       _MM_SHUFFLE(3, 2, 3, 2));

    return Matrix{// x_x, x_y, x_z, x_w
                  Vector{_mm_shuffle_ps(row1, row3, _MM_SHUFFLE(2, 0, 2, 0))},
                  // y_x, y_y, y_z, y_w
                  Vector{_mm_shuffle_ps(row1, row3, _MM_SHUFFLE(3, 1, 3, 1))},
                  // z_x, z_y, z_z, z_w
                  Vector{_mm_shuffle_ps(row2, row4, _MM_SHUFFLE(2, 0, 2, 0))},
                  // t_x, t_y, t_z, t_w
                  Vector{_mm_shuffle_ps(row2, row4, _MM_SHUFFLE(3, 1, 3, 1))}};
#endif  // USE_DOUBLE
#else
    return Matrix{matrix.getRow(0), matrix.getRow(1), matrix.getRow(2),
                  matrix.getRow(3)};
#endif  // USE_INTRINSICS
  }

  static inline Matrix inverse(const Matrix &matrix) {
    FLOAT determinant = matrix.determinant();
    if (approxEqual(determinant, 0.0)) return matrix;
  }

 private:
#if defined(USE_INTRINSICS) && defined(__AVX512F__)
#ifdef USE_DOUBLE
  std : array<__m512d, 2> m_value;
#else
  __m512 m_value;
#endif
#else
  std::array<Vector, 4> m_value{};
#endif
};

}  // namespace cgmath::internal

#endif  // CGMATH_INTERNAL_MATRIX_HPP
