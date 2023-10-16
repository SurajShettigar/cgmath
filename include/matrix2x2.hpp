// Copyright 2023 Suraj Shettigar
// SPDX-License-Identifier: Apache-2.0

#ifndef CGMATH_MATRIX2X2_HPP
#define CGMATH_MATRIX2X2_HPP

#include "internal/matrix2x2.hpp"
#include "vector2.hpp"

namespace cgmath {

class Matrix2x2 {
 public:
  // Constructors / Destructors
  /// Constructs an identity 2x2 matrix.
  explicit Matrix2x2() = default;
  /// Constructs a 2x2 matrix with x-axis and y-axis column values.
  explicit Matrix2x2(FLOAT x_x, FLOAT x_y, FLOAT y_x, FLOAT y_y)
      : m_value{x_x, x_y, y_x, y_y} {};
  /// Constructs a 2x2 matrix with 2D x-axis and y-axis array values.
  explicit Matrix2x2(const std::array<FLOAT, 2> &x_axis,
                     const std::array<FLOAT, 2> &y_axis)
      : m_value{x_axis, y_axis} {};
  /// Constructs a 2x2 matrix with 2D x-axis and y-axis vector values.
  explicit Matrix2x2(const Vector2 &x_axis, const Vector2 &y_axis)
      : m_value{x_axis.m_value, y_axis.m_value} {};
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
   * Returns the matrix row. [x_x, y_x] or [x_y, y_y]
   * @param index Row index.
   * @return 2D vector containing given row values.
   */
  [[nodiscard]] inline Vector2 operator[](size_t index) const {
    return Vector2{m_value[index]};
  }

  /**
   * Returns the matrix row. [x_x, y_x] or [x_y, y_y]
   * @param index Row index.
   * @return 2D vector containing given row values.
   */
  [[nodiscard]] inline Vector2 getRow(size_t index) const {
    return Vector2{m_value.getRow(index)};
  }

  /**
   * Returns the matrix column. [x_x, x_y] or [y_x, y_y]
   * @param index Column index.
   * @return 2D vector containing given column values.
   */
  [[nodiscard]] inline Vector2 getColumn(size_t index) const {
    return Vector2{m_value.getColumn(index)};
  }

  /// Returns x-axis (column 0) of the matrix.
  [[nodiscard]] inline Vector2 getXAxis() const {
    return Vector2{m_value.getXAxis()};
  }
  /// Returns y-axis (column 1) of the matrix.
  [[nodiscard]] inline Vector2 getYAxis() const {
    return Vector2{m_value.getYAxis()};
  }

  /// Sets the x-axis column with the given values.
  inline void setXAxis(FLOAT x_x, FLOAT x_y) { m_value.setXAxis(x_x, x_y); }
  /// Sets the x-axis column with the given 2D array values.
  inline void setXAxis(const std::array<FLOAT, 2> &x_axis) {
    m_value.setXAxis(x_axis);
  }
  /// Sets the x-axis column with the given 2D vector values.
  inline void setXAxis(const Vector2 &x_axis) {
    m_value.setXAxis(x_axis.m_value);
  }

  /// Sets the y-axis column with the given values.
  inline void setYAxis(FLOAT y_x, FLOAT y_y) { m_value.setYAxis(y_x, y_y); }
  /// Sets the y-axis column with the given 2D array values.
  inline void setYAxis(const std::array<FLOAT, 2> &y_axis) {
    m_value.setYAxis(y_axis);
  }
  /// Sets the y-axis column with the given 2D vector values.
  inline void setYAxis(const Vector2 &y_axis) {
    m_value.setYAxis(y_axis.m_value);
  }

  // Operators
  inline bool operator==(const Matrix2x2 &rhs) const {
    return m_value == rhs.m_value;
  }

  inline bool operator!=(const Matrix2x2 &rhs) const {
    return m_value != rhs.m_value;
  }

  /// Matrix-Matrix multiplication.
  inline Matrix2x2 operator*(const Matrix2x2 &rhs) const {
    return Matrix2x2{m_value * rhs.m_value};
  }

  /// Matrix-2D Column vector multiplication.
  inline Vector2 operator*(const Vector2 &vec) const {
    return Vector2{m_value * vec.m_value};
  }

  // TODO: Implement matrix2x2-scalar multiplication.

  // Type-Conversions
  /// Convert to a human-readable string value.
  explicit inline operator std::string() const {
    return static_cast<std::string>(m_value);
  }

  friend inline std::ostream &operator<<(std::ostream &out,
                                         const Matrix2x2 &mat) {
    return out << static_cast<std::string>(mat);
  }

  // Functions
  /// Returns the determinant of the 2x2 matrix.
  [[nodiscard]] inline FLOAT determinant() const {
    return m_value.determinant();
  }

  /// Transposes the given 2x2 matrix.
  static inline Matrix2x2 transpose(const Matrix2x2 &matrix) {
    return Matrix2x2{internal::Matrix2x2::transpose(matrix.m_value)};
  }

  /// Returns the inverse of the given 2x2 matrix, if it exists.
  /// Otherwise returns the original matrix as is.
  static inline Matrix2x2 inverse(const Matrix2x2 &matrix) {
    return Matrix2x2{internal::Matrix2x2::inverse(matrix.m_value)};
  }

 private:
  internal::Matrix2x2 m_value{};

  explicit Matrix2x2(const internal::Matrix2x2 &value) : m_value{value} {};
};
}  // namespace cgmath
#endif  // CGMATH_MATRIX2X2_HPP
