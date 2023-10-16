// Copyright 2023 Suraj Shettigar
// SPDX-License-Identifier: Apache-2.0

#ifndef CGMATH_MATRIX3X3_HPP
#define CGMATH_MATRIX3X3_HPP

#include "internal/matrix3x3.hpp"
#include "vector3.hpp"

namespace cgmath {

class Matrix3x3 {
 public:
  // Constructors / Destructors
  /// Constructs an identity 3x3 matrix.
  explicit Matrix3x3() : m_value{} {};
  /// Constructs a 3x3 matrix with x-axis, y-axis and z-axis column values.
  explicit Matrix3x3(FLOAT x_x, FLOAT x_y, FLOAT x_z, FLOAT y_x, FLOAT y_y,
                     FLOAT y_z, FLOAT z_x, FLOAT z_y, FLOAT z_z)
      : m_value{x_x, x_y, x_z, y_x, y_y, y_z, z_x, z_y, z_z} {};
  /// Constructs a 3x3 matrix with 3D x-axis, y-axis and z-axis array values.
  explicit Matrix3x3(const std::array<FLOAT, 3> &x_axis,
                     const std::array<FLOAT, 3> &y_axis,
                     const std::array<FLOAT, 3> &z_axis)
      : m_value{x_axis, y_axis, z_axis} {};
  /// Constructs a 3x3 matrix with 3D x-axis, y-axis and z-axis vector values.
  explicit Matrix3x3(const Vector3 &x_axis, const Vector3 &y_axis,
                     const Vector3 &z_axis)
      : m_value{x_axis.m_value, y_axis.m_value, z_axis.m_value} {};

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
   * @return 3D vector containing given row values. (Eg: [x_x, y_x, z_x])
   */
  [[nodiscard]] inline Vector3 operator[](size_t index) const {
    return Vector3{m_value[index]};
  }

  /**
   * Returns the matrix row.
   * @param index Row index.
   * @return 3D vector containing given row values. (Eg: [x_x, y_x, z_x])
   */
  [[nodiscard]] inline Vector3 getRow(size_t index) const {
    return Vector3{m_value.getRow(index)};
  }

  /**
   * Returns the matrix column.
   * @param index Column index.
   * @return 3D vector containing given column values. (Eg: [x_x, x_y, x_z])
   */
  [[nodiscard]] inline Vector3 getColumn(size_t index) const {
    return Vector3{m_value.getColumn(index)};
  }

  /// Returns x-axis (column 0) of the matrix.
  [[nodiscard]] inline Vector3 getXAxis() const {
    return Vector3{m_value.getXAxis()};
  }
  /// Returns y-axis (column 1) of the matrix.
  [[nodiscard]] inline Vector3 getYAxis() const {
    return Vector3{m_value.getYAxis()};
  }
  /// Returns y-axis (column 2) of the matrix.
  [[nodiscard]] inline Vector3 getZAxis() const {
    return Vector3{m_value.getZAxis()};
  }

  /// Sets the x-axis column with the given values.
  inline void setXAxis(FLOAT x_x, FLOAT x_y, FLOAT x_z) {
    m_value.setXAxis(x_x, x_y, x_z);
  }
  /// Sets the x-axis column with the given 3D array values.
  inline void setXAxis(const std::array<FLOAT, 3> &x_axis) {
    m_value.setXAxis(x_axis);
  }
  /// Sets the x-axis column with the given 3D vector values.
  inline void setXAxis(const Vector3 &x_axis) {
    m_value.setXAxis(x_axis.m_value);
  }

  /// Sets the y-axis column with the given values.
  inline void setYAxis(FLOAT y_x, FLOAT y_y, FLOAT y_z) {
    m_value.setYAxis(y_x, y_y, y_z);
  }
  /// Sets the y-axis column with the given 3D array values.
  inline void setYAxis(const std::array<FLOAT, 3> &y_axis) {
    m_value.setYAxis(y_axis);
  }
  /// Sets the y-axis column with the given 3D vector values.
  inline void setYAxis(const Vector3 &y_axis) {
    m_value.setYAxis(y_axis.m_value);
  }

  /// Sets the z-axis column with the given values.
  inline void setZAxis(FLOAT z_x, FLOAT z_y, FLOAT z_z) {
    m_value.setZAxis(z_x, z_y, z_z);
  }
  /// Sets the z-axis column with the given 3D array values.
  inline void setZAxis(const std::array<FLOAT, 3> &z_axis) {
    m_value.setZAxis(z_axis);
  }
  /// Sets the z-axis column with the given 3D vector values.
  inline void setZAxis(const Vector3 &z_axis) {
    m_value.setZAxis(z_axis.m_value);
  }

  // Operators
  inline bool operator==(const Matrix3x3 &rhs) const {
    return m_value == rhs.m_value;
  }

  inline bool operator!=(const Matrix3x3 &rhs) const {
    return m_value != rhs.m_value;
  }

  /// Matrix-Matrix multiplication.
  inline Matrix3x3 operator*(const Matrix3x3 &rhs) const {
    return Matrix3x3{m_value * rhs.m_value};
  }

  /// Matrix-3D Column vector multiplication.
  inline Vector3 operator*(const Vector3 &rhs) const {
    return Vector3{m_value * rhs.m_value};
  }

  friend inline Matrix3x3 operator*(const Matrix3x3 &lhs, FLOAT rhs) {
    return Matrix3x3{lhs.m_value * rhs};
  }

  friend inline Matrix3x3 operator*(FLOAT lhs, const Matrix3x3 &rhs) {
    return Matrix3x3{lhs * rhs.m_value};
  }

  // Type-Conversions
  /// Convert to a human-readable string value.
  explicit inline operator std::string() const {
    return static_cast<std::string>(m_value);
  }

  friend inline std::ostream &operator<<(std::ostream &out,
                                         const Matrix3x3 &mat) {
    return out << static_cast<std::string>(mat);
  }

  // Functions
  /// Returns the determinant of the 3x3 matrix.
  [[nodiscard]] inline FLOAT determinant() const {
    return m_value.determinant();
  }

  /// Transposes the given 3x3 matrix.
  static inline Matrix3x3 transpose(const Matrix3x3 &matrix) {
    return Matrix3x3{internal::Matrix3x3::transpose(matrix.m_value)};
  }

  /// Returns the inverse of the given 3x3 matrix, if it exists.
  /// Otherwise returns the original matrix as is.
  static inline Matrix3x3 inverse(const Matrix3x3 &matrix) {
    return Matrix3x3{internal::Matrix3x3::inverse(matrix.m_value)};
  }

 private:
  internal::Matrix3x3 m_value{};

  explicit Matrix3x3(const internal::Matrix3x3 &value) : m_value{value} {};
};
}  // namespace cgmath
#endif  // CGMATH_MATRIX3X3_HPP
