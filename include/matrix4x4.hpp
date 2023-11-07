// Copyright 2023 Suraj Shettigar
// SPDX-License-Identifier: Apache-2.0

#ifndef CGMATH_MATRIX4X4_HPP
#define CGMATH_MATRIX4X4_HPP

#include "internal/matrix4x4.hpp"
#include "vector3.hpp"
#include "vector4.hpp"

namespace cgmath {

class Matrix4x4 {
 public:
  // Constructors / Destructors

  /// Constructs an identity 4x4 matrix.
  explicit Matrix4x4() : m_value{} {};

  /// Constructs a 4x4 matrix with x, y, z and w (translation in the case of 3D
  /// affine transform matrix) column values.
  explicit Matrix4x4(FLOAT x_x, FLOAT x_y, FLOAT x_z, FLOAT x_w, FLOAT y_x,
                     FLOAT y_y, FLOAT y_z, FLOAT y_w, FLOAT z_x, FLOAT z_y,
                     FLOAT z_z, FLOAT z_w, FLOAT t_x, FLOAT t_y, FLOAT t_z,
                     FLOAT t_w)
      : m_value{x_x, x_y, x_z, x_w, y_x, y_y, y_z, y_w,
                z_x, z_y, z_z, z_w, t_x, t_y, t_z, t_w} {};

  /// Constructs a 4x4 matrix with 3D x, y, z and w (translation in the case of
  /// 3D affine transform matrix) array values.
  explicit Matrix4x4(const std::array<FLOAT, 4> &x_axis,
                     const std::array<FLOAT, 4> &y_axis,
                     const std::array<FLOAT, 4> &z_axis,
                     const std::array<FLOAT, 4> &translation)
      : m_value{x_axis, y_axis, z_axis, translation} {};

  /// Constructs a 4x4 matrix with 3D x, y, z and w (translation in the case of
  /// 3D affine transform matrix) vector values.
  explicit Matrix4x4(const Vector4 &x_axis, const Vector4 &y_axis,
                     const Vector4 &z_axis, const Vector4 &translation)
      : m_value{x_axis.m_value, y_axis.m_value, z_axis.m_value,
                translation.m_value} {};

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
  [[nodiscard]] inline Vector4 operator[](size_t index) const {
    return Vector4{m_value[index]};
  }

  /**
   * Returns the matrix row.
   * @param index Row index.
   * @return 4D vector containing given row values. (Eg: [x_x, y_x, z_x, t_x])
   */
  [[nodiscard]] inline Vector4 getRow(size_t index) const {
    return Vector4{m_value.getRow(index)};
  }

  /**
   * Returns the matrix column.
   * @param index Column index.
   * @return 4D vector containing given column values. (Eg: [x_x, x_y, x_z,
   * x_w])
   */
  [[nodiscard]] inline Vector4 getColumn(size_t index) const {
    return Vector4{m_value.getColumn(index)};
  }

  /// Returns x-axis (column 0) of the matrix.
  [[nodiscard]] inline Vector4 getXAxis() const {
    return Vector4{m_value.getXAxis()};
  }
  /// Returns y-axis (column 1) of the matrix.
  [[nodiscard]] inline Vector4 getYAxis() const {
    return Vector4{m_value.getYAxis()};
  }
  /// Returns z-axis (column 2) of the matrix.
  [[nodiscard]] inline Vector4 getZAxis() const {
    return Vector4{m_value.getZAxis()};
  }
  /// Returns the translation values (w-axis) (column 3) of the matrix.
  [[nodiscard]] inline Vector4 getTranslation() const {
    return Vector4{m_value.getTranslation()};
  }

  /// Sets the x-axis column with the given values.
  inline void setXAxis(FLOAT x_x, FLOAT x_y, FLOAT x_z, FLOAT x_w) {
    m_value.setXAxis(x_x, x_y, x_z, x_w);
  }
  /// Sets the x-axis column with the given 4D array values.
  inline void setXAxis(const std::array<FLOAT, 4> &x_axis) {
    m_value.setXAxis(x_axis);
  }
  /// Sets the x-axis column with the given 4D vector values.
  inline void setXAxis(const Vector4 &x_axis) {
    m_value.setXAxis(x_axis.m_value);
  }

  /// Sets the y-axis column with the given values.
  inline void setYAxis(FLOAT y_x, FLOAT y_y, FLOAT y_z, FLOAT y_w) {
    m_value.setYAxis(y_x, y_y, y_z, y_w);
  }
  /// Sets the y-axis column with the given 4D array values.
  inline void setYAxis(const std::array<FLOAT, 4> &y_axis) {
    m_value.setYAxis(y_axis);
  }
  /// Sets the y-axis column with the given 4D vector values.
  inline void setYAxis(const Vector4 &y_axis) {
    m_value.setYAxis(y_axis.m_value);
  }

  /// Sets the z-axis column with the given values.
  inline void setZAxis(FLOAT z_x, FLOAT z_y, FLOAT z_z, FLOAT z_w) {
    m_value.setZAxis(z_x, z_y, z_z, z_w);
  }
  /// Sets the z-axis column with the given 4D array values.
  inline void setZAxis(const std::array<FLOAT, 4> &z_axis) {
    m_value.setZAxis(z_axis);
  }
  /// Sets the z-axis column with the given 4D vector values.
  inline void setZAxis(const Vector4 &z_axis) {
    m_value.setZAxis(z_axis.m_value);
  }

  /// Sets the translation column (w-axis) with the given values.
  inline void setTranslation(FLOAT t_x, FLOAT t_y, FLOAT t_z, FLOAT t_w) {
    m_value.setTranslation(t_x, t_y, t_z, t_w);
  }
  /// Sets the translation column (w-axis) with the given 4D array values.
  inline void setTranslation(const std::array<FLOAT, 4> &translation) {
    m_value.setTranslation(translation);
  }
  /// Sets the translation column (w-axis) with the given 4D vector values.
  inline void setTranslation(const Vector4 &translation) {
    m_value.setTranslation(translation.m_value);
  }

  // Operators

  inline bool operator==(const Matrix4x4 &rhs) const {
    return m_value == rhs.m_value;
  }

  inline bool operator!=(const Matrix4x4 &rhs) const {
    return m_value != rhs.m_value;
  }

  /// Matrix-Matrix multiplication.
  friend inline Matrix4x4 operator*(const Matrix4x4 &lhs,
                                    const Matrix4x4 &rhs) {
    return Matrix4x4{lhs.m_value * rhs.m_value};
  }

  inline Vector4 operator*(const Vector4 &rhs) const {
    return Vector4{m_value * rhs.m_value};
  }

  // Type-Conversions
  /// Convert to a human-readable string value.
  explicit inline operator std::string() const {
    return static_cast<std::string>(m_value);
  }

  friend inline std::ostream &operator<<(std::ostream &out,
                                         const Matrix4x4 &mat) {
    return out << static_cast<std::string>(mat);
  }

  // Functions

  /// Returns the transpose of the given 4x4 matrix.
  static inline Matrix4x4 transpose(const Matrix4x4 &matrix) {
    return Matrix4x4{internal::Matrix4x4::transpose(matrix.m_value)};
  }

  /// Calculates the determinant of the given 4x4 matrix.
  [[nodiscard]] inline FLOAT determinant() const {
    return m_value.determinant();
  }

  /// Returns the inverse of the given 4x4 matrix.
  static inline Matrix4x4 inverse(const Matrix4x4 &matrix) {
    return Matrix4x4{internal::Matrix4x4::inverse(matrix.m_value)};
  }

  /**
   * Creates a 3D scaling matrix from the given x, y and z axis scales.
   * @param scale 3D vector containing the x, y and z axis scales.
   * @return 4x4 Affine transform matrix containing only the 3D scale component.
   */
  static inline Matrix4x4 scale(const Vector3 &scale) {
    return Matrix4x4{internal::Matrix4x4::scale(scale.m_value)};
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
    return Matrix4x4{internal::Matrix4x4::rotationX(angle)};
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
    return Matrix4x4{internal::Matrix4x4::rotationY(angle)};
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
    return Matrix4x4{internal::Matrix4x4::rotationZ(angle)};
  }

  /**
   * Creates a 3D rotation matrix from the given euler angles in degree.
   * @param angles A 3D vector containing euler angles for each axis. The angle
   * is assumed positive for counter-clockwise rotation and negative for
   * clockwise rotation.
   * @return 4x4 Affine transform matrix containing only the 3D rotational
   * component.
   */
  static inline Matrix4x4 rotation(const Vector3 &euler_angles) {
    return Matrix4x4{internal::Matrix4x4::rotation(euler_angles.m_value)};
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
  static inline Matrix4x4 rotationOverAxis(FLOAT angle, const Vector3 &axis) {
    return Matrix4x4{
        internal::Matrix4x4::rotationOverAxis(angle, axis.m_value)};
  }

  /**
   * Creates a 3D translation matrix from the given x, y and z axis values.
   * @param translate 3D vector containing the x, y and z axis translation
   * values.
   * @return 4x4 Affine transform matrix containing only the 3D translation
   * component.
   */
  static inline Matrix4x4 translation(const Vector3 &translate) {
    return Matrix4x4{internal::Matrix4x4::translation(translate.m_value)};
  }

  /**
   * Returns the inverse of the given 4x4 affine transformation matrix.
   * The matrix should have only translation, rotation and scaling components.
   * @param mat The transform matrix to invert.
   * @return 4x4 inverse affine transformation matrix
   */
  static inline Matrix4x4 transformInverse(const Matrix4x4 &matrix) {
    return Matrix4x4{internal::Matrix4x4::transformInverse(matrix.m_value)};
  }

  /**
   * Returns the inverse of the given 4x4 affine transformation matrix with unit
   * scale. The matrix should have only translation and rotation components.
   * @param mat The transform matrix to invert.
   * @return 4x4 inverse affine transformation matrix
   */
  static inline Matrix4x4 transformInverseUnitScale(const Matrix4x4 &matrix) {
    return Matrix4x4{
        internal::Matrix4x4::transformInverseUnitScale(matrix.m_value)};
  }

 private:
  internal::Matrix4x4 m_value{};

  explicit Matrix4x4(const internal::Matrix4x4 &value) : m_value{value} {};
};
}  // namespace cgmath
#endif  // CGMATH_MATRIX4X4_HPP
