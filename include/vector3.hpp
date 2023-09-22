// Copyright 2023 Suraj Shettigar
// SPDX-License-Identifier: Apache-2.0

#ifndef CGMATH_VECTOR3_HPP
#define CGMATH_VECTOR3_HPP

#include <array>
#include <iostream>
#include <string>

#include "constants.hpp"

namespace cgmath {
class Vector2;

class Vector4;

/**
 * 3D Floating-Point Vector
 */
class Vector3 {
 public:
  /// [0.0, 0.0, 0.0]
  static const Vector3 ZERO;
  /// [1.0, 1.0, 1.0]
  static const Vector3 ONE;
  /// [-1.0, 0.0, 0.0]
  static const Vector3 LEFT;
  /// [1.0, 0.0, 0.0]
  static const Vector3 RIGHT;
  /// [0.0, 1.0, 0.0]
  static const Vector3 UP;
  /// [0.0, -1.0, 0.0]
  static const Vector3 DOWN;
  /// [0.0, 0.0, 1.0]
  static const Vector3 FORWARD;
  /// [0.0, 0.0, -1.0]
  static const Vector3 BACK;

  // Constructors / Destructors
  /// Construct a 0 initialized 3D vector.
  Vector3() = default;

  ~Vector3() = default;

  /// Construct a 3D vector from 3 floating point values.
  explicit Vector3(FLOAT x, FLOAT y, FLOAT z) : m_value{x, y, z} {};

  /// Construct a 3D vector from an array of size 3.
  explicit Vector3(const std::array<FLOAT, 3> &value) : m_value{value} {};

  // Copy
  Vector3(const Vector3 &vec3) {
    if (this != &vec3) m_value = vec3.m_value;
  }

  Vector3 &operator=(const Vector3 &vec3) {
    if (this != &vec3) m_value = vec3.m_value;
    return *this;
  }

  // Type-Conversions
  /// Convert to a 2D vector.
  explicit operator Vector2() const;

  /// Convert to a 4D vector with w set to 0.
  explicit operator Vector4() const;

  /// Convert to a human-readable string value.
  explicit operator std::string() const {
    return "[" + std::to_string(m_value[0]) + ", " +
           std::to_string(m_value[1]) + ", " + std::to_string(m_value[2]) + "]";
  }

  // Getters / Setters
  [[nodiscard]] inline FLOAT getX() const { return m_value[0]; }

  [[nodiscard]] inline FLOAT getY() const { return m_value[1]; }

  [[nodiscard]] inline FLOAT getZ() const { return m_value[2]; }

  inline void setX(FLOAT x) { m_value[0] = x; }

  inline void setY(FLOAT y) { m_value[1] = y; }

  inline void setZ(FLOAT z) { m_value[2] = z; }

  inline void set(FLOAT x, FLOAT y, FLOAT z) {
    m_value[0] = x;
    m_value[1] = y;
    m_value[2] = z;
  }

  inline void set(const std::array<FLOAT, 3> &value) { m_value = value; }

  FLOAT operator[](size_t index) const { return m_value[index]; }

  FLOAT &operator[](size_t index) { return m_value[index]; }

  // Operators
  bool operator==(const Vector3 &rhs) const {
    return approxEqual(m_value[0], rhs.m_value[0]) &&
           approxEqual(m_value[1], rhs.m_value[1]) &&
           approxEqual(m_value[2], rhs.m_value[2]);
  }

  bool operator!=(const Vector3 &rhs) const {
    return approxNotEqual(m_value[0], rhs.m_value[0]) ||
           approxNotEqual(m_value[1], rhs.m_value[1]) ||
           approxNotEqual(m_value[2], rhs.m_value[2]);
  }

  Vector3 operator-() const {
    return Vector3{-m_value[0], -m_value[1], -m_value[2]};
  }

  Vector3 &operator+=(const Vector3 &rhs) {
    m_value[0] += rhs.m_value[0];
    m_value[1] += rhs.m_value[1];
    m_value[2] += rhs.m_value[2];
    return *this;
  }

  Vector3 &operator-=(const Vector3 &rhs) {
    m_value[0] -= rhs.m_value[0];
    m_value[1] -= rhs.m_value[1];
    m_value[2] -= rhs.m_value[2];
    return *this;
  }

  Vector3 &operator*=(FLOAT rhs) {
    m_value[0] *= rhs;
    m_value[1] *= rhs;
    m_value[2] *= rhs;
    return *this;
  }

  Vector3 &operator/=(FLOAT rhs) { return *this *= 1.0 / rhs; }

  friend Vector3 operator+(const Vector3 &lhs, const Vector3 &rhs) {
    return Vector3{lhs.m_value[0] + rhs.m_value[0],
                   lhs.m_value[1] + rhs.m_value[1],
                   lhs.m_value[2] + rhs.m_value[2]};
  }

  friend Vector3 operator-(const Vector3 &lhs, const Vector3 &rhs) {
    return Vector3{lhs.m_value[0] - rhs.m_value[0],
                   lhs.m_value[1] - rhs.m_value[1],
                   lhs.m_value[2] - rhs.m_value[2]};
  }

  friend Vector3 operator*(const Vector3 &lhs, FLOAT rhs) {
    return Vector3{lhs.m_value[0] * rhs, lhs.m_value[1] * rhs,
                   lhs.m_value[2] * rhs};
  }

  friend Vector3 operator*(FLOAT lhs, const Vector3 &rhs) { return rhs * lhs; }

  friend Vector3 operator/(const Vector3 &lhs, FLOAT rhs) {
    return Vector3{lhs.m_value[0] / rhs, lhs.m_value[1] / rhs,
                   lhs.m_value[2] / rhs};
  }

  friend Vector3 operator/(FLOAT lhs, const Vector3 &rhs) {
    return Vector3{lhs / rhs.m_value[0], lhs / rhs.m_value[1],
                   lhs / rhs.m_value[2]};
  }

  friend std::ostream &operator<<(std::ostream &out, const Vector3 &vec3) {
    return out << static_cast<std::string>(vec3);
  }

  // Vector specific operations
  /// Get the squared length / magnitude of the 3D vector.
  [[nodiscard]] inline FLOAT lengthSquared() const {
    return m_value[0] * m_value[0] + m_value[1] * m_value[1] +
           m_value[2] * m_value[2];
  }

  /// Get the length / magnitude of the 3D vector.
  [[nodiscard]] inline FLOAT length() const {
    return std::sqrt(lengthSquared());
  }

  /// Normalize the given 3D vector. The 3D vector will have a length of 1.
  [[nodiscard]] inline Vector3 normalized() const { return *this / length(); }

  /**
   * Get the length of the given 3D vector.
   * @param vec3 3D vector whose length / magnitude needs to be found.
   * @return The length / the magnitude.
   */
  static FLOAT length(const Vector3 &vec3) { return vec3.length(); }

  /**
   * Get a normalized vector from the given 3D vector.
   * @param vec3 3D vector to normalize.
   * @return A normalized 3D vector whose length is 1.
   */
  static Vector3 normalize(const Vector3 &vec3) { return vec3.normalized(); }

  /**
   * Find the dot product of two 3D vectors. Dot product of two normalized
   * vectors gives the cosine of the angle between them.
   * @param lhs A 3D vector.
   * @param rhs Another 3D vector.
   * @return A scalar floating point value. In case of normalized vectors, value
   * will be in the range of [-1, 1]. With -1 indicating vectors facing in the
   * opposite direction, 0 indicating orthogonal / perpendicular vectors and 1
   * indicating vectors facing in the same direction.
   */
  static FLOAT dot(const Vector3 &lhs, const Vector3 &rhs) {
    return lhs.m_value[0] * rhs.m_value[0] + lhs.m_value[1] * rhs.m_value[1] +
           lhs.m_value[2] * rhs.m_value[2];
  }

  /**
   * Find the angle between two 3D vectors. In other words, the first vector is
   * rotated by that angle to produce the second vector.
   * @param lhs First 3D vector.
   * @param rhs Other 3D vector.
   * @return Angle in radians between two 3D vectors.
   */
  static FLOAT angle(const Vector3 &lhs, const Vector3 &rhs) {
    return acos(dot(lhs, rhs) / (lhs.length() * rhs.length()));
  }

  /**
   * Find the cross product of two 3D vectors. By using two normalized, linearly
   * independent (not lying in the same plane) 3D vectors, we get another
   * normalized 3D vector which is orthogonal / perpendicular to these two
   * vectors. Cross product is not commutative. Switching LHS and RHS gives a
   * vector facing in the opposite direction.
   * @param lhs First vector.
   * @param rhs Second vector.
   * @return A 3D vector. Vector is orthogonal to LHS and RHS vector if LHS and
   * RHS are linearly independent.
   */
  static Vector3 cross(const Vector3 &lhs, const Vector3 &rhs) {
#ifndef USE_DOUBLE
    // Static cast to double to ensure there are no rounding errors when two
    // vectors are extremely similar.
    auto vx = static_cast<double>(lhs.m_value[0]);
    auto vy = static_cast<double>(lhs.m_value[1]);
    auto vz = static_cast<double>(lhs.m_value[2]);
    auto wx = static_cast<double>(rhs.m_value[0]);
    auto wy = static_cast<double>(rhs.m_value[1]);
    auto wz = static_cast<double>(rhs.m_value[2]);
    return Vector3(static_cast<FLOAT>(vy * wz - vz * wy),
                   static_cast<FLOAT>(vz * wx - vx * wz),
                   static_cast<FLOAT>(vx * wy - vy * wx));
#endif
    return Vector3{
        lhs.m_value[1] * rhs.m_value[2] - lhs.m_value[2] * rhs.m_value[1],
        lhs.m_value[2] * rhs.m_value[0] - lhs.m_value[0] * rhs.m_value[2],
        lhs.m_value[0] * rhs.m_value[1] - lhs.m_value[1] * rhs.m_value[0]};
  }

  /**
   *  Linearly interpolates between two 3D vectors.
   * @param from 3D vector to interpolate from.
   * @param to 3D vector to interpolate to.
   * @param t Scalar value within the  range [0, 1] to control interpolation.
   * @return An interpolated 3D vector.
   */
  static Vector3 lerp(const Vector3 &from, const Vector3 &to, FLOAT t) {
    return from * (1.0 - t) + to * t;
  }

  /**
   * Creates a reflected 3D vector on a plane defined by a normal vector.
   * @param vec The 3D vector to reflect.
   * @param normal The normal direction of the plane.
   * @return The reflected 3D vector.
   */
  static Vector3 reflect(const Vector3 &vec, const Vector3 &normal) {
    return vec - 2.0 * dot(vec, normal) * normal;
  }

  /**
   * Project a 3D vector onto another 3D vector.
   * @param from The 3D vector to be projected
   * @param to The 3D vector onto which the vector is projected.
   * @return 3D vector in the projected direction but with a scaled magnitude.
   */
  static Vector3 project(const Vector3 &from, const Vector3 &to) {
    return to * (dot(from, to) / to.lengthSquared());
  }

  /**
   * Projects a 3D vector onto a plane defined by the given normal.
   * @param vec The 3D vector to be projected.
   * @param planeNormal The normal direction of the plane.
   * @return Projected 3D vector.
   */
  static Vector3 projectOnPlane(const Vector3 &vec,
                                const Vector3 &planeNormal) {
    return vec - project(vec, planeNormal);
  }

  /**
   * Creates a coordinate frame / orthonormal basis around the given 3D vector.
   * Useful for surface shading and forming a basis for a camera.
   * @param normal Input 3D vector around which the coordinate frame is formed.
   * For surface shading, the normal of the surface (vector pointing out of the
   * surface) is used. And for cameras, the camera's facing direction is used.
   * @param tangent Output 3D vector perpendicular to the normal vector. In case
   * of surfaces, this vector forms a tangent to the surface.
   * @param binormal Output 3D vector perpendicular to both normal and tangent
   * vector.
   */
  static void getCoordinateFrame(const Vector3 &normal, Vector3 *tangent,
                                 Vector3 *binormal) {
    // Refer: https://graphics.pixar.com/library/OrthonormalB/paper.pdf
    auto z_sign = std::copysign(1.0, normal.m_value[2]);
    const FLOAT a = -1.0 / (z_sign + normal.m_value[2]);
    const FLOAT b = normal.m_value[0] * normal.m_value[1] * a;

    tangent->m_value[0] =
        1.0 + z_sign * normal.m_value[0] * normal.m_value[0] * a;
    tangent->m_value[1] = z_sign * b;
    tangent->m_value[2] = -z_sign * normal.m_value[0];

    binormal->m_value[0] = b;
    binormal->m_value[1] = z_sign + normal.m_value[1] * normal.m_value[1] * a;
    binormal->m_value[2] = -normal.m_value[1];
  }

 private:
  std::array<FLOAT, 3> m_value{0.0, 0.0, 0.0};
};
}  // namespace cgmath

#endif  // CGMATH_VECTOR3_HPP
