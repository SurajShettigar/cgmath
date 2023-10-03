// Copyright 2023 Suraj Shettigar
// SPDX-License-Identifier: Apache-2.0

#ifndef CGMATH_VECTOR4_HPP
#define CGMATH_VECTOR4_HPP

#include <array>
#include <iostream>
#include <string>

#include "constants.hpp"
#include "internal/vector.hpp"

namespace cgmath {
class Vector2;

class Vector3;

/**
 * 4D Floating-Point Vector
 */
class Vector4 {
 public:
  /// [0.0, 0.0, 0.0, 0.0]
  static const Vector4 ZERO;
  /// [1.0, 1.0, 1.0, 1.0]
  static const Vector4 ONE;
  /// [-1.0, 0.0, 0.0, 0.0]
  static const Vector4 LEFT;
  /// [1.0, 0.0, 0.0, 0.0]
  static const Vector4 RIGHT;
  /// [0.0, 1.0, 0.0, 0.0]
  static const Vector4 UP;
  /// [0.0, -1.0, 0.0, 0.0]
  static const Vector4 DOWN;
  /// [0.0, 0.0, 1.0, 0.0]
  static const Vector4 FORWARD;
  /// [0.0, 0.0, -1.0, 0.0]
  static const Vector4 BACK;
  /// [1.0, 1.0, 1.0, 0.0]
  static const Vector4 DIRECTION;
  /// [1.0, 1.0, 1.0, 1.0]
  static const Vector4 POSITION;

  // Constructors / Destructors
  /// Construct a 0 initialized 4D vector.
  Vector4() = default;

  ~Vector4() = default;

  /// Construct a 4D vector from 4 floating point values.
  explicit Vector4(FLOAT x, FLOAT y, FLOAT z, FLOAT w) : m_value{x, y, z, w} {};

  /// Construct a 4D vector from an array of size 4.
  explicit Vector4(const std::array<FLOAT, 4> &value) : m_value{value} {};

  // Copy
  Vector4(const Vector4 &vec4) {
    if (this != &vec4) m_value = vec4.m_value;
  }

  Vector4 &operator=(const Vector4 &vec4) {
    if (this != &vec4) m_value = vec4.m_value;
    return *this;
  }

  // Type-Conversions
  /// Convert to a 2D vector.
  explicit operator Vector2() const;

  /// Convert to a 3D vector.
  explicit operator Vector3() const;

  /// Convert to a human-readable string value.
  explicit operator std::string() const {
    return "[" + std::to_string(m_value[0]) + ", " +
           std::to_string(m_value[1]) + ", " + std::to_string(m_value[2]) +
           ", " + std::to_string(m_value[3]) + "]";
  }

  // Getters / Setters
  [[nodiscard]] inline FLOAT getX() const { return m_value[0]; }

  [[nodiscard]] inline FLOAT getY() const { return m_value[1]; }

  [[nodiscard]] inline FLOAT getZ() const { return m_value[2]; }

  [[nodiscard]] inline FLOAT getW() const { return m_value[3]; }

  inline void setX(FLOAT x) { m_value[0] = x; }

  inline void setY(FLOAT y) { m_value[1] = y; }

  inline void setZ(FLOAT z) { m_value[2] = z; }

  inline void setW(FLOAT w) { m_value[3] = w; }

  inline void set(FLOAT x, FLOAT y, FLOAT z, FLOAT w) {
    m_value[0] = x;
    m_value[1] = y;
    m_value[2] = z;
    m_value[3] = w;
  }

  inline void set(const std::array<FLOAT, 4> &value) { m_value = value; }

  FLOAT operator[](size_t index) const { return m_value[index]; }

  FLOAT &operator[](size_t index) { return m_value[index]; }

  // Operators
  bool operator==(const Vector4 &rhs) const {
    return approxEqual(m_value[0], rhs.m_value[0]) &&
           approxEqual(m_value[1], rhs.m_value[1]) &&
           approxEqual(m_value[2], rhs.m_value[2]) &&
           approxEqual(m_value[3], rhs.m_value[3]);
  }

  bool operator!=(const Vector4 &rhs) const {
    return approxNotEqual(m_value[0], rhs.m_value[0]) ||
           approxNotEqual(m_value[1], rhs.m_value[1]) ||
           approxNotEqual(m_value[2], rhs.m_value[2]) ||
           approxNotEqual(m_value[3], rhs.m_value[3]);
  }

  Vector4 operator-() const {
    return Vector4{-m_value[0], -m_value[1], -m_value[2], -m_value[3]};
  }

  Vector4 &operator+=(const Vector4 &rhs) {
    m_value[0] += rhs.m_value[0];
    m_value[1] += rhs.m_value[1];
    m_value[2] += rhs.m_value[2];
    m_value[3] += rhs.m_value[3];
    return *this;
  }

  Vector4 &operator-=(const Vector4 &rhs) {
    m_value[0] -= rhs.m_value[0];
    m_value[1] -= rhs.m_value[1];
    m_value[2] -= rhs.m_value[2];
    m_value[3] -= rhs.m_value[3];
    return *this;
  }

  Vector4 &operator*=(FLOAT rhs) {
    m_value[0] *= rhs;
    m_value[1] *= rhs;
    m_value[2] *= rhs;
    m_value[3] *= rhs;
    return *this;
  }

  Vector4 &operator/=(FLOAT rhs) { return *this *= 1.0 / rhs; }

  friend Vector4 operator+(const Vector4 &lhs, const Vector4 &rhs) {
    return Vector4{
        lhs.m_value[0] + rhs.m_value[0], lhs.m_value[1] + rhs.m_value[1],
        lhs.m_value[2] + rhs.m_value[2], lhs.m_value[3] + rhs.m_value[3]};
  }

  friend Vector4 operator-(const Vector4 &lhs, const Vector4 &rhs) {
    return Vector4{
        lhs.m_value[0] - rhs.m_value[0], lhs.m_value[1] - rhs.m_value[1],
        lhs.m_value[2] - rhs.m_value[2], lhs.m_value[3] - rhs.m_value[3]};
  }

  friend Vector4 operator*(const Vector4 &lhs, FLOAT rhs) {
    return Vector4{lhs.m_value[0] * rhs, lhs.m_value[1] * rhs,
                   lhs.m_value[2] * rhs, lhs.m_value[3] * rhs};
  }

  friend Vector4 operator*(FLOAT lhs, const Vector4 &rhs) { return rhs * lhs; }

  friend Vector4 operator/(const Vector4 &lhs, FLOAT rhs) {
    return Vector4{lhs.m_value[0] / rhs, lhs.m_value[1] / rhs,
                   lhs.m_value[2] / rhs, lhs.m_value[3] / rhs};
  }

  friend Vector4 operator/(FLOAT lhs, const Vector4 &rhs) {
    return Vector4{lhs / rhs.m_value[0], lhs / rhs.m_value[1],
                   lhs / rhs.m_value[2], lhs / rhs.m_value[3]};
  }

  friend std::ostream &operator<<(std::ostream &out, const Vector4 &vec4) {
    return out << static_cast<std::string>(vec4);
  }

  // Vector specific operations
  /// Get the squared length / magnitude of the 4D vector.
  [[nodiscard]] inline FLOAT lengthSquared() const {
    return m_value[0] * m_value[0] + m_value[1] * m_value[1] +
           m_value[2] * m_value[2] + m_value[3] * m_value[3];
  }

  /// Get the length / magnitude of the 4D vector.
  [[nodiscard]] inline FLOAT length() const {
    return std::sqrt(lengthSquared());
  }

  /// Normalize the given 4 vector. The 4D vector will have a length of 1.
  [[nodiscard]] inline Vector4 normalized() const { return *this / length(); }

  /**
   * Get the length of the given 4D vector.
   * @param vec4 4D vector whose length / magnitude needs to be found.
   * @return The length / the magnitude.
   */
  static FLOAT length(const Vector4 &vec4) { return vec4.length(); }

  /**
   * Get a normalized vector from the given 4D vector.
   * @param vec4 4D vector to normalize.
   * @return A normalized 4D vector whose length is 1.
   */
  static Vector4 normalize(const Vector4 &vec4) { return vec4.normalized(); }

  /**
   *  Linearly interpolates between two 4D vectors.
   * @param from 4D vector to interpolate from.
   * @param to 4D vector to interpolate to.
   * @param t Scalar value within the  range [0, 1] to control interpolation.
   * @return An interpolated 4D vector.
   */
  static Vector4 lerp(const Vector4 &from, const Vector4 &to, FLOAT t) {
    return from * (1.0 - t) + to * t;
  }

 private:
  std::array<FLOAT, 4> m_value{0.0, 0.0, 0.0, 0.0};
};
}  // namespace cgmath

#endif  // CGMATH_VECTOR4_HPP
