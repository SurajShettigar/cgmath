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

using Vector = internal::Vector;

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
  explicit inline operator std::string() const {
    return static_cast<std::string>(m_value);
  }

  // Getters / Setters
  [[nodiscard]] inline FLOAT getX() const { return m_value.getX(); }

  [[nodiscard]] inline FLOAT getY() const { return m_value.getY(); }

  [[nodiscard]] inline FLOAT getZ() const { return m_value.getZ(); }

  [[nodiscard]] inline FLOAT getW() const { return m_value.getW(); }

  inline void setX(FLOAT x) { m_value.setX(x); }

  inline void setY(FLOAT y) { m_value.setY(y); }

  inline void setZ(FLOAT z) { m_value.setZ(z); }

  inline void setW(FLOAT w) { m_value.setW(w); }

  inline void set(FLOAT x, FLOAT y, FLOAT z, FLOAT w) {
    m_value.set(x, y, z, w);
  }

  inline void set(const std::array<FLOAT, 4> &value) { m_value.set(value); }

  inline FLOAT operator[](size_t index) const { return m_value[index]; }

  // Operators
  inline bool operator==(const Vector4 &rhs) const {
    return m_value == rhs.m_value;
  }

  inline bool operator!=(const Vector4 &rhs) const {
    return m_value != rhs.m_value;
  }

  inline Vector4 operator-() const { return Vector4(-m_value); }

  inline Vector4 &operator+=(const Vector4 &rhs) {
    m_value += rhs.m_value;
    return *this;
  }

  inline Vector4 &operator-=(const Vector4 &rhs) {
    m_value -= rhs.m_value;
    return *this;
  }

  inline Vector4 &operator*=(FLOAT rhs) {
    m_value *= rhs;
    return *this;
  }

  inline Vector4 &operator/=(FLOAT rhs) {
    m_value /= rhs;
    return *this;
  }

  friend inline Vector4 operator+(const Vector4 &lhs, const Vector4 &rhs) {
    return Vector4{lhs.m_value + rhs.m_value};
  }

  friend inline Vector4 operator-(const Vector4 &lhs, const Vector4 &rhs) {
    return Vector4{lhs.m_value - rhs.m_value};
  }

  friend inline Vector4 operator*(const Vector4 &lhs, FLOAT rhs) {
    return Vector4{lhs.m_value * rhs};
  }

  friend inline Vector4 operator*(FLOAT lhs, const Vector4 &rhs) {
    return Vector4{lhs * rhs.m_value};
  }

  friend inline Vector4 operator/(const Vector4 &lhs, FLOAT rhs) {
    return Vector4{lhs.m_value / rhs};
  }

  friend inline Vector4 operator/(FLOAT lhs, const Vector4 &rhs) {
    return Vector4{lhs / rhs.m_value};
  }

  friend inline std::ostream &operator<<(std::ostream &out,
                                         const Vector4 &vec4) {
    return out << static_cast<std::string>(vec4);
  }

  // Vector specific operations
  /// Get the squared length / magnitude of the 4D vector.
  [[nodiscard]] inline FLOAT lengthSquared() const {
    return m_value.lengthSquared();
  }

  /// Get the length / magnitude of the 4D vector.
  [[nodiscard]] inline FLOAT length() const { return m_value.length(); }

  /// Normalize the given 4 vector. The 4D vector will have a length of 1.
  [[nodiscard]] inline Vector4 normalized() const {
    return Vector4{m_value.normalized()};
  }

  /**
   * Get the length of the given 4D vector.
   * @param vec4 4D vector whose length / magnitude needs to be found.
   * @return The length / the magnitude.
   */
  static inline FLOAT length(const Vector4 &vec4) {
    return Vector::length(vec4.m_value);
  }

  /**
   * Get a normalized vector from the given 4D vector.
   * @param vec4 4D vector to normalize.
   * @return A normalized 4D vector whose length is 1.
   */
  static inline Vector4 normalize(const Vector4 &vec4) {
    return Vector4{Vector::normalize(vec4.m_value)};
  }

  /**
   *  Linearly interpolates between two 4D vectors.
   * @param from 4D vector to interpolate from.
   * @param to 4D vector to interpolate to.
   * @param t Scalar value within the  range [0, 1] to control interpolation.
   * @return An interpolated 4D vector.
   */
  static inline Vector4 lerp(const Vector4 &from, const Vector4 &to, FLOAT t) {
    return Vector4{Vector::lerp(from.m_value, to.m_value, t)};
  }

 private:
  Vector m_value{0.0};

  explicit Vector4(const Vector &vector) : m_value{vector} {};
};
}  // namespace cgmath

#endif  // CGMATH_VECTOR4_HPP
