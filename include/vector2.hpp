// Copyright 2023 Suraj Shettigar
// SPDX-License-Identifier: Apache-2.0

#ifndef CGMATH_VECTOR2_HPP
#define CGMATH_VECTOR2_HPP

#include <array>
#include <iostream>
#include <string>

#include "constants.hpp"
#include "internal/vector.hpp"

namespace cgmath {
class Vector3;

class Vector4;

using Vector = internal::Vector;
/**
 * 2D Floating-Point Vector
 */
class Vector2 {
 public:
  /// [0.0, 0.0]
  static const Vector2 ZERO;
  /// [1.0, 1.0]
  static const Vector2 ONE;
  /// [-1.0, 0.0]
  static const Vector2 LEFT;
  /// [1.0, 0.0]
  static const Vector2 RIGHT;
  /// [0.0, 1.0]
  static const Vector2 UP;
  /// [0.0, -1.0]
  static const Vector2 DOWN;

  // Constructors / Destructors
  /// Construct a 0 initialized 2D vector.
  Vector2() = default;

  ~Vector2() = default;

  /// Construct a 2D vector from 2 floating point values.
  explicit Vector2(FLOAT x, FLOAT y) : m_value{x, y} {};

  /// Construct a 2D vector from an array of size 2.
  explicit Vector2(const std::array<FLOAT, 2> &value) : m_value{value} {};

  // Copy
  Vector2(const Vector2 &vec2) {
    if (this != &vec2) m_value = vec2.m_value;
  }

  Vector2 &operator=(const Vector2 &vec2) {
    if (this != &vec2) m_value = vec2.m_value;
    return *this;
  }

  // Type-Conversions
  /// Convert to a 3D vector.
  explicit operator Vector3() const;

  /// Convert to a 4D vector with z, w set to 0.
  explicit operator Vector4() const;

  /// Convert to a human-readable string value.
  explicit inline operator std::string() const {
    return static_cast<std::string>(m_value);
  }

  // Getters / Setters
  [[nodiscard]] inline FLOAT getX() const { return m_value.getX(); }

  [[nodiscard]] inline FLOAT getY() const { return m_value.getY(); }

  inline void setX(FLOAT x) { m_value.setX(x); }

  inline void setY(FLOAT y) { m_value.setY(y); }

  inline void set(FLOAT x, FLOAT y) { m_value.set(x, y); }

  inline void set(const std::array<FLOAT, 2> &value) { m_value.set(value); }

  inline FLOAT operator[](size_t index) const { return m_value[index]; }

  // Operators
  inline bool operator==(const Vector2 &rhs) const {
    return m_value == rhs.m_value;
  }

  inline bool operator!=(const Vector2 &rhs) const {
    return m_value != rhs.m_value;
  }

  inline Vector2 operator-() const { return Vector2{-m_value}; }

  inline Vector2 &operator+=(const Vector2 &rhs) {
    m_value += rhs.m_value;
    return *this;
  }

  inline Vector2 &operator-=(const Vector2 &rhs) {
    m_value -= rhs.m_value;
    return *this;
  }

  inline Vector2 &operator*=(FLOAT rhs) {
    m_value *= rhs;
    ;
    return *this;
  }

  inline Vector2 &operator/=(FLOAT rhs) {
    m_value /= rhs;
    return *this;
  }

  friend inline Vector2 operator+(const Vector2 &lhs, const Vector2 &rhs) {
    return Vector2{lhs.m_value + rhs.m_value};
  }

  friend inline Vector2 operator-(const Vector2 &lhs, const Vector2 &rhs) {
    return Vector2{lhs.m_value - rhs.m_value};
  }

  friend inline Vector2 operator*(const Vector2 &lhs, FLOAT rhs) {
    return Vector2{lhs.m_value * rhs};
  }

  friend inline Vector2 operator*(FLOAT lhs, const Vector2 &rhs) {
    return Vector2{lhs * rhs.m_value};
  }

  friend inline Vector2 operator/(const Vector2 &lhs, FLOAT rhs) {
    return Vector2{lhs.m_value / rhs};
  }

  friend inline Vector2 operator/(FLOAT lhs, const Vector2 &rhs) {
    return Vector2{lhs / rhs.m_value};
  }

  friend inline std::ostream &operator<<(std::ostream &out,
                                         const Vector2 &vec2) {
    return out << static_cast<std::string>(vec2);
  }

  // Vector specific operations
  /// Get the squared length / magnitude of the 2D vector.
  [[nodiscard]] inline FLOAT lengthSquared() const {
    return m_value.lengthSquared();
  }

  /// Get the length / magnitude of the 2D vector.
  [[nodiscard]] inline FLOAT length() const { return m_value.length(); }

  /// Normalize the given 2D vector. The 2D vector will have a length of 1.
  [[nodiscard]] inline Vector2 normalized() const {
    return Vector2{m_value.normalized()};
  }

  /**
   * Get the length of the given 2D vector.
   * @param vec2 2D vector whose length / magnitude needs to be found.
   * @return The length / the magnitude.
   */
  static inline FLOAT length(const Vector2 &vec2) {
    return Vector::length(vec2.m_value);
  }

  /**
   * Get a normalized vector from the given 3D vector.
   * @param vec2 2D vector to normalize.
   * @return A normalized 3D vector whose length is 1.
   */
  static inline Vector2 normalize(const Vector2 &vec2) {
    return Vector2{Vector::normalize(vec2.m_value)};
  }

  /**
   * Find the dot product of two 2D vectors. Dot product of two normalized
   * vectors gives the cosine of the angle between them.
   * @param lhs A 2D vector.
   * @param rhs Another 2D vector.
   * @return A scalar floating point value. In case of normalized vectors, value
   * will be in the range of [-1, 1]. With -1 indicating vectors facing in the
   * opposite direction, 0 indicating orthogonal / perpendicular vectors and 1
   * indicating vectors facing in the same direction.
   */
  static inline FLOAT dot(const Vector2 &lhs, const Vector2 &rhs) {
    return Vector::dot(lhs.m_value, rhs.m_value);
  }

  /**
   * Find the angle between two 2D vectors. In other words, the first vector is
   * rotated by that angle to produce the second vector.
   * @param lhs First 2D vector.
   * @param rhs Other 2D vector.
   * @return Angle in radians between two 2D vectors.
   */
  static inline FLOAT angle(const Vector2 &lhs, const Vector2 &rhs) {
    return acos(dot(lhs, rhs) / (lhs.length() * rhs.length()));
  }

  /**
   *  Linearly interpolates between two 2D vectors.
   * @param from 2D vector to interpolate from.
   * @param to 2D vector to interpolate to.
   * @param t Scalar value within the  range [0, 1] to control interpolation.
   * @return An interpolated 2D vector.
   */
  static inline Vector2 lerp(const Vector2 &from, const Vector2 &to, FLOAT t) {
    return Vector2{Vector::lerp(from.m_value, to.m_value, t)};
  }

  /**
   * Project a 2D vector onto another 2D vector.
   * @param from The 2D vector to be projected
   * @param to The 2D vector onto which the vector is projected.
   * @return 2D vector in the projected direction but with a scaled magnitude.
   */
  static inline Vector2 project(const Vector2 &from, const Vector2 &to) {
    return to * (dot(from, to) / to.lengthSquared());
  }

 private:
  Vector m_value{0.0};

  explicit Vector2(const Vector &value) : m_value{value} {};
};
}  // namespace cgmath

#endif  // CGMATH_VECTOR2_HPP
