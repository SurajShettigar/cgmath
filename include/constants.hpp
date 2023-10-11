// Copyright 2023 Suraj Shettigar
// SPDX-License-Identifier: Apache-2.0

#ifndef CGMATH_CONSTANTS_HPP
#define CGMATH_CONSTANTS_HPP

#include <cmath>
#include <limits>
#include <cstring>
#include <array>
#include <string>
#include <iostream>

namespace cgmath {
#ifdef USE_DOUBLE
typedef double FLOAT;
#else
typedef float FLOAT;
#endif  // USE_DOUBLE

/**
 * Different coordinate systems / frames. Can be either Right-Hand or Left-Hand
 * with either Y-UP or Z-UP.
 */
enum class CoordinateSystem {
  /// Right-Hand System with Y axis as up direction.
  RHS_Y_UP = 0,
  /// Left-Hand System with Y axis as up direction.
  LHS_Y_UP = 1,
  /// Right-Hand System with Z axis as up direction.
  RHS_Z_UP = 2,
  /// Left-Hand System with Z axis as up direction.
  LHS_Z_UP = 3,
};

static constexpr FLOAT MACHINE_EPSILON = std::numeric_limits<FLOAT>::epsilon();
static constexpr FLOAT EPSILON = MACHINE_EPSILON * 4.0;
static constexpr FLOAT INFINITE = std::numeric_limits<FLOAT>::infinity();

/// π
static constexpr FLOAT PI = 3.141592653589793;
/// 1 / π
static constexpr FLOAT ONE_OVER_PI = 0.318309886183791;
/// 180 / π
static constexpr FLOAT RADIANS_TO_DEGREES = 57.295779513082321;
/// π / 180
static constexpr FLOAT DEGREES_TO_RADIANS = 0.017453292519943;

/// Checks if two floating point values are equal. Because of floating point
/// imprecision,
/// == operator cannot be used directly.
inline bool approxEqual(FLOAT lhs, FLOAT rhs) {
  // Refer: https://realtimecollisiondetection.net/blog/?p=89
  return std::fabs(rhs - lhs) <=
         EPSILON * std::max(static_cast<FLOAT>(1.0),
                            std::max(std::abs(lhs), std::abs(rhs)));
}

/// Checks if two floating point values are not equal. Because of floating point
/// imprecision,
/// != operator cannot be used directly.
inline bool approxNotEqual(FLOAT lhs, FLOAT rhs) {
  return !approxEqual(lhs, rhs);
}

/**
 * Limits the given value to the specified range.
 * @param value The value to limit.
 * @param from The minimum value of the range.
 * @param to The maximum value of the range.
 * @return Clamped value.
 */
inline FLOAT clamp(FLOAT value, FLOAT from, FLOAT to) {
  return std::min(std::max(value, from), to);
}

/**
 * Transforms the given value with an input range to the output range.
 * @param value The value to remap.
 * @param from_min The minimum value of the input range.
 * @param from_max The maximum value of the input range.
 * @param to_min The minimum value of the output range.
 * @param to_max The maximum value of the output range.
 * @return Value lying within the output range.
 */
inline FLOAT remap(FLOAT value, FLOAT from_min, FLOAT from_max, FLOAT to_min,
                   FLOAT to_max) {
  return (value - from_min) * (to_max - to_min) / (from_max - from_min) +
         to_min;
}

/// Converts radian angle to degree.
inline static FLOAT degree(FLOAT radian_angle) {
  return radian_angle * RADIANS_TO_DEGREES;
}

/// Converts degree angle to radian.
inline static FLOAT radian(FLOAT degree_angle) {
  return degree_angle * DEGREES_TO_RADIANS;
}
}  // namespace cgmath

#endif  // CGMATH_CONSTANTS_HPP
