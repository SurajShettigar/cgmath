// Copyright 2023 Suraj Shettigar
// SPDX-License-Identifier: Apache-2.0

#ifndef CGMATH_CONSTANTS_HPP
#define CGMATH_CONSTANTS_HPP

#include <cmath>
#include <limits>

namespace cgmath {
#ifdef USE_DOUBLE
typedef double FLOAT;
#else
typedef float FLOAT;
#endif  // USE_DOUBLE

#ifdef USE_INTRINSICS
/**
 * x86-64 Micro-architecture levels defined based on the availability of
 * Instruction Set Extensions / features for SIMD.
 */
enum class SIMDMicroArchitectureLevel {
  /// ALL CPUs which support only upto SSE2, SSE instruction set.
  BASE = 0,
  /// ALL CPUs which support upto SSE4_2, SSSE3 instruction set.
  VERSION_2 = 1,
  /// ALL CPUs which support upto AVX2, AVX instruction set.
  VERSION_3 = 2,
  /// ALL CPUs which support AVX512 instruction set.
  VERSION_4 = 3
};
#if defined(__AVX512F__)
// 64 byte vectorization

#include <immintrin.h>

static constexpr SIMDMicroArchitectureLevel SIMD_LEVEL =
    SIMDMicroArchitectureLevel::VERSION_4;

#if USE_DOUBLE
typedef std::array<__m512d, 2> MATRIX;
typedef __m256d VECTOR;
#else
typedef __m512 MATRIX;
typedef __m128 VECTOR;
#endif

#elif defined(__AVX2__) || defined(__AVX__)
// 32 byte vectorization

#include <immintrin.h>

static constexpr SIMDMicroArchitectureLevel SIMD_LEVEL =
    SIMDMicroArchitectureLevel::VERSION_3;

#if USE_DOUBLE
typedef std::array<__m256d, 4> MATRIX;
typedef __m256d VECTOR;
#else
typedef std::array<__m128, 4> MATRIX;
typedef __m128 VECTOR;
#endif

#else
// 16 byte vectorization

#include <emmintrin.h>
static constexpr SIMDMicroArchitectureLevel SIMD_LEVEL =
    SIMDMicroArchitectureLevel::BASE;

#if USE_DOUBLE
typedef std::array<__m128d, 8> MATRIX;
typedef std::array<__m128d, 2> VECTOR;
#else
typedef std::array<__m128, 4> MATRIX;
typedef __m128 VECTOR;
#endif

#endif
#else
typdef std::array<FLOAT, 16> MATRIX;
typdef std::array<FLOAT, 4> VECTOR;
#endif  // USE_INTRINSICS

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
  return std::abs(rhs - lhs) <=
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
