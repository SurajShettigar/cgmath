// Copyright 2023 Suraj Shettigar
// SPDX-License-Identifier: Apache-2.0

#ifndef CGMATH_INTERNAL_VECTOR_HPP
#define CGMATH_INTERNAL_VECTOR_HPP

#include "../constants.hpp"
#include "simd.hpp"

namespace cgmath::internal {

class Vector {
 public:
  // Constructors / Destructors
  /// Construct a 0 initialized 4D vector.
  Vector()
      :
#ifdef USE_INTRINSICS
#ifdef USE_DOUBLE
#if defined(__AVX512F__) || defined(__AVX2__) || defined(__AVX__)
        m_value{_mm256_setzero_pd()}
#else
        m_value{_mm_setzero_pd(), _mm_setzero_pd()}
#endif  // AVX INTRINSICS
#else
        m_value{_mm_setzero_ps()}
#endif  // USE_DOUBLE
#else
        m_value{0.0, 0.0, 0.0, 0.0}
#endif  // USE_INTRINSICS
        {};

  /// Construct a 4D vector from 4 floating point values.
  explicit Vector(FLOAT x, FLOAT y, FLOAT z, FLOAT w)
      :
#ifdef USE_INTRINSICS
#ifdef USE_DOUBLE
#if defined(__AVX512F__) || defined(__AVX2__) || defined(__AVX__)
        m_value{_mm256_set_pd(w, z, y, x)}
#else
        m_value{_mm_set_pd(y, x), _mm_set_pd(w, z)}
#endif  // AVX INTRINSICS
#else
        m_value{_mm_set_ps(w, z, y, x)}
#endif  // USE_DOUBLE
#else
        m_value{x, y, z, w}
#endif  // USE_INTRINSICS
        {};

  /// Construct a 4D vector from an array of size 4.
  explicit Vector(const std::array<FLOAT, 4> &value)
      :
#ifdef USE_INTRINSICS
#ifdef USE_DOUBLE
#if defined(__AVX512F__) || defined(__AVX2__) || defined(__AVX__)
        m_value{_mm256_load_pd(value.data())}
#else
        m_value{_mm_set_pd(value[1], value[0]), _mm_set_pd(value[3], value[2])}
#endif  // AVX INTRINSICS
#else
        m_value{_mm_load_ps(value.data())}
#endif  // USE_DOUBLE
#else
        m_value{value}
#endif  // USE_INTRINSICS
        {};

  ~Vector() = default;

  // Copy
  Vector(const Vector &vec) {
    if (this != &vec) m_value = vec.m_value;
  }

  Vector &operator=(const Vector &vec) {
    if (this != &vec) m_value = vec.m_value;
    return *this;
  }

  // Getters / Setters
//  [[nodiscard]] inline FLOAT getX() const { return m_value[0]; }
//
//  [[nodiscard]] inline FLOAT getY() const { return m_value[1]; }
//
//  [[nodiscard]] inline FLOAT getZ() const { return m_value[2]; }
//
//  [[nodiscard]] inline FLOAT getW() const { return m_value[3]; }

  inline void setX(FLOAT x) {
#ifdef USE_INTRINSICS
#ifdef USE_DOUBLE
#if defined(__AVX512F__) || defined(__AVX2__) || defined(__AVX__)
    __m128d low = _mm256_castpd256_pd128(m_value);
    __m128d setX = _mm_set_sd(x);
    low = _mm_move_sd(low, setX);
    m_value = _mm256_insertf128_pd(m_value, low, 0);
#else
    __m128d setX = _mm_set_sd(x);
    m_value[0] = _mm_move_sd(m_value[0], setX);
#endif  // AVX INTRINSICS
#else
    __m128 setX = _mm_set_ss(x);
    m_value = _mm_move_ss(m_value, setX);
#endif  // USE_DOUBLE
#else
    m_value[0] = x;
#endif  // USE_INTRINSICS
  }

  inline void setY(FLOAT y) {
#ifdef USE_INTRINSICS
#ifdef USE_DOUBLE
#if defined(__AVX512F__) || defined(__AVX2__) || defined(__AVX__)
    __m128d low = _mm256_castpd256_pd128(m_value);
    low = _mm_permute_pd(low, _MM_SHUFFLE2(0, 1));
    __m128d setY = _mm_set_sd(y);
    low = _mm_move_sd(low, setY);
    low = _mm_permute_pd(low, _MM_SHUFFLE2(0, 1));
    m_value = _mm256_insertf128_pd(m_value, low, 0);
#else
    __m128d setY = _mm_set_sd(y);
    m_value[0] = _mm_shuffle_pd(m_value[0], m_value[0], _MM_SHUFFLE2(0, 1));
    m_value[0] = _mm_move_sd(m_value[0], setY);
    m_value[0] = _mm_shuffle_pd(m_value[0], m_value[0], _MM_SHUFFLE2(0, 1));
#endif  // AVX INTRINSICS
#else
    __m128 setY = _mm_set_ss(y);
    m_value = _mm_shuffle_ps(m_value, m_value, _MM_SHUFFLE(3, 2, 0, 1));
    m_value = _mm_move_sd(m_value, setY);
    m_value = _mm_shuffle_ps(m_value, m_value, _MM_SHUFFLE(3, 2, 0, 1));
#endif  // USE_DOUBLE
#else
    m_value[1] = y;
#endif  // USE_INTRINSICS
  }

  inline void setZ(FLOAT z) {
#ifdef USE_INTRINSICS
#ifdef USE_DOUBLE
#if defined(__AVX512F__) || defined(__AVX2__) || defined(__AVX__)
    __m128d high = _mm256_extractf128_pd(m_value, 1);
    __m128d setZ = _mm_set_sd(z);
    high = _mm_move_sd(high, setZ);
    m_value = _mm256_insertf128_pd(m_value, high, 1);
#else
    __m128d setZ = _mm_set_sd(z);
    m_value[1] = _mm_move_sd(m_value[1], setZ);
#endif  // AVX INTRINSICS
#else
    __m128 setZ = _mm_set_ss(z);
    m_value = _mm_shuffle_ps(m_value, m_value, _MM_SHUFFLE(3, 0, 1, 2));
    m_value = _mm_move_sd(m_value, setZ);
    m_value = _mm_shuffle_ps(m_value, m_value, _MM_SHUFFLE(3, 0, 1, 2));
#endif  // USE_DOUBLE
#else
    m_value[2] = z;
#endif  // USE_INTRINSICS
  }

  inline void setW(FLOAT w) {
#ifdef USE_INTRINSICS
#ifdef USE_DOUBLE
#if defined(__AVX512F__) || defined(__AVX2__) || defined(__AVX__)
    __m128d high = _mm256_extractf128_pd(m_value, 1);
    high = _mm_permute_pd(high, _MM_SHUFFLE2(0, 1));
    __m128d setW = _mm_set_sd(w);
    high = _mm_move_sd(high, setW);
    high = _mm_permute_pd(high, _MM_SHUFFLE2(0, 1));
    m_value = _mm256_insertf128_pd(m_value, high, 1);
#else
    __m128d setW = _mm_set_sd(w);
    m_value[1] = _mm_shuffle_pd(m_value[1], m_value[1], _MM_SHUFFLE2(0, 1));
    m_value[1] = _mm_move_sd(m_value[1], setW);
    m_value[1] = _mm_shuffle_pd(m_value[1], m_value[1], _MM_SHUFFLE2(0, 1));
#endif  // AVX INTRINSICS
#else
    __m128 setW = _mm_set_ss(w);
    m_value = _mm_shuffle_ps(m_value, m_value, _MM_SHUFFLE(0, 2, 1, 3));
    m_value = _mm_move_sd(m_value, setW);
    m_value = _mm_shuffle_ps(m_value, m_value, _MM_SHUFFLE(0, 2, 1, 3));
#endif  // USE_DOUBLE
#else
    m_value[3] = w;
#endif  // USE_INTRINSICS
  }

  inline void set(FLOAT x, FLOAT y, FLOAT z, FLOAT w) {
#ifdef USE_INTRINSICS
#ifdef USE_DOUBLE
#if defined(__AVX512F__) || defined(__AVX2__) || defined(__AVX__)
    m_value = _mm256_set_pd(w, z, y, x);
#else
    m_value[0] = _mm_set_pd(y, x);
    m_value[1] = _mm_set_pd(w, z);
#endif  // AVX INTRINSICS
#else
    m_value = _mm_set_ps(w, z, y, x);
#endif  // USE_DOUBLE
#else
    m_value[0] = x;
    m_value[1] = y;
    m_value[2] = z;
    m_value[3] = w;
#endif  // USE_INTRINSICS
  }

//  inline void set(const std::array<FLOAT, 4> &value) { m_value = value; }
//
//  FLOAT operator[](size_t index) const { return m_value[index]; }
//
//  FLOAT &operator[](size_t index) { return m_value[index]; }

  // Type-Conversions
  /// Convert to a human-readable string value.
//  explicit inline operator std::string() const {
//    return "[" + std::to_string(m_value[0]) + ", " +
//           std::to_string(m_value[1]) + ", " + std::to_string(m_value[2]) +
//           ", " + std::to_string(m_value[3]) + "]";
//  }

 private:
#ifdef USE_INTRINSICS
/**<---USE_INTRINSICS BEGIN--->**/
#if USE_DOUBLE
#if defined(__AVX512F__) || defined(__AVX2__) || defined(__AVX__)
  __m256d m_value{};
#else
  std::array<__m128d, 2> m_value{};
#endif
#else
  __m128 m_value{};
#endif
/**<---USE_INTRINSICS END--->**/
#else
  /**<---SCALAR BEGIN--->**/
  std::array<FLOAT, 4> m_value{};
  /**<---SCALAR END--->**/
#endif
};

}  // namespace cgmath::internal

#endif  // CGMATH_INTERNAL_VECTOR_HPP
