// Copyright 2023 Suraj Shettigar
// SPDX-License-Identifier: Apache-2.0

#ifndef CGMATH_INTERNAL_VECTOR_HPP
#define CGMATH_INTERNAL_VECTOR_HPP

#include "../constants.hpp"
#include "simd.hpp"

namespace cgmath::internal {

class Matrix2x2;
class Matrix3x3;
class Matrix4x4;

class Vector {
  friend class Matrix2x2;
  friend class Matrix3x3;
  friend class Matrix4x4;

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

  /// Construct a 4D vector with each dimension set to the given value
  explicit Vector(FLOAT value)
      :
#ifdef USE_INTRINSICS
#ifdef USE_DOUBLE
#if defined(__AVX512F__) || defined(__AVX2__) || defined(__AVX__)
        m_value{_mm256_set1_pd(value)}
#else
        m_value{_mm_set1_pd(value), _mm_set1_pd(value)}
#endif  // AVX INTRINSICS
#else
        m_value{_mm_set1_ps(value)}
#endif  // USE_DOUBLE
#else
        m_value{value, value, value, value}
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

  /// Construct a 3D vector from 3 floating point values.
  explicit Vector(FLOAT x, FLOAT y, FLOAT z) : Vector{x, y, z, 0.0} {};

  /// Construct a 2D vector from 2 floating point values.
  explicit Vector(FLOAT x, FLOAT y) : Vector{x, y, 0.0, 0.0} {};

  /// Construct a 4D vector from an array of size 4.
  explicit Vector(const std::array<FLOAT, 4> &value)
      :
#ifdef USE_INTRINSICS
#ifdef USE_DOUBLE
#if defined(__AVX512F__) || defined(__AVX2__) || defined(__AVX__)
        m_value{_mm256_load_pd(value.data())}
#else
        m_value{_mm_load_pd(value.data()),
                _mm_load_pd(reinterpret_cast<const FLOAT *>(
                    reinterpret_cast<const uint8_t *>(value.data()) +
                    (2 * sizeof(FLOAT))))}
#endif  // AVX INTRINSICS
#else
        m_value{_mm_load_ps(value.data())}
#endif  // USE_DOUBLE
#else
        m_value{value}
#endif  // USE_INTRINSICS
        {};

  /// Construct a 3D vector from an array of size 3.
  explicit Vector(const std::array<FLOAT, 3> &value)
      : Vector{std::array<FLOAT, 4>{value[0], value[1], value[2], 0.0}} {};

  /// Construct a 2D vector from an array of size 2.
  explicit Vector(const std::array<FLOAT, 2> &value)
      : Vector{std::array<FLOAT, 4>{value[0], value[1], 0.0, 0.0}} {};

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
  inline void get(FLOAT *data) const {
#ifdef USE_INTRINSICS
#ifdef USE_DOUBLE
#if defined(__AVX512F__) || defined(__AVX2__) || defined(__AVX__)
    _mm256_store_pd(data, m_value);
#else
    _mm_store_pd(data, m_value[0]);
    _mm_store_pd(reinterpret_cast<FLOAT *>(reinterpret_cast<uint8_t *>(data) +
                                           2 * sizeof(FLOAT)),
                 m_value[1]);
#endif  // AVX INTRINSICS
#else
    _mm_store_ps(data, m_value);
#endif  // USE_DOUBLE
#else
    std::memcpy(data, m_value.data(), m_value.size());
#endif  // USE_INTRINSICS
  }

  inline FLOAT operator[](size_t index) const {
#ifdef USE_INTRINSICS
    std::array<FLOAT, 4> val{};
    get(val.data());
    return val[index];
#else
    return m_value[index];
#endif  // USE_INTRINSICS
  }

  [[nodiscard]] inline FLOAT getX() const { return (*this)[0]; }

  [[nodiscard]] inline FLOAT getY() const { return (*this)[1]; }

  [[nodiscard]] inline FLOAT getZ() const { return (*this)[2]; }

  [[nodiscard]] inline FLOAT getW() const { return (*this)[3]; }

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
    /// SSE2 compatible method
    //    m_value = _mm_shuffle_ps(m_value, m_value, _MM_SHUFFLE(3, 2, 0, 1));
    //    m_value = _mm_move_sd(m_value, setY);
    //    m_value = _mm_shuffle_ps(m_value, m_value, _MM_SHUFFLE(3, 2, 0, 1));

    /// SSE4_1 compatible method
    m_value = _mm_insert_ps(m_value, setY, 0x10);
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

    /// SSE2 compatible method
    //    m_value = _mm_shuffle_ps(m_value, m_value, _MM_SHUFFLE(3, 0, 1, 2));
    //    m_value = _mm_move_sd(m_value, setZ);
    //    m_value = _mm_shuffle_ps(m_value, m_value, _MM_SHUFFLE(3, 0, 1, 2));

    /// SSE4_1 compatible method
    m_value = _mm_insert_ps(m_value, setZ, 0x20);
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

    /// SSE2 compatible method
    //    m_value = _mm_shuffle_ps(m_value, m_value, _MM_SHUFFLE(0, 2, 1, 3));
    //    m_value = _mm_move_sd(m_value, setW);
    //    m_value = _mm_shuffle_ps(m_value, m_value, _MM_SHUFFLE(0, 2, 1, 3));

    /// SSE4_1 compatible method
    m_value = _mm_insert_ps(m_value, setW, 0x30);
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

  inline void set(FLOAT x, FLOAT y, FLOAT z) { set(x, y, z, 0.0); }

  inline void set(FLOAT x, FLOAT y) { set(x, y, 0.0, 0.0); }

  inline void set(const std::array<FLOAT, 4> &value) {
#ifdef USE_INTRINSICS
#ifdef USE_DOUBLE
#if defined(__AVX512F__) || defined(__AVX2__) || defined(__AVX__)
    m_value = _mm256_load_pd(value.data());
#else
    m_value[0] = _mm_load_pd(value.data());
    m_value[1] = _mm_load_pd(reinterpret_cast<const FLOAT *>(
        reinterpret_cast<const uint8_t *>(value.data()) + (2 * sizeof(FLOAT))));
#endif  // AVX INTRINSICS
#else
    m_value = _mm_load_ps(value.data());
#endif  // USE_DOUBLE
#else
    m_value = value;
#endif  // USE_INTRINSICS
  }

  inline void set(const std::array<FLOAT, 3> &value) {
    set(std::array<FLOAT, 4>{value[0], value[1], value[2], 0.0});
  }

  inline void set(const std::array<FLOAT, 2> &value) {
    set(std::array<FLOAT, 4>{value[0], value[1], 0.0, 0.0});
  }
  // Type-Conversions
  /// Convert to a human-readable string value.
  explicit inline operator std::string() const {
#ifdef USE_INTRINSICS
    std::array<FLOAT, 4> val{};
    get(val.data());
    return "[" + std::to_string(val[0]) + ", " + std::to_string(val[1]) + ", " +
           std::to_string(val[2]) + ", " + std::to_string(val[3]) + "]";
#else
    return "[" + std::to_string(m_value[0]) + ", " +
           std::to_string(m_value[1]) + ", " + std::to_string(m_value[2]) +
           ", " + std::to_string(m_value[3]) + "]";
#endif  // USE_INTRINSICS
  }

  // Operators
  inline friend Vector operator+(const Vector &lhs, const Vector &rhs) {
#ifdef USE_INTRINSICS
#ifdef USE_DOUBLE
#if defined(__AVX512F__) || defined(__AVX2__) || defined(__AVX__)
    return Vector(_mm256_add_pd(lhs.m_value, rhs.m_value));
#else
    return Vector({_mm_add_pd(lhs.m_value[0], rhs.m_value[0]),
                   _mm_add_pd(lhs.m_value[1], rhs.m_value[1])});
#endif  // AVX INTRINSICS
#else
    return Vector(_mm_add_ps(lhs.m_value, rhs.m_value));
#endif  // USE_DOUBLE
#else
    return Vector{
        lhs.m_value[0] + rhs.m_value[0], lhs.m_value[1] + rhs.m_value[1],
        lhs.m_value[2] + rhs.m_value[2], lhs.m_value[3] + rhs.m_value[3]};
#endif  // USE_INTRINSICS
  }

  inline friend Vector operator-(const Vector &lhs, const Vector &rhs) {
#ifdef USE_INTRINSICS
#ifdef USE_DOUBLE
#if defined(__AVX512F__) || defined(__AVX2__) || defined(__AVX__)
    return Vector(_mm256_sub_pd(lhs.m_value, rhs.m_value));
#else
    return Vector({_mm_sub_pd(lhs.m_value[0], rhs.m_value[0]),
                   _mm_sub_pd(lhs.m_value[1], rhs.m_value[1])});
#endif  // AVX INTRINSICS
#else
    return Vector(_mm_sub_ps(lhs.m_value, rhs.m_value));
#endif  // USE_DOUBLE
#else
    return Vector{
        lhs.m_value[0] - rhs.m_value[0], lhs.m_value[1] - rhs.m_value[1],
        lhs.m_value[2] - rhs.m_value[2], lhs.m_value[3] - rhs.m_value[3]};
#endif  // USE_INTRINSICS
  }

  inline friend Vector operator*(const Vector &lhs, FLOAT rhs) {
#ifdef USE_INTRINSICS
#ifdef USE_DOUBLE
#if defined(__AVX512F__) || defined(__AVX2__) || defined(__AVX__)
    return Vector(_mm256_mul_pd(lhs.m_value, Vector(rhs).m_value));
#else
    auto right = Vector(rhs);
    return Vector({_mm_mul_pd(lhs.m_value[0], right.m_value[0]),
                   _mm_mul_pd(lhs.m_value[1], right.m_value[1])});
#endif  // AVX INTRINSICS
#else
    return Vector(_mm_mul_ps(lhs.m_value, Vector(rhs).m_value));
#endif  // USE_DOUBLE
#else
    return Vector{lhs.m_value[0] * rhs, lhs.m_value[1] * rhs,
                  lhs.m_value[2] * rhs, lhs.m_value[3] * rhs};
#endif  // USE_INTRINSICS
  }

  inline friend Vector operator*(FLOAT lhs, const Vector &rhs) {
    return rhs * lhs;
  }

  inline friend Vector operator/(const Vector &lhs, FLOAT rhs) {
#ifdef USE_INTRINSICS
#ifdef USE_DOUBLE
#if defined(__AVX512F__) || defined(__AVX2__) || defined(__AVX__)
    return Vector(_mm256_div_pd(lhs.m_value, Vector(rhs).m_value));
#else
    auto right = Vector(rhs);
    return Vector({_mm_div_pd(lhs.m_value[0], right.m_value[0]),
                   _mm_div_pd(lhs.m_value[1], right.m_value[1])});
#endif  // AVX INTRINSICS
#else
    return Vector(_mm_div_ps(lhs.m_value, Vector(rhs).m_value));
#endif  // USE_DOUBLE
#else
    return Vector{lhs.m_value[0] / rhs, lhs.m_value[1] / rhs,
                  lhs.m_value[2] / rhs, lhs.m_value[3] / rhs};
#endif  // USE_INTRINSICS
  }

  inline friend Vector operator/(FLOAT lhs, const Vector &rhs) {
#ifdef USE_INTRINSICS
#ifdef USE_DOUBLE
#if defined(__AVX512F__) || defined(__AVX2__) || defined(__AVX__)
    return Vector(_mm256_div_pd(Vector(lhs).m_value, rhs.m_value));
#else
    return Vector({_mm_div_pd(Vector(lhs).m_value[0], rhs.m_value[0]),
                   _mm_div_pd(Vector(lhs).m_value[1], rhs.m_value[1])});
#endif  // AVX INTRINSICS
#else
    return Vector(_mm_div_ps(Vector(lhs).m_value, rhs.m_value));
#endif  // USE_DOUBLE
#else
    return Vector{lhs / rhs.m_value[0], lhs / rhs.m_value[1],
                  lhs / rhs.m_value[2], lhs / rhs.m_value[3]};
#endif  // USE_INTRINSICS
  }

  inline bool operator<(const Vector &rhs) const {
#ifdef USE_INTRINSICS
#ifdef USE_DOUBLE
#if defined(__AVX512F__) || defined(__AVX2__) || defined(__AVX__)
    __m256d cmp = _mm256_cmp_pd(m_value, rhs.m_value, 0x01);
    return _mm256_movemask_pd(cmp) == 0b1111;
#else
    __m128d cmp1 = _mm_cmplt_pd(m_value[0], rhs.m_value[0]);
    __m128d cmp2 = _mm_cmplt_pd(m_value[1], rhs.m_value[1]);
    return (_mm_movemask_pd(cmp1) & _mm_movemask_pd(cmp2)) == 0b11;
#endif  // AVX INTRINSICS
#else
    __m128 cmp = _mm_cmplt_ps(m_value, rhs.m_value);
    return _mm_movemask_ps(cmp) == 0b1111;
#endif  // USE_DOUBLE
#else
    return m_value[0] < rhs.m_value[0] && m_value[1] < rhs.m_value[1] &&
           m_value[2] < rhs.m_value[2] && m_value[3] < rhs.m_value[3];
#endif  // USE_INTRINSICS
  }

  inline bool operator<=(const Vector &rhs) const {
#ifdef USE_INTRINSICS
#ifdef USE_DOUBLE
#if defined(__AVX512F__) || defined(__AVX2__) || defined(__AVX__)
    __m256d cmp = _mm256_cmp_pd(m_value, rhs.m_value, 0x02);
    return _mm256_movemask_pd(cmp) == 0b1111;
#else
    __m128d cmp1 = _mm_cmple_pd(m_value[0], rhs.m_value[0]);
    __m128d cmp2 = _mm_cmple_pd(m_value[1], rhs.m_value[1]);
    return (_mm_movemask_pd(cmp1) & _mm_movemask_pd(cmp2)) == 0b11;
#endif  // AVX INTRINSICS
#else
    __m128 cmp = _mm_cmple_ps(m_value, rhs.m_value);
    return _mm_movemask_ps(cmp) == 0b1111;
#endif  // USE_DOUBLE
#else
    return m_value[0] <= rhs.m_value[0] && m_value[1] <= rhs.m_value[1] &&
           m_value[2] <= rhs.m_value[2] && m_value[3] <= rhs.m_value[3];
#endif  // USE_INTRINSICS
  }

  /// Checks if two floating point values are equal. Because of floating point
  /// imprecision,
  /// == operator cannot be used directly.
  inline bool operator==(const Vector &rhs) const {
    // Refer: https://realtimecollisiondetection.net/blog/?p=89
    return Vector::abs(*this - rhs) <=
           EPSILON * Vector::max(Vector(1.0), Vector::max(Vector::abs(*this),
                                                          Vector::abs(rhs)));
  }

  inline bool operator!=(const Vector &rhs) const { return !(*this == rhs); }

  inline bool operator>(const Vector &rhs) const { return rhs < *this; }

  inline bool operator>=(const Vector &rhs) const { return rhs <= *this; }

  inline Vector operator-() const {
#ifdef USE_INTRINSICS
#ifdef USE_DOUBLE
#if defined(__AVX512F__) || defined(__AVX2__) || defined(__AVX__)
    return Vector(_mm256_xor_pd(m_value, _mm256_set1_pd(-0.0)));
#else
    __m128d mask = _mm_set1_pd(-0.0);
    return Vector({_mm_xor_pd(m_value[0], mask), _mm_xor_pd(m_value[1], mask)});
#endif  // AVX INTRINSICS
#else
    return Vector(_mm_xor_ps(m_value, _mm_set1_ps(-0.0f)));
#endif  // USE_DOUBLE
#else
    return Vector(-m_value[0], -m_value[1], -m_value[2], -m_value[3]);
#endif  // USE_INTRINSICS
  }

  inline Vector &operator+=(const Vector &rhs) {
#ifdef USE_INTRINSICS
    *this = *this + rhs;
#else
    m_value[0] += rhs.m_value[0];
    m_value[1] += rhs.m_value[1];
    m_value[2] += rhs.m_value[2];
    m_value[3] += rhs.m_value[3];
#endif  // USE_INTRINSICS
    return *this;
  }

  inline Vector &operator-=(const Vector &rhs) {
#ifdef USE_INTRINSICS
    *this = *this - rhs;
#else
    m_value[0] -= rhs.m_value[0];
    m_value[1] -= rhs.m_value[1];
    m_value[2] -= rhs.m_value[2];
    m_value[3] -= rhs.m_value[3];
#endif  // USE_INTRINSICS
    return *this;
  }

  inline Vector &operator*=(FLOAT rhs) {
#ifdef USE_INTRINSICS
    *this = *this * rhs;
#else
    m_value[0] *= rhs;
    m_value[1] *= rhs;
    m_value[2] *= rhs;
    m_value[3] *= rhs;
#endif  // USE_INTRINSICS
    return *this;
  }

  inline Vector &operator/=(FLOAT rhs) {
    return *this *= static_cast<float>(1.0) / rhs;
  }

  friend inline std::ostream &operator<<(std::ostream &out, const Vector &vec) {
    return out << static_cast<std::string>(vec);
  }

  // Vector specific operations
  /// Get the absolute values of each vector dimension.
  [[nodiscard]] inline Vector abs() const {
#ifdef USE_INTRINSICS
#ifdef USE_DOUBLE
#if defined(__AVX512F__) || defined(__AVX2__) || defined(__AVX__)
    __m256d mask = _mm256_castsi256_pd(
        _mm256_set1_epi64x(static_cast<int64_t>(0x8000000000000000)));
    return Vector(_mm256_andnot_pd(mask, m_value));
#else
    __m128d mask = _mm_set1_pd(-0.0);
    return Vector(
        {_mm_andnot_pd(mask, m_value[0]), _mm_andnot_pd(mask, m_value[1])});
#endif  // AVX INTRINSICS
#else
    return Vector(_mm_andnot_ps(_mm_set1_ps(-0.0f), m_value));
#endif  // USE_DOUBLE
#else
    return Vector(std::fabs(m_value[0]), std::fabs(m_value[1]),
                  std::fabs(m_value[2]), std::fabs(m_value[3]));
#endif  // USE_INTRINSICS
  }

  /// Get the squared length / magnitude of the 4D vector.
  [[nodiscard]] inline FLOAT lengthSquared() const {
    return Vector::dot(*this, *this);
  }

  /// Get the length / magnitude of the 4D vector.
  [[nodiscard]] inline FLOAT length() const {
    return std::sqrt(lengthSquared());
  }

  /// Normalize the given 4 vector. The 4D vector will have a length of 1.
  [[nodiscard]] inline Vector normalized() const { return *this / length(); }

  /// Get the absolute values of each vector dimension.
  static inline Vector abs(const Vector &vec) { return vec.abs(); }

  /// Returns a vector with the maximum values among two given vectors.
  static inline Vector max(const Vector &lhs, const Vector &rhs) {
#ifdef USE_INTRINSICS
#ifdef USE_DOUBLE
#if defined(__AVX512F__) || defined(__AVX2__) || defined(__AVX__)
    return Vector(_mm256_max_pd(lhs.m_value, rhs.m_value));
#else
    return Vector({_mm_max_pd(lhs.m_value[0], rhs.m_value[0]),
                   _mm_max_pd(lhs.m_value[1], rhs.m_value[1])});
#endif  // AVX INTRINSICS
#else
    return Vector(_mm_max_ps(lhs.m_value, rhs.m_value));
#endif  // USE_DOUBLE
#else
    return Vector(std::max(lhs.m_value[0], rhs.m_value[0]),
                  std::max(lhs.m_value[1], rhs.m_value[1]),
                  std::max(lhs.m_value[2], rhs.m_value[2]),
                  std::max(lhs.m_value[3], rhs.m_value[3]));
#endif  // USE_INTRINSICS
  }

  /**
   * Get the length of the given 4D vector.
   * @param vec 4D vector whose length / magnitude needs to be found.
   * @return The length / the magnitude.
   */
  static inline FLOAT length(const Vector &vec) { return vec.length(); }

  /**
   * Get a normalized vector from the given 4D vector.
   * @param vec 4D vector to normalize.
   * @return A normalized 4D vector whose length is 1.
   */
  static inline Vector normalize(const Vector &vec) { return vec.normalized(); }

  /**
   *  Linearly interpolates between two 4D vectors.
   * @param from 4D vector to interpolate from.
   * @param to 4D vector to interpolate to.
   * @param t Scalar value within the  range [0, 1] to control interpolation.
   * @return An interpolated 4D vector.
   */
  static inline Vector lerp(const Vector &from, const Vector &to, FLOAT t) {
    return from * (static_cast<FLOAT>(1.0) - t) + to * t;
  }

  /**
   * Find the dot product of two 4D vectors.
   * @param lhs A 3D vector.
   * @param rhs Another 3D vector.
   * @return A scalar floating point value. In case of normalized vectors, value
   * will be in the range of [-1, 1].
   */
  static inline FLOAT dot(const Vector &lhs, const Vector &rhs) {
#ifdef USE_INTRINSICS
#ifdef USE_DOUBLE
#if defined(__AVX512F__) || defined(__AVX2__) || defined(__AVX__)
    __m256d dp = _mm256_mul_pd(lhs.m_value, rhs.m_value);
    dp = _mm256_hadd_pd(dp, dp);
    __m128d high = _mm256_extractf128_pd(dp, 0x01);
    high = _mm_add_pd(_mm256_castpd256_pd128(dp), high);
    return _mm_cvtsd_f64(high);
#else
    __m128d a = _mm_mul_pd(lhs.m_value[0], rhs.m_value[0]);
    __m128d b = _mm_mul_pd(lhs.m_value[1], rhs.m_value[1]);
    b = _mm_add_pd(a, b);
    b = _mm_hadd_pd(b, b);
    return _mm_cvtsd_f64(b);
#endif  // AVX INTRINSICS
#else
    __m128 len = _mm_dp_ps(lhs.m_value, rhs.m_value, 0xFF);
    return _mm_cvtss_f32(len);
#endif  // USE_DOUBLE
#else
    return lhs.m_value[0] * rhs.m_value[0] + lhs.m_value[1] * rhs.m_value[1] +
           lhs.m_value[2] * rhs.m_value[2] + lhs.m_value[3] * rhs.m_value[3];
#endif  // USE_INTRINSICS
  }

  /**
   * Find the cross product of two 3D vectors.
   * @param lhs A 3D vector.
   * @param rhs Another 3D vector.
   * @return A 3D vector. Vector is orthogonal to LHS and RHS vector if LHS and
   * RHS are linearly independent.
   */
  static inline Vector cross(const Vector &lhs, const Vector &rhs) {
#ifdef USE_INTRINSICS
#ifdef USE_DOUBLE
#if defined(__AVX512F__) || defined(__AVX2__) || defined(__AVX__)
    __m256d l1 = _mm256_permute4x64_pd(lhs.m_value, _MM_SHUFFLE(3, 0, 2, 1));
    __m256d r1 = _mm256_permute4x64_pd(rhs.m_value, _MM_SHUFFLE(3, 1, 0, 2));
    __m256d l2 = _mm256_permute4x64_pd(lhs.m_value, _MM_SHUFFLE(3, 1, 0, 2));
    __m256d r2 = _mm256_permute4x64_pd(rhs.m_value, _MM_SHUFFLE(3, 0, 2, 1));
    l1 = _mm256_mul_pd(l1, r1);
    l2 = _mm256_mul_pd(l2, r2);
    return Vector(_mm256_sub_pd(l1, l2));
#else
    __m128d l11 = _mm_shuffle_pd(lhs.m_value[0], lhs.m_value[1], 0b01);
    __m128d l12 = _mm_shuffle_pd(lhs.m_value[0], lhs.m_value[1], 0b10);
    __m128d r11 = _mm_shuffle_pd(rhs.m_value[1], rhs.m_value[0], 0b00);
    __m128d r12 = _mm_shuffle_pd(rhs.m_value[0], rhs.m_value[1], 0b11);

    __m128d l21 = _mm_shuffle_pd(lhs.m_value[1], lhs.m_value[0], 0b00);
    __m128d l22 = _mm_shuffle_pd(lhs.m_value[0], lhs.m_value[1], 0b11);
    __m128d r21 = _mm_shuffle_pd(rhs.m_value[0], rhs.m_value[1], 0b01);
    __m128d r22 = _mm_shuffle_pd(rhs.m_value[0], rhs.m_value[1], 0b10);

    l11 = _mm_mul_pd(l11, r11);
    l12 = _mm_mul_pd(l12, r12);
    l21 = _mm_mul_pd(l21, r21);
    l22 = _mm_mul_pd(l22, r22);

    return Vector({_mm_sub_pd(l11, l21), _mm_sub_pd(l12, l22)});
#endif  // AVX INTRINSICS
#else
    __m128 l1 =
        _mm_shuffle_ps(lhs.m_value, lhs.m_value, _MM_SHUFFLE(3, 0, 2, 1));
    __m128 r1 =
        _mm_shuffle_ps(rhs.m_value, rhs.m_value, _MM_SHUFFLE(3, 1, 0, 2));
    __m128 l2 =
        _mm_shuffle_ps(lhs.m_value, lhs.m_value, _MM_SHUFFLE(3, 1, 0, 2));
    __m128 r2 =
        _mm_shuffle_ps(rhs.m_value, rhs.m_value, _MM_SHUFFLE(3, 0, 2, 1));
    l1 = _mm_mul_ps(l1, r1);
    l2 = _mm_mul_ps(l2, r2);
    return Vector(_mm_sub_ps(l1, l2));
#endif  // USE_DOUBLE
#else
#ifndef USE_DOUBLE
    // Static cast to double to ensure there are no rounding errors when two
    // vectors are extremely similar.
    auto vx = static_cast<double>(lhs.m_value[0]);
    auto vy = static_cast<double>(lhs.m_value[1]);
    auto vz = static_cast<double>(lhs.m_value[2]);
    auto wx = static_cast<double>(rhs.m_value[0]);
    auto wy = static_cast<double>(rhs.m_value[1]);
    auto wz = static_cast<double>(rhs.m_value[2]);
    return Vector(static_cast<FLOAT>(vy * wz - vz * wy),
                  static_cast<FLOAT>(vz * wx - vx * wz),
                  static_cast<FLOAT>(vx * wy - vy * wx), 0.0f);
#endif
    return Vector{
        lhs.m_value[1] * rhs.m_value[2] - lhs.m_value[2] * rhs.m_value[1],
        lhs.m_value[2] * rhs.m_value[0] - lhs.m_value[0] * rhs.m_value[2],
        lhs.m_value[0] * rhs.m_value[1] - lhs.m_value[1] * rhs.m_value[0], 0.0};
#endif  // USE_INTRINSICS
  }

 private:
#ifdef USE_INTRINSICS
/**<---USE_INTRINSICS BEGIN--->**/
#if USE_DOUBLE
#if defined(__AVX512F__) || defined(__AVX2__) || defined(__AVX__)
  __m256d m_value{};
  explicit Vector(__m256d value) : m_value{value} {};
#else
  std::array<__m128d, 2> m_value{};
  explicit Vector(const std::array<__m128d, 2> &value) : m_value{value} {};
#endif
#else
  __m128 m_value{};
  explicit Vector(__m128 value) : m_value{value} {};
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
