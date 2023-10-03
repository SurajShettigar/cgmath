// Copyright 2023 Suraj Shettigar
// SPDX-License-Identifier: Apache-2.0

#ifndef CGMATH_INTERNAL_SIMD_HPP
#define CGMATH_INTERNAL_SIMD_HPP

#ifdef USE_INTRINSICS
namespace cgmath::internal {

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
#elif defined(__AVX2__) || defined(__AVX__)
// 32 byte vectorization
#include <immintrin.h>
static constexpr SIMDMicroArchitectureLevel SIMD_LEVEL =
    SIMDMicroArchitectureLevel::VERSION_3;
#else
// 16 byte vectorization
#include <emmintrin.h>
static constexpr SIMDMicroArchitectureLevel SIMD_LEVEL =
    SIMDMicroArchitectureLevel::BASE;
#endif

}  // namespace cgmath::internal
#endif  // USE_INTRINSICS

#endif  // CGMATH_INTERNAL_SIMD_HPP
