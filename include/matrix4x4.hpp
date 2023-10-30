// Copyright 2023 Suraj Shettigar
// SPDX-License-Identifier: Apache-2.0

#ifndef CGMATH_MATRIX4X4_HPP
#define CGMATH_MATRIX4X4_HPP

#include "constants.hpp"
#include "internal/matrix4x4.hpp"

namespace cgmath {

class Matrix4x4 {
 private:
  internal::Matrix4x4 m_value{};
};
}  // namespace cgmath
#endif  // CGMATH_MATRIX4X4_HPP
