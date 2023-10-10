// Copyright 2023 Suraj Shettigar
// SPDX-License-Identifier: Apache-2.0

#ifndef CGMATH_MATRIX4X4_HPP
#define CGMATH_MATRIX4X4_HPP

#include "constants.hpp"
#include "internal/matrix.hpp"

namespace cgmath {

class Matrix4x4 {

 private:
    internal::Matrix m_value {};
};
}
#endif  // CGMATH_MATRIX4X4_HPP
