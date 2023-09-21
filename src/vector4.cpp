// Copyright 2023 Suraj Shettigar
// SPDX-License-Identifier: Apache-2.0

#include "vector4.hpp"

#include "vector2.hpp"
#include "vector3.hpp"

using namespace cgmath;

const Vector4 Vector4::ZERO{0.0, 0.0, 0.0, 0.0};
const Vector4 Vector4::ONE{1.0, 1.0, 1.0, 1.0};
const Vector4 Vector4::LEFT{-1.0, 0.0, 0.0, 0.0};
const Vector4 Vector4::RIGHT{1.0, 0.0, 0.0, 0.0};
const Vector4 Vector4::UP{0.0, 1.0, 0.0, 0.0};
const Vector4 Vector4::DOWN{0.0, -1.0, 0.0, 0.0};
const Vector4 Vector4::FORWARD{0.0, 0.0, 1.0, 0.0};
const Vector4 Vector4::BACK{0.0, 0.0, -1.0, 0.0};
const Vector4 Vector4::DIRECTION{1.0, 1.0, 1.0, 0.0};
const Vector4 Vector4::POSITION{1.0, 1.0, 1.0, 1.0};

Vector4::operator Vector2() const {
    return Vector2{m_value[0], m_value[1]};
}

Vector4::operator Vector3() const {
    return Vector3{m_value[0], m_value[1], m_value[2]};
}
