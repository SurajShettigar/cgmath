#include "vector3.hpp"

#include "vector2.hpp"
#include "vector4.hpp"

using namespace cgmath;

const Vector3 Vector3::ZERO{0.0, 0.0, 0.0};
const Vector3 Vector3::ONE{1.0, 1.0, 1.0};
const Vector3 Vector3::RIGHT{1.0, 0.0, 0.0};
const Vector3 Vector3::LEFT{-1.0, 0.0, 0.0};
const Vector3 Vector3::UP{0.0, 1.0, 0.0};
const Vector3 Vector3::DOWN{0.0, -1.0, 0.0};
const Vector3 Vector3::FORWARD{0.0, 0.0, 1.0};
const Vector3 Vector3::BACK{0.0, 0.0, -1.0};

Vector3::operator Vector2() const {
    return Vector2{m_value[0], m_value[1]};
}

Vector3::operator Vector4() const {
    return Vector4{m_value[0], m_value[1], m_value[2], 0.0};
}