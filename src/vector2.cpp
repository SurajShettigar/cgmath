#include "vector2.hpp"

#include "vector3.hpp"
#include "vector4.hpp"

using namespace cgmath;

const Vector2 Vector2::ZERO{0.0, 0.0};
const Vector2 Vector2::ONE{1.0, 1.0};
const Vector2 Vector2::LEFT{-1.0, 0.0};
const Vector2 Vector2::RIGHT{1.0, 0.0};
const Vector2 Vector2::UP{0.0, 1.0};
const Vector2 Vector2::DOWN{0.0, -1.0};

Vector2::operator Vector3() const {
  return Vector3{m_value[0], m_value[1], 0.0};
}

Vector2::operator Vector4() const {
  return Vector4{m_value[0], m_value[1], 0.0, 0.0};
}