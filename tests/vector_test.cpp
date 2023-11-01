// Copyright 2023 Suraj Shettigar
// SPDX-License-Identifier: Apache-2.0

#include <gtest/gtest.h>

#include <vector2.hpp>
#include <vector3.hpp>
#include <vector4.hpp>

#ifdef USE_DOUBLE
#define ASSERT_EQ_FLT(val1, val2) ASSERT_DOUBLE_EQ(val1, val2)
#else
#define ASSERT_EQ_FLT(val1, val2) ASSERT_FLOAT_EQ(val1, val2)
#endif

#ifdef BENCHMARK_TEST
#include <benchmark/benchmark.h>
static constexpr size_t NUM_ITERATIONS = 1000000000;
static constexpr size_t NUM_REPEATS = 5;
static constexpr benchmark::TimeUnit TIME_UNIT =
    benchmark::TimeUnit::kNanosecond;
#endif

using namespace cgmath;

class VectorTest : public ::testing::Test {
 protected:
  Vector2 vec2_def_init;
  Vector2 vec2_scalar_init;
  Vector2 vec2_arr_init;
  Vector2 vec2_set;
  Vector3 vec3_def_init;
  Vector3 vec3_scalar_init;
  Vector3 vec3_arr_init;
  Vector3 vec3_set;
  Vector4 vec4_def_init;
  Vector4 vec4_scalar_init;
  Vector4 vec4_arr_init;
  Vector4 vec4_set;

  void SetUp() override {
    vec2_def_init = Vector2{};
    vec2_scalar_init = Vector2{1.0, 0.5};
    vec2_arr_init = Vector2{std::array<FLOAT, 2>{1.0, 0.5}};
    vec2_set = Vector2{};

    vec3_def_init = Vector3{};
    vec3_scalar_init = Vector3{1.0, 0.5, 0.25};
    vec3_arr_init = Vector3{std::array<FLOAT, 3>{1.0, 0.5, 0.25}};
    vec3_set = Vector3{};

    vec4_def_init = Vector4{};
    vec4_scalar_init = Vector4{1.0, 0.5, 0.25, 0.125};
    vec4_arr_init = Vector4{std::array<FLOAT, 4>{1.0, 0.5, 0.25, 0.125}};
    vec4_set = Vector4{};
  }
};

TEST_F(VectorTest, Vecto2Getters) {
  ASSERT_EQ_FLT(vec2_def_init[0], 0.0);
  ASSERT_EQ_FLT(vec2_def_init[1], 0.0);

  ASSERT_EQ_FLT(vec2_scalar_init[0], 1.0);
  ASSERT_EQ_FLT(vec2_scalar_init[1], 0.5);

  ASSERT_EQ_FLT(vec2_arr_init[0], 1.0);
  ASSERT_EQ_FLT(vec2_arr_init[1], 0.5);

  ASSERT_EQ_FLT(vec2_def_init.getX(), 0.0);
  ASSERT_EQ_FLT(vec2_def_init.getY(), 0.0);

  ASSERT_EQ_FLT(vec2_scalar_init.getX(), 1.0);
  ASSERT_EQ_FLT(vec2_scalar_init.getY(), 0.5);

  ASSERT_EQ_FLT(vec2_arr_init.getX(), 1.0);
  ASSERT_EQ_FLT(vec2_arr_init.getY(), 0.5);
}

TEST_F(VectorTest, Vector2Setters) {
  vec2_set.set(0.25, 0.125);
  ASSERT_EQ_FLT(vec2_set[0], 0.25);
  ASSERT_EQ_FLT(vec2_set[1], 0.125);

  vec2_set.set(std::array<FLOAT, 2>{0.0625, 0.03125});
  ASSERT_EQ_FLT(vec2_set[0], 0.0625);
  ASSERT_EQ_FLT(vec2_set[1], 0.03125);

  vec2_set.setX(0.015625);
  vec2_set.setY(0.0078125);
  ASSERT_EQ_FLT(vec2_set[0], 0.015625);
  ASSERT_EQ_FLT(vec2_set[1], 0.0078125);
}

TEST_F(VectorTest, Vector2Operators) {
  ASSERT_EQ(vec2_scalar_init, vec2_arr_init);
  ASSERT_NE(vec2_scalar_init, vec2_def_init);

  Vector2 val = Vector2{-1.0, -0.5};
  ASSERT_EQ(-vec2_scalar_init, val);

  val += vec2_scalar_init;
  ASSERT_EQ(val, vec2_def_init);

  val -= vec2_scalar_init;
  ASSERT_EQ(val, -vec2_scalar_init);

  val *= 2.0;
  Vector2 newVal = Vector2{-2.0, -1.0};
  ASSERT_EQ(val, newVal);

  val /= 2.0;
  ASSERT_EQ(val, -vec2_scalar_init);

  val = vec2_scalar_init + vec2_scalar_init;
  newVal = Vector2{2.0, 1.0};
  ASSERT_EQ(val, newVal);

  val = vec2_scalar_init - vec2_scalar_init;
  newVal = Vector2{0.0, 0.0};
  ASSERT_EQ(val, newVal);

  val = vec2_scalar_init * 2.0;
  newVal = Vector2{2.0, 1.0};
  ASSERT_EQ(val, newVal);

  val = 2.0 * vec2_scalar_init;
  newVal = Vector2{2.0, 1.0};
  ASSERT_EQ(val, newVal);

  val = vec2_scalar_init / 2.0;
  newVal = Vector2{0.5, 0.25};
  ASSERT_EQ(val, newVal);

  val = 1.0 / vec2_scalar_init;
  newVal = Vector2{1.0, 2.0};
  ASSERT_EQ(val, newVal);
}

TEST_F(VectorTest, Vector2FunctionsLength) {
  ASSERT_EQ_FLT(vec2_scalar_init.lengthSquared(), 1.25);
  ASSERT_EQ_FLT(vec2_scalar_init.length(), 1.118033988749895);
  ASSERT_EQ_FLT(Vector2::length(vec2_scalar_init), 1.118033988749895);
}

TEST_F(VectorTest, Vector2FunctionsNormalize) {
  Vector2 vec1 = Vector2{0.894427190999916, 0.447213595499958};
  ASSERT_EQ(vec2_scalar_init.normalized(), vec1);
  ASSERT_EQ(Vector2::normalize(vec2_scalar_init), vec1);
}

TEST_F(VectorTest, Vector2FunctionsDotProduct) {
  Vector2 vec1 = Vector2{-0.5, 0.25};
  Vector2 vec2 = Vector2{1.0, 0.5};
  ASSERT_EQ_FLT(Vector2::dot(vec2, vec1), -0.375);
  vec1 = vec2_scalar_init.normalized();
  ASSERT_EQ_FLT(Vector2::dot(vec1, vec1), 1.0);
  ASSERT_EQ_FLT(Vector2::dot(vec1, -vec1), -1.0);
  vec1 = Vector2{1.0, 0.0};
  vec2 = Vector2{0.0, 1.0};
  ASSERT_EQ_FLT(Vector2::dot(vec1, vec2), 0.0);
}

TEST_F(VectorTest, Vector2FunctionsAngle) {
  Vector2 vec1 = Vector2{1.0, 0.0};
  Vector2 vec2 = Vector2{0.0, 1.0};
  ASSERT_EQ_FLT(Vector2::angle(vec1, vec2), 90.0);
  vec1 = Vector2{1.0, 0.0};
  vec2 = Vector2{-1.0, 0.0};
  ASSERT_EQ_FLT(Vector2::angle(vec1, vec2), 180.0);
  vec1 = Vector2{1.0, 0.0};
  vec2 = Vector2{1.0, 0.001};

#ifdef USE_DOUBLE
  ASSERT_EQ_FLT(Vector2::angle(vec1, vec2), 0.057295760416576934);
#else
  ASSERT_EQ_FLT(Vector2::angle(vec1, vec2), 0.05595291);
#endif
}

TEST_F(VectorTest, Vector2FunctionsLerp) {
  ASSERT_EQ(Vector2::lerp(vec2_def_init, vec2_scalar_init, 0.00001),
            vec2_scalar_init / 100000);
  ASSERT_EQ(Vector2::lerp(vec2_def_init, vec2_scalar_init, 0.5),
            vec2_scalar_init / 2.0);
  ASSERT_EQ(Vector2::lerp(vec2_def_init, vec2_scalar_init, 0.75),
            vec2_scalar_init / 1.333333333333333);
  ASSERT_EQ(Vector2::lerp(vec2_def_init, vec2_scalar_init, 1.0),
            vec2_scalar_init);

  Vector2 vec1 = Vector2{0.00025, 1000000.000025};
  Vector2 vec2 = Vector2{1.00333, 45084500.0000006533};
  Vector2 vec3 = Vector2{0.501790, 23042250.00001282665};
  ASSERT_EQ(Vector2::lerp(vec1, vec2, 0.5), vec3);
}

TEST_F(VectorTest, Vector2FunctionsProject) {
  Vector2 vec1 = Vector2{-4.0, 2.0};
  Vector2 vec2 = Vector2{3.0, 0.5};
  Vector2 vec3 = Vector2{-3.567567567567568, -0.594594594594595};
  ASSERT_EQ(Vector2::project(vec1, vec2), vec3);
}

#ifdef BENCHMARK_TEST
static void Vector2BenchmarkFunction(benchmark::State &state) {
  for (auto s : state) {
    Vector2 a{10.0, 5.0};
    Vector2 b{2.0, -0.54};
    a += b;
    a -= b;
    b *= 2.0;
    b /= 2.0;
    FLOAT dp = Vector2::dot(a, b);
    Vector2 c = Vector2::normalize(a);
    c = c * dp;
    Vector2 d = Vector2::project(c, b);
    benchmark::DoNotOptimize(d);
  }
}

TEST_F(VectorTest, Vector2Benchmark) {
  benchmark::internal::Benchmark *b = benchmark::RegisterBenchmark(
      "VECTOR2_BENCHMARK", Vector2BenchmarkFunction);
  b->Iterations(NUM_ITERATIONS);
  b->Repetitions(NUM_REPEATS);
  b->ReportAggregatesOnly(true);
  b->Unit(TIME_UNIT);
}
#endif

TEST_F(VectorTest, Vector3Getters) {
  ASSERT_EQ_FLT(vec3_def_init[0], 0.0);
  ASSERT_EQ_FLT(vec3_def_init[1], 0.0);
  ASSERT_EQ_FLT(vec3_def_init[2], 0.0);

  ASSERT_EQ_FLT(vec3_scalar_init[0], 1.0);
  ASSERT_EQ_FLT(vec3_scalar_init[1], 0.5);
  ASSERT_EQ_FLT(vec3_scalar_init[2], 0.25);

  ASSERT_EQ_FLT(vec3_arr_init[0], 1.0);
  ASSERT_EQ_FLT(vec3_arr_init[1], 0.5);
  ASSERT_EQ_FLT(vec3_arr_init[2], 0.25);

  ASSERT_EQ_FLT(vec3_def_init.getX(), 0.0);
  ASSERT_EQ_FLT(vec3_def_init.getY(), 0.0);
  ASSERT_EQ_FLT(vec3_def_init.getZ(), 0.0);

  ASSERT_EQ_FLT(vec3_scalar_init.getX(), 1.0);
  ASSERT_EQ_FLT(vec3_scalar_init.getY(), 0.5);
  ASSERT_EQ_FLT(vec3_scalar_init.getZ(), 0.25);

  ASSERT_EQ_FLT(vec3_arr_init.getX(), 1.0);
  ASSERT_EQ_FLT(vec3_arr_init.getY(), 0.5);
  ASSERT_EQ_FLT(vec3_arr_init.getZ(), 0.25);
}

TEST_F(VectorTest, Vector3Setters) {
  vec3_set.set(0.25, 0.125, 0.0625);
  ASSERT_EQ_FLT(vec3_set[0], 0.25);
  ASSERT_EQ_FLT(vec3_set[1], 0.125);
  ASSERT_EQ_FLT(vec3_set[2], 0.0625);

  vec3_set.set(std::array<FLOAT, 3>{0.03125, 0.015625, 0.0078125});
  ASSERT_EQ_FLT(vec3_set[0], 0.03125);
  ASSERT_EQ_FLT(vec3_set[1], 0.015625);
  ASSERT_EQ_FLT(vec3_set[2], 0.0078125);

  vec3_set.setX(0.00390625);
  vec3_set.setY(0.001953125);
  vec3_set.setZ(0.0009765625);
  ASSERT_EQ_FLT(vec3_set[0], 0.00390625);
  ASSERT_EQ_FLT(vec3_set[1], 0.001953125);
  ASSERT_EQ_FLT(vec3_set[2], 0.0009765625);
}

TEST_F(VectorTest, Vector3Operators) {
  ASSERT_EQ(vec3_scalar_init, vec3_arr_init);
  ASSERT_NE(vec3_scalar_init, vec3_def_init);

  Vector3 val = Vector3{-1.0, -0.5, -0.25};
  ASSERT_EQ(-vec3_scalar_init, val);

  val += vec3_scalar_init;
  ASSERT_EQ(val, vec3_def_init);

  val -= vec3_scalar_init;
  ASSERT_EQ(val, -vec3_scalar_init);

  val *= 2.0;
  Vector3 newVal = Vector3{-2.0, -1.0, -0.5};
  ASSERT_EQ(val, newVal);

  val /= 2.0;
  ASSERT_EQ(val, -vec3_scalar_init);

  val = vec3_scalar_init + vec3_scalar_init;
  newVal = Vector3{2.0, 1.0, 0.5};
  ASSERT_EQ(val, newVal);

  val = vec3_scalar_init - vec3_scalar_init;
  newVal = Vector3{0.0, 0.0, 0.0};
  ASSERT_EQ(val, newVal);

  val = vec3_scalar_init * 2.0;
  newVal = Vector3{2.0, 1.0, 0.5};
  ASSERT_EQ(val, newVal);

  val = 2.0 * vec3_scalar_init;
  newVal = Vector3{2.0, 1.0, 0.5};
  ASSERT_EQ(val, newVal);

  val = vec3_scalar_init / 2.0;
  newVal = Vector3{0.5, 0.25, 0.125};
  ASSERT_EQ(val, newVal);

  val = 1.0 / vec3_scalar_init;
  newVal = Vector3{1.0, 2.0, 4.0};
  ASSERT_EQ(val, newVal);
}

TEST_F(VectorTest, Vector3FunctionsLength) {
  ASSERT_EQ_FLT(vec3_scalar_init.lengthSquared(), 1.3125);
  ASSERT_EQ_FLT(vec3_scalar_init.length(), 1.1456439237389600);
  ASSERT_EQ_FLT(Vector3::length(vec3_scalar_init), 1.1456439237389600);
}

TEST_F(VectorTest, Vector3FunctionsNormalize) {
  Vector3 vec1 =
      Vector3{0.872871560943970, 0.436435780471985, 0.218217890235992};
  ASSERT_EQ(vec3_scalar_init.normalized(), vec1);
  ASSERT_EQ(Vector3::normalize(vec3_scalar_init), vec1);
}

TEST_F(VectorTest, Vector3FunctionsDotProduct) {
  Vector3 vec1 = Vector3{-0.5, 0.25, -0.75};
  Vector3 vec2 = Vector3{1.0, 0.5, 0.25};
  ASSERT_EQ_FLT(Vector3::dot(vec1, vec2), -0.5625);
  vec1 = vec3_scalar_init.normalized();
  ASSERT_EQ_FLT(Vector3::dot(vec1, vec1), 1.0);
  ASSERT_EQ_FLT(Vector3::dot(vec1, -vec1), -1.0);
  vec1 = Vector3{1.0, 0.0, 0.0};
  vec2 = Vector3{0.0, 1.0, 0.0};
  ASSERT_EQ_FLT(Vector3::dot(vec1, vec2), 0.0);
}

TEST_F(VectorTest, Vector3FunctionsAngle) {
  Vector3 vec1 = Vector3{1.0, 0.0, 0.0};
  Vector3 vec2 = Vector3{0.0, 1.0, 0.0};
  ASSERT_EQ_FLT(Vector3::angle(vec1, vec2), 90.0);
  vec1 = Vector3{1.0, 0.0, 0.0};
  vec2 = Vector3{-1.0, 0.0, 0.0};
  ASSERT_EQ_FLT(Vector3::angle(vec1, vec2), 180.0);
  vec1 = Vector3{1.0, 0.0, 0.0};
  vec2 = Vector3{1.0, 0.001, 0.0};

#ifdef USE_DOUBLE
  ASSERT_EQ_FLT(Vector3::angle(vec1, vec2), 0.057295760416576934);
#else
  ASSERT_EQ_FLT(Vector3::angle(vec1, vec2), 0.05595291);
#endif
}

TEST_F(VectorTest, Vector3FunctionsCrossProduct) {
  Vector3 vec1 = Vector3{1.0, 0.0, 0.0};
  Vector3 vec2 = Vector3{0.0, 1.0, 0.0};
  Vector3 vec3 = Vector3{0.0, 0.0, 1.0};
  ASSERT_EQ(Vector3::cross(vec1, vec2), vec3);
  vec1 = Vector3{0.1, 0.24, 0.45};
  vec2 = Vector3{0.10001, 0.24, 0.45};
  vec3 = Vector3{0.0, 0.0000045, -0.0000024};
  ASSERT_EQ(Vector3::cross(vec1, vec2), vec3);
  vec1 = Vector3{1.0, 0.0, 0.0};
  vec2 = Vector3{-1.0, 0.0, 0.0};
  vec3 = Vector3{0.0, 0.0, 0.0};
  ASSERT_EQ(Vector3::cross(vec1, vec2), vec3);
}

TEST_F(VectorTest, Vector3FunctionsLerp) {
  ASSERT_EQ(Vector3::lerp(vec3_def_init, vec3_scalar_init, 0.00001),
            vec3_scalar_init / 100000);
  ASSERT_EQ(Vector3::lerp(vec3_def_init, vec3_scalar_init, 0.5),
            vec3_scalar_init / 2.0);
  ASSERT_EQ(Vector3::lerp(vec3_def_init, vec3_scalar_init, 0.75),
            vec3_scalar_init / 1.333333333333333);
  ASSERT_EQ(Vector3::lerp(vec3_def_init, vec3_scalar_init, 1.0),
            vec3_scalar_init);

  Vector3 vec1 = Vector3{0.00025, 0.0012, 1000000.000025};
  Vector3 vec2 = Vector3{1.00333, 0.0013, 45084500.0000006533};
  Vector3 vec3 = Vector3{0.501790, 0.001250, 23042250.00001282665};
  ASSERT_EQ(Vector3::lerp(vec1, vec2, 0.5), vec3);
}

TEST_F(VectorTest, Vector3FunctionsReflect) {
  Vector3 vec1 = Vector3{-1.0, 0.00001, 0.0};
  Vector3 vec2 = Vector3{-1.0, -0.00001, 0.0};
  Vector3 vec3 = Vector3{0.0, 1.0, 0.0};
  ASSERT_EQ(Vector3::reflect(vec1, vec3), vec2);
}

TEST_F(VectorTest, Vector3FunctionsProject) {
  Vector3 vec1 = Vector3{-4.0, 2.0, 5.0};
  Vector3 vec2 = Vector3{3.0, 0.5, 1.5};
  Vector3 vec3 =
      Vector3{-0.91304347826087, -0.152173913043478, -0.456521739130435};
  ASSERT_EQ(Vector3::project(vec1, vec2), vec3);
}

TEST_F(VectorTest, Vector3FunctionsProjectOnPlane) {
  Vector3 vec1 = Vector3{-4.0, 2.0, 5.0};
  Vector3 vec2 = Vector3{0.0, 1.0, 0.0};
  Vector3 vec3 = Vector3{-4.0, 0.0, 5.0};
  ASSERT_EQ(Vector3::projectOnPlane(vec1, vec2), vec3);
}

TEST_F(VectorTest, Vector3FunctionsCoordinateFrame) {
  Vector3 vec1 = Vector3{0.0, 0.0, 1.0};
  Vector3 vec2, vec3;
  Vector3::getCoordinateFrame(vec1, &vec2, &vec3);
  ASSERT_EQ(Vector3::cross(vec2, vec3), vec1);

  vec1 = Vector3{0.000385273149999, 0.000384600159999, -0.999999851823647};
  Vector3::getCoordinateFrame(vec1, &vec2, &vec3);
  ASSERT_EQ(Vector3::cross(vec2, vec3), vec1);
}

#ifdef BENCHMARK_TEST
static void Vector3BenchmarkFunction(benchmark::State &state) {
  for (auto s : state) {
    Vector3 a{10.0, 5.0, -1.250};
    Vector3 b{2.0, -0.54, 2.25};
    a += b;
    a -= b;
    b *= 2.0;
    b /= 2.0;
    FLOAT dp = Vector3::dot(a, b);
    Vector3 c = Vector3::cross(a, b);
    c = Vector3::normalize(c);
    c = c * dp;
    Vector3 d = Vector3::project(c, b);
    benchmark::DoNotOptimize(d);
  }
}

TEST_F(VectorTest, Vector3Benchmark) {
  benchmark::internal::Benchmark *b = benchmark::RegisterBenchmark(
      "VECTOR3_BENCHMARK", Vector3BenchmarkFunction);
  b->Iterations(NUM_ITERATIONS);
  b->Repetitions(NUM_REPEATS);
  b->ReportAggregatesOnly(true);
  b->Unit(TIME_UNIT);
}
#endif

TEST_F(VectorTest, Vector4Getters) {
  ASSERT_EQ_FLT(vec4_def_init[0], 0.0);
  ASSERT_EQ_FLT(vec4_def_init[1], 0.0);
  ASSERT_EQ_FLT(vec4_def_init[2], 0.0);
  ASSERT_EQ_FLT(vec4_def_init[3], 0.0);

  ASSERT_EQ_FLT(vec4_scalar_init[0], 1.0);
  ASSERT_EQ_FLT(vec4_scalar_init[1], 0.5);
  ASSERT_EQ_FLT(vec4_scalar_init[2], 0.25);
  ASSERT_EQ_FLT(vec4_scalar_init[3], 0.125);

  ASSERT_EQ_FLT(vec4_arr_init[0], 1.0);
  ASSERT_EQ_FLT(vec4_arr_init[1], 0.5);
  ASSERT_EQ_FLT(vec4_arr_init[2], 0.25);
  ASSERT_EQ_FLT(vec4_arr_init[3], 0.125);

  ASSERT_EQ_FLT(vec4_def_init.getX(), 0.0);
  ASSERT_EQ_FLT(vec4_def_init.getY(), 0.0);
  ASSERT_EQ_FLT(vec4_def_init.getZ(), 0.0);
  ASSERT_EQ_FLT(vec4_def_init.getW(), 0.0);

  ASSERT_EQ_FLT(vec4_scalar_init.getX(), 1.0);
  ASSERT_EQ_FLT(vec4_scalar_init.getY(), 0.5);
  ASSERT_EQ_FLT(vec4_scalar_init.getZ(), 0.25);
  ASSERT_EQ_FLT(vec4_scalar_init.getW(), 0.125);

  ASSERT_EQ_FLT(vec4_arr_init.getX(), 1.0);
  ASSERT_EQ_FLT(vec4_arr_init.getY(), 0.5);
  ASSERT_EQ_FLT(vec4_arr_init.getZ(), 0.25);
  ASSERT_EQ_FLT(vec4_arr_init.getW(), 0.125);
}

TEST_F(VectorTest, Vector4Setters) {
  vec4_set.set(0.25, 0.125, 0.0625, 0.03125);
  ASSERT_EQ_FLT(vec4_set[0], 0.25);
  ASSERT_EQ_FLT(vec4_set[1], 0.125);
  ASSERT_EQ_FLT(vec4_set[2], 0.0625);
  ASSERT_EQ_FLT(vec4_set[3], 0.03125);

  vec4_set.set(
      std::array<FLOAT, 4>{0.015625, 0.0078125, 0.00390625, 0.001953125});
  ASSERT_EQ_FLT(vec4_set[0], 0.015625);
  ASSERT_EQ_FLT(vec4_set[1], 0.0078125);
  ASSERT_EQ_FLT(vec4_set[2], 0.00390625);
  ASSERT_EQ_FLT(vec4_set[3], 0.001953125);

  vec4_set.setX(0.0009765625);
  vec4_set.setY(0.00048828125);
  vec4_set.setZ(0.001953125);
  vec4_set.setW(0.000244140625);
  ASSERT_EQ_FLT(vec4_set[0], 0.0009765625);
  ASSERT_EQ_FLT(vec4_set[1], 0.00048828125);
  ASSERT_EQ_FLT(vec4_set[2], 0.001953125);
  ASSERT_EQ_FLT(vec4_set[3], 0.000244140625);
}

TEST_F(VectorTest, Vector4Operators) {
  ASSERT_EQ(vec4_scalar_init, vec4_arr_init);
  ASSERT_NE(vec4_scalar_init, vec4_def_init);

  Vector4 val = Vector4{-1.0, -0.5, -0.25, -0.125};
  ASSERT_EQ(-vec4_scalar_init, val);

  val += vec4_scalar_init;
  ASSERT_EQ(val, vec4_def_init);

  val -= vec4_scalar_init;
  ASSERT_EQ(val, -vec4_scalar_init);

  val *= 2.0;
  Vector4 newVal = Vector4{-2.0, -1.0, -0.5, -0.25};
  ASSERT_EQ(val, newVal);

  val /= 2.0;
  ASSERT_EQ(val, -vec4_scalar_init);

  val = vec4_scalar_init + vec4_scalar_init;
  newVal = Vector4{2.0, 1.0, 0.5, 0.25};
  ASSERT_EQ(val, newVal);

  val = vec4_scalar_init - vec4_scalar_init;
  newVal = Vector4{0.0, 0.0, 0.0, 0.0};
  ASSERT_EQ(val, newVal);

  val = vec4_scalar_init * 2.0;
  newVal = Vector4{2.0, 1.0, 0.5, 0.25};
  ASSERT_EQ(val, newVal);

  val = 2.0 * vec4_scalar_init;
  newVal = Vector4{2.0, 1.0, 0.5, 0.25};
  ASSERT_EQ(val, newVal);

  val = vec4_scalar_init / 2.0;
  newVal = Vector4{0.5, 0.25, 0.125, 0.0625};
  ASSERT_EQ(val, newVal);

  val = 1.0 / vec4_scalar_init;
  newVal = Vector4{1.0, 2.0, 4.0, 8.0};
  ASSERT_EQ(val, newVal);
}

TEST_F(VectorTest, Vector4FunctionsLength) {
  ASSERT_EQ_FLT(vec4_scalar_init.lengthSquared(), 1.328125);
  ASSERT_EQ_FLT(vec4_scalar_init.length(), 1.152443057161611);
  ASSERT_EQ_FLT(Vector4::length(vec4_scalar_init), 1.152443057161611);
}

TEST_F(VectorTest, Vector4FunctionsNormalize) {
  Vector4 vec1 = Vector4{0.867721831274625, 0.433860915637312,
                         0.216930457818656, 0.108465228909328};
  ASSERT_EQ(vec4_scalar_init.normalized(), vec1);
  ASSERT_EQ(Vector4::normalize(vec4_scalar_init), vec1);
}

TEST_F(VectorTest, Vector4FunctionsLerp) {
  ASSERT_EQ(Vector4::lerp(vec4_def_init, vec4_scalar_init, 0.00001),
            vec4_scalar_init / 100000);
  ASSERT_EQ(Vector4::lerp(vec4_def_init, vec4_scalar_init, 0.5),
            vec4_scalar_init / 2.0);
  ASSERT_EQ(Vector4::lerp(vec4_def_init, vec4_scalar_init, 0.75),
            vec4_scalar_init / 1.333333333333333);
  ASSERT_EQ(Vector4::lerp(vec4_def_init, vec4_scalar_init, 1.0),
            vec4_scalar_init);

  Vector4 vec1 = Vector4{0.00025, 0.0012, 1000000.000025, -10.0024232};
  Vector4 vec2 =
      Vector4{1.00333, 0.0013, 45084500.0000006533, -1.333333333333333};
  Vector4 vec3 =
      Vector4{0.501790, 0.001250, 23042250.00001282665, -5.667878266666667};
  ASSERT_EQ(Vector4::lerp(vec1, vec2, 0.5), vec3);
}

#ifdef BENCHMARK_TEST
static void Vector4BenchmarkFunction(benchmark::State &state) {
  for (auto s : state) {
    Vector4 a{10.0, 5.0, -1.25, -3.0};
    Vector4 b{2.0, -0.54, 2.25, 1.34432};
    a += b;
    a -= b;
    b *= 2.0;
    b /= 2.0;
    benchmark::DoNotOptimize(a);
    benchmark::DoNotOptimize(b);
    benchmark::ClobberMemory();
  }
}

TEST_F(VectorTest, Vector4Benchmark) {
  benchmark::internal::Benchmark *b = benchmark::RegisterBenchmark(
      "VECTOR4_BENCHMARK", Vector4BenchmarkFunction);
  b->Iterations(NUM_ITERATIONS);
  b->Repetitions(NUM_REPEATS);
  b->ReportAggregatesOnly(true);
  b->Unit(TIME_UNIT);
}
#endif

TEST_F(VectorTest, VectorConversions) {
  vec2_set = static_cast<Vector2>(vec4_scalar_init);
  ASSERT_EQ_FLT(vec2_set[0], vec4_scalar_init[0]);
  ASSERT_EQ_FLT(vec2_set[1], vec4_scalar_init[1]);

  std::string val = static_cast<std::string>(vec2_def_init);
  ASSERT_EQ(val, "[0.000000, 0.000000, 0.000000, 0.000000]");

  vec3_set = static_cast<Vector3>(vec2_scalar_init);
  ASSERT_EQ_FLT(vec3_set[0], vec2_scalar_init[0]);
  ASSERT_EQ_FLT(vec3_set[1], vec2_scalar_init[1]);
  ASSERT_EQ_FLT(vec3_set[2], 0.0);

  vec3_set = static_cast<Vector3>(vec4_scalar_init);
  ASSERT_EQ_FLT(vec3_set[0], vec4_scalar_init[0]);
  ASSERT_EQ_FLT(vec3_set[1], vec4_scalar_init[1]);
  ASSERT_EQ_FLT(vec3_set[2], vec4_scalar_init[2]);

  val = static_cast<std::string>(vec3_def_init);
  ASSERT_EQ(val, "[0.000000, 0.000000, 0.000000, 0.000000]");

  vec4_set = static_cast<Vector4>(vec2_scalar_init);
  ASSERT_EQ_FLT(vec4_set[0], vec2_scalar_init[0]);
  ASSERT_EQ_FLT(vec4_set[1], vec2_scalar_init[1]);
  ASSERT_EQ_FLT(vec4_set[2], 0.0);
  ASSERT_EQ_FLT(vec4_set[3], 0.0);

  vec4_set = static_cast<Vector4>(vec3_scalar_init);
  ASSERT_EQ_FLT(vec4_set[0], vec3_scalar_init[0]);
  ASSERT_EQ_FLT(vec4_set[1], vec3_scalar_init[1]);
  ASSERT_EQ_FLT(vec4_set[2], vec3_scalar_init[2]);
  ASSERT_EQ_FLT(vec4_set[3], 0.0);

  val = static_cast<std::string>(vec4_def_init);
  ASSERT_EQ(val, "[0.000000, 0.000000, 0.000000, 0.000000]");
}

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  int testResult = RUN_ALL_TESTS();
#ifdef BENCHMARK_TEST
  if (testResult == 0) {
    benchmark::Initialize(&argc, argv);

    std::string archValue{};
#ifdef USE_DOUBE
    archValue = "DOUBLE_";
#else
    archValue = "FLOAT_";
#endif
#ifdef USE_INTRINSICS
#if defined(__AVX512F__) || defined(__AVX2__) || defined(__AVX__)
    archValue += "AVX";
#else
    archValue += "SSE";
#endif  // AVX INTRINSICS
#else
    archValue += "NONE";
#endif  // USE_INTRINSICS

    benchmark::AddCustomContext("SIMD Architecture", archValue);
    benchmark::RunSpecifiedBenchmarks();
    benchmark::Shutdown();
  }
#endif
  return testResult;
}