// Copyright 2023 Suraj Shettigar
// SPDX-License-Identifier: Apache-2.0

#include <gtest/gtest.h>

#include <matrix2x2.hpp>
#include <matrix3x3.hpp>
#include <matrix4x4.hpp>

using namespace cgmath;

#ifdef USE_DOUBLE
#define ASSERT_EQ_FLT(val1, val2) ASSERT_DOUBLE_EQ(val1, val2)
#else
#define ASSERT_EQ_FLT(val1, val2) ASSERT_FLOAT_EQ(val1, val2)
#endif

class MatrixTest : public ::testing::Test {
 protected:
  Matrix2x2 mat2_def_init;
  Matrix2x2 mat2_scalar_init;
  Matrix2x2 mat2_arr_init;
  Matrix2x2 mat2_vec_init;

  Matrix3x3 mat3_def_init;
  Matrix3x3 mat3_scalar_init;
  Matrix3x3 mat3_arr_init;
  Matrix3x3 mat3_vec_init;

  Matrix4x4 mat4_def_init;
  Matrix4x4 mat4_scalar_init;
  Matrix4x4 mat4_arr_init;
  Matrix4x4 mat4_vec_init;

  void SetUp() override {
    mat2_def_init = Matrix2x2{};
    mat2_scalar_init = Matrix2x2{1.0, 0.5, -0.25, -0.125};
    std::array<FLOAT, 2> a{-1.0, 0.5};
    std::array<FLOAT, 2> b{0.25, -0.125};
    mat2_arr_init = Matrix2x2{a, b};
    mat2_vec_init =
        Matrix2x2{Vector2{10000.12345, -10.825}, Vector2{12.565, 2540.9832235}};

    mat3_def_init = Matrix3x3{};
    mat3_scalar_init =
        Matrix3x3{1.0, 0.5, -0.25, -0.125, 0.8, 1.5, 0.45, 0.35, -2.0};
    std::array<FLOAT, 3> x{-1.0, 0.5, 0.25};
    std::array<FLOAT, 3> y{0.25, -0.125, -0.65};
    std::array<FLOAT, 3> z{1.25, -2.125, 3.65};
    mat3_arr_init = Matrix3x3{x, y, z};
    mat3_vec_init = Matrix3x3{Vector3{10000.12345, -10.825, 2512.483},
                              Vector3{12.565, 2540.9832235, 5543.12},
                              Vector3{-3221.87, 1233.9832235, 4322.12}};

    mat4_def_init = Matrix4x4{};
    mat4_scalar_init =
        Matrix4x4{1.0,  0.5,  -0.25, -0.125, 0.8,    1.5, 0.45, 0.35,
                  -2.0, 1.45, 3.0,   0.5,    -0.155, 3.1, 0.1,  1.0};
    std::array<FLOAT, 4> x4{-1.0, 0.5, 0.25, 0.0};
    std::array<FLOAT, 4> y4{0.25, -0.125, -0.65, 0.0};
    std::array<FLOAT, 4> z4{1.25, -2.125, 3.65, 0.0};
    std::array<FLOAT, 4> t4{0.25, 2.0, -1.5, 1.0};
    mat4_arr_init = Matrix4x4{x4, y4, z4, t4};
    mat4_vec_init = Matrix4x4{Vector4{10000.12345, -10.825, 2512.483, 0.0},
                              Vector4{12.565, 2540.9832235, 5543.12, 0.0},
                              Vector4{-3221.87, 1233.9832235, 4322.12, 0.0},
                              Vector4{-124.12, 22.4728, 5.1138287, 1.0}};
  }
};

TEST_F(MatrixTest, Matrix2x2Getters) {
  Vector2 vec1{1.0, 0.0};
  Vector2 vec2{0.0, 1.0};
  ASSERT_EQ(mat2_def_init[0], vec1);
  ASSERT_EQ(mat2_def_init[1], vec2);

  ASSERT_EQ(mat2_def_init.getRow(0), vec1);
  ASSERT_EQ(mat2_def_init.getRow(1), vec2);
  ASSERT_EQ(mat2_def_init.getColumn(0), vec1);
  ASSERT_EQ(mat2_def_init.getColumn(1), vec2);
  ASSERT_EQ(mat2_def_init.getXAxis(), vec1);
  ASSERT_EQ(mat2_def_init.getYAxis(), vec2);

  vec1 = Vector2{1.0, -0.25};
  vec2 = Vector2{0.5, -0.125};
  ASSERT_EQ(mat2_scalar_init[0], vec1);
  ASSERT_EQ(mat2_scalar_init[1], vec2);
  ASSERT_EQ(mat2_scalar_init.getRow(0), vec1);
  ASSERT_EQ(mat2_scalar_init.getRow(1), vec2);
  vec1 = Vector2{1.0, 0.5};
  vec2 = Vector2{-0.25, -0.125};
  ASSERT_EQ(mat2_scalar_init.getColumn(0), vec1);
  ASSERT_EQ(mat2_scalar_init.getColumn(1), vec2);
  ASSERT_EQ(mat2_scalar_init.getXAxis(), vec1);
  ASSERT_EQ(mat2_scalar_init.getYAxis(), vec2);

  vec1 = Vector2{-1.0, 0.25};
  vec2 = Vector2{0.5, -0.125};
  ASSERT_EQ(mat2_arr_init[0], vec1);
  ASSERT_EQ(mat2_arr_init[1], vec2);
  ASSERT_EQ(mat2_arr_init.getRow(0), vec1);
  ASSERT_EQ(mat2_arr_init.getRow(1), vec2);
  vec1 = Vector2{-1.0, 0.5};
  vec2 = Vector2{0.25, -0.125};
  ASSERT_EQ(mat2_arr_init.getColumn(0), vec1);
  ASSERT_EQ(mat2_arr_init.getColumn(1), vec2);
  ASSERT_EQ(mat2_arr_init.getXAxis(), vec1);
  ASSERT_EQ(mat2_arr_init.getYAxis(), vec2);

  vec1 = Vector2{10000.12345, 12.565};
  vec2 = Vector2{-10.825, 2540.9832235};
  ASSERT_EQ(mat2_vec_init[0], vec1);
  ASSERT_EQ(mat2_vec_init[1], vec2);
  ASSERT_EQ(mat2_vec_init.getRow(0), vec1);
  ASSERT_EQ(mat2_vec_init.getRow(1), vec2);
  vec1 = Vector2{10000.12345, -10.825};
  vec2 = Vector2{12.565, 2540.9832235};
  ASSERT_EQ(mat2_vec_init.getColumn(0), vec1);
  ASSERT_EQ(mat2_vec_init.getColumn(1), vec2);
  ASSERT_EQ(mat2_vec_init.getXAxis(), vec1);
  ASSERT_EQ(mat2_vec_init.getYAxis(), vec2);
}

TEST_F(MatrixTest, Matrix2x2Operators) {
  ASSERT_EQ(mat2_def_init, mat2_def_init);
  ASSERT_EQ(mat2_scalar_init, mat2_scalar_init);

  ASSERT_NE(mat2_def_init, mat2_scalar_init);
  ASSERT_NE(mat2_scalar_init, mat2_arr_init);

  Matrix2x2 mat1 = mat2_scalar_init * mat2_def_init;
  ASSERT_EQ(mat2_scalar_init, mat1);
  mat1 = mat2_arr_init * mat2_scalar_init;
  Matrix2x2 mat2{-0.875, 0.4375, 0.21875, -0.109375};
  ASSERT_EQ(mat1, mat2);
  mat1 = mat2_vec_init * mat2_arr_init;
  mat2 =
      Matrix2x2{-9993.840950, 1281.316611750, 2498.46023750, -320.32915293750};
  ASSERT_EQ(mat1, mat2);

  Vector2 vec1{10.25, -0.025530};
  Vector2 vec2 = mat2_def_init * vec1;
  ASSERT_EQ(vec1, vec2);

  vec2 = mat2_scalar_init * vec1;
  Vector2 vec3 = Vector2{10.2563825, 5.128191250};
  ASSERT_EQ(vec2, vec3);

  mat1 = mat2_def_init * 2.0;
  mat2 = Matrix2x2{2.0, 0.0, 0.0, 2.0};
  ASSERT_EQ(mat1, mat2);

  mat1 = 10.0 * mat2_scalar_init;
  mat2 = Matrix2x2{10.0, 5.0, -2.5, -1.25};
  ASSERT_EQ(mat1, mat2);

  std::string matStr = "[1.000000, 0.000000\n0.000000, 1.000000]";
  ASSERT_EQ(static_cast<std::string>(mat2_def_init), matStr);
}

TEST_F(MatrixTest, Matrix2x2Setters) {
  Matrix2x2 mat1 = mat2_def_init;
  mat1.setXAxis(2.0, 4.0);
  mat1.setYAxis(-0.5, 10.4);
  Matrix2x2 mat2{2.0, 4.0, -0.5, 10.4};
  ASSERT_EQ(mat1, mat2);

  std::array<cgmath::FLOAT, 2> x_axis{7.200045, 0.00019929};
  std::array<cgmath::FLOAT, 2> y_axis{-0.1, -100.25};
  mat1.setXAxis(x_axis);
  mat1.setYAxis(y_axis);
  mat2 = Matrix2x2{x_axis, y_axis};
  ASSERT_EQ(mat1, mat2);

  Vector2 x_vec{7.200045, 0.00019929};
  Vector2 y_vec{-0.1, -100.25};
  mat1.setXAxis(x_vec);
  mat1.setYAxis(y_vec);
  mat2 = Matrix2x2{Vector2{7.200045, 0.00019929}, Vector2{-0.1, -100.25}};
  ASSERT_EQ(mat1, mat2);
}

TEST_F(MatrixTest, Matrix2x2FunctionsDeterminant) {
  ASSERT_EQ_FLT(mat2_def_init.determinant(), 1.0);
  ASSERT_EQ_FLT(mat2_scalar_init.determinant(), 0.0);
  ASSERT_EQ_FLT(mat2_arr_init.determinant(), 0.0);
  ASSERT_EQ_FLT(mat2_vec_init.determinant(), 25410281.93550394);

  Matrix2x2 mat{4.0, -0.5, 2.0, 3.0};
  ASSERT_EQ_FLT(mat.determinant(), 13.0);
}

TEST_F(MatrixTest, Matrix2x2FunctionsTranspose) {
  ASSERT_EQ(Matrix2x2::transpose(mat2_def_init), mat2_def_init);
  Matrix2x2 mat{1.0, -0.25, 0.5, -0.125};
  ASSERT_EQ(Matrix2x2::transpose(mat2_scalar_init), mat);
  mat = Matrix2x2{-1.0, 0.25, 0.5, -0.125};
  ASSERT_EQ(Matrix2x2::transpose(mat2_arr_init), mat);
  mat = Matrix2x2{10000.12345, 12.565, -10.825, 2540.9832235};
  ASSERT_EQ(Matrix2x2::transpose(mat2_vec_init), mat);
}

TEST_F(MatrixTest, Matrix2x2FunctionsInverse) {
  ASSERT_EQ(Matrix2x2::inverse(mat2_def_init), mat2_def_init);
  ASSERT_EQ(Matrix2x2::inverse(mat2_scalar_init), mat2_scalar_init);
  ASSERT_EQ(Matrix2x2::inverse(mat2_arr_init), mat2_arr_init);
  Matrix2x2 mat{9.999823024197417e-5, 4.260086537991149e-7,
                -4.944848715922289e-7, 3.935463398392111e-4};
  ASSERT_EQ(Matrix2x2::inverse(mat2_vec_init), mat);

  mat = Matrix2x2{4.0, -0.5, 2.0, 3.0};
  Matrix2x2 mat2{0.230769230769231, 0.038461538461539, -0.153846153846154,
                 0.307692307692308};
  ASSERT_EQ(Matrix2x2::inverse(mat), mat2);
}

TEST_F(MatrixTest, Matrix3x3Getters) {
  Vector3 vec1{1.0, 0.0, 0.0};
  Vector3 vec2{0.0, 1.0, 0.0};
  Vector3 vec3{0.0, 0.0, 1.0};
  ASSERT_EQ(mat3_def_init[0], vec1);
  ASSERT_EQ(mat3_def_init[1], vec2);
  ASSERT_EQ(mat3_def_init[2], vec3);
  ;
  ASSERT_EQ(mat3_def_init.getRow(0), vec1);
  ASSERT_EQ(mat3_def_init.getRow(1), vec2);
  ASSERT_EQ(mat3_def_init.getRow(2), vec3);
  ASSERT_EQ(mat3_def_init.getColumn(0), vec1);
  ASSERT_EQ(mat3_def_init.getColumn(1), vec2);
  ASSERT_EQ(mat3_def_init.getColumn(2), vec3);
  ASSERT_EQ(mat3_def_init.getXAxis(), vec1);
  ASSERT_EQ(mat3_def_init.getYAxis(), vec2);
  ASSERT_EQ(mat3_def_init.getZAxis(), vec3);

  vec1 = Vector3{1.0, -0.125, 0.45};
  vec2 = Vector3{0.5, 0.8, 0.35};
  vec3 = Vector3{-0.25, 1.5, -2.0};
  ASSERT_EQ(mat3_scalar_init[0], vec1);
  ASSERT_EQ(mat3_scalar_init[1], vec2);
  ASSERT_EQ(mat3_scalar_init[2], vec3);
  ASSERT_EQ(mat3_scalar_init.getRow(0), vec1);
  ASSERT_EQ(mat3_scalar_init.getRow(1), vec2);
  ASSERT_EQ(mat3_scalar_init.getRow(2), vec3);
  vec1 = Vector3{1.0, 0.5, -0.25};
  vec2 = Vector3{-0.125, 0.8, 1.5};
  vec3 = Vector3{0.45, 0.35, -2.0};
  ASSERT_EQ(mat3_scalar_init.getColumn(0), vec1);
  ASSERT_EQ(mat3_scalar_init.getColumn(1), vec2);
  ASSERT_EQ(mat3_scalar_init.getColumn(2), vec3);
  ASSERT_EQ(mat3_scalar_init.getXAxis(), vec1);
  ASSERT_EQ(mat3_scalar_init.getYAxis(), vec2);
  ASSERT_EQ(mat3_scalar_init.getZAxis(), vec3);

  vec1 = Vector3{-1.0, 0.25, 1.25};
  vec2 = Vector3{0.5, -0.125, -2.125};
  vec3 = Vector3{0.25, -0.65, 3.65};
  ASSERT_EQ(mat3_arr_init[0], vec1);
  ASSERT_EQ(mat3_arr_init[1], vec2);
  ASSERT_EQ(mat3_arr_init[2], vec3);
  ASSERT_EQ(mat3_arr_init.getRow(0), vec1);
  ASSERT_EQ(mat3_arr_init.getRow(1), vec2);
  ASSERT_EQ(mat3_arr_init.getRow(2), vec3);
  vec1 = Vector3{-1.0, 0.5, 0.25};
  vec2 = Vector3{0.25, -0.125, -0.65};
  vec3 = Vector3{1.25, -2.125, 3.65};
  ASSERT_EQ(mat3_arr_init.getColumn(0), vec1);
  ASSERT_EQ(mat3_arr_init.getColumn(1), vec2);
  ASSERT_EQ(mat3_arr_init.getColumn(2), vec3);
  ASSERT_EQ(mat3_arr_init.getXAxis(), vec1);
  ASSERT_EQ(mat3_arr_init.getYAxis(), vec2);
  ASSERT_EQ(mat3_arr_init.getZAxis(), vec3);

  vec1 = Vector3{10000.12345, 12.565, -3221.87};
  vec2 = Vector3{-10.825, 2540.9832235, 1233.9832235};
  vec3 = Vector3{2512.483, 5543.12, 4322.12};
  ASSERT_EQ(mat3_vec_init[0], vec1);
  ASSERT_EQ(mat3_vec_init[1], vec2);
  ASSERT_EQ(mat3_vec_init[2], vec3);
  ASSERT_EQ(mat3_vec_init.getRow(0), vec1);
  ASSERT_EQ(mat3_vec_init.getRow(1), vec2);
  ASSERT_EQ(mat3_vec_init.getRow(2), vec3);
  vec1 = Vector3{10000.12345, -10.825, 2512.483};
  vec2 = Vector3{12.565, 2540.9832235, 5543.12};
  vec3 = Vector3{-3221.87, 1233.9832235, 4322.12};
  ASSERT_EQ(mat3_vec_init.getColumn(0), vec1);
  ASSERT_EQ(mat3_vec_init.getColumn(1), vec2);
  ASSERT_EQ(mat3_vec_init.getColumn(2), vec3);
  ASSERT_EQ(mat3_vec_init.getXAxis(), vec1);
  ASSERT_EQ(mat3_vec_init.getYAxis(), vec2);
  ASSERT_EQ(mat3_vec_init.getZAxis(), vec3);
}

TEST_F(MatrixTest, Matrix3x3Operators) {
  ASSERT_EQ(mat3_def_init, mat3_def_init);
  ASSERT_EQ(mat3_scalar_init, mat3_scalar_init);

  ASSERT_NE(mat3_def_init, mat3_scalar_init);
  ASSERT_NE(mat3_scalar_init, mat3_arr_init);

  Matrix3x3 mat1 = mat3_scalar_init * mat3_def_init;
  ASSERT_EQ(mat3_scalar_init, mat1);
  mat1 = mat3_arr_init * mat3_scalar_init;
  Matrix3x3 mat2{-1.1875, 0.96875, -0.9875, 2.2,   -3.35,
                 4.92375, -2.8625, 4.43125, -7.415};
  ASSERT_EQ(mat1, mat2);
  mat1 = mat3_vec_init * mat3_arr_init;
  mat2 = Matrix3x3{-10799.308450, 1589.812417625,   1339.607,
                   4592.6757375,  -1122.4182482125, -2874.14725,
                   713.6281875,   -909.0818341625,  7137.21175};
  ASSERT_EQ(mat1, mat2);

  Vector3 vec1{10.25, -0.025530, 3.456};
  Vector3 vec2 = mat3_def_init * vec1;
  ASSERT_EQ(vec1, vec2);

  vec2 = mat3_scalar_init * vec1;
  Vector3 vec3 = Vector3{11.80839125, 6.314176, -9.512795};
  ASSERT_EQ(vec2, vec3);

  mat1 = mat3_def_init * 2.0;
  mat2 = Matrix3x3{2.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 2.0};
  ASSERT_EQ(mat1, mat2);

  mat1 = 10.0 * mat3_scalar_init;
  mat2 = Matrix3x3{10.0, 5.0, -2.5, -1.25, 8.0, 15.0, 4.5, 3.5, -20.0};
  ASSERT_EQ(mat1, mat2);

  std::string matStr =
      "[1.000000, 0.000000, 0.000000\n0.000000, 1.000000, 0.000000\n0.000000, "
      "0.000000, 1.000000]";
  ASSERT_EQ(static_cast<std::string>(mat3_def_init), matStr);
}

TEST_F(MatrixTest, Matrix3x3Setters) {
  Matrix3x3 mat1 = mat3_def_init;
  mat1.setXAxis(2.0, 4.0, 10.0);
  mat1.setYAxis(-0.5, 10.4, -7.8);
  mat1.setZAxis(2.5, 0.789, 1.1);
  Matrix3x3 mat2{2.0, 4.0, 10.0, -0.5, 10.4, -7.8, 2.5, 0.789, 1.1};
  ASSERT_EQ(mat1, mat2);

  std::array<cgmath::FLOAT, 3> x_axis{7.200045, 0.00019929, 0.1926345623};
  std::array<cgmath::FLOAT, 3> y_axis{-0.1, -100.25, 123.456};
  std::array<cgmath::FLOAT, 3> z_axis{2.1, 10.25, 13.46};
  mat1.setXAxis(x_axis);
  mat1.setYAxis(y_axis);
  mat1.setZAxis(z_axis);
  mat2 = Matrix3x3{x_axis, y_axis, z_axis};
  ASSERT_EQ(mat1, mat2);

  Vector3 x_vec{7.200045, 0.00019929, 0.1926345623};
  Vector3 y_vec{-0.1, -100.25, 123.456};
  Vector3 z_vec{2.1, 10.25, 13.46};
  mat1.setXAxis(x_vec);
  mat1.setYAxis(y_vec);
  mat1.setZAxis(z_vec);
  mat2 = Matrix3x3{x_vec, y_vec, z_vec};
  ASSERT_EQ(mat1, mat2);
}

TEST_F(MatrixTest, Matrix3x3FunctionsDeterminant) {
  ASSERT_EQ_FLT(mat3_def_init.determinant(), 1.0);
  ASSERT_EQ_FLT(mat3_scalar_init.determinant(), -1.8115625);
  ASSERT_EQ_FLT(mat3_arr_init.determinant(), 0.88125);
  ASSERT_EQ_FLT(mat3_vec_init.determinant(), 62225543347.1125499578925);
}

TEST_F(MatrixTest, Matrix3x3FunctionsTranspose) {
  ASSERT_EQ(Matrix3x3::transpose(mat3_def_init), mat3_def_init);
  Matrix3x3 mat{1.0, -0.125, 0.45, 0.5, 0.8, 0.35, -0.25, 1.5, -2.0};
  ASSERT_EQ(Matrix3x3::transpose(mat3_scalar_init), mat);
  mat = Matrix3x3{-1.0, 0.25, 1.25, 0.5, -0.125, -2.125, 0.25, -0.65, 3.65};
  ASSERT_EQ(Matrix3x3::transpose(mat3_arr_init), mat);
  mat = Matrix3x3{10000.12345,  12.565,   -3221.87, -10.825, 2540.9832235,
                  1233.9832235, 2512.483, 5543.12,  4322.12};
  ASSERT_EQ(Matrix3x3::transpose(mat3_vec_init), mat);
}

TEST_F(MatrixTest, Matrix3x3FunctionsInverse) {
  ASSERT_EQ(Matrix3x3::inverse(mat3_def_init), mat3_def_init);
  Matrix3x3 mat{1.173020527859238,   -0.5037088149042608, -0.5244091771606003,
                -0.2346041055718475, 1.041918233569087,   0.8107641883732965,
                0.2228739002932551,  0.0690012075211316,  -0.4761083318958082};
  ASSERT_EQ(Matrix3x3::inverse(mat3_scalar_init), mat);
  mat = Matrix3x3{-2.085106382978723, -2.673758865248227, -0.3333333333333333,
                  -1.957446808510638, -4.49645390070922,  -0.6666666666666667,
                  -0.425531914893617, -1.702127659574468, 0.0};
  ASSERT_EQ(Matrix3x3::inverse(mat3_arr_init), mat);
  mat = Matrix3x3{
      0.000066569403837897021, 0.000050576477938862198, -0.00010356167386729587,
      -0.00028788048297582672, 0.00082468749019457542,  -0.00089031468379841395,
      0.00013181439930137358,  -0.00019775010655732485, 0.00040835773492178042};
  ASSERT_EQ(Matrix3x3::inverse(mat3_vec_init), mat);
}

TEST_F(MatrixTest, Matrix3x3FunctionsScale) {
  Matrix3x3 mat = Matrix3x3::scale(Vector2{0.5, 0.5});
  Matrix3x3 res{0.5, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 1.0};
  ASSERT_EQ(mat, res);
}

TEST_F(MatrixTest, Matrix3x3FunctionsRotation) {
  Matrix3x3 mat = Matrix3x3::rotation(45.0);
  Matrix3x3 res{0.70710678118654752,
                0.70710678118654752,
                0.0,
                -0.70710678118654752,
                0.70710678118654752,
                0.0,
                0.0,
                0.0,
                1.0};
  ASSERT_EQ(mat, res);
  mat = Matrix3x3::rotation(0.00001);
  res = Matrix3x3{0.9999999999999848,
                  0.0000001745329252,
                  0.0,
                  -0.0000001745329252,
                  0.9999999999999848,
                  0.0,
                  0.0,
                  0.0,
                  1.0};
  ASSERT_EQ(mat, res);
  mat = Matrix3x3::rotation(-10.5);
  res = Matrix3x3{0.9832549075639546,
                  -0.1822355254921475,
                  0.0,
                  0.1822355254921475,
                  0.9832549075639546,
                  0.0,
                  0.0,
                  0.0,
                  1.0};
  ASSERT_EQ(mat, res);
}

TEST_F(MatrixTest, Matrix3x3FunctionsTranslation) {
  Matrix3x3 mat = Matrix3x3::translation(Vector2{1.0, 0.5});
  Matrix3x3 res{1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.5, 1.0};
  ASSERT_EQ(mat, res);
}

TEST_F(MatrixTest, Matrix3x3FunctionsTransformInverse) {
  Matrix3x3 mat = Matrix3x3::translation(Vector2{1.0, 0.5}) *
                  Matrix3x3::rotation(45.0) *
                  Matrix3x3::scale(Vector2{0.5, 0.5});
  Matrix3x3 res{0.3535533905932738,
                0.3535533905932738,
                0.0,
                -0.3535533905932738,
                0.3535533905932738,
                0.0,
                1.0,
                0.5,
                1.0};
  ASSERT_EQ(mat, res);

  mat = Matrix3x3::transformInverse(mat);
  res = Matrix3x3{1.4142135623730949,  -1.4142135623730949, 0.0,
                  1.4142135623730949,  1.4142135623730949,  0.0,
                  -2.1213203435596424, 0.7071067811865475,  1.0};
  ASSERT_EQ(mat, res);
  mat = Matrix3x3::inverse(Matrix3x3::scale(Vector2{0.5, 0.5})) *
        Matrix3x3::inverse(Matrix3x3::rotation(45.0)) *
        Matrix3x3::inverse(Matrix3x3::translation(Vector2{1.0, 0.5}));
  ASSERT_EQ(mat, res);
}

TEST_F(MatrixTest, Matrix3x3FunctionsTransformInverseUnitScale) {
  Matrix3x3 mat =
      Matrix3x3::translation(Vector2{1.0, 0.5}) * Matrix3x3::rotation(45.0);
  Matrix3x3 res{0.70710678118654752,
                0.70710678118654752,
                0.0,
                -0.70710678118654752,
                0.70710678118654752,
                0.0,
                1.0,
                0.5,
                1.0};
  ASSERT_EQ(mat, res);

  mat = Matrix3x3::transformInverseUnitScale(mat);
  res = Matrix3x3{0.7071067811865475288016887242097,
                  -0.7071067811865475288016887242097,
                  0.0,
                  0.7071067811865475288016887242097,
                  0.7071067811865475288016887242097,
                  0.0,
                  -1.0606601717798213,
                  0.3535533905932738,
                  1.0};
  ASSERT_EQ(mat, res);
  mat = Matrix3x3::inverse(Matrix3x3::rotation(45.0)) *
        Matrix3x3::inverse(Matrix3x3::translation(Vector2{1.0, 0.5}));
  ASSERT_EQ(mat, res);
}

TEST_F(MatrixTest, Matrix4x4Getters) {
  Vector4 vec1{1.0, 0.0, 0.0, 0.0};
  Vector4 vec2{0.0, 1.0, 0.0, 0.0};
  Vector4 vec3{0.0, 0.0, 1.0, 0.0};
  Vector4 vec4{0.0, 0.0, 0.0, 1.0};
  ASSERT_EQ(mat4_def_init[0], vec1);
  ASSERT_EQ(mat4_def_init[1], vec2);
  ASSERT_EQ(mat4_def_init[2], vec3);
  ASSERT_EQ(mat4_def_init[3], vec4);
  ;
  ASSERT_EQ(mat4_def_init.getRow(0), vec1);
  ASSERT_EQ(mat4_def_init.getRow(1), vec2);
  ASSERT_EQ(mat4_def_init.getRow(2), vec3);
  ASSERT_EQ(mat4_def_init.getRow(3), vec4);
  ASSERT_EQ(mat4_def_init.getColumn(0), vec1);
  ASSERT_EQ(mat4_def_init.getColumn(1), vec2);
  ASSERT_EQ(mat4_def_init.getColumn(2), vec3);
  ASSERT_EQ(mat4_def_init.getColumn(3), vec4);
  ASSERT_EQ(mat4_def_init.getXAxis(), vec1);
  ASSERT_EQ(mat4_def_init.getYAxis(), vec2);
  ASSERT_EQ(mat4_def_init.getZAxis(), vec3);
  ASSERT_EQ(mat4_def_init.getTranslation(), vec4);

  vec1 = Vector4{1.0, 0.8, -2.0, -0.155};
  vec2 = Vector4{0.5, 1.5, 1.45, 3.1};
  vec3 = Vector4{-0.25, 0.45, 3.0, 0.1};
  vec4 = Vector4{-0.125, 0.35, 0.5, 1.0};
  ASSERT_EQ(mat4_scalar_init[0], vec1);
  ASSERT_EQ(mat4_scalar_init[1], vec2);
  ASSERT_EQ(mat4_scalar_init[2], vec3);
  ASSERT_EQ(mat4_scalar_init[3], vec4);
  ASSERT_EQ(mat4_scalar_init.getRow(0), vec1);
  ASSERT_EQ(mat4_scalar_init.getRow(1), vec2);
  ASSERT_EQ(mat4_scalar_init.getRow(2), vec3);
  ASSERT_EQ(mat4_scalar_init.getRow(3), vec4);
  vec1 = Vector4{1.0, 0.5, -0.25, -0.125};
  vec2 = Vector4{0.8, 1.5, 0.45, 0.35};
  vec3 = Vector4{-2.0, 1.45, 3.0, 0.5};
  vec4 = Vector4{-0.155, 3.1, 0.1, 1.0};
  ASSERT_EQ(mat4_scalar_init.getColumn(0), vec1);
  ASSERT_EQ(mat4_scalar_init.getColumn(1), vec2);
  ASSERT_EQ(mat4_scalar_init.getColumn(2), vec3);
  ASSERT_EQ(mat4_scalar_init.getColumn(3), vec4);
  ASSERT_EQ(mat4_scalar_init.getXAxis(), vec1);
  ASSERT_EQ(mat4_scalar_init.getYAxis(), vec2);
  ASSERT_EQ(mat4_scalar_init.getZAxis(), vec3);
  ASSERT_EQ(mat4_scalar_init.getTranslation(), vec4);

  vec1 = Vector4{-1.0, 0.25, 1.25, 0.25};
  vec2 = Vector4{0.5, -0.125, -2.125, 2.0};
  vec3 = Vector4{0.25, -0.65, 3.65, -1.5};
  vec4 = Vector4{0.0, 0.0, 0.0, 1.0};
  ASSERT_EQ(mat4_arr_init[0], vec1);
  ASSERT_EQ(mat4_arr_init[1], vec2);
  ASSERT_EQ(mat4_arr_init[2], vec3);
  ASSERT_EQ(mat4_arr_init[3], vec4);
  ASSERT_EQ(mat4_arr_init.getRow(0), vec1);
  ASSERT_EQ(mat4_arr_init.getRow(1), vec2);
  ASSERT_EQ(mat4_arr_init.getRow(2), vec3);
  ASSERT_EQ(mat4_arr_init.getRow(3), vec4);
  vec1 = Vector4{-1.0, 0.5, 0.25, 0.0};
  vec2 = Vector4{0.25, -0.125, -0.65, 0.0};
  vec3 = Vector4{1.25, -2.125, 3.65, 0.0};
  vec4 = Vector4{0.25, 2.0, -1.5, 1.0};
  ASSERT_EQ(mat4_arr_init.getColumn(0), vec1);
  ASSERT_EQ(mat4_arr_init.getColumn(1), vec2);
  ASSERT_EQ(mat4_arr_init.getColumn(2), vec3);
  ASSERT_EQ(mat4_arr_init.getColumn(3), vec4);
  ASSERT_EQ(mat4_arr_init.getXAxis(), vec1);
  ASSERT_EQ(mat4_arr_init.getYAxis(), vec2);
  ASSERT_EQ(mat4_arr_init.getZAxis(), vec3);
  ASSERT_EQ(mat4_arr_init.getTranslation(), vec4);

  vec1 = Vector4{10000.12345, 12.565, -3221.87, -124.12};
  vec2 = Vector4{-10.825, 2540.9832235, 1233.9832235, 22.4728};
  vec3 = Vector4{2512.483, 5543.12, 4322.12, 5.1138287};
  vec4 = Vector4{0.0, 0.0, 0.0, 1.0};
  ASSERT_EQ(mat4_vec_init[0], vec1);
  ASSERT_EQ(mat4_vec_init[1], vec2);
  ASSERT_EQ(mat4_vec_init[2], vec3);
  ASSERT_EQ(mat4_vec_init[3], vec4);
  ASSERT_EQ(mat4_vec_init.getRow(0), vec1);
  ASSERT_EQ(mat4_vec_init.getRow(1), vec2);
  ASSERT_EQ(mat4_vec_init.getRow(2), vec3);
  ASSERT_EQ(mat4_vec_init.getRow(3), vec4);
  vec1 = Vector4{10000.12345, -10.825, 2512.483, 0.0};
  vec2 = Vector4{12.565, 2540.9832235, 5543.12, 0.0};
  vec3 = Vector4{-3221.87, 1233.9832235, 4322.12, 0.0};
  vec4 = Vector4{-124.12, 22.4728, 5.1138287, 1.0};
  ASSERT_EQ(mat4_vec_init.getColumn(0), vec1);
  ASSERT_EQ(mat4_vec_init.getColumn(1), vec2);
  ASSERT_EQ(mat4_vec_init.getColumn(2), vec3);
  ASSERT_EQ(mat4_vec_init.getColumn(3), vec4);
  ASSERT_EQ(mat4_vec_init.getXAxis(), vec1);
  ASSERT_EQ(mat4_vec_init.getYAxis(), vec2);
  ASSERT_EQ(mat4_vec_init.getZAxis(), vec3);
  ASSERT_EQ(mat4_vec_init.getTranslation(), vec4);
}

TEST_F(MatrixTest, Matrix4x4Operators) {
  ASSERT_EQ(mat4_def_init, mat4_def_init);
  ASSERT_EQ(mat4_scalar_init, mat4_scalar_init);

  ASSERT_NE(mat4_def_init, mat4_scalar_init);
  ASSERT_NE(mat4_scalar_init, mat4_arr_init);

  Matrix4x4 mat1 = mat4_scalar_init * mat4_def_init;
  ASSERT_EQ(mat4_scalar_init, mat1);
  mat1 = mat4_arr_init * mat4_scalar_init;
  Matrix4x4 mat2{-1.21875, 0.71875, -0.8,     -0.125,   0.225,  -0.04375,
                 0.3425,   0.35,    6.2375,   -6.55625, 8.7575, 0.5,
                 1.305,    1.3225,  -3.18875, 1.0};
  ASSERT_EQ(mat1, mat2);
  mat1 = mat4_vec_init * mat4_arr_init;
  mat2 = Matrix4x4{-10799.308450, 1589.812417625,   1339.607,     0.0,
                   4592.6757375,  -1122.4182482125, -2874.14725,  0.0,
                   713.6281875,   -909.0818341625,  7137.21175,   0.0,
                   7233.8458625,  3250.75816175,    5236.2945787, 1.0};
  ASSERT_EQ(mat1, mat2);

  Vector4 vec1{10.25, -0.025530, 3.456, 1.0};
  Vector4 vec2 = mat4_def_init * vec1;
  ASSERT_EQ(vec1, vec2);

  vec2 = mat4_scalar_init * vec1;
  Vector4 vec3 = Vector4{3.162576, 13.197905, 7.8940115, 1.4378145};
  ASSERT_EQ(vec2, vec3);

  std::string matStr =
      "[1.000000, 0.800000, -2.000000, -0.155000\n0.500000, 1.500000, "
      "1.450000, 3.100000\n-0.250000, 0.450000, 3.000000, 0.100000\n-0.125000, "
      "0.350000, 0.500000, 1.000000]";
  ASSERT_EQ(static_cast<std::string>(mat4_scalar_init), matStr);
}

TEST_F(MatrixTest, Matrix4x4Setters) {
  Matrix4x4 mat1 = mat4_def_init;
  mat1.setXAxis(2.0, 4.0, 10.0, 1.0);
  mat1.setYAxis(-0.5, 10.4, -7.8, 2.1);
  mat1.setZAxis(2.5, 0.789, 1.1, 0.15);
  mat1.setTranslation(2.15, 0.1789, 1.12, 1.0);
  Matrix4x4 mat2{2.0, 4.0,   10.0, 1.0,  -0.5, 10.4,   -7.8, 2.1,
                 2.5, 0.789, 1.1,  0.15, 2.15, 0.1789, 1.12, 1.0};
  ASSERT_EQ(mat1, mat2);

  std::array<cgmath::FLOAT, 4> x_axis{7.200045, 0.00019929, 0.1926345623,
                                      0.783251667};
  std::array<cgmath::FLOAT, 4> y_axis{-0.1, -100.25, 123.456, 0.00000004588};
  std::array<cgmath::FLOAT, 4> z_axis{2.1, 10.25, 13.46, 1.0002465};
  std::array<cgmath::FLOAT, 4> w_axis{2.1, 10.25, 13.46, 0.12000043};
  mat1.setXAxis(x_axis);
  mat1.setYAxis(y_axis);
  mat1.setZAxis(z_axis);
  mat1.setTranslation(w_axis);
  mat2 = Matrix4x4{x_axis, y_axis, z_axis, w_axis};
  ASSERT_EQ(mat1, mat2);

  Vector4 x_vec{7.20045, 0.019929, 0.1935623, 10.38378110123};
  Vector4 y_vec{-0.041, -10.2115, 12.0000456, 1.0};
  Vector4 z_vec{3.4, 11.14, 13.46733416, -5.2};
  Vector4 w_vec{0.21, 1.25, 1.46, -1.05};
  mat1.setXAxis(x_vec);
  mat1.setYAxis(y_vec);
  mat1.setZAxis(z_vec);
  mat1.setTranslation(w_vec);
  mat2 = Matrix4x4{x_vec, y_vec, z_vec, w_vec};
  ASSERT_EQ(mat1, mat2);
}

TEST_F(MatrixTest, Matrix4x4FunctionsDeterminant) {
  ASSERT_EQ_FLT(mat4_def_init.determinant(), 1.0);
  ASSERT_EQ_FLT(mat4_scalar_init.determinant(), -1.8725859375);
  ASSERT_EQ_FLT(mat4_arr_init.determinant(), 0.88125);
  ASSERT_EQ_FLT(mat4_vec_init.determinant(), 62225543347.1125499578925);
}

TEST_F(MatrixTest, Matrix4x4FunctionsTranspose) {
  ASSERT_EQ(Matrix4x4::transpose(mat4_def_init), mat4_def_init);

  Matrix4x4 mat{1.0,   0.8,  -2.0, -0.155, 0.5,    1.5,  1.45, 3.1,
                -0.25, 0.45, 3.0,  0.1,    -0.125, 0.35, 0.5,  1.0};
  ASSERT_EQ(Matrix4x4::transpose(mat4_scalar_init), mat);

  mat = Matrix4x4{-1.0, 0.25,  1.25, 0.25, 0.5, -0.125, -2.125, 2.0,
                  0.25, -0.65, 3.65, -1.5, 0.0, 0.0,    0.0,    1.0};
  ASSERT_EQ(Matrix4x4::transpose(mat4_arr_init), mat);

  mat =
      Matrix4x4{10000.12345,  12.565,  -3221.87, -124.12, -10.825, 2540.9832235,
                1233.9832235, 22.4728, 2512.483, 5543.12, 4322.12, 5.1138287,
                0.0,          0.0,     0.0,      1.0};
  ASSERT_EQ(Matrix4x4::transpose(mat4_vec_init), mat);
}

TEST_F(MatrixTest, Matrix4x4FunctionsInverse) {
  ASSERT_EQ(Matrix4x4::inverse(mat4_def_init), mat4_def_init);

  Matrix4x4 mat{-0.6759369354710857, 1.3854504341005712, -0.2493210007885152,
                -0.4447392684748280, 1.7718145445594536, -1.3010083816246751,
                0.3256692992227493,  0.5139951020271934, -0.3804431538939718,
                0.8587890242019934,  0.1875414596292727, -0.4419022825220806,
                -5.5593509977429274, 4.1619918979018820, -1.0669737286756699,
                -0.6181291746456897};
  ASSERT_EQ(Matrix4x4::inverse(mat4_scalar_init), mat);

  mat = Matrix4x4{-2.085106382978723,
                  -2.673758865248227,
                  -0.3333333333333333,
                  0.0,
                  -1.957446808510638,
                  -4.49645390070922,
                  -0.6666666666666667,
                  0.0,
                  -0.425531914893617,
                  -1.702127659574468,
                  0.0,
                  0.0,
                  3.7978723404255319,
                  7.1081560283687943,
                  1.4166666666666667,
                  1.0};
  ASSERT_EQ(Matrix4x4::inverse(mat4_arr_init), mat);

  mat = Matrix4x4{0.000066569403837897021, 0.000050576477938862198,
                  -0.00010356167386729587, 0.0,
                  -0.00028788048297582672, 0.00082468749019457542,
                  -0.00089031468379841395, 0.0,
                  0.00013181439930137358,  -0.00019775010655732485,
                  0.00040835773492178042,  0.0,
                  0.0140579986639583,      -0.0112442244175322,
                  0.0050655173609462,      1.0};
  Matrix4x4 tmp = Matrix4x4::inverse(mat4_vec_init);
  ASSERT_EQ(tmp, mat);
}