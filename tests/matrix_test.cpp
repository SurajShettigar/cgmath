// Copyright 2023 Suraj Shettigar
// SPDX-License-Identifier: Apache-2.0

#include <gtest/gtest.h>

#include <matrix2x2.hpp>
#include <matrix3x3.hpp>

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
  }
};

TEST_F(MatrixTest, Matrix2x2Getters) {
  Vector2 vec1{1.0, 0.0};
  Vector2 vec2{0.0, 1.0};
  ASSERT_EQ(mat2_def_init[0], vec1);
  ASSERT_EQ(mat2_def_init[1], vec2);
  ;
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
  ASSERT_EQ(mat2_arr_init, mat2_arr_init);
  ASSERT_EQ(mat2_vec_init, mat2_vec_init);

  ASSERT_NE(mat2_def_init, mat2_scalar_init);
  ASSERT_NE(mat2_scalar_init, mat2_arr_init);
  ASSERT_NE(mat2_arr_init, mat2_vec_init);
  ASSERT_NE(mat2_vec_init, mat2_def_init);

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

  // TODO: Test matrix2x2 scalar multiplication.
  // TODO: Test matrix2x2 string conversion.
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
  ASSERT_EQ(mat3_arr_init, mat3_arr_init);
  ASSERT_EQ(mat3_vec_init, mat3_vec_init);

  ASSERT_NE(mat3_def_init, mat3_scalar_init);
  ASSERT_NE(mat3_scalar_init, mat3_arr_init);
  ASSERT_NE(mat3_arr_init, mat3_vec_init);
  ASSERT_NE(mat3_vec_init, mat3_def_init);

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

  // TODO: Test matrix3x3 scalar multiplication.
  // TODO: Test matrix3x3 string conversion.
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
  mat = Matrix3x3{0.000066569403837897021,
                  0.000050576477938862198,
                  -0.00010356167386729587,
                  -0.00028788048297582672,
                  0.00082468749019457542,
                  -0.00089031468379841395,
                  0.00013181439930137358,
                  -0.00019775010655732485,
                  0.00040835773492178042};
  ASSERT_EQ(Matrix3x3::inverse(mat3_vec_init), mat);
}