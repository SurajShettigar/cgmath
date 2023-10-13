// Copyright 2023 Suraj Shettigar
// SPDX-License-Identifier: Apache-2.0

#include <gtest/gtest.h>

#include <matrix4x4.hpp>

using namespace cgmath;
using namespace cgmath::internal;

class MatrixTest : public ::testing::Test {
 protected:
  Matrix2x2 mat2_def_init;
  Matrix2x2 mat2_scalar_init;
  Matrix2x2 mat2_arr_init;
  Matrix2x2 mat2_vec_init;

  void SetUp() override {
    mat2_def_init = Matrix2x2{};
    mat2_scalar_init = Matrix2x2{1.0, 0.5, -0.25, -0.125};
    std::array<FLOAT, 2> a{-1.0, 0.5};
    std::array<FLOAT, 2> b{0.25, -0.125};
    mat2_arr_init = Matrix2x2{a, b};
    mat2_vec_init =
        Matrix2x2{Vector{10000.12345, 12.565, -10.825, 2540.9832235}};
  }
};

TEST_F(MatrixTest, Matrix2x2Getters) {
  Vector vec1 = Vector{1.0, 0.0, 0.0, 0.0};
  Vector vec2 = Vector{0.0, 1.0, 0.0, 0.0};
  ASSERT_EQ(mat2_def_init[0], vec1);
  ASSERT_EQ(mat2_def_init[1], vec2);
  ;
  ASSERT_EQ(mat2_def_init.getRow(0), vec1);
  ASSERT_EQ(mat2_def_init.getRow(1), vec2);
  ASSERT_EQ(mat2_def_init.getColumn(0), vec1);
  ASSERT_EQ(mat2_def_init.getColumn(1), vec2);
  ASSERT_EQ(mat2_def_init.getXAxis(), vec1);
  ASSERT_EQ(mat2_def_init.getYAxis(), vec2);

  vec1 = Vector{1.0, -0.25, 0.0, 0.0};
  vec2 = Vector{0.5, -0.125, 0.0, 0.0};
  ASSERT_EQ(mat2_scalar_init[0], vec1);
  ASSERT_EQ(mat2_scalar_init[1], vec2);
  ASSERT_EQ(mat2_scalar_init.getRow(0), vec1);
  ASSERT_EQ(mat2_scalar_init.getRow(1), vec2);
  vec1 = Vector{1.0, 0.5, 0.0, 0.0};
  vec2 = Vector{-0.25, -0.125, 0.0, 0.0};
  ASSERT_EQ(mat2_scalar_init.getColumn(0), vec1);
  ASSERT_EQ(mat2_scalar_init.getColumn(1), vec2);
  ASSERT_EQ(mat2_scalar_init.getXAxis(), vec1);
  ASSERT_EQ(mat2_scalar_init.getYAxis(), vec2);

  vec1 = Vector{-1.0, 0.25, 0.0, 0.0};
  vec2 = Vector{0.5, -0.125, 0.0, 0.0};
  ASSERT_EQ(mat2_arr_init[0], vec1);
  ASSERT_EQ(mat2_arr_init[1], vec2);
  ASSERT_EQ(mat2_arr_init.getRow(0), vec1);
  ASSERT_EQ(mat2_arr_init.getRow(1), vec2);
  vec1 = Vector{-1.0, 0.5, 0.0, 0.0};
  vec2 = Vector{0.25, -0.125, 0.0, 0.0};
  ASSERT_EQ(mat2_arr_init.getColumn(0), vec1);
  ASSERT_EQ(mat2_arr_init.getColumn(1), vec2);
  ASSERT_EQ(mat2_arr_init.getXAxis(), vec1);
  ASSERT_EQ(mat2_arr_init.getYAxis(), vec2);

  vec1 = Vector{10000.12345, 12.565, 0.0, 0.0};
  vec2 = Vector{-10.825, 2540.9832235, 0.0, 0.0};
  ASSERT_EQ(mat2_vec_init[0], vec1);
  ASSERT_EQ(mat2_vec_init[1], vec2);
  ASSERT_EQ(mat2_vec_init.getRow(0), vec1);
  ASSERT_EQ(mat2_vec_init.getRow(1), vec2);
  vec1 = Vector{10000.12345, -10.825, 0.0, 0.0};
  vec2 = Vector{12.565, 2540.9832235, 0.0, 0.0};
  ASSERT_EQ(mat2_vec_init.getColumn(0), vec1);
  ASSERT_EQ(mat2_vec_init.getColumn(1), vec2);
  ASSERT_EQ(mat2_vec_init.getXAxis(), vec1);
  ASSERT_EQ(mat2_vec_init.getYAxis(), vec2);
}

TEST_F(MatrixTest, Matrix2x2Setters) {
  Matrix2x2 mat1 = mat2_def_init;
  mat1.setXAxis(2.0, 4.0);
  mat1.setYAxis(-0.5, 10.4);
  Matrix2x2 mat2 {2.0, 4.0, -0.5, 10.4};
  ASSERT_EQ(mat1, mat2);

  std::array<FLOAT, 2> x_axis{7.200045, 0.00019929};
  std::array<FLOAT, 2> y_axis{-0.1, -100.25};
  mat1.setXAxis(x_axis);
  mat1.setYAxis(y_axis);
  mat2 = Matrix2x2{x_axis, y_axis};
  ASSERT_EQ(mat1, mat2);

  Vector x_vec{7.200045, 0.00019929};
  Vector y_vec{-0.1, -100.25};
  mat1.setXAxis(x_vec);
  mat1.setYAxis(y_vec);
  mat2 = Matrix2x2{Vector{7.200045, -0.1, 0.00019929, -100.25}};
  ASSERT_EQ(mat1, mat2);
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
  Matrix2x2 mat2 {-0.875, 0.4375, 0.21875, -0.109375};
  ASSERT_EQ(mat1, mat2);
  mat1 = mat2_vec_init * mat2_arr_init;
  mat2 = Matrix2x2{-9993.840950, 1281.316611750, 2498.46023750, -320.32915293750};
  ASSERT_EQ(mat1, mat2);

  Vector vec1 {10.25, -0.025530};
  Vector vec2 = mat2_def_init * vec1;
  ASSERT_EQ(vec1, vec2);

  vec2 = mat2_scalar_init * vec1;
  Vector vec3 = Vector {10.2563825, 5.128191250};
  ASSERT_EQ(vec2, vec3);
}