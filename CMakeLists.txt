cmake_minimum_required(VERSION 3.26)
set(CMAKE_CXX_STANDARD 17)

# Library setup
project(cgmath
        VERSION 0.1.0
        DESCRIPTION "A simple, fast and easy to use math library targeted towards computer graphics applications."
        HOMEPAGE_URL "https://github.com/SurajShettigar/cgmath"
        LANGUAGES C CXX)
add_library(cgmath STATIC
        src/vector2.cpp
        src/vector3.cpp
        src/vector4.cpp
        src/matrix2x2.cpp
        src/matrix3x3.cpp
        src/matrix4x4.cpp)
set_target_properties(cgmath PROPERTIES VERSION ${PROJECT_VERSION})
target_include_directories(cgmath PUBLIC include PRIVATE internal)

# Library install setup
include(GNUInstallDirs)
install(TARGETS cgmath LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# Library configuration options
option(INCLUDE_TESTS "Include test files and test dependencies?" ON)
option(INCLUDE_BENCHMARKING_TESTS "Include benchmarking test cases?" ON)
option(DOUBLE_PRECISION "Uses double precision floating point values." ON)
option(SIMD_INTRINSICS "Use SIMD intrinsics for vector math. Much faster than scalar calculations." ON)
set(SIMD_ARCH "AVX2" CACHE STRING "Target Architecture used to compile SIMD intrinsics.")
set_property(CACHE SIMD_ARCH PROPERTY STRINGS "SSE4_2" "AVX" "AVX2" "AVX512")

# Cmake includes
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake ${PROJECT_SOURCE_DIR}/ext)
# Library dependencies
include(Packages)
# Compiler setup
include(CompilerSetup)

if (SIMD_INTRINSICS)
    message("Math Library using SIMD intrinsics")
    target_compile_definitions(cgmath PUBLIC USE_INTRINSICS)
endif ()
if (DOUBLE_PRECISION)
    message("Math Library using double precision")
    target_compile_definitions(cgmath PUBLIC USE_DOUBLE)
endif ()
if (INCLUDE_TESTS)
    include(CTest)
    enable_testing()
    add_executable(cgmath_test tests/vector_test.cpp
            tests/matrix_test.cpp)
    target_link_libraries(cgmath_test
            PUBLIC GTest::gtest
            PRIVATE cgmath)
    if(INCLUDE_BENCHMARKING_TESTS)
        target_link_libraries(cgmath_test
                PUBLIC benchmark::benchmark)
        target_compile_definitions(cgmath PUBLIC BENCHMARK_TEST)
    endif ()
    add_test(NAME vector_test COMMAND cgmath_test)
endif ()


