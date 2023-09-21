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
        src/vector4.cpp)
set_target_properties(cgmath PROPERTIES VERSION ${PROJECT_VERSION}
        PUBLIC_HEADER include/vector.hpp)
target_include_directories(cgmath PUBLIC include)

# Library install setup
include(GNUInstallDirs)
install(TARGETS cgmath LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# Library configuration options
option(INCLUDE_TESTS "Include test files and test dependencies?" OFF)
option(USE_SIMD_INTRINSICS "Use SIMD intrinsics for vector math. Much faster than scalar calculations." ON)
option(DOUBLE_PRECISION "Uses double precision floating point values." OFF)

# Library dependencies
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake ${PROJECT_SOURCE_DIR}/ext)
include(Packages)

if (USE_SIMD_INTRINSICS)
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
    add_executable(cgmath_test tests/vector_test.cpp)
    target_link_libraries(cgmath_test
            PUBLIC GTest::gtest GTest::gtest_main
            PRIVATE cgmath)
    add_test(NAME vector_test COMMAND cgmath_test)
endif ()


