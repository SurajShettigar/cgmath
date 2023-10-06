if (INCLUDE_TESTS)
    find_package(GTest 1.14.0 EXACT)
    if (NOT GTest_FOUND)
        include(FetchContent)
        message("Downloading Google Test...")
        FetchContent_Declare(GTest
                GIT_REPOSITORY https://github.com/google/googletest.git
                GIT_TAG v1.14.0
                FIND_PACKAGE_ARGS NAMES GTest
        )
        FetchContent_MakeAvailable(GTest)
        find_package(GTest 1.14.0 EXACT REQUIRED)
    endif ()
    if ($CACHE{INCLUDE_BENCHMARKING_TESTS})
        message("Downloading Google Benchmark...")
        FetchContent_Declare(GBenchmark
                GIT_REPOSITORY https://github.com/google/benchmark.git
                GIT_TAG v1.8.3
        )
        FetchContent_MakeAvailable(GBenchmark)
        set(BENCHMARK_ENABLE_TESTING OFF)
        set(BENCHMARK_ENABLE_GTEST_TESTS OFF)
        set(BENCHMARK_USE_BUNDLED_GTEST OFF)
    endif ()
endif ()