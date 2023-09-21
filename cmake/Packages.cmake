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
endif ()