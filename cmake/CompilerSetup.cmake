list(APPEND COMPILE_OPTIONS_LIST)

# Find compiler option suitable for the selected architecture.
set(ARCH_OPTION $CACHE{SIMD_ARCH})
if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    string(PREPEND ARCH_OPTION "/arch:")
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    if(${ARCH_OPTION} STREQUAL "AVX512")
        set(ARCH_OPTION "-march=x86-64-v4")
    elseif (${ARCH_OPTION} MATCHES "AVX")
        set(ARCH_OPTION "-march=x86-64-v3")
    else()
        set(ARCH_OPTION "-march=x86-64")
    endif ()
endif ()
list(APPEND COMPILE_OPTIONS_LIST ${ARCH_OPTION})

# Set Compiler options
list(JOIN COMPILE_OPTIONS_LIST " " COMPILE_OPTIONS)
message("Compiler Options: " ${COMPILE_OPTIONS})
add_compile_options(${COMPILE_OPTIONS})