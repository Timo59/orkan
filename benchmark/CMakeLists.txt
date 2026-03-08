# benchmark/CMakeLists.txt - Mixed state benchmark executables

# Benchmarks require Release build for meaningful results
if(NOT CMAKE_BUILD_TYPE STREQUAL "Release")
    message(STATUS "Skipping benchmark (requires -DCMAKE_BUILD_TYPE=Release)")
    return()
endif()

set(BENCH_SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")
set(EXTERN_DIR "${CMAKE_SOURCE_DIR}/extern")

#=======================================================================================================================
# Find external frameworks (auto-detect)
#=======================================================================================================================

# QuEST - check extern/ first, then system paths
# QuEST v4 requires the root directory as include path (headers use quest/include/...)
find_path(QuEST_INCLUDE_DIR quest/include/quest.h
    HINTS
        ${EXTERN_DIR}/QuEST
        ${QuEST_DIR}
        $ENV{QuEST_DIR}
        /usr/local
        /opt/quest
)
find_library(QuEST_LIBRARY
    NAMES QuEST quest
    HINTS
        ${EXTERN_DIR}/QuEST/build
        ${QuEST_DIR}/lib
        $ENV{QuEST_DIR}/lib
        /usr/local/lib
        /opt/quest/lib
)
if(QuEST_INCLUDE_DIR AND QuEST_LIBRARY)
    set(QuEST_FOUND TRUE)
endif()

# OpenMP - for parallelization in Qulacs
# macOS requires explicit hints for Homebrew's libomp
if(APPLE)
    execute_process(
        COMMAND brew --prefix libomp
        OUTPUT_VARIABLE LIBOMP_PREFIX
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
    )
    if(LIBOMP_PREFIX)
        set(OpenMP_ROOT ${LIBOMP_PREFIX})
        set(OpenMP_C_FLAGS "-Xpreprocessor -fopenmp -I${LIBOMP_PREFIX}/include")
        set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I${LIBOMP_PREFIX}/include")
        set(OpenMP_C_LIB_NAMES "omp")
        set(OpenMP_CXX_LIB_NAMES "omp")
        set(OpenMP_omp_LIBRARY "${LIBOMP_PREFIX}/lib/libomp.dylib")
    endif()
endif()
find_package(OpenMP)

# Qulacs - build csim_static from the submodule sources
# Only the csim sub-library is needed; it depends on Eigen but not Boost.
# We build it in-tree so `cmake --build` handles everything.
set(QULACS_CSIM_DIR "${EXTERN_DIR}/qulacs/src/csim")
if(EXISTS "${QULACS_CSIM_DIR}/update_ops_dm.hpp")
    # Collect csim sources
    file(GLOB QULACS_CSIM_SRC "${QULACS_CSIM_DIR}/*.cpp")
    add_library(csim_static STATIC ${QULACS_CSIM_SRC})

    # Eigen - check Homebrew, then system paths
    find_path(EIGEN3_INCLUDE_DIR Eigen/Core
        HINTS
            /opt/homebrew/include/eigen3
            /usr/local/include/eigen3
            /usr/include/eigen3
    )
    if(EIGEN3_INCLUDE_DIR)
        target_include_directories(csim_static SYSTEM PRIVATE ${EIGEN3_INCLUDE_DIR})
    endif()

    # Qulacs headers reference each other via <csim/...>
    target_include_directories(csim_static PRIVATE "${EXTERN_DIR}/qulacs/src")
    # Suppress warnings in third-party code
    target_compile_options(csim_static PRIVATE -w)
    set_property(TARGET csim_static PROPERTY CXX_STANDARD 17)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(csim_static PRIVATE OpenMP::OpenMP_CXX)
    endif()

    set(QULACS_INCLUDE_DIR "${EXTERN_DIR}/qulacs/src")
    set(QULACS_FOUND TRUE)
endif()

#=======================================================================================================================
# bench_mixed - Main mixed state benchmark
#=======================================================================================================================

set(BENCH_MIXED_SRC
        "${BENCH_SRC_DIR}/bench_mixed.c"
        "${BENCH_SRC_DIR}/bench_baselines.c"
        "${BENCH_SRC_DIR}/bench_main.c"
)

if(QuEST_FOUND)
    list(APPEND BENCH_MIXED_SRC "${BENCH_SRC_DIR}/bench_quest.c")
endif()

if(QULACS_FOUND)
    list(APPEND BENCH_MIXED_SRC "${BENCH_SRC_DIR}/bench_qulacs.cpp")
endif()

add_executable(bench_mixed ${BENCH_MIXED_SRC})

# Output binary to cmake-build-*/benchmark/
set_target_properties(bench_mixed PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/benchmark"
)

target_include_directories(bench_mixed
        PRIVATE
            "${CMAKE_CURRENT_SOURCE_DIR}/include"
)

# Link libq — inherits compiler flags, BLAS, OpenMP, and include dirs via PUBLIC deps
target_link_libraries(bench_mixed
        PRIVATE
            q
)

target_compile_definitions(bench_mixed
        PRIVATE
            _GNU_SOURCE
)

if(QuEST_FOUND)
    target_compile_definitions(bench_mixed PRIVATE WITH_QUEST)
    target_include_directories(bench_mixed PRIVATE
        ${QuEST_INCLUDE_DIR}
        ${EXTERN_DIR}/QuEST/build  # For generated config.h
    )
    target_link_libraries(bench_mixed PRIVATE ${QuEST_LIBRARY})
endif()

if(QULACS_FOUND)
    target_compile_definitions(bench_mixed PRIVATE WITH_QULACS)
    target_include_directories(bench_mixed PRIVATE
        ${QULACS_INCLUDE_DIR}
    )
    target_link_libraries(bench_mixed PRIVATE csim_static)
    # Qulacs requires C++14 minimum (we use C++17 for consistency)
    set_property(TARGET bench_mixed PROPERTY CXX_STANDARD 17)
    set_property(TARGET bench_mixed PROPERTY CXX_STANDARD_REQUIRED ON)
    # Enable OpenMP if available (Qulacs uses it internally)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(bench_mixed PRIVATE OpenMP::OpenMP_CXX)
    endif()
endif()

#=======================================================================================================================
# bench target - Run tests then benchmark in one command
#=======================================================================================================================

if(BUILD_TESTS)
    add_custom_target(bench
        COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure
        COMMAND $<TARGET_FILE:bench_mixed>
        DEPENDS test_state test_gate bench_mixed
        COMMENT "Running tests, then benchmark..."
    )
endif()

#=======================================================================================================================
# Status summary
#=======================================================================================================================

message(STATUS "Benchmark configuration:")
message(STATUS "  QuEST:     ${QuEST_FOUND}")
message(STATUS "  Qulacs:    ${QULACS_FOUND}")
message(STATUS "  OpenMP:    ${OpenMP_CXX_FOUND}")
