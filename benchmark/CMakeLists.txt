# benchmark/CMakeLists.txt - Mixed state benchmark executables

set(SRC_DIR "${CMAKE_SOURCE_DIR}/src")
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

# Qulacs - check extern/ first (requires Boost)
# Qulacs provides dm_* functions for density matrix gate operations
# CMP0167 (CMake 3.30+): Use Boost's own BoostConfig.cmake instead of removed FindBoost module
if(POLICY CMP0167)
    cmake_policy(SET CMP0167 NEW)
endif()
find_package(Boost 1.71.0 QUIET CONFIG)
if(Boost_FOUND)
    find_path(QULACS_INCLUDE_DIR csim/update_ops_dm.hpp
        HINTS
            ${EXTERN_DIR}/qulacs/src
            ${QULACS_DIR}/include
            $ENV{QULACS_DIR}/include
            /usr/local/include/qulacs
    )
    find_library(QULACS_CSIM_LIBRARY
        NAMES csim_static csim
        HINTS
            ${EXTERN_DIR}/qulacs/build/lib
            ${EXTERN_DIR}/qulacs/lib
            ${QULACS_DIR}/lib
            $ENV{QULACS_DIR}/lib
            /usr/local/lib
    )
    if(QULACS_INCLUDE_DIR AND QULACS_CSIM_LIBRARY)
        set(QULACS_FOUND TRUE)
    endif()
endif()

#=======================================================================================================================
# bench_mixed - Main mixed state benchmark
#=======================================================================================================================

set(BENCH_MIXED_SRC
        "${SRC_DIR}/state.c"
        "${SRC_DIR}/state_pure.c"
        "${SRC_DIR}/state_packed.c"
        "${SRC_DIR}/state_tiled.c"
        "${SRC_DIR}/gate.c"
        "${SRC_DIR}/gate_pure.c"
        "${SRC_DIR}/gate_packed.c"
        "${SRC_DIR}/gate_tiled.c"
        "${BENCH_SRC_DIR}/bench_mixed.c"
)

if(QuEST_FOUND)
    list(APPEND BENCH_MIXED_SRC "${BENCH_SRC_DIR}/bench_quest.c")
endif()

if(QULACS_FOUND)
    list(APPEND BENCH_MIXED_SRC "${BENCH_SRC_DIR}/bench_qulacs.cpp")
endif()

add_executable(bench_mixed ${BENCH_MIXED_SRC})

# Output binary to cmake-build-debug/benchmark/
set_target_properties(bench_mixed PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/benchmark"
)

target_include_directories(bench_mixed
        PRIVATE
            "${CMAKE_SOURCE_DIR}/include"
            "${CMAKE_CURRENT_SOURCE_DIR}/include"
)

target_link_libraries(bench_mixed
        PRIVATE
            blas_compiler_flags
            omp_compiler_flags
            ${PROJECT_COMPILER_FLAGS}
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
    target_link_libraries(bench_mixed PRIVATE ${QULACS_CSIM_LIBRARY})
    # Qulacs requires C++14 minimum (we use C++17 for consistency)
    set_property(TARGET bench_mixed PROPERTY CXX_STANDARD 17)
    set_property(TARGET bench_mixed PROPERTY CXX_STANDARD_REQUIRED ON)
    # Enable OpenMP if available (Qulacs uses it internally)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(bench_mixed PRIVATE OpenMP::OpenMP_CXX)
    endif()
endif()

#=======================================================================================================================
# Status summary
#=======================================================================================================================

message(STATUS "Benchmark configuration:")
message(STATUS "  QuEST:     ${QuEST_FOUND}")
message(STATUS "  Qulacs:    ${QULACS_FOUND}")
message(STATUS "  OpenMP:    ${OpenMP_CXX_FOUND}")
