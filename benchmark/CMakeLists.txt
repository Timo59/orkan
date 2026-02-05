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

# OpenMP - for Eigen3 parallelization in Quantum++
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

# Quantum++ (requires Eigen3) - check extern/ first
find_package(Eigen3 QUIET CONFIG HINTS /opt/homebrew/opt/eigen/share/eigen3/cmake)
if(Eigen3_FOUND)
    find_path(QPP_INCLUDE_DIR qpp/qpp.hpp
        HINTS
            ${EXTERN_DIR}/qpp/include
            ${QPP_DIR}/include
            $ENV{QPP_DIR}/include
            /usr/local/include
            /opt/qpp/include
    )
    if(QPP_INCLUDE_DIR)
        set(QPP_FOUND TRUE)
    endif()
endif()

# Qulacs - check extern/ first (requires Boost)
# Qulacs provides dm_* functions for density matrix gate operations
find_package(Boost 1.71.0 QUIET)
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
        "${SRC_DIR}/qhipster.c"
        "${SRC_DIR}/mhipster.c"
        "${BENCH_SRC_DIR}/bench_mixed.c"
)

if(QuEST_FOUND)
    list(APPEND BENCH_MIXED_SRC "${BENCH_SRC_DIR}/bench_quest.c")
endif()

if(QPP_FOUND)
    list(APPEND BENCH_MIXED_SRC "${BENCH_SRC_DIR}/bench_qpp.cpp")
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

if(QPP_FOUND)
    target_compile_definitions(bench_mixed PRIVATE WITH_QPP)
    target_include_directories(bench_mixed PRIVATE
        ${QPP_INCLUDE_DIR}
        ${EXTERN_DIR}/qpp/qasmtools/include  # qasmtools dependency
    )
    target_link_libraries(bench_mixed PRIVATE Eigen3::Eigen)
    set_property(TARGET bench_mixed PROPERTY CXX_STANDARD 17)
    set_property(TARGET bench_mixed PROPERTY CXX_STANDARD_REQUIRED ON)
    # Enable OpenMP for Eigen3 parallelization in Quantum++
    if(OpenMP_CXX_FOUND)
        target_link_libraries(bench_mixed PRIVATE OpenMP::OpenMP_CXX)
    endif()
endif()

if(QULACS_FOUND)
    target_compile_definitions(bench_mixed PRIVATE WITH_QULACS)
    # Only include the csim headers directory, NOT the Qulacs Eigen directory
    # (to avoid conflicts with Quantum++'s system Eigen)
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
message(STATUS "  Quantum++: ${QPP_FOUND}")
message(STATUS "  Qulacs:    ${QULACS_FOUND}")
message(STATUS "  OpenMP:    ${OpenMP_CXX_FOUND}")
