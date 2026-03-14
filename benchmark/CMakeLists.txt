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
    # Suppress warnings in third-party code; match qlib's Release vectorisation flags
    target_compile_options(csim_static PRIVATE -w $<$<CONFIG:Release>:-march=native>)
    set_property(TARGET csim_static PROPERTY CXX_STANDARD 17)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(csim_static PRIVATE OpenMP::OpenMP_CXX)
    endif()

    set(QULACS_INCLUDE_DIR "${EXTERN_DIR}/qulacs/src")
    set(QULACS_FOUND TRUE)
endif()

# Qiskit-Aer DM — header-only kernel vendored to extern/aer-dm/
# Vendored from commit: c5841316078b2f141418459faf12e67070888fd7
set(AER_DM_INCLUDE_DIR "${EXTERN_DIR}/aer-dm/include")
if(EXISTS "${AER_DM_INCLUDE_DIR}/simulators/density_matrix/densitymatrix.hpp")
    set(AER_DM_FOUND TRUE)
endif()

# nlohmann/json — required by Aer's framework/json.hpp
# Check Homebrew first (macOS), then system paths
if(APPLE)
    execute_process(
        COMMAND brew --prefix nlohmann-json
        OUTPUT_VARIABLE NLOHMANN_JSON_PREFIX
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
    )
endif()
find_path(NLOHMANN_JSON_INCLUDE_DIR nlohmann/json.hpp
    HINTS
        "${NLOHMANN_JSON_PREFIX}/include"
        /opt/homebrew/include
        /usr/local/include
        /usr/include
)

if(AER_DM_FOUND)
    add_library(aer_dm_static INTERFACE)
    # SYSTEM suppresses vendored-header warnings without leaking -w to consumer TUs
    target_include_directories(aer_dm_static SYSTEM INTERFACE "${AER_DM_INCLUDE_DIR}")
    if(NLOHMANN_JSON_INCLUDE_DIR)
        target_include_directories(aer_dm_static SYSTEM INTERFACE "${NLOHMANN_JSON_INCLUDE_DIR}")
    else()
        message(WARNING "nlohmann/json.hpp not found — Aer DM headers will fail to compile. Install with: brew install nlohmann-json")
    endif()
    # -march=native is omitted: bench_gate receives it transitively through libq
    # CXX_STANDARD cannot be set on INTERFACE targets; set it on consuming targets instead
endif()

#=======================================================================================================================
# bench_gate - Main mixed state benchmark
#=======================================================================================================================

set(BENCH_MIXED_SRC
        "${BENCH_SRC_DIR}/bench_util.c"
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

if(AER_DM_FOUND)
    list(APPEND BENCH_MIXED_SRC "${BENCH_SRC_DIR}/bench_aer_dm.cpp")
endif()

add_executable(bench_gate ${BENCH_MIXED_SRC})

# Output binary to cmake-build-*/benchmark/
set_target_properties(bench_gate PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/benchmark"
    EXCLUDE_FROM_ALL TRUE
)

target_include_directories(bench_gate
        PRIVATE
            "${CMAKE_CURRENT_SOURCE_DIR}/include"
)

# Link libq — inherits compiler flags, BLAS, OpenMP, and include dirs via PUBLIC deps
target_link_libraries(bench_gate
        PRIVATE
            q
)

if(QuEST_FOUND)
    target_compile_definitions(bench_gate PRIVATE WITH_QUEST)
    target_include_directories(bench_gate PRIVATE
        ${QuEST_INCLUDE_DIR}
        ${EXTERN_DIR}/QuEST/build  # For generated config.h
    )
    target_link_libraries(bench_gate PRIVATE ${QuEST_LIBRARY})
endif()

if(QULACS_FOUND)
    target_compile_definitions(bench_gate PRIVATE WITH_QULACS)
    target_include_directories(bench_gate PRIVATE
        ${QULACS_INCLUDE_DIR}
    )
    target_link_libraries(bench_gate PRIVATE csim_static)
    # Qulacs requires C++14 minimum (we use C++17 for consistency)
    set_property(TARGET bench_gate PROPERTY CXX_STANDARD 17)
    set_property(TARGET bench_gate PROPERTY CXX_STANDARD_REQUIRED ON)
    # Enable OpenMP if available (Qulacs uses it internally)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(bench_gate PRIVATE OpenMP::OpenMP_CXX)
    endif()
endif()

if(AER_DM_FOUND)
    target_compile_definitions(bench_gate PRIVATE WITH_AER_DM)
    target_link_libraries(bench_gate PRIVATE aer_dm_static)
    # Required if Qulacs is absent (Qulacs block already sets CXX_STANDARD 17 when present)
    set_property(TARGET bench_gate PROPERTY CXX_STANDARD 17)
    set_property(TARGET bench_gate PROPERTY CXX_STANDARD_REQUIRED ON)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(bench_gate PRIVATE OpenMP::OpenMP_CXX)
    endif()
endif()

#=======================================================================================================================
# Status summary
#=======================================================================================================================

message(STATUS "Benchmark configuration:")
message(STATUS "  QuEST:     ${QuEST_FOUND}")
message(STATUS "  Qulacs:    ${QULACS_FOUND}")
message(STATUS "  Aer DM:    ${AER_DM_FOUND}")
message(STATUS "  OpenMP:    ${OpenMP_CXX_FOUND}")

#=======================================================================================================================
# bench_aer_only — Qiskit-Aer data-collection tool (temporary; delete after data collection)
#=======================================================================================================================
if(AER_DM_FOUND)
    # Same source list as bench_gate — no new source file needed.
    # BENCH_AER_ONLY_MODE=1 gates out all non-Aer method calls in bench_main.c.
    add_executable(bench_aer_only ${BENCH_MIXED_SRC})
    set_target_properties(bench_aer_only PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/benchmark"
        EXCLUDE_FROM_ALL TRUE
        CXX_STANDARD 17
        CXX_STANDARD_REQUIRED ON
    )
    target_include_directories(bench_aer_only PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}/include"
    )
    # WITH_AER_DM enables bench_aer_dm() calls; BENCH_AER_ONLY_MODE skips all other methods.
    target_compile_definitions(bench_aer_only PRIVATE WITH_AER_DM BENCH_AER_ONLY_MODE=1)
    # q is required: bench_mixed.c and bench_baselines.c reference qlib/BLAS symbols unconditionally.
    target_link_libraries(bench_aer_only PRIVATE q aer_dm_static)
    # QuEST and Qulacs omitted: without their WITH_* defines, their .c/.cpp TUs are empty.
    if(OpenMP_CXX_FOUND)
        target_link_libraries(bench_aer_only PRIVATE OpenMP::OpenMP_CXX)
    endif()
    message(STATUS "  bench_aer_only: enabled (temporary data-collection tool)")
endif()

#=======================================================================================================================
# bench_circuit - Circuit-level benchmark
#=======================================================================================================================

set(BENCH_CIRCUIT_SRC
        "${BENCH_SRC_DIR}/bench_circuit_main.c"
        "${BENCH_SRC_DIR}/bench_circuit.c"
        "${BENCH_SRC_DIR}/bench_circuit_qlib.c"
)

if(QuEST_FOUND)
    list(APPEND BENCH_CIRCUIT_SRC "${BENCH_SRC_DIR}/bench_circuit_quest.c")
endif()

if(QULACS_FOUND)
    list(APPEND BENCH_CIRCUIT_SRC "${BENCH_SRC_DIR}/bench_circuit_qulacs.cpp")
endif()

if(AER_DM_FOUND)
    list(APPEND BENCH_CIRCUIT_SRC "${BENCH_SRC_DIR}/bench_circuit_aer_dm.cpp")
endif()

add_executable(bench_circuit ${BENCH_CIRCUIT_SRC})

set_target_properties(bench_circuit PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/benchmark"
    EXCLUDE_FROM_ALL TRUE
)

target_include_directories(bench_circuit
        PRIVATE
            "${CMAKE_CURRENT_SOURCE_DIR}/include"
)

# q target provides qlib gates, state, and transitive BLAS/OpenMP dependencies.
target_link_libraries(bench_circuit
        PRIVATE
            q
)

if(QuEST_FOUND)
    target_compile_definitions(bench_circuit PRIVATE WITH_QUEST)
    target_include_directories(bench_circuit PRIVATE
        ${QuEST_INCLUDE_DIR}
        ${EXTERN_DIR}/QuEST/build  # For generated config.h
    )
    target_link_libraries(bench_circuit PRIVATE ${QuEST_LIBRARY})
endif()

if(QULACS_FOUND)
    target_compile_definitions(bench_circuit PRIVATE WITH_QULACS)
    target_include_directories(bench_circuit PRIVATE
        ${QULACS_INCLUDE_DIR}
    )
    target_link_libraries(bench_circuit PRIVATE csim_static)
    set_property(TARGET bench_circuit PROPERTY CXX_STANDARD 17)
    set_property(TARGET bench_circuit PROPERTY CXX_STANDARD_REQUIRED ON)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(bench_circuit PRIVATE OpenMP::OpenMP_CXX)
    endif()
endif()

if(AER_DM_FOUND)
    target_compile_definitions(bench_circuit PRIVATE WITH_AER_DM)
    target_link_libraries(bench_circuit PRIVATE aer_dm_static)
    set_property(TARGET bench_circuit PROPERTY CXX_STANDARD 17)
    set_property(TARGET bench_circuit PROPERTY CXX_STANDARD_REQUIRED ON)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(bench_circuit PRIVATE OpenMP::OpenMP_CXX)
    endif()
endif()

message(STATUS "  bench_circuit: enabled (circuit-level benchmarks)")
