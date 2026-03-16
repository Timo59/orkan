# benchmark/CMakeLists.txt — Gate benchmark executable

# Benchmarks require Release build for meaningful results
if(NOT CMAKE_BUILD_TYPE STREQUAL "Release")
    message(STATUS "Skipping benchmark (requires -DCMAKE_BUILD_TYPE=Release)")
    return()
endif()

set(BENCH_SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")
set(EXTERN_DIR "${CMAKE_SOURCE_DIR}/extern")

# ======================================================================
# Find external frameworks (auto-detect)
# ======================================================================

# QuEST — check extern/ first, then system paths
find_path(QuEST_INCLUDE_DIR quest/include/quest.h
    HINTS ${EXTERN_DIR}/QuEST ${QuEST_DIR} $ENV{QuEST_DIR}
          /usr/local /opt/quest)
find_library(QuEST_LIBRARY NAMES QuEST quest
    HINTS ${EXTERN_DIR}/QuEST/build ${QuEST_DIR}/lib $ENV{QuEST_DIR}/lib
          /usr/local/lib /opt/quest/lib)
if(QuEST_INCLUDE_DIR AND QuEST_LIBRARY)
    set(QuEST_FOUND TRUE)
endif()

# OpenMP — macOS needs Homebrew libomp hints
if(APPLE)
    execute_process(COMMAND brew --prefix libomp
        OUTPUT_VARIABLE LIBOMP_PREFIX OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
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

# Qulacs — build csim_static from submodule sources
set(QULACS_CSIM_DIR "${EXTERN_DIR}/qulacs/src/csim")
if(EXISTS "${QULACS_CSIM_DIR}/update_ops_dm.hpp")
    file(GLOB QULACS_CSIM_SRC "${QULACS_CSIM_DIR}/*.cpp")
    add_library(csim_static STATIC ${QULACS_CSIM_SRC})
    find_path(EIGEN3_INCLUDE_DIR Eigen/Core
        HINTS /opt/homebrew/include/eigen3 /usr/local/include/eigen3 /usr/include/eigen3)
    if(EIGEN3_INCLUDE_DIR)
        target_include_directories(csim_static SYSTEM PRIVATE ${EIGEN3_INCLUDE_DIR})
    endif()
    target_include_directories(csim_static PRIVATE "${EXTERN_DIR}/qulacs/src")
    target_compile_options(csim_static PRIVATE -w $<$<CONFIG:Release>:-march=native>)
    set_property(TARGET csim_static PROPERTY CXX_STANDARD 17)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(csim_static PRIVATE OpenMP::OpenMP_CXX)
    endif()
    set(QULACS_INCLUDE_DIR "${EXTERN_DIR}/qulacs/src")
    set(QULACS_FOUND TRUE)
endif()

# Qiskit-Aer DM — header-only kernel vendored to extern/aer-dm/
set(AER_DM_INCLUDE_DIR "${EXTERN_DIR}/aer-dm/include")
if(EXISTS "${AER_DM_INCLUDE_DIR}/simulators/density_matrix/densitymatrix.hpp")
    set(AER_DM_FOUND TRUE)
endif()

if(APPLE)
    execute_process(COMMAND brew --prefix nlohmann-json
        OUTPUT_VARIABLE NLOHMANN_JSON_PREFIX OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
endif()
find_path(NLOHMANN_JSON_INCLUDE_DIR nlohmann/json.hpp
    HINTS "${NLOHMANN_JSON_PREFIX}/include" /opt/homebrew/include /usr/local/include /usr/include)

if(AER_DM_FOUND)
    add_library(aer_dm_static INTERFACE)
    target_include_directories(aer_dm_static SYSTEM INTERFACE "${AER_DM_INCLUDE_DIR}")
    if(NLOHMANN_JSON_INCLUDE_DIR)
        target_include_directories(aer_dm_static SYSTEM INTERFACE "${NLOHMANN_JSON_INCLUDE_DIR}")
    else()
        message(WARNING "nlohmann/json.hpp not found — Aer DM headers will fail to compile.")
    endif()
endif()

# ======================================================================
# bench_gate — Single-gate throughput benchmark
# ======================================================================

set(BENCH_GATE_SRC
    "${BENCH_SRC_DIR}/bench_main.c"
    "${BENCH_SRC_DIR}/bench_run.c"
    "${BENCH_SRC_DIR}/bench_result.c"
    "${BENCH_SRC_DIR}/bench_util.c"
    "${BENCH_SRC_DIR}/bench_qlib.c"
    "${BENCH_SRC_DIR}/bench_blas.c"
    "${BENCH_SRC_DIR}/bench_console.c"
    "${BENCH_SRC_DIR}/bench_csv.c"
    "${BENCH_SRC_DIR}/bench_pgfplots.c"
)

if(QuEST_FOUND)
    list(APPEND BENCH_GATE_SRC "${BENCH_SRC_DIR}/bench_quest.c")
endif()
if(QULACS_FOUND)
    list(APPEND BENCH_GATE_SRC "${BENCH_SRC_DIR}/bench_qulacs.cpp")
endif()
if(AER_DM_FOUND)
    list(APPEND BENCH_GATE_SRC "${BENCH_SRC_DIR}/bench_aer.cpp")
endif()

add_executable(bench_gate ${BENCH_GATE_SRC})
set_target_properties(bench_gate PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/benchmark"
    EXCLUDE_FROM_ALL TRUE)
target_include_directories(bench_gate PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(bench_gate PRIVATE q)

if(QuEST_FOUND)
    target_compile_definitions(bench_gate PRIVATE HAVE_QUEST)
    target_include_directories(bench_gate PRIVATE ${QuEST_INCLUDE_DIR} ${EXTERN_DIR}/QuEST/build)
    target_link_libraries(bench_gate PRIVATE ${QuEST_LIBRARY})
endif()

if(QULACS_FOUND)
    target_compile_definitions(bench_gate PRIVATE HAVE_QULACS)
    target_include_directories(bench_gate PRIVATE ${QULACS_INCLUDE_DIR})
    target_link_libraries(bench_gate PRIVATE csim_static)
    set_property(TARGET bench_gate PROPERTY CXX_STANDARD 17)
    set_property(TARGET bench_gate PROPERTY CXX_STANDARD_REQUIRED ON)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(bench_gate PRIVATE OpenMP::OpenMP_CXX)
    endif()
endif()

if(AER_DM_FOUND)
    target_compile_definitions(bench_gate PRIVATE HAVE_AER)
    target_link_libraries(bench_gate PRIVATE aer_dm_static)
    set_property(TARGET bench_gate PROPERTY CXX_STANDARD 17)
    set_property(TARGET bench_gate PROPERTY CXX_STANDARD_REQUIRED ON)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(bench_gate PRIVATE OpenMP::OpenMP_CXX)
    endif()
endif()

# ======================================================================
# Status summary
# ======================================================================

message(STATUS "Benchmark configuration:")
message(STATUS "  QuEST:     ${QuEST_FOUND}")
message(STATUS "  Qulacs:    ${QULACS_FOUND}")
message(STATUS "  Aer DM:    ${AER_DM_FOUND}")
message(STATUS "  OpenMP:    ${OpenMP_CXX_FOUND}")
