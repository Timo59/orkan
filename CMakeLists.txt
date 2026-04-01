cmake_minimum_required(VERSION 3.27)

# Set the name of the project and store it to CMAKE_PROJECT_NAME. Set PROJECT_SOURCE_DIR and PROJECT_BINARY_DIR
project(QSim
        VERSION     1.0.0   # Sets PROJECT_VERSION and <PROJECT_NAME>_VERSION
        DESCRIPTION "Quantum Computer simulator" # Sets PROJECT_DESCRIPTION and <PROJECT_NAME>_DESCRIPTION
        LANGUAGES   C CXX   # Programming languages to build the project
)

# Default to Debug if no build type specified (single-config generators only)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type (Debug, Release, RelWithDebInfo)" FORCE)
    message(STATUS "No build type specified, defaulting to Release")
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# Guard: catch mismatches between build directory name and build type.
# Prevents accidental Release builds inside cmake-build-debug (and vice versa).
get_filename_component(_BUILD_DIR_NAME "${CMAKE_BINARY_DIR}" NAME)
if(_BUILD_DIR_NAME MATCHES "[Dd]ebug" AND NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(FATAL_ERROR
        "Build directory '${_BUILD_DIR_NAME}' implies Debug "
        "but CMAKE_BUILD_TYPE='${CMAKE_BUILD_TYPE}'.\n"
        "Run: cmake --preset debug")
endif()
if(_BUILD_DIR_NAME MATCHES "[Rr]elease" AND NOT CMAKE_BUILD_TYPE STREQUAL "Release")
    message(FATAL_ERROR
        "Build directory '${_BUILD_DIR_NAME}' implies Release "
        "but CMAKE_BUILD_TYPE='${CMAKE_BUILD_TYPE}'.\n"
        "Run: cmake --preset release")
endif()
unset(_BUILD_DIR_NAME)

set(CMAKE_C_STANDARD 17)
set(CMAKE_CXX_STANDARD 17)

# Lowercase project name for target and install paths
string(TOLOWER "${PROJECT_NAME}" PROJECT_NAME_LOWER)

# Add CMake modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# Options
option(BUILD_TESTS "Build tests" OFF)
option(BUILD_BENCHMARKS "Build benchmarks" OFF)
option(ENABLE_ASAN "Enable AddressSanitizer for memory debugging" OFF)
option(ENABLE_OPENMP "Enable OpenMP parallelization" ON)

# Tile dimension configuration (2^LOG_TILE_DIM x 2^LOG_TILE_DIM tiles)
# Debug uses small tiles (8x8) so multi-tile code paths are exercised with few qubits.
# Release uses production tiles (32x32) for cache-optimal performance.
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(LOG_TILE_DIM 3 CACHE STRING "Log2 of tile dimension (default: 3 for Debug, 5 for Release)")
    set(OMP_THRESHOLD 4 CACHE STRING "Min dim for OpenMP (debug override; Release uses per-file defaults)")
    message(STATUS "OMP_THRESHOLD: ${OMP_THRESHOLD} (debug override)")
else()
    set(LOG_TILE_DIM 5 CACHE STRING "Log2 of tile dimension (default: 3 for Debug, 5 for Release)")
endif()

if(LOG_TILE_DIM LESS 3 OR LOG_TILE_DIM GREATER 8)
    message(FATAL_ERROR "LOG_TILE_DIM must be between 3 and 8")
endif()

math(EXPR TILE_DIM "1 << ${LOG_TILE_DIM}")
message(STATUS "Tile dimension: ${TILE_DIM}x${TILE_DIM}")

if(ENABLE_ASAN)
    message(STATUS "AddressSanitizer enabled")
    set(SANITIZE_FLAGS "-fsanitize=address,undefined -fno-omit-frame-pointer")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${SANITIZE_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SANITIZE_FLAGS}")
    # Optional: add for linker flags too (needed on some Linux/macOS setups)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${SANITIZE_FLAGS}")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${SANITIZE_FLAGS}")
endif()

# Platform detection and compiler setup
include(cmake/PlatformConfig.cmake)

# Find dependencies (BLAS only needed for tests/benchmarks, not the library itself)
if(BUILD_TESTS OR BUILD_BENCHMARKS)
    include(cmake/Dependencies.cmake)
endif()

# Create compiler flags interface library
include(cmake/CompilerFlags.cmake)

# Build main library
add_subdirectory(src)

if (BUILD_TESTS AND EXISTS "${CMAKE_SOURCE_DIR}/test")
    enable_testing()
    add_subdirectory(test)

    # Run tests then install (usage: cmake --build . --target verified_install)
    add_custom_target(verified_install
            COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure
            COMMAND ${CMAKE_COMMAND} --install ${CMAKE_BINARY_DIR}
            DEPENDS ${PROJECT_NAME_LOWER} test_state test_gate test_channel
            COMMENT "Running tests before install..."
    )
endif ()

if(BUILD_BENCHMARKS AND EXISTS "${CMAKE_SOURCE_DIR}/benchmark")
    add_subdirectory(benchmark)
endif()

if(EXISTS "${CMAKE_SOURCE_DIR}/profile")
    add_subdirectory(profile)
endif()
