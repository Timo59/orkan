cmake_minimum_required(VERSION 3.27)

# Set the name of the project and store it to CMAKE_PROJECT_NAME. Set PROJECT_SOURCE_DIR and PROJECT_BINARY_DIR
project(QSim
        VERSION     1.0.0   # Sets PROJECT_VERSION and <PROJECT_NAME>_VERSION
        DESCRIPTION "Quantum Computer simulator" # Sets PROJECT_DESCRIPTION and <PROJECT_NAME>_DESCRIPTION
        LANGUAGES   C CXX   # Programming languages to build the project
)

set(CMAKE_C_STANDARD 17)
set(CMAKE_CXX_STANDARD 17)

# Add CMake modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# Options
option(BUILD_EXAMPLES "Build examples" OFF)
option(BUILD_TESTS "Build tests" ON)
option(BUILD_BENCHMARKS "Build benchmarks" OFF)
option(ENABLE_ASAN "Enable AddressSanitizer for memory debugging" OFF)
option(ENABLE_OPENMP "Enable OpenMP parallelization" ON)

# Tile dimension configuration (2^LOG_TILE_DIM x 2^LOG_TILE_DIM tiles)
set(LOG_TILE_DIM 5 CACHE STRING "Log2 of tile dimension (default: 5 = 32x32)")

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

# Find dependencies
include(cmake/Dependencies.cmake)

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
            DEPENDS test_state test_gate
            COMMENT "Running tests before install..."
    )
endif ()

if(BUILD_EXAMPLES AND EXISTS "${CMAKE_SOURCE_DIR}/examples")
    add_subdirectory(examples)
endif()

if(EXISTS "${CMAKE_SOURCE_DIR}/benchmark")
    add_subdirectory(benchmark)
endif()
