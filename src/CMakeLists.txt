# src/CMakeLists.txt - Main library target

# Collect state module source files
set(STATE_SOURCES
        "${CMAKE_CURRENT_SOURCE_DIR}/state.c"
        "${CMAKE_CURRENT_SOURCE_DIR}/state_pure.c"
        "${CMAKE_CURRENT_SOURCE_DIR}/state_packed.c"
        "${CMAKE_CURRENT_SOURCE_DIR}/state_tiled.c"
)

# Collect all source files
set(Q_SOURCES
        ${STATE_SOURCES}
        "${CMAKE_CURRENT_SOURCE_DIR}/gate.c"
        "${CMAKE_CURRENT_SOURCE_DIR}/gate_pure.c"
        "${CMAKE_CURRENT_SOURCE_DIR}/gate_packed.c"
        "${CMAKE_CURRENT_SOURCE_DIR}/gate_tiled.c"
)

# Include directory
set(INCLUDE_DIR "${CMAKE_SOURCE_DIR}/include")

# Create library target
add_library(q SHARED ${Q_SOURCES})

set_target_properties(q PROPERTIES
        VERSION         ${PROJECT_VERSION}
        SOVERSION       ${PROJECT_VERSION_MAJOR}
        PUBLIC_HEADER   "${INCLUDE_DIR}/qlib.h"
)
target_include_directories(q
        PUBLIC
            $<BUILD_INTERFACE:${INCLUDE_DIR}>
            $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDE_DIR}>
        PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}
)
target_compile_definitions(q
        PRIVATE
            LOG_TILE_DIM=${LOG_TILE_DIM}
)
target_link_libraries(q
        PUBLIC
            ${PROJECT_COMPILER_FLAGS}
            blas_compiler_flags
            omp_compiler_flags
)