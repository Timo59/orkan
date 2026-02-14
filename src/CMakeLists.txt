# src/CMakeLists.txt - Main library target

# Collect state module source files
set(STATE_SOURCES
        "${CMAKE_CURRENT_SOURCE_DIR}/state/state.c"
        "${CMAKE_CURRENT_SOURCE_DIR}/state/state_pure.c"
        "${CMAKE_CURRENT_SOURCE_DIR}/state/state_packed.c"
        "${CMAKE_CURRENT_SOURCE_DIR}/state/state_tiled.c"
)

# Collect packed gate sources (one file per gate)
file(GLOB GATE_PACKED_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/gate/packed/gate_packed_*.c")

# Collect tiled gate sources
file(GLOB GATE_TILED_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/gate/tiled/gate_tiled_*.c")

# Collect all source files
set(Q_SOURCES
        ${STATE_SOURCES}
        "${CMAKE_CURRENT_SOURCE_DIR}/gate/gate.c"
        "${CMAKE_CURRENT_SOURCE_DIR}/gate/gate_pure.c"
        ${GATE_PACKED_SOURCES}
        ${GATE_TILED_SOURCES}
)

# Include directory
set(INCLUDE_DIR "${CMAKE_SOURCE_DIR}/include")

# Create library target
add_library(q SHARED ${Q_SOURCES})

set_target_properties(q PROPERTIES
        VERSION         ${PROJECT_VERSION}
        SOVERSION       ${PROJECT_VERSION_MAJOR}
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

# Install rules
include(GNUInstallDirs)

install(TARGETS q
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

install(FILES
    ${INCLUDE_DIR}/qlib.h
    ${INCLUDE_DIR}/q_types.h
    ${INCLUDE_DIR}/state.h
    ${INCLUDE_DIR}/gate.h
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)