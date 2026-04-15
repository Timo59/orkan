# src/CMakeLists.txt - Main library target

# Collect state module source files
set(STATE_SOURCES
        "${CMAKE_CURRENT_SOURCE_DIR}/state/state.c"
        "${CMAKE_CURRENT_SOURCE_DIR}/state/state_pure.c"
        "${CMAKE_CURRENT_SOURCE_DIR}/state/state_packed.c"
        "${CMAKE_CURRENT_SOURCE_DIR}/state/state_tiled.c"
)

# Collect channel module source files
set(CHANNEL_SOURCES
        "${CMAKE_CURRENT_SOURCE_DIR}/channel/channel.c"
        "${CMAKE_CURRENT_SOURCE_DIR}/channel/channel_packed.c"
        "${CMAKE_CURRENT_SOURCE_DIR}/channel/channel_tiled.c"
)

# Collect measurement module source files
set(MEAS_SOURCES
        "${CMAKE_CURRENT_SOURCE_DIR}/meas/meas.c"
        "${CMAKE_CURRENT_SOURCE_DIR}/meas/meas_pure.c"
        "${CMAKE_CURRENT_SOURCE_DIR}/meas/meas_packed.c"
        "${CMAKE_CURRENT_SOURCE_DIR}/meas/meas_tiled.c"
)

# Collect circuit module source files
set(CIRC_SOURCES
        "${CMAKE_CURRENT_SOURCE_DIR}/circ/circ.c"
        "${CMAKE_CURRENT_SOURCE_DIR}/circ/circ_pure.c"
        "${CMAKE_CURRENT_SOURCE_DIR}/circ/circ_packed.c"
        "${CMAKE_CURRENT_SOURCE_DIR}/circ/circ_tiled.c"
)

# Collect packed gate sources (one file per gate)
set(GATE_PACKED_SOURCES
        "${CMAKE_CURRENT_SOURCE_DIR}/gate/packed/gate_packed_1q.c"
        "${CMAKE_CURRENT_SOURCE_DIR}/gate/packed/gate_packed_ccx.c"
        "${CMAKE_CURRENT_SOURCE_DIR}/gate/packed/gate_packed_cx.c"
        "${CMAKE_CURRENT_SOURCE_DIR}/gate/packed/gate_packed_cy.c"
        "${CMAKE_CURRENT_SOURCE_DIR}/gate/packed/gate_packed_cz.c"
        "${CMAKE_CURRENT_SOURCE_DIR}/gate/packed/gate_packed_swap.c"
)

# Collect tiled gate sources (one file per gate)
set(GATE_TILED_SOURCES
        "${CMAKE_CURRENT_SOURCE_DIR}/gate/tiled/gate_tiled_1q.c"
        "${CMAKE_CURRENT_SOURCE_DIR}/gate/tiled/gate_tiled_2q.c"
        "${CMAKE_CURRENT_SOURCE_DIR}/gate/tiled/gate_tiled_ccx.c"
        "${CMAKE_CURRENT_SOURCE_DIR}/gate/tiled/gate_tiled_cx.c"
        "${CMAKE_CURRENT_SOURCE_DIR}/gate/tiled/gate_tiled_cy.c"
        "${CMAKE_CURRENT_SOURCE_DIR}/gate/tiled/gate_tiled_cz.c"
        "${CMAKE_CURRENT_SOURCE_DIR}/gate/tiled/gate_tiled_h.c"
        "${CMAKE_CURRENT_SOURCE_DIR}/gate/tiled/gate_tiled_hy.c"
        "${CMAKE_CURRENT_SOURCE_DIR}/gate/tiled/gate_tiled_p.c"
        "${CMAKE_CURRENT_SOURCE_DIR}/gate/tiled/gate_tiled_rx.c"
        "${CMAKE_CURRENT_SOURCE_DIR}/gate/tiled/gate_tiled_ry.c"
        "${CMAKE_CURRENT_SOURCE_DIR}/gate/tiled/gate_tiled_rz.c"
        "${CMAKE_CURRENT_SOURCE_DIR}/gate/tiled/gate_tiled_s.c"
        "${CMAKE_CURRENT_SOURCE_DIR}/gate/tiled/gate_tiled_sdg.c"
        "${CMAKE_CURRENT_SOURCE_DIR}/gate/tiled/gate_tiled_swap.c"
        "${CMAKE_CURRENT_SOURCE_DIR}/gate/tiled/gate_tiled_t.c"
        "${CMAKE_CURRENT_SOURCE_DIR}/gate/tiled/gate_tiled_tdg.c"
        "${CMAKE_CURRENT_SOURCE_DIR}/gate/tiled/gate_tiled_x.c"
        "${CMAKE_CURRENT_SOURCE_DIR}/gate/tiled/gate_tiled_y.c"
        "${CMAKE_CURRENT_SOURCE_DIR}/gate/tiled/gate_tiled_z.c"
)

# Collect all source files
set(Q_SOURCES
        ${STATE_SOURCES}
        ${CHANNEL_SOURCES}
        ${MEAS_SOURCES}
        ${CIRC_SOURCES}
        "${CMAKE_CURRENT_SOURCE_DIR}/gate/gate.c"
        "${CMAKE_CURRENT_SOURCE_DIR}/gate/gate_pure.c"
        ${GATE_PACKED_SOURCES}
        ${GATE_TILED_SOURCES}
)

# Include directory
set(INCLUDE_DIR "${CMAKE_SOURCE_DIR}/include")

include(GNUInstallDirs)

# Create library target
add_library(${PROJECT_NAME_LOWER} SHARED ${Q_SOURCES})

set_target_properties(${PROJECT_NAME_LOWER} PROPERTIES
        VERSION         ${PROJECT_VERSION}
        SOVERSION       ${PROJECT_VERSION_MAJOR}
)
target_include_directories(${PROJECT_NAME_LOWER}
        PUBLIC
            $<BUILD_INTERFACE:${INCLUDE_DIR}>
            $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
        PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}
            "${CMAKE_CURRENT_SOURCE_DIR}/internal"
            "${CMAKE_CURRENT_SOURCE_DIR}/gate/tiled"
)
target_compile_definitions(${PROJECT_NAME_LOWER}
        PUBLIC
            LOG_TILE_DIM=${LOG_TILE_DIM}
)
if(DEFINED OMP_THRESHOLD)
    target_compile_definitions(${PROJECT_NAME_LOWER} PRIVATE OMP_THRESHOLD=${OMP_THRESHOLD})
endif()
target_link_libraries(${PROJECT_NAME_LOWER}
        PRIVATE
            ${PROJECT_COMPILER_FLAGS}
)
if(ENABLE_OPENMP AND TARGET OpenMP::OpenMP_C)
    target_link_libraries(${PROJECT_NAME_LOWER} PUBLIC OpenMP::OpenMP_C)
    # FindOpenMP does not always propagate the include directory on macOS
    if(APPLE AND LIBOMP_PREFIX)
        target_include_directories(${PROJECT_NAME_LOWER} PUBLIC "${LIBOMP_PREFIX}/include")
    endif()
endif()

# Install rules
include(CMakePackageConfigHelpers)

# Install library and register export set
install(TARGETS ${PROJECT_NAME_LOWER}
    EXPORT ${PROJECT_NAME}Targets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

# Install public headers
install(FILES
    ${INCLUDE_DIR}/qlib.h
    ${INCLUDE_DIR}/q_types.h
    ${INCLUDE_DIR}/state.h
    ${INCLUDE_DIR}/gate.h
    ${INCLUDE_DIR}/channel.h
    ${INCLUDE_DIR}/meas.h
    ${INCLUDE_DIR}/circ.h
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME_LOWER}
)

# Install CMake targets file (provides Orkan::orkan imported target)
install(EXPORT ${PROJECT_NAME}Targets
    NAMESPACE ${PROJECT_NAME}::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME_LOWER}
)

# Generate and install package config file
configure_package_config_file(
    "${CMAKE_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME_LOWER}
)

# Generate and install version file
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME_LOWER}
)

