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

# Collect packed gate sources (one file per gate)
file(GLOB GATE_PACKED_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/gate/packed/gate_packed_*.c")

# Collect tiled gate sources
file(GLOB GATE_TILED_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/gate/tiled/gate_tiled_*.c")

# Collect all source files
set(Q_SOURCES
        ${STATE_SOURCES}
        ${CHANNEL_SOURCES}
        "${CMAKE_CURRENT_SOURCE_DIR}/gate/gate.c"
        "${CMAKE_CURRENT_SOURCE_DIR}/gate/gate_pure.c"
        ${GATE_PACKED_SOURCES}
        ${GATE_TILED_SOURCES}
        "${CMAKE_CURRENT_SOURCE_DIR}/circ/pauli.c"
        "${CMAKE_CURRENT_SOURCE_DIR}/circ/pauli_pure.c"
)

# Include directory
set(INCLUDE_DIR "${CMAKE_SOURCE_DIR}/include")

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
            omp_compiler_flags
)

# For the installed/exported target, link OpenMP directly (imported targets
# are resolved by find_dependency in the package config file)
if(ENABLE_OPENMP AND TARGET OpenMP::OpenMP_C)
    target_link_libraries(${PROJECT_NAME_LOWER} PUBLIC OpenMP::OpenMP_C)
endif()

# Install rules
include(GNUInstallDirs)
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
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME_LOWER}
)

# Install CMake targets file (provides QSim::qsim imported target)
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

