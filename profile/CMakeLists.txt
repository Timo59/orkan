# profile/CMakeLists.txt - Micro-profiling executables for gate functions
#
# Profiling builds use Release optimization (-O3 -march=native) to measure
# production-representative performance, but can also be built in Debug
# to exercise the small-tile (8x8) code paths.

set(PROF_SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")

# --- profile_x_tiled ---
add_executable(profile_x_tiled "${PROF_SRC_DIR}/profile_x_tiled.c")

set_target_properties(profile_x_tiled PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/profile"
)

target_include_directories(profile_x_tiled PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
)

# Link libq — inherits compiler flags, BLAS, OpenMP, and include dirs via PUBLIC deps
target_link_libraries(profile_x_tiled PRIVATE ${PROJECT_NAME_LOWER})

# --- profile_h_tiled ---
add_executable(profile_h_tiled "${PROF_SRC_DIR}/profile_h_tiled.c")

set_target_properties(profile_h_tiled PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/profile"
)

target_include_directories(profile_h_tiled PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
)

target_link_libraries(profile_h_tiled PRIVATE ${PROJECT_NAME_LOWER})

# --- profile_channel_tiled ---
add_executable(profile_channel_tiled "${PROF_SRC_DIR}/profile_channel_tiled.c")

set_target_properties(profile_channel_tiled PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/profile"
)

target_include_directories(profile_channel_tiled PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
)

target_link_libraries(profile_channel_tiled PRIVATE ${PROJECT_NAME_LOWER})
