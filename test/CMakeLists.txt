# Initialize submodule if not present
if (NOT EXISTS "${CMAKE_SOURCE_DIR}/extern/Unity/src/unity.h")
    find_package(Git QUIET)
    if(GIT_FOUND)
        message(STATUS "Initializing Unity submodule...")
        execute_process(
                COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive extern/Unity
                WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
                RESULT_VARIABLE GIT_SUBMOD_RESULT
        )
        if(NOT GIT_SUBMOD_RESULT EQUAL "0")
            message(FATAL_ERROR "git submodule update --init failed with ${GIT_SUBMOD_RESULT}")
        endif()
    else()
        message(FATAL_ERROR "Git not found. Cannot initialize Unity submodule")
    endif()
endif()

# Find unity source directory
find_path(UNITY_SRC
        NAMES unity.h
        PATHS "${CMAKE_SOURCE_DIR}/extern/Unity/src"
        NO_DEFAULT_PATH
)
if (NOT UNITY_SRC)
    message(FATAL_ERROR "Unity source directory not found.")
endif ()

# Set source directories
set(SRC_DIR "${CMAKE_SOURCE_DIR}/src")
set(TEST_SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")

#=======================================================================================================================
# test_state
#=======================================================================================================================

# State module source files
set(STATE_SOURCES
        "${SRC_DIR}/state.c"
        "${SRC_DIR}/state_pure.c"
        "${SRC_DIR}/state_packed.c"
        "${SRC_DIR}/state_tiled.c"
)

# Source files for test_state
set(TEST_STATE_SRC
        ${STATE_SOURCES}
        "${TEST_SRC_DIR}/test_state.c"
        "${TEST_SRC_DIR}/test_state_pure.c"
        "${TEST_SRC_DIR}/test_state_packed.c"
        "${TEST_SRC_DIR}/test_state_tiled.c"
        "${UNITY_SRC}/unity.c"
)

add_executable(test_state ${TEST_STATE_SRC})
target_compile_definitions(test_state
        PRIVATE
            UNITY_INCLUDE_DOUBLE
            UNITY_SUPPORT_64
            LOG_TILE_DIM=${LOG_TILE_DIM}
)
target_include_directories(test_state
        PRIVATE
            ${UNITY_SRC}
            "${CMAKE_SOURCE_DIR}/include"
            "${CMAKE_CURRENT_SOURCE_DIR}/include"
)
target_link_libraries(test_state
        PRIVATE
            blas_compiler_flags
            ${PROJECT_COMPILER_FLAGS}
)
add_test(NAME test_state COMMAND test_state)

#=======================================================================================================================
# test_gate
#=======================================================================================================================

# Source files required for test_gate
set(TEST_GATE_SRC
        ${STATE_SOURCES}
        "${SRC_DIR}/gate.c"
        "${SRC_DIR}/qhipster.c"
        "${SRC_DIR}/mhipster.c"
        "${TEST_SRC_DIR}/test_gate.c"
        "${TEST_SRC_DIR}/test_qhipster.c"
        "${TEST_SRC_DIR}/test_mhipster.c"
        "${TEST_SRC_DIR}/test_pure_states.c"
        "${TEST_SRC_DIR}/test_mixed_states.c"
        "${TEST_SRC_DIR}/gatemat.c"
        "${TEST_SRC_DIR}/linalg.c"
        "${UNITY_SRC}/unity.c"
)

add_executable(test_gate ${TEST_GATE_SRC})
target_compile_definitions(test_gate
        PRIVATE
            UNITY_INCLUDE_DOUBLE
            UNITY_SUPPORT_64
            LOG_TILE_DIM=${LOG_TILE_DIM}
)
target_include_directories(test_gate
        PRIVATE
        ${UNITY_SRC}
            "${CMAKE_SOURCE_DIR}/include"
            "${CMAKE_CURRENT_SOURCE_DIR}/include"
)
target_link_libraries(test_gate
        PRIVATE
            blas_compiler_flags
            omp_compiler_flags
            ${PROJECT_COMPILER_FLAGS}
)

