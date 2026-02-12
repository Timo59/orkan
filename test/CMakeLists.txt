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
set(STATE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/state")
set(GATE_DIR  "${CMAKE_CURRENT_SOURCE_DIR}/gate")
set(UTILITY_DIR "${CMAKE_CURRENT_SOURCE_DIR}/utility")

#=======================================================================================================================
# test_state
#=======================================================================================================================

set(TEST_STATE_SRC
        "${STATE_DIR}/test_state.c"
        "${STATE_DIR}/test_state_pure.c"
        "${STATE_DIR}/test_state_packed.c"
        "${STATE_DIR}/test_state_tiled.c"
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
            "${CMAKE_CURRENT_SOURCE_DIR}/include"
)
target_link_libraries(test_state
        PRIVATE
            q
)
add_test(NAME test_state COMMAND test_state)

#=======================================================================================================================
# test_gate
#=======================================================================================================================

set(TEST_GATE_SRC
        "${GATE_DIR}/test_gate.c"
        "${GATE_DIR}/test_gate_pure_1q.c"
        "${GATE_DIR}/test_gate_pure_2q.c"
        "${GATE_DIR}/test_gate_packed_1q.c"
        "${GATE_DIR}/test_gate_packed_2q.c"
        "${GATE_DIR}/test_gate_packed_3q.c"
        "${GATE_DIR}/test_gate_tiled_1q.c"
        "${GATE_DIR}/test_gate_tiled_2q.c"
        "${UTILITY_DIR}/test_pure_states.c"
        "${UTILITY_DIR}/test_mixed_states.c"
        "${UTILITY_DIR}/test_mixed_utils.c"
        "${UTILITY_DIR}/gatemat.c"
        "${UTILITY_DIR}/linalg.c"
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
            "${CMAKE_CURRENT_SOURCE_DIR}/include"
)
target_link_libraries(test_gate
        PRIVATE
            q
)
add_test(NAME test_gate COMMAND test_gate)
