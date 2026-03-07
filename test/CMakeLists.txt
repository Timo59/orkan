include(FetchContent)
FetchContent_Declare(
    unity
    GIT_REPOSITORY https://github.com/ThrowTheSwitch/Unity.git
    GIT_TAG        v2.6.1
)
FetchContent_MakeAvailable(unity)
set(UNITY_SRC "${unity_SOURCE_DIR}/src")

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
        "${GATE_DIR}/test_gate_pure_3q.c"
        "${GATE_DIR}/test_gate_packed_1q.c"
        "${GATE_DIR}/test_gate_packed_2q.c"
        "${GATE_DIR}/test_gate_packed_3q.c"
        "${GATE_DIR}/test_gate_tiled_1q.c"
        "${GATE_DIR}/test_gate_tiled_2q.c"
        "${GATE_DIR}/test_gate_tiled_2q_mat.c"
        "${GATE_DIR}/test_gate_tiled_3q.c"
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

#=======================================================================================================================
# test_circuit
#=======================================================================================================================

set(CIRCUIT_DIR "${CMAKE_CURRENT_SOURCE_DIR}/circuit")

set(TEST_CIRCUIT_SRC
        "${CIRCUIT_DIR}/test_circuit.c"
        "${CIRCUIT_DIR}/test_pauli.c"
        "${UTILITY_DIR}/linalg.c"
        "${UTILITY_DIR}/gatemat.c"
        "${UTILITY_DIR}/test_pure_states.c"
        "${UNITY_SRC}/unity.c"
)

add_executable(test_circuit ${TEST_CIRCUIT_SRC})
target_compile_definitions(test_circuit
        PRIVATE
            UNITY_INCLUDE_DOUBLE
            UNITY_SUPPORT_64
            LOG_TILE_DIM=${LOG_TILE_DIM}
)
if(DEFINED OMP_THRESHOLD)
    target_compile_definitions(test_circuit PRIVATE OMP_THRESHOLD=${OMP_THRESHOLD})
endif()
target_include_directories(test_circuit
        PRIVATE
            ${UNITY_SRC}
            "${CMAKE_CURRENT_SOURCE_DIR}/include"
)
target_link_libraries(test_circuit
        PRIVATE
            q
)
add_test(NAME test_circuit COMMAND test_circuit)
