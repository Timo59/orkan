/**
 * @file test_state.c
 * @brief Test runner for all state module tests
 *
 * This test runner executes comprehensive unit tests for the state module,
 * covering all three state types:
 * - PURE: Pure quantum states (statevectors)
 * - MIXED_PACKED: Mixed states with packed lower-triangular storage
 * - MIXED_TILED: Mixed states with tiled/blocked storage for cache efficiency
 */

#include "test.h"
#include "state.h"

/*
 * =====================================================================================================================
 * Test function declarations from test_state_pure.c
 * =====================================================================================================================
 */

void test_pure_len(void);
void test_pure_init_null_data(void);
void test_pure_init_ownership(void);
void test_pure_free_double_free(void);
void test_pure_plus_values(void);
void test_pure_get_set_roundtrip(void);
void test_pure_cp_independence(void);
void test_pure_single_qubit(void);
void test_pure_four_qubits(void);
void test_pure_reinit(void);

/*
 * =====================================================================================================================
 * Test function declarations from test_state_packed.c
 * =====================================================================================================================
 */

void test_packed_len(void);
void test_packed_init_null_data(void);
void test_packed_init_ownership(void);
void test_packed_free_double_free(void);
void test_packed_plus_values(void);
void test_packed_get_lower(void);
void test_packed_get_upper(void);
void test_packed_set_lower(void);
void test_packed_set_upper(void);
void test_packed_cp_independence(void);
void test_packed_hermitian_symmetry(void);
void test_packed_single_qubit(void);
void test_packed_zero_qubit(void);

/*
 * =====================================================================================================================
 * Test function declarations from test_state_tiled.c
 * =====================================================================================================================
 */

void test_tiled_len_small(void);
void test_tiled_len_8_qubits(void);
void test_tiled_init_null_data(void);
void test_tiled_init_ownership(void);
void test_tiled_free_double_free(void);
void test_tiled_plus_values(void);
void test_tiled_get_within_tile(void);
void test_tiled_get_cross_tile(void);
void test_tiled_get_upper(void);
void test_tiled_set_roundtrip(void);
void test_tiled_cp_independence(void);
void test_tiled_hermitian_symmetry(void);
void test_tiled_single_qubit(void);
void test_tiled_tile_boundaries(void);
void test_tiled_set_upper(void);
void test_tiled_6_qubits(void);

/*
 * =====================================================================================================================
 * Dispatch layer tests (implemented in this file)
 * =====================================================================================================================
 */

/*
 * Test: Dispatch API roundtrip - test all three state types through dispatch layer
 */
void test_dispatch_roundtrip(void) {
    // Test PURE state through dispatch API
    state_t pure_state = {.type = PURE, .data = NULL, .qubits = 0};
    state_plus(&pure_state, 2);

    // Verify plus state: all elements = 1/sqrt(4) = 0.5
    TEST_ASSERT_COMPLEX_WITHIN(0.5 + 0.0*I, state_get(&pure_state, 0, 0), PRECISION);

    // Set a value and verify roundtrip
    state_set(&pure_state, 1, 0, 3.0 + 1.0*I);
    TEST_ASSERT_COMPLEX_WITHIN(3.0 + 1.0*I, state_get(&pure_state, 1, 0), PRECISION);

    state_free(&pure_state);

    // Test MIXED_PACKED state through dispatch API
    state_t packed_state = {.type = MIXED_PACKED, .data = NULL, .qubits = 0};
    state_plus(&packed_state, 2);

    // Verify plus state: all elements = 1/4 = 0.25
    TEST_ASSERT_COMPLEX_WITHIN(0.25 + 0.0*I, state_get(&packed_state, 1, 0), PRECISION);

    // Set a value and verify roundtrip
    state_set(&packed_state, 1, 0, 3.0 + 1.0*I);
    TEST_ASSERT_COMPLEX_WITHIN(3.0 + 1.0*I, state_get(&packed_state, 1, 0), PRECISION);

    state_free(&packed_state);

    // Test MIXED_TILED state through dispatch API
    state_t tiled_state = {.type = MIXED_TILED, .data = NULL, .qubits = 0};
    state_plus(&tiled_state, 2);

    // Verify plus state: all elements = 1/4 = 0.25
    TEST_ASSERT_COMPLEX_WITHIN(0.25 + 0.0*I, state_get(&tiled_state, 1, 0), PRECISION);

    // Set a value and verify roundtrip
    state_set(&tiled_state, 1, 0, 3.0 + 1.0*I);
    TEST_ASSERT_COMPLEX_WITHIN(3.0 + 1.0*I, state_get(&tiled_state, 1, 0), PRECISION);

    state_free(&tiled_state);
}

/*
 * Test: state_print() smoke test - verify it doesn't crash
 */
void test_dispatch_print(void) {
    // PURE state
    state_t pure_state = {.type = PURE, .data = NULL, .qubits = 0};
    state_plus(&pure_state, 1);
    state_print(&pure_state);
    state_free(&pure_state);

    // MIXED_PACKED state
    state_t packed_state = {.type = MIXED_PACKED, .data = NULL, .qubits = 0};
    state_plus(&packed_state, 1);
    state_print(&packed_state);
    state_free(&packed_state);

    // MIXED_TILED state
    state_t tiled_state = {.type = MIXED_TILED, .data = NULL, .qubits = 0};
    state_plus(&tiled_state, 1);
    state_print(&tiled_state);
    state_free(&tiled_state);
}

/*
 * =====================================================================================================================
 * Unity: setUp and tearDown
 * =====================================================================================================================
 */

void setUp(void) {}

void tearDown(void) {}

/*
 * =====================================================================================================================
 * Main test runner
 * =====================================================================================================================
 */

int main(void) {
    UNITY_BEGIN();

    /*
     * PURE state tests (10 tests)
     */
    printf("\n========== PURE State Tests ==========\n");
    RUN_TEST(test_pure_len);
    RUN_TEST(test_pure_init_null_data);
    RUN_TEST(test_pure_init_ownership);
    RUN_TEST(test_pure_free_double_free);
    RUN_TEST(test_pure_plus_values);
    RUN_TEST(test_pure_get_set_roundtrip);
    RUN_TEST(test_pure_cp_independence);
    RUN_TEST(test_pure_single_qubit);
    RUN_TEST(test_pure_four_qubits);
    RUN_TEST(test_pure_reinit);

    /*
     * MIXED_PACKED state tests (13 tests)
     */
    printf("\n========== MIXED_PACKED State Tests ==========\n");
    RUN_TEST(test_packed_len);
    RUN_TEST(test_packed_init_null_data);
    RUN_TEST(test_packed_init_ownership);
    RUN_TEST(test_packed_free_double_free);
    RUN_TEST(test_packed_plus_values);
    RUN_TEST(test_packed_get_lower);
    RUN_TEST(test_packed_get_upper);
    RUN_TEST(test_packed_set_lower);
    RUN_TEST(test_packed_set_upper);
    RUN_TEST(test_packed_cp_independence);
    RUN_TEST(test_packed_hermitian_symmetry);
    RUN_TEST(test_packed_single_qubit);
    RUN_TEST(test_packed_zero_qubit);

    /*
     * MIXED_TILED state tests (16 tests)
     */
    printf("\n========== MIXED_TILED State Tests ==========\n");
    RUN_TEST(test_tiled_len_small);
    RUN_TEST(test_tiled_len_8_qubits);
    RUN_TEST(test_tiled_init_null_data);
    RUN_TEST(test_tiled_init_ownership);
    RUN_TEST(test_tiled_free_double_free);
    RUN_TEST(test_tiled_plus_values);
    RUN_TEST(test_tiled_get_within_tile);
    RUN_TEST(test_tiled_get_cross_tile);
    RUN_TEST(test_tiled_get_upper);
    RUN_TEST(test_tiled_set_roundtrip);
    RUN_TEST(test_tiled_cp_independence);
    RUN_TEST(test_tiled_hermitian_symmetry);
    RUN_TEST(test_tiled_single_qubit);
    RUN_TEST(test_tiled_tile_boundaries);
    RUN_TEST(test_tiled_set_upper);
    RUN_TEST(test_tiled_6_qubits);

    /*
     * Dispatch layer tests (2 tests)
     */
    printf("\n========== Dispatch Layer Tests ==========\n");
    RUN_TEST(test_dispatch_roundtrip);
    RUN_TEST(test_dispatch_print);

    return UNITY_END();
}
