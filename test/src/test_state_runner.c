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
     * PURE state tests (9 tests)
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

    /*
     * MIXED_PACKED state tests (12 tests)
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

    /*
     * MIXED_TILED state tests (14 tests)
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

    return UNITY_END();
}
