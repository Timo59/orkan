/**
 * @file test_meas.c
 * @brief Test runner for measurement module tests
 *
 * Tests the mean() function across all three state representations:
 * PURE, MIXED_PACKED, MIXED_TILED.
 */

#include "test.h"
#include "meas.h"

/*
 * =====================================================================================================================
 * Test function declarations from test_meas_pure.c
 * =====================================================================================================================
 */

void test_mean_pure_uniform(void);
void test_mean_pure_identity(void);
void test_mean_pure_basis(void);
void test_mean_pure_basis_k(void);

/*
 * =====================================================================================================================
 * Test function declarations from test_meas_packed.c
 * =====================================================================================================================
 */

void test_mean_packed_uniform(void);
void test_mean_packed_identity(void);
void test_mean_packed_pure_embed(void);

/*
 * =====================================================================================================================
 * Test function declarations from test_meas_tiled.c
 * =====================================================================================================================
 */

void test_mean_tiled_uniform(void);
void test_mean_tiled_identity(void);
void test_mean_tiled_pure_embed(void);
void test_mean_tiled_small(void);
void test_mean_tiled_multi(void);

/*
 * =====================================================================================================================
 * Cross-representation consistency test
 * =====================================================================================================================
 */

/*
 * Same state (|+>^n) and observable (obs[x] = x) must produce equal results
 * across all three representations.
 */
void test_mean_consistency(void) {
    for (qubit_t n = 1; n <= MAXQUBITS; ++n) {
        const idx_t dim = POW2(n, idx_t);
        double obs[dim];
        for (idx_t x = 0; x < dim; ++x)
            obs[x] = (double)x;

        state_t pure = {.type = PURE, .data = NULL, .qubits = 0};
        state_plus(&pure, n);
        double val_pure = mean(&pure, obs);
        state_free(&pure);

        state_t packed = {.type = MIXED_PACKED, .data = NULL, .qubits = 0};
        state_plus(&packed, n);
        double val_packed = mean(&packed, obs);
        state_free(&packed);

        state_t tiled = {.type = MIXED_TILED, .data = NULL, .qubits = 0};
        state_plus(&tiled, n);
        double val_tiled = mean(&tiled, obs);
        state_free(&tiled);

        TEST_ASSERT_DOUBLE_WITHIN(PRECISION, val_pure, val_packed);
        TEST_ASSERT_DOUBLE_WITHIN(PRECISION, val_pure, val_tiled);
    }
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
     * PURE state tests
     */
    printf("\n========== PURE Mean Tests ==========\n");
    RUN_TEST(test_mean_pure_uniform);
    RUN_TEST(test_mean_pure_identity);
    RUN_TEST(test_mean_pure_basis);
    RUN_TEST(test_mean_pure_basis_k);

    /*
     * MIXED_PACKED state tests
     */
    printf("\n========== MIXED_PACKED Mean Tests ==========\n");
    RUN_TEST(test_mean_packed_uniform);
    RUN_TEST(test_mean_packed_identity);
    RUN_TEST(test_mean_packed_pure_embed);

    /*
     * MIXED_TILED state tests
     */
    printf("\n========== MIXED_TILED Mean Tests ==========\n");
    RUN_TEST(test_mean_tiled_uniform);
    RUN_TEST(test_mean_tiled_identity);
    RUN_TEST(test_mean_tiled_pure_embed);
    RUN_TEST(test_mean_tiled_small);
    RUN_TEST(test_mean_tiled_multi);

    /*
     * Cross-representation consistency
     */
    printf("\n========== Consistency Tests ==========\n");
    RUN_TEST(test_mean_consistency);

    return UNITY_END();
}
