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
void test_mean_pure_imaginary(void);
void test_mean_pure_negative_obs(void);
void test_mean_pure_mixed_sign_obs(void);

void test_matel_diag_pure_self(void);
void test_matel_diag_pure_orthogonal(void);
void test_matel_diag_pure_conjugate_symmetry(void);
void test_matel_diag_pure_conjugate_symmetry_complex(void);
void test_matel_diag_pure_known(void);
void test_matel_diag_pure_complex_2q(void);

/*
 * =====================================================================================================================
 * Test function declarations from test_meas_packed.c
 * =====================================================================================================================
 */

void test_mean_packed_uniform(void);
void test_mean_packed_identity(void);
void test_mean_packed_pure_embed(void);
void test_mean_packed_basis_k(void);
void test_mean_packed_negative_obs(void);
void test_mean_packed_mixed_sign_obs(void);

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
void test_mean_tiled_basis_k(void);
void test_mean_tiled_negative_obs(void);
void test_mean_tiled_mixed_sign_obs(void);

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
    RUN_TEST(test_mean_pure_imaginary);
    RUN_TEST(test_mean_pure_negative_obs);
    RUN_TEST(test_mean_pure_mixed_sign_obs);

    /*
     * PURE matel_diag tests
     */
    printf("\n========== PURE Matel Diag Tests ==========\n");
    RUN_TEST(test_matel_diag_pure_self);
    RUN_TEST(test_matel_diag_pure_orthogonal);
    RUN_TEST(test_matel_diag_pure_conjugate_symmetry);
    RUN_TEST(test_matel_diag_pure_conjugate_symmetry_complex);
    RUN_TEST(test_matel_diag_pure_known);
    RUN_TEST(test_matel_diag_pure_complex_2q);

    /*
     * MIXED_PACKED state tests
     */
    printf("\n========== MIXED_PACKED Mean Tests ==========\n");
    RUN_TEST(test_mean_packed_uniform);
    RUN_TEST(test_mean_packed_identity);
    RUN_TEST(test_mean_packed_pure_embed);
    RUN_TEST(test_mean_packed_basis_k);
    RUN_TEST(test_mean_packed_negative_obs);
    RUN_TEST(test_mean_packed_mixed_sign_obs);

    /*
     * MIXED_TILED state tests
     */
    printf("\n========== MIXED_TILED Mean Tests ==========\n");
    RUN_TEST(test_mean_tiled_uniform);
    RUN_TEST(test_mean_tiled_identity);
    RUN_TEST(test_mean_tiled_pure_embed);
    RUN_TEST(test_mean_tiled_small);
    RUN_TEST(test_mean_tiled_multi);
    RUN_TEST(test_mean_tiled_basis_k);
    RUN_TEST(test_mean_tiled_negative_obs);
    RUN_TEST(test_mean_tiled_mixed_sign_obs);

    /*
     * Cross-representation consistency
     */
    printf("\n========== Consistency Tests ==========\n");
    RUN_TEST(test_mean_consistency);

    return UNITY_END();
}
