/**
 * @file test_circ.c
 * @brief Test runner for circuit module tests
 *
 * Tests the exp_diag() function across all three state representations:
 * PURE, MIXED_PACKED, MIXED_TILED.
 */

#include "test.h"
#include "test_circ.h"
#include "circ.h"

/*
 * =====================================================================================================================
 * Test function declarations from test_circ_pure.c
 * =====================================================================================================================
 */

void test_exp_diag_pure_basis_k(void);
void test_exp_diag_pure_plus(void);
void test_exp_diag_pure_zero_obs(void);
void test_exp_diag_pure_global_phase(void);
void test_exp_diag_pure_unitarity(void);

/*
 * =====================================================================================================================
 * Test function declarations from test_circ_packed.c
 * =====================================================================================================================
 */

void test_exp_diag_packed_basis_k(void);
void test_exp_diag_packed_plus(void);
void test_exp_diag_packed_double_apply(void);
void test_exp_diag_packed_zero_obs(void);
void test_exp_diag_packed_identity_obs(void);

/*
 * =====================================================================================================================
 * Test function declarations from test_circ_tiled.c
 * =====================================================================================================================
 */

void test_exp_diag_tiled_basis_k(void);
void test_exp_diag_tiled_plus(void);
void test_exp_diag_tiled_double_apply(void);
void test_exp_diag_tiled_zero_obs(void);
void test_exp_diag_tiled_identity_obs(void);
void test_exp_diag_tiled_small(void);
void test_exp_diag_tiled_multi(void);

/*
 * =====================================================================================================================
 * Cross-representation consistency test
 * =====================================================================================================================
 */

/*
 * Same state (|+>^n) and observable (diag[x] = x) must produce equivalent
 * results across all three representations.  Since |+>^n is pure, we verify
 *   rho_packed = rho_tiled = |psi><psi|
 * element-wise after evolution.
 */
void test_exp_diag_consistency(void) {
    const double t = 0.77;
    for (qubit_t n = 1; n <= MAXQUBITS; ++n) {
        const idx_t dim = POW2(n, idx_t);
        double diag[dim];
        fill_ramp_diag(diag, dim);

        state_t pure   = {.type = PURE,         .data = NULL, .qubits = 0};
        state_t packed = {.type = MIXED_PACKED,  .data = NULL, .qubits = 0};
        state_t tiled  = {.type = MIXED_TILED,   .data = NULL, .qubits = 0};

        state_plus(&pure, n);
        state_plus(&packed, n);
        state_plus(&tiled, n);

        exp_diag(&pure, diag, t);
        exp_diag(&packed, diag, t);
        exp_diag(&tiled, diag, t);

        /* Compare packed and tiled against each other element-wise */
        for (idx_t r = 0; r < dim; ++r)
            for (idx_t c = 0; c <= r; ++c)
                TEST_ASSERT_COMPLEX_WITHIN(state_get(&packed, r, c),
                                           state_get(&tiled, r, c), PRECISION);

        /* Compare pure |psi><psi| outer product against packed */
        for (idx_t r = 0; r < dim; ++r) {
            for (idx_t c = 0; c <= r; ++c) {
                cplx_t expected = state_get(&pure, r, 0) * conj(state_get(&pure, c, 0));
                TEST_ASSERT_COMPLEX_WITHIN(expected,
                                           state_get(&packed, r, c), PRECISION);
            }
        }

        state_free(&pure);
        state_free(&packed);
        state_free(&tiled);
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
    printf("\n========== PURE exp_diag Tests ==========\n");
    RUN_TEST(test_exp_diag_pure_basis_k);
    RUN_TEST(test_exp_diag_pure_plus);
    RUN_TEST(test_exp_diag_pure_zero_obs);
    RUN_TEST(test_exp_diag_pure_global_phase);
    RUN_TEST(test_exp_diag_pure_unitarity);

    /*
     * MIXED_PACKED state tests
     */
    printf("\n========== MIXED_PACKED exp_diag Tests ==========\n");
    RUN_TEST(test_exp_diag_packed_basis_k);
    RUN_TEST(test_exp_diag_packed_plus);
    RUN_TEST(test_exp_diag_packed_double_apply);
    RUN_TEST(test_exp_diag_packed_zero_obs);
    RUN_TEST(test_exp_diag_packed_identity_obs);

    /*
     * MIXED_TILED state tests
     */
    printf("\n========== MIXED_TILED exp_diag Tests ==========\n");
    RUN_TEST(test_exp_diag_tiled_basis_k);
    RUN_TEST(test_exp_diag_tiled_plus);
    RUN_TEST(test_exp_diag_tiled_double_apply);
    RUN_TEST(test_exp_diag_tiled_zero_obs);
    RUN_TEST(test_exp_diag_tiled_identity_obs);
    RUN_TEST(test_exp_diag_tiled_small);
    RUN_TEST(test_exp_diag_tiled_multi);

    /*
     * Cross-representation consistency
     */
    printf("\n========== Consistency Tests ==========\n");
    RUN_TEST(test_exp_diag_consistency);

    return UNITY_END();
}
