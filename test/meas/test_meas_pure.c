/**
 * @file test_meas_pure.c
 * @brief Tests for mean() with PURE state representation
 */

#include "test.h"
#include "meas.h"


/*
 * =====================================================================================================================
 * Helper: fill observable obs[x] = x for x = 0 .. dim-1
 * =====================================================================================================================
 */

static void fill_identity_obs(double *obs, idx_t dim) {
    for (idx_t x = 0; x < dim; ++x)
        obs[x] = (double)x;
}

static void fill_uniform_obs(double *obs, idx_t dim) {
    for (idx_t x = 0; x < dim; ++x)
        obs[x] = 1.0;
}

/*
 * =====================================================================================================================
 * Tests
 * =====================================================================================================================
 */

/*
 * |+>^n with obs[x] = 1:  sum |psi|^2 = 1
 */
void test_mean_pure_uniform(void) {
    for (qubit_t n = 1; n <= MAXQUBITS; ++n) {
        state_t s = {.type = PURE, .data = NULL, .qubits = 0};
        state_plus(&s, n);

        const idx_t dim = POW2(n, idx_t);
        double obs[dim];
        fill_uniform_obs(obs, dim);

        TEST_ASSERT_DOUBLE_WITHIN(PRECISION, 1.0, mean(&s, obs));
        state_free(&s);
    }
}

/*
 * |+>^n with obs[x] = x:  expected = (2^n - 1) / 2
 */
void test_mean_pure_identity(void) {
    for (qubit_t n = 1; n <= MAXQUBITS; ++n) {
        state_t s = {.type = PURE, .data = NULL, .qubits = 0};
        state_plus(&s, n);

        const idx_t dim = POW2(n, idx_t);
        double obs[dim];
        fill_identity_obs(obs, dim);

        double expected = (double)(dim - 1) / 2.0;
        TEST_ASSERT_DOUBLE_WITHIN(PRECISION, expected, mean(&s, obs));
        state_free(&s);
    }
}

/*
 * |0> with obs[x] = x:  expected = 0
 */
void test_mean_pure_basis(void) {
    for (qubit_t n = 1; n <= MAXQUBITS; ++n) {
        state_t s = {.type = PURE, .data = NULL, .qubits = 0};
        state_init(&s, n, NULL);
        state_set(&s, 0, 0, 1.0 + 0.0 * I);

        const idx_t dim = POW2(n, idx_t);
        double obs[dim];
        fill_identity_obs(obs, dim);

        TEST_ASSERT_DOUBLE_WITHIN(PRECISION, 0.0, mean(&s, obs));
        state_free(&s);
    }
}

/*
 * |k> with obs[x] = x:  expected = k
 */
void test_mean_pure_basis_k(void) {
    for (qubit_t n = 2; n <= MAXQUBITS; ++n) {
        const idx_t dim = POW2(n, idx_t);
        double obs[dim];
        fill_identity_obs(obs, dim);

        for (idx_t k = 0; k < dim; ++k) {
            state_t s = {.type = PURE, .data = NULL, .qubits = 0};
            state_init(&s, n, NULL);
            state_set(&s, k, 0, 1.0 + 0.0 * I);

            TEST_ASSERT_DOUBLE_WITHIN(PRECISION, (double)k, mean(&s, obs));
            state_free(&s);
        }
    }
}
