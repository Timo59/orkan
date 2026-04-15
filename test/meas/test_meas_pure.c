/**
 * @file test_meas_pure.c
 * @brief Tests for mean() with PURE state representation
 */

#include "test.h"
#include "test_meas.h"
#include "meas.h"


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

/*
 * Purely imaginary amplitudes: |psi> = (i/sqrt(2), i/sqrt(2))
 * Verifies that cimag(a)^2 contributes to |psi|^2.
 */
void test_mean_pure_imaginary(void) {
    state_t s = {.type = PURE, .data = NULL, .qubits = 0};
    state_init(&s, 1, NULL);
    state_set(&s, 0, 0, 0.0 + INVSQRT2 * I);
    state_set(&s, 1, 0, 0.0 + INVSQRT2 * I);

    double obs[2] = {0.0, 1.0};
    TEST_ASSERT_DOUBLE_WITHIN(PRECISION, 0.5, mean(&s, obs));
    state_free(&s);
}

/*
 * |+>^n with obs[x] = -1:  expected = -1
 */
void test_mean_pure_negative_obs(void) {
    for (qubit_t n = 1; n <= MAXQUBITS; ++n) {
        state_t s = {.type = PURE, .data = NULL, .qubits = 0};
        state_plus(&s, n);

        const idx_t dim = POW2(n, idx_t);
        double obs[dim];
        for (idx_t x = 0; x < dim; ++x)
            obs[x] = -1.0;

        TEST_ASSERT_DOUBLE_WITHIN(PRECISION, -1.0, mean(&s, obs));
        state_free(&s);
    }
}

/*
 * |+>^2 with obs = {-3, +1, -1, +3}:  expected = 0 (exact cancellation)
 */
void test_mean_pure_mixed_sign_obs(void) {
    state_t s = {.type = PURE, .data = NULL, .qubits = 0};
    state_plus(&s, 2);

    double obs[4] = {-3.0, 1.0, -1.0, 3.0};
    TEST_ASSERT_DOUBLE_WITHIN(PRECISION, 0.0, mean(&s, obs));
    state_free(&s);
}
