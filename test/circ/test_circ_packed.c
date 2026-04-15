/**
 * @file test_circ_packed.c
 * @brief Tests for exp_diag() with MIXED_PACKED state representation
 */

#include "test.h"
#include "test_circ.h"
#include "circ.h"

#include <math.h>

/*
 * =====================================================================================================================
 * Tests
 * =====================================================================================================================
 */

/*
 * |k><k| with diag[x] = x, t = 1.0:
 *   Diagonal density matrix — all off-diagonal elements are zero,
 *   so the state must be unchanged.
 */
void test_exp_diag_packed_basis_k(void) {
    for (qubit_t n = 1; n <= MAXQUBITS; ++n) {
        const idx_t dim = POW2(n, idx_t);
        double diag[dim];
        fill_ramp_diag(diag, dim);

        for (idx_t k = 0; k < dim; ++k) {
            state_t s = {.type = MIXED_PACKED, .data = NULL, .qubits = 0};
            state_init(&s, n, NULL);
            state_set(&s, k, k, 1.0 + 0.0 * I);

            exp_diag(&s, diag, 1.0);

            for (idx_t r = 0; r < dim; ++r) {
                for (idx_t c = 0; c <= r; ++c) {
                    cplx_t expected = (r == k && c == k) ? 1.0 : 0.0;
                    TEST_ASSERT_COMPLEX_WITHIN(expected,
                                               state_get(&s, r, c), PRECISION);
                }
            }
            state_free(&s);
        }
    }
}

/*
 * |+><+| with diag[x] = x, t = pi/4:
 *   rho[r][c] = (1/dim) * e^{-i (r-c) pi/4}
 */
void test_exp_diag_packed_plus(void) {
    const double t = M_PI / 4.0;
    for (qubit_t n = 1; n <= MAXQUBITS; ++n) {
        const idx_t dim = POW2(n, idx_t);
        double diag[dim];
        fill_ramp_diag(diag, dim);

        state_t s = {.type = MIXED_PACKED, .data = NULL, .qubits = 0};
        state_plus(&s, n);
        exp_diag(&s, diag, t);

        const double val = 1.0 / (double)dim;
        for (idx_t r = 0; r < dim; ++r) {
            for (idx_t c = 0; c <= r; ++c) {
                double angle = ((double)r - (double)c) * t;
                cplx_t expected = CMPLX(val * cos(angle), -val * sin(angle));
                TEST_ASSERT_COMPLEX_WITHIN(expected,
                                           state_get(&s, r, c), PRECISION);
            }
        }
        state_free(&s);
    }
}

/*
 * Apply exp_diag twice:  exp(-iDt) exp(-iDt) = exp(-i D 2t).
 * First application creates complex off-diagonals from the real |+><+| state;
 * second exercises the full complex * complex multiplication path.
 */
void test_exp_diag_packed_double_apply(void) {
    const double t = M_PI / 4.0;
    for (qubit_t n = 1; n <= MAXQUBITS; ++n) {
        const idx_t dim = POW2(n, idx_t);
        double diag[dim];
        fill_ramp_diag(diag, dim);

        state_t s = {.type = MIXED_PACKED, .data = NULL, .qubits = 0};
        state_plus(&s, n);
        exp_diag(&s, diag, t);
        exp_diag(&s, diag, t);

        const double val = 1.0 / (double)dim;
        for (idx_t r = 0; r < dim; ++r) {
            for (idx_t c = 0; c <= r; ++c) {
                double angle = ((double)r - (double)c) * 2.0 * t;
                cplx_t expected = CMPLX(val * cos(angle), -val * sin(angle));
                TEST_ASSERT_COMPLEX_WITHIN(expected,
                                           state_get(&s, r, c), PRECISION);
            }
        }
        state_free(&s);
    }
}

/*
 * diag[x] = 0:  state unchanged.
 */
void test_exp_diag_packed_zero_obs(void) {
    for (qubit_t n = 1; n <= MAXQUBITS; ++n) {
        const idx_t dim = POW2(n, idx_t);
        double diag[dim];
        fill_zero_diag(diag, dim);

        state_t s = {.type = MIXED_PACKED, .data = NULL, .qubits = 0};
        state_plus(&s, n);
        state_t ref = state_cp(&s);

        exp_diag(&s, diag, 1.0);

        for (idx_t r = 0; r < dim; ++r)
            for (idx_t c = 0; c <= r; ++c)
                TEST_ASSERT_COMPLEX_WITHIN(state_get(&ref, r, c),
                                           state_get(&s, r, c), PRECISION);
        state_free(&s);
        state_free(&ref);
    }
}

/*
 * diag[x] = 1:  all phase differences are zero, state unchanged.
 */
void test_exp_diag_packed_identity_obs(void) {
    for (qubit_t n = 1; n <= MAXQUBITS; ++n) {
        const idx_t dim = POW2(n, idx_t);
        double diag[dim];
        fill_uniform_diag(diag, dim);

        state_t s = {.type = MIXED_PACKED, .data = NULL, .qubits = 0};
        state_plus(&s, n);
        state_t ref = state_cp(&s);

        exp_diag(&s, diag, 2.71);

        for (idx_t r = 0; r < dim; ++r)
            for (idx_t c = 0; c <= r; ++c)
                TEST_ASSERT_COMPLEX_WITHIN(state_get(&ref, r, c),
                                           state_get(&s, r, c), PRECISION);
        state_free(&s);
        state_free(&ref);
    }
}
