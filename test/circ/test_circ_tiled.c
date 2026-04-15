/**
 * @file test_circ_tiled.c
 * @brief Tests for exp_diag() with MIXED_TILED state representation
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
 * |k><k| with diag[x] = x:  diagonal state unchanged.
 */
void test_exp_diag_tiled_basis_k(void) {
    for (qubit_t n = 1; n <= MAXQUBITS; ++n) {
        const idx_t dim = POW2(n, idx_t);
        double diag[dim];
        fill_ramp_diag(diag, dim);

        for (idx_t k = 0; k < dim; ++k) {
            state_t s = {.type = MIXED_TILED, .data = NULL, .qubits = 0};
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
 * |+><+| with diag[x] = x, t = pi/4.
 */
void test_exp_diag_tiled_plus(void) {
    const double t = M_PI / 4.0;
    for (qubit_t n = 1; n <= MAXQUBITS; ++n) {
        const idx_t dim = POW2(n, idx_t);
        double diag[dim];
        fill_ramp_diag(diag, dim);

        state_t s = {.type = MIXED_TILED, .data = NULL, .qubits = 0};
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
 * Double application: exercises complex * complex path (see packed equivalent).
 */
void test_exp_diag_tiled_double_apply(void) {
    const double t = M_PI / 4.0;
    for (qubit_t n = 1; n <= MAXQUBITS; ++n) {
        const idx_t dim = POW2(n, idx_t);
        double diag[dim];
        fill_ramp_diag(diag, dim);

        state_t s = {.type = MIXED_TILED, .data = NULL, .qubits = 0};
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
void test_exp_diag_tiled_zero_obs(void) {
    for (qubit_t n = 1; n <= MAXQUBITS; ++n) {
        const idx_t dim = POW2(n, idx_t);
        double diag[dim];
        fill_zero_diag(diag, dim);

        state_t s = {.type = MIXED_TILED, .data = NULL, .qubits = 0};
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
 * diag[x] = 1:  uniform diagonal, all phase diffs zero, state unchanged.
 */
void test_exp_diag_tiled_identity_obs(void) {
    for (qubit_t n = 1; n <= MAXQUBITS; ++n) {
        const idx_t dim = POW2(n, idx_t);
        double diag[dim];
        fill_uniform_diag(diag, dim);

        state_t s = {.type = MIXED_TILED, .data = NULL, .qubits = 0};
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

/*
 * n = 1 (single tile with padding).
 */
void test_exp_diag_tiled_small(void) {
    state_t s = {.type = MIXED_TILED, .data = NULL, .qubits = 0};
    state_plus(&s, 1);

    double diag[2] = {0.0, 1.0};
    const double t = M_PI / 2.0;
    exp_diag(&s, diag, t);

    /* rho[0][0] = 0.5 (unchanged), rho[1][1] = 0.5 (unchanged) */
    TEST_ASSERT_COMPLEX_WITHIN(CMPLX(0.5, 0.0), state_get(&s, 0, 0), PRECISION);
    TEST_ASSERT_COMPLEX_WITHIN(CMPLX(0.5, 0.0), state_get(&s, 1, 1), PRECISION);

    /* rho[1][0] = 0.5 * e^{-i(1-0)*pi/2} = 0.5 * (-i) = -0.5i */
    TEST_ASSERT_COMPLEX_WITHIN(CMPLX(0.0, -0.5), state_get(&s, 1, 0), PRECISION);

    state_free(&s);
}

/*
 * n = LOG_TILE_DIM + 1:  exercises inter-tile traversal.
 */
void test_exp_diag_tiled_multi(void) {
    const qubit_t n = LOG_TILE_DIM + 1;
    const idx_t dim = POW2(n, idx_t);

    state_t s = {.type = MIXED_TILED, .data = NULL, .qubits = 0};
    state_plus(&s, n);

    double diag[2 * TILE_DIM];
    fill_ramp_diag(diag, dim);

    const double t = M_PI / 4.0;
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
