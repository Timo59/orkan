/**
 * @file test_circ_pure.c
 * @brief Tests for exp_diag() with PURE state representation
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
 * |k> with diag[x] = x, t = 1.0:
 *   amplitude at k should be e^{-ik} = cos(k) - i sin(k), all others zero.
 */
void test_exp_diag_pure_basis_k(void) {
    for (qubit_t n = 1; n <= MAXQUBITS; ++n) {
        const idx_t dim = POW2(n, idx_t);
        double diag[dim];
        fill_ramp_diag(diag, dim);

        for (idx_t k = 0; k < dim; ++k) {
            state_t s = {.type = PURE, .data = NULL, .qubits = 0};
            state_init(&s, n, NULL);
            state_set(&s, k, 0, 1.0 + 0.0 * I);

            exp_diag(&s, diag, 1.0);

            for (idx_t x = 0; x < dim; ++x) {
                cplx_t expected = (x == k)
                    ? CMPLX(cos((double)k), -sin((double)k))
                    : 0.0;
                TEST_ASSERT_COMPLEX_WITHIN(expected, state_get(&s, x, 0), PRECISION);
            }
            state_free(&s);
        }
    }
}

/*
 * |+>^n with diag[x] = x, t = pi/4:
 *   each amplitude = (1/sqrt(dim)) * e^{-i x pi/4}
 */
void test_exp_diag_pure_plus(void) {
    const double t = M_PI / 4.0;
    for (qubit_t n = 1; n <= MAXQUBITS; ++n) {
        const idx_t dim = POW2(n, idx_t);
        double diag[dim];
        fill_ramp_diag(diag, dim);

        state_t s = {.type = PURE, .data = NULL, .qubits = 0};
        state_plus(&s, n);
        exp_diag(&s, diag, t);

        const double amp = 1.0 / sqrt((double)dim);
        for (idx_t x = 0; x < dim; ++x) {
            double angle = (double)x * t;
            cplx_t expected = CMPLX(amp * cos(angle), -amp * sin(angle));
            TEST_ASSERT_COMPLEX_WITHIN(expected, state_get(&s, x, 0), PRECISION);
        }
        state_free(&s);
    }
}

/*
 * diag[x] = 0 for all x:  state must be unchanged.
 */
void test_exp_diag_pure_zero_obs(void) {
    for (qubit_t n = 1; n <= MAXQUBITS; ++n) {
        const idx_t dim = POW2(n, idx_t);
        double diag[dim];
        fill_zero_diag(diag, dim);

        state_t s = {.type = PURE, .data = NULL, .qubits = 0};
        state_plus(&s, n);
        state_t ref = state_cp(&s);

        exp_diag(&s, diag, 1.0);

        for (idx_t x = 0; x < dim; ++x)
            TEST_ASSERT_COMPLEX_WITHIN(state_get(&ref, x, 0),
                                       state_get(&s, x, 0), PRECISION);
        state_free(&s);
        state_free(&ref);
    }
}

/*
 * diag[x] = 1 for all x, t = pi/3:
 *   all amplitudes acquire global phase e^{-i pi/3}.
 */
void test_exp_diag_pure_global_phase(void) {
    const double t = M_PI / 3.0;
    for (qubit_t n = 1; n <= MAXQUBITS; ++n) {
        const idx_t dim = POW2(n, idx_t);
        double diag[dim];
        fill_uniform_diag(diag, dim);

        state_t s = {.type = PURE, .data = NULL, .qubits = 0};
        state_plus(&s, n);
        state_t ref = state_cp(&s);

        exp_diag(&s, diag, t);

        const cplx_t phase = CMPLX(cos(t), -sin(t));
        for (idx_t x = 0; x < dim; ++x) {
            cplx_t expected = state_get(&ref, x, 0) * phase;
            TEST_ASSERT_COMPLEX_WITHIN(expected, state_get(&s, x, 0), PRECISION);
        }
        state_free(&s);
        state_free(&ref);
    }
}

/*
 * Apply exp_diag twice: exp(-iDt) exp(-iDt) = exp(-iD 2t).
 * First application creates complex amplitudes from the real |+> state;
 * second exercises the full complex * complex multiplication path.
 */
void test_exp_diag_pure_double_apply(void) {
    const double t = M_PI / 4.0;
    for (qubit_t n = 1; n <= MAXQUBITS; ++n) {
        const idx_t dim = POW2(n, idx_t);
        double diag[dim];
        fill_ramp_diag(diag, dim);

        state_t s = {.type = PURE, .data = NULL, .qubits = 0};
        state_plus(&s, n);
        exp_diag(&s, diag, t);
        exp_diag(&s, diag, t);

        const double amp = 1.0 / sqrt((double)dim);
        for (idx_t x = 0; x < dim; ++x) {
            double angle = (double)x * 2.0 * t;
            cplx_t expected = CMPLX(amp * cos(angle), -amp * sin(angle));
            TEST_ASSERT_COMPLEX_WITHIN(expected, state_get(&s, x, 0), PRECISION);
        }
        state_free(&s);
    }
}

/*
 * Negative t: exp(-iDt) exp(+iDt) = I, state must return to original.
 */
void test_exp_diag_pure_negative_t(void) {
    const double t = 0.77;
    for (qubit_t n = 1; n <= MAXQUBITS; ++n) {
        const idx_t dim = POW2(n, idx_t);
        double diag[dim];
        fill_ramp_diag(diag, dim);

        state_t s = {.type = PURE, .data = NULL, .qubits = 0};
        state_plus(&s, n);
        state_t ref = state_cp(&s);

        exp_diag(&s, diag, t);
        exp_diag(&s, diag, -t);

        for (idx_t x = 0; x < dim; ++x)
            TEST_ASSERT_COMPLEX_WITHIN(state_get(&ref, x, 0),
                                       state_get(&s, x, 0), PRECISION);
        state_free(&s);
        state_free(&ref);
    }
}

/*
 * Mixed-sign diagonal: diag[x] = +x (even) or -x (odd).
 * Stresses all sign-quadrant combinations in the trig product formula.
 */
void test_exp_diag_pure_mixed_sign_diag(void) {
    const double t = M_PI / 4.0;
    for (qubit_t n = 1; n <= MAXQUBITS; ++n) {
        const idx_t dim = POW2(n, idx_t);
        double diag[dim];
        fill_alternating_diag(diag, dim);

        state_t s = {.type = PURE, .data = NULL, .qubits = 0};
        state_plus(&s, n);
        exp_diag(&s, diag, t);

        const double amp = 1.0 / sqrt((double)dim);
        for (idx_t x = 0; x < dim; ++x) {
            double angle = diag[x] * t;
            cplx_t expected = CMPLX(amp * cos(angle), -amp * sin(angle));
            TEST_ASSERT_COMPLEX_WITHIN(expected, state_get(&s, x, 0), PRECISION);
        }
        state_free(&s);
    }
}

/*
 * t = 0 with non-trivial diagonal: state must be unchanged.
 * Distinct from zero_obs (which uses diag = 0).
 */
void test_exp_diag_pure_zero_t(void) {
    for (qubit_t n = 1; n <= MAXQUBITS; ++n) {
        const idx_t dim = POW2(n, idx_t);
        double diag[dim];
        fill_ramp_diag(diag, dim);

        state_t s = {.type = PURE, .data = NULL, .qubits = 0};
        state_plus(&s, n);
        state_t ref = state_cp(&s);

        exp_diag(&s, diag, 0.0);

        for (idx_t x = 0; x < dim; ++x)
            TEST_ASSERT_COMPLEX_WITHIN(state_get(&ref, x, 0),
                                       state_get(&s, x, 0), PRECISION);
        state_free(&s);
        state_free(&ref);
    }
}

/*
 * Unitarity check: sum |psi[x]|^2 must remain 1 after evolution.
 */
void test_exp_diag_pure_unitarity(void) {
    for (qubit_t n = 1; n <= MAXQUBITS; ++n) {
        const idx_t dim = POW2(n, idx_t);
        double diag[dim];
        fill_ramp_diag(diag, dim);

        state_t s = {.type = PURE, .data = NULL, .qubits = 0};
        state_plus(&s, n);
        exp_diag(&s, diag, 1.23);

        double norm = 0.0;
        for (idx_t x = 0; x < dim; ++x) {
            cplx_t a = state_get(&s, x, 0);
            norm += creal(a) * creal(a) + cimag(a) * cimag(a);
        }
        TEST_ASSERT_DOUBLE_WITHIN(PRECISION, 1.0, norm);
        state_free(&s);
    }
}
