/**
 * @file test_meas_pure.c
 * @brief Tests for mean() with PURE state representation
 */

#include "test.h"
#include "test_meas.h"
#include "meas.h"
#include "gate.h"


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


/*
 * =====================================================================================================================
 * matel_diag tests
 * =====================================================================================================================
 */

/*
 * matel_diag(s, s, obs) == mean(s, obs)  (diagonal = expectation value)
 */
void test_matel_diag_pure_self(void) {
    for (qubit_t n = 1; n <= MAXQUBITS; ++n) {
        state_t s = {.type = PURE, .data = NULL, .qubits = 0};
        state_plus(&s, n);

        const idx_t dim = POW2(n, idx_t);
        double obs[dim];
        fill_identity_obs(obs, dim);

        cplx_t result = matel_diag(&s, &s, obs);
        double expected = mean(&s, obs);

        TEST_ASSERT_DOUBLE_WITHIN(PRECISION, expected, creal(result));
        TEST_ASSERT_DOUBLE_WITHIN(PRECISION, 0.0,      cimag(result));
        state_free(&s);
    }
}

/*
 * <0|H|1> = 0  for any diagonal H  (orthogonal basis states)
 */
void test_matel_diag_pure_orthogonal(void) {
    for (qubit_t n = 1; n <= MAXQUBITS; ++n) {
        state_t bra = {.type = PURE, .data = NULL, .qubits = 0};
        state_t ket = {.type = PURE, .data = NULL, .qubits = 0};
        state_init(&bra, n, NULL);
        state_init(&ket, n, NULL);
        state_set(&bra, 0, 0, 1.0 + 0.0 * I);
        state_set(&ket, 1, 0, 1.0 + 0.0 * I);

        const idx_t dim = POW2(n, idx_t);
        double obs[dim];
        fill_identity_obs(obs, dim);

        cplx_t result = matel_diag(&bra, &ket, obs);
        TEST_ASSERT_COMPLEX_WITHIN(0.0 + 0.0 * I, result, PRECISION);

        state_free(&bra);
        state_free(&ket);
    }
}

/*
 * matel_diag(a, b, obs) == conj(matel_diag(b, a, obs))
 */
void test_matel_diag_pure_conjugate_symmetry(void) {
    for (qubit_t n = 1; n <= MAXQUBITS; ++n) {
        state_t a = {.type = PURE, .data = NULL, .qubits = 0};
        state_t b = {.type = PURE, .data = NULL, .qubits = 0};
        state_plus(&a, n);
        state_init(&b, n, NULL);
        state_set(&b, 0, 0, 1.0 + 0.0 * I);

        const idx_t dim = POW2(n, idx_t);
        double obs[dim];
        fill_identity_obs(obs, dim);

        cplx_t ab = matel_diag(&a, &b, obs);
        cplx_t ba = matel_diag(&b, &a, obs);
        TEST_ASSERT_COMPLEX_WITHIN(conj(ba), ab, PRECISION);

        state_free(&a);
        state_free(&b);
    }
}

/*
 * Conjugate-symmetry test on genuinely complex states.
 *
 * The earlier conjugate_symmetry test uses |+>^n and |0>, both real-valued, so
 * matel_diag(a,b) and matel_diag(b,a) are real and the relation conj(ba)==ab
 * collapses to ba==ab — a sign flip in the imaginary-part accumulation
 * (sum_im in matel_diag_pure) would still pass.  Here a is purely imaginary
 * on |1>, making the result purely imaginary and non-degenerate:
 *
 *   a = (1/sqrt2)|0> + (i/sqrt2)|1>,  b = |+> = (1/sqrt2)|0> + (1/sqrt2)|1>
 *   matel(a,b,obs={0,1}) = conj(i/sqrt2)*(1/sqrt2) = -i/2
 *   matel(b,a,obs={0,1}) = conj(1/sqrt2)*(i/sqrt2) =  i/2
 */
void test_matel_diag_pure_conjugate_symmetry_complex(void) {
    state_t a = {.type = PURE, .data = NULL, .qubits = 0};
    state_t b = {.type = PURE, .data = NULL, .qubits = 0};
    state_init(&a, 1, NULL);
    state_set(&a, 0, 0, INVSQRT2 + 0.0 * I);
    state_set(&a, 1, 0, 0.0 + INVSQRT2 * I);
    state_plus(&b, 1);

    double obs[2] = {0.0, 1.0};
    cplx_t ab = matel_diag(&a, &b, obs);
    cplx_t ba = matel_diag(&b, &a, obs);

    TEST_ASSERT_COMPLEX_WITHIN(0.0 - 0.5 * I, ab, PRECISION);
    TEST_ASSERT_COMPLEX_WITHIN(0.0 + 0.5 * I, ba, PRECISION);
    TEST_ASSERT_COMPLEX_WITHIN(conj(ba), ab, PRECISION);
    /* Non-degeneracy: ab != ba.  Without this, a sign flip in sum_im that
       made matel_diag plain symmetric (instead of conjugate-symmetric)
       would still satisfy the conj relation above. */
    TEST_ASSERT_TRUE(cabs(ab - ba) > 0.5);

    state_free(&a);
    state_free(&b);
}

/*
 * Known value: bra = (1/sqrt2, i/sqrt2), ket = (1/sqrt2, 1/sqrt2), obs = {0, 1}
 * Expected: 0 * conj(1/sqrt2) * (1/sqrt2) + 1 * conj(i/sqrt2) * (1/sqrt2)
 *         = (-i/sqrt2) * (1/sqrt2) = -i/2
 */
void test_matel_diag_pure_known(void) {
    state_t bra = {.type = PURE, .data = NULL, .qubits = 0};
    state_t ket = {.type = PURE, .data = NULL, .qubits = 0};
    state_init(&bra, 1, NULL);
    state_init(&ket, 1, NULL);
    state_set(&bra, 0, 0, INVSQRT2 + 0.0 * I);
    state_set(&bra, 1, 0, 0.0 + INVSQRT2 * I);
    state_set(&ket, 0, 0, INVSQRT2 + 0.0 * I);
    state_set(&ket, 1, 0, INVSQRT2 + 0.0 * I);

    double obs[2] = {0.0, 1.0};
    cplx_t result = matel_diag(&bra, &ket, obs);
    TEST_ASSERT_COMPLEX_WITHIN(0.0 - 0.5 * I, result, PRECISION);

    state_free(&bra);
    state_free(&ket);
}



/*
 * Known-value matel_diag with complex amplitudes at n=2.
 * Breaks the circular matel_diag == mean dependency.
 *
 * bra = |+>^2 = (1/2, 1/2, 1/2, 1/2)
 * ket = i*|+>^2 = (i/2, i/2, i/2, i/2)
 * obs = {0, 1, 2, 3}
 *
 * matel = sum_x obs[x] * conj(1/2) * (i/2) = (i/4) * (0+1+2+3) = 6i/4 = 1.5i
 */
void test_matel_diag_pure_complex_2q(void) {
    state_t bra = {.type = PURE, .data = NULL, .qubits = 0};
    state_t ket = {.type = PURE, .data = NULL, .qubits = 0};
    state_plus(&bra, 2);

    /* ket = i * |+>^2 */
    state_init(&ket, 2, NULL);
    const idx_t dim = POW2(2, idx_t);
    for (idx_t x = 0; x < dim; ++x)
        state_set(&ket, x, 0, 0.0 + 0.5 * I);

    double obs[4] = {0.0, 1.0, 2.0, 3.0};
    cplx_t result = matel_diag(&bra, &ket, obs);
    TEST_ASSERT_COMPLEX_WITHIN(0.0 + 1.5 * I, result, PRECISION);

    state_free(&bra);
    state_free(&ket);
}
