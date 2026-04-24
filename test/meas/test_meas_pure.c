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
 * qop_t wrappers for sample_matrix tests
 * =====================================================================================================================
 */

static void op_identity(state_t *s) { (void)s; }
static void op_x0(state_t *s) { x(s, 0); }
static void op_h0(state_t *s) { h(s, 0); }
static void op_z0(state_t *s) { z(s, 0); }


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
 * =====================================================================================================================
 * sample_matrix tests
 * =====================================================================================================================
 */

/*
 * L=1, sizes={1}, single identity op.
 * K=1, S = [<iota|H|iota>] = [mean(iota, obs)]
 */
void test_sample_matrix_pure_trivial(void) {
    for (qubit_t n = 1; n <= MAXQUBITS; ++n) {
        state_t s = {.type = PURE, .data = NULL, .qubits = 0};
        state_plus(&s, n);

        const idx_t dim = POW2(n, idx_t);
        double obs[dim];
        fill_identity_obs(obs, dim);

        qop_t layer0[] = {op_identity};
        qop_t *ops[] = {layer0};
        idx_t sizes[] = {1};

        cplx_t out[1];
        sample_matrix(&s, obs, ops, sizes, 1, out);

        double expected = mean(&s, obs);
        TEST_ASSERT_DOUBLE_WITHIN(PRECISION, expected, creal(out[0]));
        TEST_ASSERT_DOUBLE_WITHIN(PRECISION, 0.0,      cimag(out[0]));

        state_free(&s);
    }
}

/*
 * 1 qubit, |iota> = |0>, obs = {0, 1}, L=1, ops = {I, X}.
 * K=2.  V_0 = I -> |0>,  V_1 = X -> |1>
 * S = [[0, 0], [0, 1]]
 * Packed (col-major lower): {0, 0, 1}
 */
void test_sample_matrix_pure_two_ops(void) {
    state_t s = {.type = PURE, .data = NULL, .qubits = 0};
    state_init(&s, 1, NULL);
    state_set(&s, 0, 0, 1.0 + 0.0 * I);

    double obs[2] = {0.0, 1.0};

    qop_t layer0[] = {op_identity, op_x0};
    qop_t *ops[] = {layer0};
    idx_t sizes[] = {2};

    cplx_t out[3];
    sample_matrix(&s, obs, ops, sizes, 1, out);

    cplx_t expected[3] = {0.0, 0.0, 1.0};
    TEST_ASSERT_EQUAL_COMPLEX_ARRAY_TOL(expected, out, 3, PRECISION);

    state_free(&s);
}

/*
 * Known-value test for sample matrix off-diagonals.
 * 1 qubit, |iota> = |+> = (1/sqrt2)(|0> + |1>), obs = {0, 1}, L=1, ops = {I, H}.
 *
 * V_0 = I:  |phi_0> = |+> = (1/sqrt2, 1/sqrt2)
 * V_1 = H:  |phi_1> = H|+> = |0> = (1, 0)
 *
 * S[0,0] = <+|H|+> = 0*(1/2) + 1*(1/2) = 0.5
 * S[1,0] = <0|H|+> = 0*1*(1/sqrt2) + 1*0*(1/sqrt2) = 0
 * S[1,1] = <0|H|0> = 0*1 + 1*0 = 0
 *
 * Packed (col-major lower): {0.5, 0, 0}
 */
void test_sample_matrix_pure_known_values(void) {
    state_t s = {.type = PURE, .data = NULL, .qubits = 0};
    state_plus(&s, 1);

    double obs[2] = {0.0, 1.0};

    qop_t layer0[] = {op_identity, op_h0};
    qop_t *ops[] = {layer0};
    idx_t sizes[] = {2};

    cplx_t out[3];
    sample_matrix(&s, obs, ops, sizes, 1, out);

    cplx_t expected[3] = {0.5, 0.0, 0.0};
    TEST_ASSERT_EQUAL_COMPLEX_ARRAY_TOL(expected, out, 3, PRECISION);

    state_free(&s);
}

/*
 * Multi-layer with non-trivial ops to verify index decomposition.
 * L=2, sizes={2, 2}, K=4. k_1 fastest.
 *
 * 1 qubit, |iota> = |0>, obs = {0, 1}.
 * ops[0] = {I, X},  ops[1] = {I, X}.
 *
 * Index decomposition (k_1 least significant):
 *   k=0: k_1=0, k_2=0 -> V=I*I -> |0>     k=1: k_1=1, k_2=0 -> V=I*X -> |1>
 *   k=2: k_1=0, k_2=1 -> V=X*I -> |1>     k=3: k_1=1, k_2=1 -> V=X*X -> |0>
 *
 * S (full, obs={0,1}):
 *        |0> |1> |1> |0>
 *   |0> [ 0   0   0   0 ]
 *   |1> [ 0   1   1   0 ]
 *   |1> [ 0   1   1   0 ]
 *   |0> [ 0   0   0   0 ]
 *
 * Packed lower triangle (col-major): {0, 0, 0, 0, 1, 1, 0, 1, 0, 0}
 */
void test_sample_matrix_pure_multi_layer(void) {
    state_t s = {.type = PURE, .data = NULL, .qubits = 0};
    state_init(&s, 1, NULL);
    state_set(&s, 0, 0, 1.0 + 0.0 * I);

    double obs[2] = {0.0, 1.0};

    qop_t layer0[] = {op_identity, op_x0};
    qop_t layer1[] = {op_identity, op_x0};
    qop_t *ops[] = {layer0, layer1};
    idx_t sizes[] = {2, 2};

    const idx_t K = 4;
    const idx_t packed_len = K * (K + 1) / 2;  /* 10 */
    cplx_t out[10];
    sample_matrix(&s, obs, ops, sizes, 2, out);

    cplx_t expected[10] = {0, 0, 0, 0, 1, 1, 0, 1, 0, 0};
    TEST_ASSERT_EQUAL_COMPLEX_ARRAY_TOL(expected, out, packed_len, PRECISION);

    state_free(&s);
}

/*
 * Asymmetric multi-layer test — discriminates k_1-fastest from k_L-fastest.
 *
 * 1 qubit, |iota> = |0>, obs = {0, 1}.
 * L=2, ops[0] = {I, X, H} (size 3), ops[1] = {I, Z} (size 2).
 * K = 3 * 2 = 6.  k_1 fastest (layer 0).
 *
 * k=0: (0,0) -> I*I|0> = |0>,  mean=0     k=1: (1,0) -> I*X|0> = |1>,  mean=1
 * k=2: (2,0) -> I*H|0> = |+>,  mean=0.5   k=3: (0,1) -> Z*I|0> = |0>,  mean=0
 * k=4: (1,1) -> Z*X|0> = -|1>, mean=1     k=5: (2,1) -> Z*H|0> = |->,  mean=0.5
 *
 * If decomposition were wrong (k_2 fastest):
 *   k=1 would decode as (k_1=0, k_2=1) -> Z*I|0> = |0>, mean=0 (not 1).
 */
void test_sample_matrix_pure_asymmetric_layers(void) {
    state_t s = {.type = PURE, .data = NULL, .qubits = 0};
    state_init(&s, 1, NULL);
    state_set(&s, 0, 0, 1.0 + 0.0 * I);

    double obs[2] = {0.0, 1.0};

    qop_t layer0[] = {op_identity, op_x0, op_h0};
    qop_t layer1[] = {op_identity, op_z0};
    qop_t *ops[] = {layer0, layer1};
    idx_t sizes[] = {3, 2};

    const idx_t K = 6;
    cplx_t out[21];  /* K*(K+1)/2 = 21 */
    sample_matrix(&s, obs, ops, sizes, 2, out);

    /* col-major lower triangle — diagonal at positions:
     * col 0: 6 entries, pos 0..5,   diag at 0
     * col 1: 5 entries, pos 6..10,  diag at 6
     * col 2: 4 entries, pos 11..14, diag at 11
     * col 3: 3 entries, pos 15..17, diag at 15
     * col 4: 2 entries, pos 18..19, diag at 18
     * col 5: 1 entry,   pos 20,     diag at 20
     */
    const double diag_expected[6] = {0.0, 1.0, 0.5, 0.0, 1.0, 0.5};
    const idx_t  diag_pos[6]      = {0, 6, 11, 15, 18, 20};
    for (int d = 0; d < 6; ++d) {
        TEST_ASSERT_DOUBLE_WITHIN(PRECISION, diag_expected[d], creal(out[diag_pos[d]]));
        TEST_ASSERT_DOUBLE_WITHIN(PRECISION, 0.0, cimag(out[diag_pos[d]]));
    }

    state_free(&s);
}

/*
 * State mutation safety: sample_matrix must not modify the input state.
 */
void test_sample_matrix_pure_state_unchanged(void) {
    state_t s = {.type = PURE, .data = NULL, .qubits = 0};
    state_plus(&s, 2);

    const idx_t dim = 4;  /* 2^2 */
    double obs[4];
    fill_identity_obs(obs, dim);

    /* snapshot amplitudes before call */
    cplx_t ref[4];
    for (idx_t x = 0; x < dim; ++x)
        ref[x] = s.data[x];

    qop_t layer0[] = {op_identity, op_x0};
    qop_t *ops[] = {layer0};
    idx_t sizes[] = {2};

    cplx_t out[3];
    sample_matrix(&s, obs, ops, sizes, 1, out);

    /* verify state unchanged */
    TEST_ASSERT_EQUAL_COMPLEX_ARRAY_TOL(ref, s.data, dim, PRECISION);
    TEST_ASSERT_EQUAL_INT(PURE, s.type);
    TEST_ASSERT_EQUAL_INT(2, s.qubits);

    state_free(&s);
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
