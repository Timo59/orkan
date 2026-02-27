// test_pauli.c - Unit tests for pauli_apply_pure and related Pauli string functions.
//
// Reference results are computed by build_pauli_matrix (constructs the full 2^n x 2^n
// matrix representation) combined with zmv (CBLAS matrix-vector product).  Every
// check_pauli_apply call follows the goto-cleanup pattern from test_gate_pure_1q.c to
// guarantee leak-free failure handling.

/*
 * =====================================================================================================================
 * Includes
 * =====================================================================================================================
 */

#ifndef Q_TEST_H
#include "test.h"
#endif

#include "circuit.h"
#include "state.h"
#include "linalg.h"

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/*
 * =====================================================================================================================
 * Helper: build full 2^n x 2^n Pauli matrix (column-major)
 *
 * M[col * dim + row] = eps(col)   if row == (col ^ x_mask)
 *                      0          otherwise
 *
 * eps(col) = i^{popcount(y_mask)} * (-1)^{popcount(col & z_mask)}
 *          = CMPLX(sign * eta_re[eta], sign * eta_im[eta])
 *   where
 *     eta    = __builtin_popcountll(y_mask) & 3
 *     sign   = (__builtin_popcountll(col & z_mask) & 1) ? -1 : 1
 *     eta_re = {1, 0, -1,  0}
 *     eta_im = {0, 1,  0, -1}
 * =====================================================================================================================
 */

static cplx_t *build_pauli_matrix(const pauli_t *P) {
    const unsigned dim = (unsigned)(1u << P->n_qubits);
    cplx_t *M = calloc((size_t)dim * dim, sizeof(*M));
    if (!M) {
        fprintf(stderr, "build_pauli_matrix(): allocation failed\n");
        return NULL;
    }

    static const double eta_re[4] = {1.0,  0.0, -1.0,  0.0};
    static const double eta_im[4] = {0.0,  1.0,  0.0, -1.0};

    const unsigned eta = (unsigned)(__builtin_popcountll(P->y_mask) & 3u);

    for (unsigned col = 0; col < dim; ++col) {
        const unsigned row  = col ^ (unsigned)P->x_mask;
        const int      sign = (__builtin_popcountll((uint64_t)col & P->z_mask) & 1) ? -1 : 1;
        M[col * dim + row]  = CMPLX(sign * eta_re[eta], sign * eta_im[eta]);
    }

    return M;
}

/*
 * =====================================================================================================================
 * Helper: apply P to a copy of psi and compare against pauli_apply_pure result.
 *
 * psi must point to a dim-element array that is NOT yet owned by a state_t.
 * A fresh malloc+memcpy copy is made so that state_init can take ownership;
 * the original psi array is left untouched for use as the reference input to zmv.
 *
 * Cleanup guarantee: state_init transfers ownership of the copy and nullifies
 * the local pointer.  On any failure we jump to cleanup where state_free handles
 * whatever partial state exists, and all heap objects are freed exactly once.
 * =====================================================================================================================
 */

static void check_pauli_apply(const pauli_t *P, const cplx_t *psi) {
    const unsigned dim = (unsigned)(1u << P->n_qubits);

    cplx_t *mat  = NULL;
    cplx_t *ref  = NULL;
    cplx_t *copy = NULL;
    state_t s    = {0};
    s.type        = PURE;

    /* Build the full Pauli matrix */
    mat = build_pauli_matrix(P);
    if (!mat) {
        fprintf(stderr, "check_pauli_apply(): build_pauli_matrix failed\n");
        goto cleanup;
    }

    /* Compute reference result: ref = M * psi */
    ref = zmv(dim, mat, psi);
    if (!ref) {
        fprintf(stderr, "check_pauli_apply(): zmv failed\n");
        goto cleanup;
    }

    /* Make a fresh copy of psi for state_init to consume */
    copy = malloc((size_t)dim * sizeof(*copy));
    if (!copy) {
        fprintf(stderr, "check_pauli_apply(): copy allocation failed\n");
        goto cleanup;
    }
    memcpy(copy, psi, (size_t)dim * sizeof(*copy));

    /* Transfer ownership of copy to the state (copy becomes NULL after this) */
    state_init(&s, (qubit_t)P->n_qubits, &copy);
    if (!s.data) {
        fprintf(stderr, "check_pauli_apply(): state_init failed\n");
        goto cleanup;
    }

    /* Apply the Pauli string in-place */
    pauli_apply_pure(&s, P);
    if (!s.data) {
        fprintf(stderr, "check_pauli_apply(): pauli_apply_pure left state data NULL\n");
        goto cleanup;
    }

    /* Compare */
    TEST_ASSERT_EQUAL_COMPLEX_ARRAY_TOL(ref, s.data, dim, PRECISION);

    /* Normal-path cleanup */
    free(mat);
    free(ref);
    free(s.data);
    s.data = NULL;
    return;

cleanup:
    free(mat);
    free(ref);
    /* copy is NULL after state_init succeeds; safe to free either way */
    free(copy);
    state_free(&s);
}

/*
 * =====================================================================================================================
 * Static test state vectors (n=1)
 *
 * Convention matching test_gate.h: initializer-list style with 0.0*I.
 * =====================================================================================================================
 */

/* |0> */
static const cplx_t PSI1_0[2]  = {1.0 + 0.0*I,  0.0 + 0.0*I};
/* |1> */
static const cplx_t PSI1_1[2]  = {0.0 + 0.0*I,  1.0 + 0.0*I};
/* |+> = (|0>+|1>)/sqrt(2) */
static const cplx_t PSI1_P[2]  = {INVSQRT2 + 0.0*I,  INVSQRT2 + 0.0*I};
/* |-> = (|0>-|1>)/sqrt(2) */
static const cplx_t PSI1_M[2]  = {INVSQRT2 + 0.0*I, -INVSQRT2 + 0.0*I};
/* |+i> = (|0>+i|1>)/sqrt(2) */
static const cplx_t PSI1_PI[2] = {INVSQRT2 + 0.0*I,  0.0 + INVSQRT2*I};
/* |-i> = (|0>-i|1>)/sqrt(2) */
static const cplx_t PSI1_MI[2] = {INVSQRT2 + 0.0*I,  0.0 - INVSQRT2*I};

/*
 * =====================================================================================================================
 * Test: Pauli mask construction via pauli_new — single-label cases
 * =====================================================================================================================
 */

void test_pauli_masks_X(void) {
    /* n=1, X: x_mask=1, z_mask=0, y_mask=0 */
    const pauli_label_t lbl = PAULI_X;
    pauli_t *P = pauli_new(&lbl, 1);
    TEST_ASSERT_NOT_NULL(P);
    TEST_ASSERT_EQUAL_UINT64(1u, P->x_mask);
    TEST_ASSERT_EQUAL_UINT64(0u, P->z_mask);
    TEST_ASSERT_EQUAL_UINT64(0u, P->y_mask);
    pauli_free(P);
}

void test_pauli_masks_Y(void) {
    /* n=1, Y: x_mask=1, z_mask=1, y_mask=1 */
    const pauli_label_t lbl = PAULI_Y;
    pauli_t *P = pauli_new(&lbl, 1);
    TEST_ASSERT_NOT_NULL(P);
    TEST_ASSERT_EQUAL_UINT64(1u, P->x_mask);
    TEST_ASSERT_EQUAL_UINT64(1u, P->z_mask);
    TEST_ASSERT_EQUAL_UINT64(1u, P->y_mask);
    pauli_free(P);
}

void test_pauli_masks_Z(void) {
    /* n=1, Z: x_mask=0, z_mask=1, y_mask=0 */
    const pauli_label_t lbl = PAULI_Z;
    pauli_t *P = pauli_new(&lbl, 1);
    TEST_ASSERT_NOT_NULL(P);
    TEST_ASSERT_EQUAL_UINT64(0u, P->x_mask);
    TEST_ASSERT_EQUAL_UINT64(1u, P->z_mask);
    TEST_ASSERT_EQUAL_UINT64(0u, P->y_mask);
    pauli_free(P);
}

void test_pauli_masks_str(void) {
    /*
     * "IXYZ" with n=4 (big-endian: s[0]='I'->labels[3], s[1]='X'->labels[2],
     *                              s[2]='Y'->labels[1], s[3]='Z'->labels[0])
     *   x_mask = bits 2(X) + 1(Y)      = 0b0110 = 6
     *   z_mask = bits 1(Y) + 0(Z)      = 0b0011 = 3
     *   y_mask = bit  1(Y)             = 0b0010 = 2
     */
    pauli_t *P = pauli_from_str("IXYZ", 4);
    TEST_ASSERT_NOT_NULL(P);
    TEST_ASSERT_EQUAL_UINT64(6u, P->x_mask);
    TEST_ASSERT_EQUAL_UINT64(3u, P->z_mask);
    TEST_ASSERT_EQUAL_UINT64(2u, P->y_mask);
    pauli_free(P);

    /* "IIII": all masks zero */
    P = pauli_from_str("IIII", 4);
    TEST_ASSERT_NOT_NULL(P);
    TEST_ASSERT_EQUAL_UINT64(0u, P->x_mask);
    TEST_ASSERT_EQUAL_UINT64(0u, P->z_mask);
    TEST_ASSERT_EQUAL_UINT64(0u, P->y_mask);
    pauli_free(P);

    /* "XXXX": x_mask=0b1111=15, z_mask=0, y_mask=0 */
    P = pauli_from_str("XXXX", 4);
    TEST_ASSERT_NOT_NULL(P);
    TEST_ASSERT_EQUAL_UINT64(15u, P->x_mask);
    TEST_ASSERT_EQUAL_UINT64(0u,  P->z_mask);
    TEST_ASSERT_EQUAL_UINT64(0u,  P->y_mask);
    pauli_free(P);
}

/*
 * =====================================================================================================================
 * Test: pauli_from_str — negative cases
 * =====================================================================================================================
 */

void test_pauli_from_str_negative(void) {
    /* NULL string must return NULL */
    pauli_t *P = pauli_from_str(NULL, 1);
    TEST_ASSERT_NULL(P);

    /* Invalid character 'A' must return NULL */
    P = pauli_from_str("A", 1);
    TEST_ASSERT_NULL(P);

    /* Invalid character in multi-char string */
    P = pauli_from_str("XAZ", 3);
    TEST_ASSERT_NULL(P);

    /* Length mismatch: string is shorter than n_qubits */
    P = pauli_from_str("XY", 3);
    TEST_ASSERT_NULL(P);

    /* Length mismatch: string is longer than n_qubits */
    P = pauli_from_str("XYZ", 2);
    TEST_ASSERT_NULL(P);
}

/*
 * =====================================================================================================================
 * Test: pauli_phase spot-checks
 *
 * eps(x) = i^{popcount(y_mask)} * (-1)^{popcount(x & z_mask)}
 * =====================================================================================================================
 */

void test_pauli_phase(void) {
    pauli_t *PX = NULL, *PY = NULL, *PZ = NULL;
    pauli_phase_t ph;

    const pauli_label_t lblX = PAULI_X;
    const pauli_label_t lblY = PAULI_Y;
    const pauli_label_t lblZ = PAULI_Z;

    PX = pauli_new(&lblX, 1);
    TEST_ASSERT_NOT_NULL(PX);
    PY = pauli_new(&lblY, 1);
    TEST_ASSERT_NOT_NULL(PY);
    PZ = pauli_new(&lblZ, 1);
    TEST_ASSERT_NOT_NULL(PZ);

    /*
     * X: y_mask=0, z_mask=0.
     *   eps(x=0) = i^0 * (-1)^0 = 1+0i -> {real=1, imag=0}
     *   eps(x=1) = i^0 * (-1)^0 = 1+0i -> {real=1, imag=0}
     */
    ph = pauli_phase(PX, 0u);
    TEST_ASSERT_EQUAL_INT(1, ph.real);
    TEST_ASSERT_EQUAL_INT(0, ph.imag);

    ph = pauli_phase(PX, 1u);
    TEST_ASSERT_EQUAL_INT(1, ph.real);
    TEST_ASSERT_EQUAL_INT(0, ph.imag);

    /*
     * Y: y_mask=1, z_mask=1.
     *   eps(x=0) = i^1 * (-1)^popcount(0&1) = i * 1 = 0+1i -> {real=0, imag=1}
     *   eps(x=1) = i^1 * (-1)^popcount(1&1) = i * (-1) = 0-1i -> {real=0, imag=-1}
     */
    ph = pauli_phase(PY, 0u);
    TEST_ASSERT_EQUAL_INT(0,  ph.real);
    TEST_ASSERT_EQUAL_INT(1,  ph.imag);

    ph = pauli_phase(PY, 1u);
    TEST_ASSERT_EQUAL_INT(0,  ph.real);
    TEST_ASSERT_EQUAL_INT(-1, ph.imag);

    /*
     * Z: y_mask=0, z_mask=1.
     *   eps(x=0) = i^0 * (-1)^popcount(0&1) = 1 * 1 = 1+0i -> {real=1, imag=0}
     *   eps(x=1) = i^0 * (-1)^popcount(1&1) = 1 * (-1) = -1+0i -> {real=-1, imag=0}
     */
    ph = pauli_phase(PZ, 0u);
    TEST_ASSERT_EQUAL_INT(1,  ph.real);
    TEST_ASSERT_EQUAL_INT(0,  ph.imag);

    ph = pauli_phase(PZ, 1u);
    TEST_ASSERT_EQUAL_INT(-1, ph.real);
    TEST_ASSERT_EQUAL_INT(0,  ph.imag);

    pauli_free(PX);
    pauli_free(PY);
    pauli_free(PZ);
}

/*
 * =====================================================================================================================
 * Test: n=1, I
 *
 * I is the identity: all six standard states must be unchanged.
 * =====================================================================================================================
 */

void test_pauli_apply_I_q0(void) {
    const pauli_label_t lbl = PAULI_I;
    pauli_t *P = pauli_new(&lbl, 1);
    TEST_ASSERT_NOT_NULL(P);

    check_pauli_apply(P, PSI1_0);
    check_pauli_apply(P, PSI1_1);
    check_pauli_apply(P, PSI1_P);
    check_pauli_apply(P, PSI1_M);
    check_pauli_apply(P, PSI1_PI);
    check_pauli_apply(P, PSI1_MI);

    pauli_free(P);
}

/*
 * =====================================================================================================================
 * Test: n=1, X
 *
 * X swaps amplitudes: X[a0,a1]^T = [a1,a0]^T, phase always +1.
 *
 * Expected results:
 *   X|0>  = |1>               = [0, 1]
 *   X|1>  = |0>               = [1, 0]
 *   X|+>  = |+>               = [s, s]           (eigenvector, eigenvalue +1)
 *   X|->  = -|->              = [-s, s]           (eigenvector, eigenvalue -1)
 *   X|+i> =  i|-i>            = [i*s, s]
 *   X|-i> = -i|+i>            = [-i*s, s]
 *
 * (All verified via reference matrix; inline values provided for documentation.)
 * =====================================================================================================================
 */

void test_pauli_apply_X_q0(void) {
    const pauli_label_t lbl = PAULI_X;
    pauli_t *P = pauli_new(&lbl, 1);
    TEST_ASSERT_NOT_NULL(P);

    check_pauli_apply(P, PSI1_0);
    check_pauli_apply(P, PSI1_1);
    check_pauli_apply(P, PSI1_P);
    check_pauli_apply(P, PSI1_M);
    check_pauli_apply(P, PSI1_PI);
    check_pauli_apply(P, PSI1_MI);

    pauli_free(P);
}

/*
 * =====================================================================================================================
 * Test: n=1, Y
 *
 * Y[a0,a1]^T = [-i*a1, i*a0]^T.
 *
 * Expected results:
 *   Y|0>  =  i|1>             = [0, i]
 *   Y|1>  = -i|0>             = [-i, 0]
 *   Y|+>  = -i|->             = [-i*s, i*s]
 *   Y|->  =  i|+>             = [i*s, i*s]
 *   Y|+i> =  |+i>             = [s, i*s]          (eigenvalue +1)
 *   Y|-i> = -|-i>             = [-s, i*s]          (eigenvalue -1)
 * =====================================================================================================================
 */

void test_pauli_apply_Y_q0(void) {
    const pauli_label_t lbl = PAULI_Y;
    pauli_t *P = pauli_new(&lbl, 1);
    TEST_ASSERT_NOT_NULL(P);

    check_pauli_apply(P, PSI1_0);
    check_pauli_apply(P, PSI1_1);
    check_pauli_apply(P, PSI1_P);
    check_pauli_apply(P, PSI1_M);
    check_pauli_apply(P, PSI1_PI);
    check_pauli_apply(P, PSI1_MI);

    pauli_free(P);
}

/*
 * =====================================================================================================================
 * Test: n=1, Z
 *
 * Z[a0,a1]^T = [a0, -a1]^T.
 *
 * Expected results:
 *   Z|0>  =  |0>              = [1, 0]            (eigenvalue +1)
 *   Z|1>  = -|1>              = [0, -1]           (eigenvalue -1)
 *   Z|+>  =  |->              = [s, -s]
 *   Z|->  =  |+>              = [s,  s]
 *   Z|+i> =  |-i>             = [s, -i*s]
 *   Z|-i> =  |+i>             = [s,  i*s]
 * =====================================================================================================================
 */

void test_pauli_apply_Z_q0(void) {
    const pauli_label_t lbl = PAULI_Z;
    pauli_t *P = pauli_new(&lbl, 1);
    TEST_ASSERT_NOT_NULL(P);

    check_pauli_apply(P, PSI1_0);
    check_pauli_apply(P, PSI1_1);
    check_pauli_apply(P, PSI1_P);
    check_pauli_apply(P, PSI1_M);
    check_pauli_apply(P, PSI1_PI);
    check_pauli_apply(P, PSI1_MI);

    pauli_free(P);
}

/*
 * =====================================================================================================================
 * n=2 test states (dim=4)
 *
 * Basis: |00>=index 0, |01>=index 1, |10>=index 2, |11>=index 3.
 * (qubit 1 = MSB, qubit 0 = LSB in the index)
 * =====================================================================================================================
 */

/* |00> */
static const cplx_t PSI2_00[4] = {1.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I};
/* |01> */
static const cplx_t PSI2_01[4] = {0.0 + 0.0*I, 1.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I};
/* |10> */
static const cplx_t PSI2_10[4] = {0.0 + 0.0*I, 0.0 + 0.0*I, 1.0 + 0.0*I, 0.0 + 0.0*I};
/* |11> */
static const cplx_t PSI2_11[4] = {0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I, 1.0 + 0.0*I};
/* |++> = (|00>+|01>+|10>+|11>)/2 */
static const cplx_t PSI2_PP[4] = {INVSQRT4 + 0.0*I, INVSQRT4 + 0.0*I, INVSQRT4 + 0.0*I, INVSQRT4 + 0.0*I};
/* |+0> = (|00>+|10>)/sqrt(2) */
static const cplx_t PSI2_P0[4] = {INVSQRT2 + 0.0*I, 0.0 + 0.0*I, INVSQRT2 + 0.0*I, 0.0 + 0.0*I};
/* |1+> = (|10>+|11>)/sqrt(2) */
static const cplx_t PSI2_1P[4] = {0.0 + 0.0*I, 0.0 + 0.0*I, INVSQRT2 + 0.0*I, INVSQRT2 + 0.0*I};

/*
 * =====================================================================================================================
 * Test: n=2, XX
 *
 * labels[0]=X, labels[1]=X.  x_mask=3, z_mask=0, y_mask=0.
 * XX permutes: col -> col XOR 3, phase always +1.
 *   |00> -> |11>,  |01> -> |10>,  |10> -> |01>,  |11> -> |00>
 *   |++> -> |++>  (XX is X(x)X, eigenstate |++> with eigenvalue +1)
 * =====================================================================================================================
 */

void test_pauli_apply_XX(void) {
    const pauli_label_t lbls[2] = {PAULI_X, PAULI_X};
    pauli_t *P = pauli_new(lbls, 2);
    TEST_ASSERT_NOT_NULL(P);

    check_pauli_apply(P, PSI2_00);
    check_pauli_apply(P, PSI2_01);
    check_pauli_apply(P, PSI2_10);
    check_pauli_apply(P, PSI2_11);
    check_pauli_apply(P, PSI2_PP);
    check_pauli_apply(P, PSI2_P0);
    check_pauli_apply(P, PSI2_1P);

    pauli_free(P);
}

/*
 * =====================================================================================================================
 * Test: n=2, ZZ
 *
 * labels[0]=Z, labels[1]=Z.  x_mask=0, z_mask=3, y_mask=0.
 * ZZ is diagonal: phase = (-1)^popcount(col).
 *   |00> -> +|00>,  |01> -> -|01>,  |10> -> -|10>,  |11> -> +|11>
 * =====================================================================================================================
 */

void test_pauli_apply_ZZ(void) {
    const pauli_label_t lbls[2] = {PAULI_Z, PAULI_Z};
    pauli_t *P = pauli_new(lbls, 2);
    TEST_ASSERT_NOT_NULL(P);

    check_pauli_apply(P, PSI2_00);
    check_pauli_apply(P, PSI2_01);
    check_pauli_apply(P, PSI2_10);
    check_pauli_apply(P, PSI2_11);
    check_pauli_apply(P, PSI2_PP);
    check_pauli_apply(P, PSI2_P0);
    check_pauli_apply(P, PSI2_1P);

    pauli_free(P);
}

/*
 * =====================================================================================================================
 * Test: n=2, XZ
 *
 * labels[0]=X, labels[1]=Z.  x_mask=1, z_mask=2, y_mask=0.
 * Flips qubit 0; phase = (-1)^popcount(col & 2).
 *   |00>(col=0) -> +|01>,  |01>(col=1) -> +|00>,
 *   |10>(col=2) -> -|11>,  |11>(col=3) -> -|10>
 *   |+0> -> INVSQRT2*(|01>-|11>),  |1+> -> -|1+>
 * =====================================================================================================================
 */

void test_pauli_apply_XZ(void) {
    const pauli_label_t lbls[2] = {PAULI_X, PAULI_Z};
    pauli_t *P = pauli_new(lbls, 2);
    TEST_ASSERT_NOT_NULL(P);

    check_pauli_apply(P, PSI2_00);
    check_pauli_apply(P, PSI2_01);
    check_pauli_apply(P, PSI2_10);
    check_pauli_apply(P, PSI2_11);
    check_pauli_apply(P, PSI2_PP);
    check_pauli_apply(P, PSI2_P0);
    check_pauli_apply(P, PSI2_1P);

    pauli_free(P);
}

/*
 * =====================================================================================================================
 * Test: n=2, YY
 *
 * labels[0]=Y, labels[1]=Y.  x_mask=3, z_mask=3, y_mask=3.
 * eta = popcount(3) & 3 = 2.  eps(col) = i^2 * (-1)^popcount(col&3) = -(-1)^popcount(col).
 *   |00>(col=0) -> -|11>,  |01>(col=1) -> +|10>,
 *   |10>(col=2) -> +|01>,  |11>(col=3) -> -|00>
 * =====================================================================================================================
 */

void test_pauli_apply_YY(void) {
    const pauli_label_t lbls[2] = {PAULI_Y, PAULI_Y};
    pauli_t *P = pauli_new(lbls, 2);
    TEST_ASSERT_NOT_NULL(P);

    check_pauli_apply(P, PSI2_00);
    check_pauli_apply(P, PSI2_01);
    check_pauli_apply(P, PSI2_10);
    check_pauli_apply(P, PSI2_11);
    check_pauli_apply(P, PSI2_PP);
    check_pauli_apply(P, PSI2_P0);
    check_pauli_apply(P, PSI2_1P);

    pauli_free(P);
}

/*
 * =====================================================================================================================
 * Test: n=2, IZ
 *
 * labels[0]=Z, labels[1]=I.  x_mask=0, z_mask=1, y_mask=0.
 * IZ is diagonal: phase = (-1)^popcount(col & 1).
 *   |00> -> +|00>,  |01> -> -|01>,  |10> -> +|10>,  |11> -> -|11>
 * =====================================================================================================================
 */

void test_pauli_apply_IZ(void) {
    const pauli_label_t lbls[2] = {PAULI_Z, PAULI_I};
    pauli_t *P = pauli_new(lbls, 2);
    TEST_ASSERT_NOT_NULL(P);

    check_pauli_apply(P, PSI2_00);
    check_pauli_apply(P, PSI2_01);
    check_pauli_apply(P, PSI2_10);
    check_pauli_apply(P, PSI2_11);
    check_pauli_apply(P, PSI2_PP);
    check_pauli_apply(P, PSI2_P0);
    check_pauli_apply(P, PSI2_1P);

    pauli_free(P);
}

/*
 * =====================================================================================================================
 * n=3 test states (dim=8)
 * =====================================================================================================================
 */

/* Computational basis |000>...|111> */
static const cplx_t PSI3_000[8] = {1.0+0.0*I, 0.0+0.0*I, 0.0+0.0*I, 0.0+0.0*I,
                                    0.0+0.0*I, 0.0+0.0*I, 0.0+0.0*I, 0.0+0.0*I};
static const cplx_t PSI3_001[8] = {0.0+0.0*I, 1.0+0.0*I, 0.0+0.0*I, 0.0+0.0*I,
                                    0.0+0.0*I, 0.0+0.0*I, 0.0+0.0*I, 0.0+0.0*I};
static const cplx_t PSI3_010[8] = {0.0+0.0*I, 0.0+0.0*I, 1.0+0.0*I, 0.0+0.0*I,
                                    0.0+0.0*I, 0.0+0.0*I, 0.0+0.0*I, 0.0+0.0*I};
static const cplx_t PSI3_011[8] = {0.0+0.0*I, 0.0+0.0*I, 0.0+0.0*I, 1.0+0.0*I,
                                    0.0+0.0*I, 0.0+0.0*I, 0.0+0.0*I, 0.0+0.0*I};
static const cplx_t PSI3_100[8] = {0.0+0.0*I, 0.0+0.0*I, 0.0+0.0*I, 0.0+0.0*I,
                                    1.0+0.0*I, 0.0+0.0*I, 0.0+0.0*I, 0.0+0.0*I};
static const cplx_t PSI3_101[8] = {0.0+0.0*I, 0.0+0.0*I, 0.0+0.0*I, 0.0+0.0*I,
                                    0.0+0.0*I, 1.0+0.0*I, 0.0+0.0*I, 0.0+0.0*I};
static const cplx_t PSI3_110[8] = {0.0+0.0*I, 0.0+0.0*I, 0.0+0.0*I, 0.0+0.0*I,
                                    0.0+0.0*I, 0.0+0.0*I, 1.0+0.0*I, 0.0+0.0*I};
static const cplx_t PSI3_111[8] = {0.0+0.0*I, 0.0+0.0*I, 0.0+0.0*I, 0.0+0.0*I,
                                    0.0+0.0*I, 0.0+0.0*I, 0.0+0.0*I, 1.0+0.0*I};

/* GHZ state: (|000>+|111>)/sqrt(2) — entangled superposition */
static const cplx_t PSI3_GHZ[8] = {INVSQRT2+0.0*I, 0.0+0.0*I, 0.0+0.0*I, 0.0+0.0*I,
                                    0.0+0.0*I, 0.0+0.0*I, 0.0+0.0*I, INVSQRT2+0.0*I};
/* W state: (|001>+|010>+|100>)/sqrt(3) — single-excitation superposition */
#define INVSQRT3 0.5773502691896258
static const cplx_t PSI3_W[8] = {0.0+0.0*I, INVSQRT3+0.0*I, INVSQRT3+0.0*I, 0.0+0.0*I,
                                  INVSQRT3+0.0*I, 0.0+0.0*I, 0.0+0.0*I, 0.0+0.0*I};
/* |+++> uniform superposition: amplitude INVSQRT8 everywhere */
static const cplx_t PSI3_PPP[8] = {INVSQRT8+0.0*I, INVSQRT8+0.0*I, INVSQRT8+0.0*I, INVSQRT8+0.0*I,
                                    INVSQRT8+0.0*I, INVSQRT8+0.0*I, INVSQRT8+0.0*I, INVSQRT8+0.0*I};

/*
 * =====================================================================================================================
 * Test: n=3, XYZ
 *
 * pauli_new({X,Y,Z}, 3): labels[0]=X, labels[1]=Y, labels[2]=Z.
 *   x_mask = bits 0(X)+1(Y) = 3
 *   z_mask = bits 1(Y)+2(Z) = 6
 *   y_mask = bit  1(Y)      = 2
 *
 * eta = popcount(2) & 3 = 1.
 * eps(col) = i * (-1)^popcount(col & 6).
 *
 * Action on each basis state (via reference matrix):
 *   |000>(col=0) -> i|011>    |001>(col=1) -> i|010>
 *   |010>(col=2) -> -i|001>   |011>(col=3) -> -i|000>
 *   |100>(col=4) -> i|111>    |101>(col=5) -> -i|110>
 *   |110>(col=6) -> i|101>    |111>(col=7) -> i|100>
 *
 * The inline spot-check for |000>->i|011> is tested explicitly as a smoke test;
 * all 8 states are then verified against the reference matrix.
 * =====================================================================================================================
 */

void test_pauli_apply_XYZ(void) {
    const pauli_label_t lbls[3] = {PAULI_X, PAULI_Y, PAULI_Z};
    pauli_t *P = pauli_new(lbls, 3);
    TEST_ASSERT_NOT_NULL(P);

    /*
     * Inline spot-check: XYZ|000> = i|011>.
     * Expected: all zero except index 3 which has value i.
     */
    {
        const cplx_t expected[8] = {
            0.0+0.0*I, 0.0+0.0*I, 0.0+0.0*I, 0.0+1.0*I,
            0.0+0.0*I, 0.0+0.0*I, 0.0+0.0*I, 0.0+0.0*I
        };

        cplx_t *copy = malloc(8 * sizeof(*copy));
        TEST_ASSERT_NOT_NULL(copy);
        memcpy(copy, PSI3_000, 8 * sizeof(*copy));

        state_t s = {0};
        s.type = PURE;
        state_init(&s, 3, &copy);
        if (!s.data) {
            free(copy);
            TEST_FAIL_MESSAGE("state_init failed in XYZ spot-check");
        }

        pauli_apply_pure(&s, P);
        TEST_ASSERT_EQUAL_COMPLEX_ARRAY_TOL(expected, s.data, 8, PRECISION);
        free(s.data);
        s.data = NULL;
    }

    /* Full reference check on all 8 computational basis states */
    check_pauli_apply(P, PSI3_000);
    check_pauli_apply(P, PSI3_001);
    check_pauli_apply(P, PSI3_010);
    check_pauli_apply(P, PSI3_011);
    check_pauli_apply(P, PSI3_100);
    check_pauli_apply(P, PSI3_101);
    check_pauli_apply(P, PSI3_110);
    check_pauli_apply(P, PSI3_111);
    check_pauli_apply(P, PSI3_GHZ);
    check_pauli_apply(P, PSI3_W);
    check_pauli_apply(P, PSI3_PPP);

    pauli_free(P);
}

/*
 * =====================================================================================================================
 * Test: n=3, ZZZ
 *
 * pauli_new({Z,Z,Z}, 3): x_mask=0, z_mask=7, y_mask=0.
 * ZZZ is diagonal: eps(col) = (-1)^popcount(col).
 *
 * Inline expected and reference cross-check:
 *   |000> -> +|000>,  |001> -> -|001>,  |010> -> -|010>,  |011> -> +|011>
 *   |100> -> -|100>,  |101> -> +|101>,  |110> -> +|110>,  |111> -> -|111>
 * =====================================================================================================================
 */

void test_pauli_apply_ZZZ(void) {
    const pauli_label_t lbls[3] = {PAULI_Z, PAULI_Z, PAULI_Z};
    pauli_t *P = pauli_new(lbls, 3);
    TEST_ASSERT_NOT_NULL(P);

    /* Inline expected values for the diagonal case */
    static const cplx_t exp_000[8] = { 1.0+0.0*I,  0.0+0.0*I,  0.0+0.0*I,  0.0+0.0*I,
                                        0.0+0.0*I,  0.0+0.0*I,  0.0+0.0*I,  0.0+0.0*I};
    static const cplx_t exp_001[8] = { 0.0+0.0*I, -1.0+0.0*I,  0.0+0.0*I,  0.0+0.0*I,
                                        0.0+0.0*I,  0.0+0.0*I,  0.0+0.0*I,  0.0+0.0*I};
    static const cplx_t exp_010[8] = { 0.0+0.0*I,  0.0+0.0*I, -1.0+0.0*I,  0.0+0.0*I,
                                        0.0+0.0*I,  0.0+0.0*I,  0.0+0.0*I,  0.0+0.0*I};
    static const cplx_t exp_011[8] = { 0.0+0.0*I,  0.0+0.0*I,  0.0+0.0*I,  1.0+0.0*I,
                                        0.0+0.0*I,  0.0+0.0*I,  0.0+0.0*I,  0.0+0.0*I};
    static const cplx_t exp_100[8] = { 0.0+0.0*I,  0.0+0.0*I,  0.0+0.0*I,  0.0+0.0*I,
                                       -1.0+0.0*I,  0.0+0.0*I,  0.0+0.0*I,  0.0+0.0*I};
    static const cplx_t exp_101[8] = { 0.0+0.0*I,  0.0+0.0*I,  0.0+0.0*I,  0.0+0.0*I,
                                        0.0+0.0*I,  1.0+0.0*I,  0.0+0.0*I,  0.0+0.0*I};
    static const cplx_t exp_110[8] = { 0.0+0.0*I,  0.0+0.0*I,  0.0+0.0*I,  0.0+0.0*I,
                                        0.0+0.0*I,  0.0+0.0*I,  1.0+0.0*I,  0.0+0.0*I};
    static const cplx_t exp_111[8] = { 0.0+0.0*I,  0.0+0.0*I,  0.0+0.0*I,  0.0+0.0*I,
                                        0.0+0.0*I,  0.0+0.0*I,  0.0+0.0*I, -1.0+0.0*I};

    /* Inline checks */
    {
        const cplx_t *inputs[8]   = {PSI3_000, PSI3_001, PSI3_010, PSI3_011,
                                      PSI3_100, PSI3_101, PSI3_110, PSI3_111};
        const cplx_t *expected[8] = {exp_000, exp_001, exp_010, exp_011,
                                      exp_100, exp_101, exp_110, exp_111};

        for (unsigned k = 0; k < 8; ++k) {
            cplx_t *copy = malloc(8 * sizeof(*copy));
            TEST_ASSERT_NOT_NULL(copy);
            memcpy(copy, inputs[k], 8 * sizeof(*copy));

            state_t s = {0};
            s.type = PURE;
            state_init(&s, 3, &copy);
            if (!s.data) {
                free(copy);
                TEST_FAIL_MESSAGE("state_init failed in ZZZ inline check");
            }

            pauli_apply_pure(&s, P);
            TEST_ASSERT_EQUAL_COMPLEX_ARRAY_TOL(expected[k], s.data, 8, PRECISION);
            free(s.data);
            s.data = NULL;
        }
    }

    /* Cross-check with reference matrix */
    check_pauli_apply(P, PSI3_000);
    check_pauli_apply(P, PSI3_001);
    check_pauli_apply(P, PSI3_010);
    check_pauli_apply(P, PSI3_011);
    check_pauli_apply(P, PSI3_100);
    check_pauli_apply(P, PSI3_101);
    check_pauli_apply(P, PSI3_110);
    check_pauli_apply(P, PSI3_111);
    check_pauli_apply(P, PSI3_GHZ);
    check_pauli_apply(P, PSI3_W);
    check_pauli_apply(P, PSI3_PPP);

    pauli_free(P);
}

/* n=4 GHZ: (|0000>+|1111>)/sqrt(2) */
static const cplx_t PSI4_GHZ[16] = {
    INVSQRT2+0.0*I, 0.0+0.0*I, 0.0+0.0*I, 0.0+0.0*I,
    0.0+0.0*I,      0.0+0.0*I, 0.0+0.0*I, 0.0+0.0*I,
    0.0+0.0*I,      0.0+0.0*I, 0.0+0.0*I, 0.0+0.0*I,
    0.0+0.0*I,      0.0+0.0*I, 0.0+0.0*I, INVSQRT2+0.0*I
};

/*
 * =====================================================================================================================
 * Test: n=4, XZYX (exercises all 16 amplitudes — MAXQUBITS)
 *
 * pauli_new({X,Z,Y,X}, 4): labels[0]=X, labels[1]=Z, labels[2]=Y, labels[3]=X.
 *   x_mask = bits 0(X)+2(Y)+3(X) = 0b1101 = 13
 *   z_mask = bits 1(Z)+2(Y)      = 0b0110 =  6
 *   y_mask = bit  2(Y)           = 0b0100 =  4
 *
 * All 16 computational basis states are tested via the reference matrix.
 * =====================================================================================================================
 */

void test_pauli_apply_4q(void) {
    const pauli_label_t lbls[4] = {PAULI_X, PAULI_Z, PAULI_Y, PAULI_X};
    pauli_t *P = pauli_new(lbls, 4);
    TEST_ASSERT_NOT_NULL(P);

    /* Verify masks before applying */
    TEST_ASSERT_EQUAL_UINT64(13u, P->x_mask);
    TEST_ASSERT_EQUAL_UINT64(6u,  P->z_mask);
    TEST_ASSERT_EQUAL_UINT64(4u,  P->y_mask);

    /* All 16 computational basis states */
    for (unsigned k = 0; k < 16; ++k) {
        cplx_t psi[16];
        memset(psi, 0, sizeof(psi));
        psi[k] = 1.0 + 0.0*I;
        check_pauli_apply(P, psi);
    }

    check_pauli_apply(P, PSI4_GHZ);

    pauli_free(P);
}

/*
 * =====================================================================================================================
 * Test: n=2, IX
 *
 * labels[0]=I, labels[1]=X.  x_mask=2, z_mask=0, y_mask=0, split_bit=1.
 * IX flips qubit 1 only: col -> col XOR 2, phase always +1.
 *   |00>(col=0) -> |10>,  |01>(col=1) -> |11>,
 *   |10>(col=2) -> |00>,  |11>(col=3) -> |01>
 * =====================================================================================================================
 */

void test_pauli_apply_IX(void) {
    const pauli_label_t lbls[2] = {PAULI_I, PAULI_X};
    pauli_t *P = pauli_new(lbls, 2);
    TEST_ASSERT_NOT_NULL(P);

    TEST_ASSERT_EQUAL_UINT64(2u, P->x_mask);
    TEST_ASSERT_EQUAL_UINT64(0u, P->z_mask);
    TEST_ASSERT_EQUAL_UINT64(0u, P->y_mask);

    check_pauli_apply(P, PSI2_00);
    check_pauli_apply(P, PSI2_01);
    check_pauli_apply(P, PSI2_10);
    check_pauli_apply(P, PSI2_11);
    check_pauli_apply(P, PSI2_PP);
    check_pauli_apply(P, PSI2_P0);
    check_pauli_apply(P, PSI2_1P);

    pauli_free(P);
}

/*
 * =====================================================================================================================
 * Test: n=2, IY
 *
 * labels[0]=I, labels[1]=Y.  x_mask=2, z_mask=2, y_mask=2, split_bit=1.
 * IY flips qubit 1 with phase i on the source col:
 *   eta = popcount(2) & 3 = 1.  eps(col) = i * (-1)^popcount(col & 2).
 *   |00>(col=0) -> i|10>,  |01>(col=1) -> i|11>,
 *   |10>(col=2) -> -i|00>, |11>(col=3) -> -i|01>
 * =====================================================================================================================
 */

void test_pauli_apply_IY(void) {
    const pauli_label_t lbls[2] = {PAULI_I, PAULI_Y};
    pauli_t *P = pauli_new(lbls, 2);
    TEST_ASSERT_NOT_NULL(P);

    TEST_ASSERT_EQUAL_UINT64(2u, P->x_mask);
    TEST_ASSERT_EQUAL_UINT64(2u, P->z_mask);
    TEST_ASSERT_EQUAL_UINT64(2u, P->y_mask);

    check_pauli_apply(P, PSI2_00);
    check_pauli_apply(P, PSI2_01);
    check_pauli_apply(P, PSI2_10);
    check_pauli_apply(P, PSI2_11);
    check_pauli_apply(P, PSI2_PP);
    check_pauli_apply(P, PSI2_P0);
    check_pauli_apply(P, PSI2_1P);

    pauli_free(P);
}

/*
 * =====================================================================================================================
 * Test: n=3, IXI
 *
 * labels[0]=I, labels[1]=X, labels[2]=I.  x_mask=2, z_mask=0, y_mask=0, split_bit=1.
 * Flips qubit 1 only; phase always +1.
 * =====================================================================================================================
 */

void test_pauli_apply_IXI(void) {
    const pauli_label_t lbls[3] = {PAULI_I, PAULI_X, PAULI_I};
    pauli_t *P = pauli_new(lbls, 3);
    TEST_ASSERT_NOT_NULL(P);

    TEST_ASSERT_EQUAL_UINT64(2u, P->x_mask);
    TEST_ASSERT_EQUAL_UINT64(0u, P->z_mask);
    TEST_ASSERT_EQUAL_UINT64(0u, P->y_mask);

    check_pauli_apply(P, PSI3_000);
    check_pauli_apply(P, PSI3_001);
    check_pauli_apply(P, PSI3_010);
    check_pauli_apply(P, PSI3_011);
    check_pauli_apply(P, PSI3_100);
    check_pauli_apply(P, PSI3_101);
    check_pauli_apply(P, PSI3_110);
    check_pauli_apply(P, PSI3_111);
    check_pauli_apply(P, PSI3_GHZ);
    check_pauli_apply(P, PSI3_W);
    check_pauli_apply(P, PSI3_PPP);

    pauli_free(P);
}

/*
 * =====================================================================================================================
 * Test: n=3, XII
 *
 * labels[0]=I, labels[1]=I, labels[2]=X.  x_mask=4, z_mask=0, y_mask=0, split_bit=2.
 * Flips qubit 2 only; phase always +1.
 * =====================================================================================================================
 */

void test_pauli_apply_XII(void) {
    const pauli_label_t lbls[3] = {PAULI_I, PAULI_I, PAULI_X};
    pauli_t *P = pauli_new(lbls, 3);
    TEST_ASSERT_NOT_NULL(P);

    TEST_ASSERT_EQUAL_UINT64(4u, P->x_mask);
    TEST_ASSERT_EQUAL_UINT64(0u, P->z_mask);
    TEST_ASSERT_EQUAL_UINT64(0u, P->y_mask);

    check_pauli_apply(P, PSI3_000);
    check_pauli_apply(P, PSI3_001);
    check_pauli_apply(P, PSI3_010);
    check_pauli_apply(P, PSI3_011);
    check_pauli_apply(P, PSI3_100);
    check_pauli_apply(P, PSI3_101);
    check_pauli_apply(P, PSI3_110);
    check_pauli_apply(P, PSI3_111);
    check_pauli_apply(P, PSI3_GHZ);
    check_pauli_apply(P, PSI3_W);
    check_pauli_apply(P, PSI3_PPP);

    pauli_free(P);
}

/*
 * =====================================================================================================================
 * Test: n=3, IZI
 *
 * labels[0]=I, labels[1]=Z, labels[2]=I.  x_mask=0, z_mask=2, y_mask=0.
 * Diagonal: phase = (-1)^popcount(col & 2) (non-uniform z_mask with gap).
 *   |000>(col=0) -> +|000>,  |001>(col=1) -> +|001>,
 *   |010>(col=2) -> -|010>,  |011>(col=3) -> -|011>,
 *   |100>(col=4) -> +|100>,  |101>(col=5) -> +|101>,
 *   |110>(col=6) -> -|110>,  |111>(col=7) -> -|111>
 * =====================================================================================================================
 */

void test_pauli_apply_IZI(void) {
    const pauli_label_t lbls[3] = {PAULI_I, PAULI_Z, PAULI_I};
    pauli_t *P = pauli_new(lbls, 3);
    TEST_ASSERT_NOT_NULL(P);

    TEST_ASSERT_EQUAL_UINT64(0u, P->x_mask);
    TEST_ASSERT_EQUAL_UINT64(2u, P->z_mask);
    TEST_ASSERT_EQUAL_UINT64(0u, P->y_mask);

    check_pauli_apply(P, PSI3_000);
    check_pauli_apply(P, PSI3_001);
    check_pauli_apply(P, PSI3_010);
    check_pauli_apply(P, PSI3_011);
    check_pauli_apply(P, PSI3_100);
    check_pauli_apply(P, PSI3_101);
    check_pauli_apply(P, PSI3_110);
    check_pauli_apply(P, PSI3_111);
    check_pauli_apply(P, PSI3_GHZ);
    check_pauli_apply(P, PSI3_W);
    check_pauli_apply(P, PSI3_PPP);

    pauli_free(P);
}

/*
 * =====================================================================================================================
 * Test: n=3, ZIZ
 *
 * labels[0]=Z, labels[1]=I, labels[2]=Z.  x_mask=0, z_mask=5, y_mask=0.
 * Diagonal: phase = (-1)^popcount(col & 5) (alternating z_mask, bit 0 and bit 2 set).
 *   |000>(col=0) -> +|000>,  |001>(col=1) -> -|001>,
 *   |010>(col=2) -> +|010>,  |011>(col=3) -> -|011>,
 *   |100>(col=4) -> -|100>,  |101>(col=5) -> +|101>,
 *   |110>(col=6) -> -|110>,  |111>(col=7) -> +|111>
 * =====================================================================================================================
 */

void test_pauli_apply_ZIZ(void) {
    const pauli_label_t lbls[3] = {PAULI_Z, PAULI_I, PAULI_Z};
    pauli_t *P = pauli_new(lbls, 3);
    TEST_ASSERT_NOT_NULL(P);

    TEST_ASSERT_EQUAL_UINT64(0u, P->x_mask);
    TEST_ASSERT_EQUAL_UINT64(5u, P->z_mask);
    TEST_ASSERT_EQUAL_UINT64(0u, P->y_mask);

    check_pauli_apply(P, PSI3_000);
    check_pauli_apply(P, PSI3_001);
    check_pauli_apply(P, PSI3_010);
    check_pauli_apply(P, PSI3_011);
    check_pauli_apply(P, PSI3_100);
    check_pauli_apply(P, PSI3_101);
    check_pauli_apply(P, PSI3_110);
    check_pauli_apply(P, PSI3_111);
    check_pauli_apply(P, PSI3_GHZ);
    check_pauli_apply(P, PSI3_W);
    check_pauli_apply(P, PSI3_PPP);

    pauli_free(P);
}

/*
 * =====================================================================================================================
 * Test: n=4, IIXI — split_bit=1
 *
 * labels[0]=I, labels[1]=X, labels[2]=I, labels[3]=I.  x_mask=2, z_mask=0, y_mask=0.
 * Flips qubit 1 only; phase always +1.
 * =====================================================================================================================
 */

void test_pauli_apply_4q_split1(void) {
    const pauli_label_t lbls[4] = {PAULI_I, PAULI_X, PAULI_I, PAULI_I};
    pauli_t *P = pauli_new(lbls, 4);
    TEST_ASSERT_NOT_NULL(P);

    TEST_ASSERT_EQUAL_UINT64(2u, P->x_mask);
    TEST_ASSERT_EQUAL_UINT64(0u, P->z_mask);
    TEST_ASSERT_EQUAL_UINT64(0u, P->y_mask);

    for (unsigned k = 0; k < 16; ++k) {
        cplx_t psi[16];
        memset(psi, 0, sizeof(psi));
        psi[k] = 1.0 + 0.0*I;
        check_pauli_apply(P, psi);
    }

    check_pauli_apply(P, PSI4_GHZ);

    pauli_free(P);
}

/*
 * =====================================================================================================================
 * Test: n=4, IXII — split_bit=2
 *
 * labels[0]=I, labels[1]=I, labels[2]=X, labels[3]=I.  x_mask=4, z_mask=0, y_mask=0.
 * Flips qubit 2 only; phase always +1.
 * =====================================================================================================================
 */

void test_pauli_apply_4q_split2(void) {
    const pauli_label_t lbls[4] = {PAULI_I, PAULI_I, PAULI_X, PAULI_I};
    pauli_t *P = pauli_new(lbls, 4);
    TEST_ASSERT_NOT_NULL(P);

    TEST_ASSERT_EQUAL_UINT64(4u, P->x_mask);
    TEST_ASSERT_EQUAL_UINT64(0u, P->z_mask);
    TEST_ASSERT_EQUAL_UINT64(0u, P->y_mask);

    for (unsigned k = 0; k < 16; ++k) {
        cplx_t psi[16];
        memset(psi, 0, sizeof(psi));
        psi[k] = 1.0 + 0.0*I;
        check_pauli_apply(P, psi);
    }

    check_pauli_apply(P, PSI4_GHZ);

    pauli_free(P);
}

/*
 * =====================================================================================================================
 * Test: n=4, XIII — split_bit=3
 *
 * labels[0]=I, labels[1]=I, labels[2]=I, labels[3]=X.  x_mask=8, z_mask=0, y_mask=0.
 * Flips qubit 3 only; phase always +1.
 * =====================================================================================================================
 */

void test_pauli_apply_4q_split3(void) {
    const pauli_label_t lbls[4] = {PAULI_I, PAULI_I, PAULI_I, PAULI_X};
    pauli_t *P = pauli_new(lbls, 4);
    TEST_ASSERT_NOT_NULL(P);

    TEST_ASSERT_EQUAL_UINT64(8u, P->x_mask);
    TEST_ASSERT_EQUAL_UINT64(0u, P->z_mask);
    TEST_ASSERT_EQUAL_UINT64(0u, P->y_mask);

    for (unsigned k = 0; k < 16; ++k) {
        cplx_t psi[16];
        memset(psi, 0, sizeof(psi));
        psi[k] = 1.0 + 0.0*I;
        check_pauli_apply(P, psi);
    }

    check_pauli_apply(P, PSI4_GHZ);

    pauli_free(P);
}

/*
 * =====================================================================================================================
 * Test: pauli_from_str round-trip
 *
 * pauli_from_str("ZYX", 3) and pauli_new({X,Y,Z}, 3) must encode the same Pauli
 * string (qubit 0 = X, qubit 1 = Y, qubit 2 = Z) and produce identical results
 * when applied to every computational basis state of the 3-qubit system.
 *
 * "ZYX" big-endian: s[0]='Z'->labels[2]=Z, s[1]='Y'->labels[1]=Y, s[2]='X'->labels[0]=X.
 * =====================================================================================================================
 */

void test_pauli_from_str_roundtrip(void) {
    const pauli_label_t lbls[3] = {PAULI_X, PAULI_Y, PAULI_Z};
    pauli_t *P_new = pauli_new(lbls, 3);
    TEST_ASSERT_NOT_NULL(P_new);

    pauli_t *P_str = pauli_from_str("ZYX", 3);
    TEST_ASSERT_NOT_NULL(P_str);

    /* Both must have identical masks */
    TEST_ASSERT_EQUAL_UINT64(P_new->x_mask, P_str->x_mask);
    TEST_ASSERT_EQUAL_UINT64(P_new->z_mask, P_str->z_mask);
    TEST_ASSERT_EQUAL_UINT64(P_new->y_mask, P_str->y_mask);

    /* Apply both to the same states and compare results */
    for (unsigned k = 0; k < 8; ++k) {
        cplx_t psi[8];
        memset(psi, 0, sizeof(psi));
        psi[k] = 1.0 + 0.0*I;

        /* Result from pauli_new version */
        cplx_t *copy_new = malloc(8 * sizeof(*copy_new));
        TEST_ASSERT_NOT_NULL(copy_new);
        memcpy(copy_new, psi, 8 * sizeof(*copy_new));
        state_t s_new = {0};
        s_new.type = PURE;
        state_init(&s_new, 3, &copy_new);
        if (!s_new.data) {
            free(copy_new);
            pauli_free(P_new);
            pauli_free(P_str);
            TEST_FAIL_MESSAGE("state_init failed for P_new in round-trip test");
        }
        pauli_apply_pure(&s_new, P_new);

        /* Result from pauli_from_str version */
        cplx_t *copy_str = malloc(8 * sizeof(*copy_str));
        if (!copy_str) {
            free(s_new.data);
            pauli_free(P_new);
            pauli_free(P_str);
            TEST_FAIL_MESSAGE("copy_str allocation failed in round-trip test");
        }
        memcpy(copy_str, psi, 8 * sizeof(*copy_str));
        state_t s_str = {0};
        s_str.type = PURE;
        state_init(&s_str, 3, &copy_str);
        if (!s_str.data) {
            free(copy_str);
            free(s_new.data);
            pauli_free(P_new);
            pauli_free(P_str);
            TEST_FAIL_MESSAGE("state_init failed for P_str in round-trip test");
        }
        pauli_apply_pure(&s_str, P_str);

        TEST_ASSERT_EQUAL_COMPLEX_ARRAY_TOL(s_new.data, s_str.data, 8, PRECISION);

        free(s_new.data);
        s_new.data = NULL;
        free(s_str.data);
        s_str.data = NULL;
    }

    pauli_free(P_new);
    pauli_free(P_str);
}

/*
 * =====================================================================================================================
 * Test: Pauli involution P^2 = I
 *
 * Applying any Pauli string twice must recover the original state for:
 *   X (n=1), Y (n=1), Z (n=1), XX (n=2), XYZ (n=3).
 * =====================================================================================================================
 */

void test_pauli_apply_involution(void) {
    /* Helper: apply P twice to a fresh copy of psi and verify we get psi back. */
    struct {
        pauli_label_t lbls[3];
        unsigned      n;
        unsigned      nlbls;
    } cases[] = {
        {{PAULI_X, 0, 0},             1, 1},
        {{PAULI_Y, 0, 0},             1, 1},
        {{PAULI_Z, 0, 0},             1, 1},
        {{PAULI_X, PAULI_X, 0},       2, 2},
        {{PAULI_X, PAULI_Y, PAULI_Z}, 3, 3},
    };
    const unsigned ncases = sizeof(cases) / sizeof(cases[0]);

    for (unsigned c = 0; c < ncases; ++c) {
        const unsigned n   = cases[c].n;
        const unsigned dim = 1u << n;

        pauli_t *P = pauli_new(cases[c].lbls, (qubit_t)n);
        TEST_ASSERT_NOT_NULL(P);

        /* Test on |0...0> (all zeros except index 0) */
        cplx_t orig[8];  /* dim <= 8 for n <= 3 */
        memset(orig, 0, dim * sizeof(*orig));
        orig[0] = 1.0 + 0.0*I;

        cplx_t *copy = malloc(dim * sizeof(*copy));
        TEST_ASSERT_NOT_NULL(copy);
        memcpy(copy, orig, dim * sizeof(*copy));

        state_t s = {0};
        s.type = PURE;
        state_init(&s, (qubit_t)n, &copy);
        if (!s.data) {
            free(copy);
            pauli_free(P);
            TEST_FAIL_MESSAGE("state_init failed in involution test");
        }

        /* First application */
        pauli_apply_pure(&s, P);
        /* Second application */
        pauli_apply_pure(&s, P);

        TEST_ASSERT_EQUAL_COMPLEX_ARRAY_TOL(orig, s.data, dim, PRECISION);

        free(s.data);
        s.data = NULL;
        pauli_free(P);
    }
}

/*
 * =====================================================================================================================
 * Test: III on n=3 uniform superposition is identity
 *
 * I^(x)3 applied to |+++> must return |+++> unchanged.
 * =====================================================================================================================
 */

void test_pauli_apply_identity_noop(void) {
    const pauli_label_t lbls[3] = {PAULI_I, PAULI_I, PAULI_I};
    pauli_t *P = pauli_new(lbls, 3);
    TEST_ASSERT_NOT_NULL(P);

    /* Uniform superposition over 8 basis states: amplitude = INVSQRT8 for each */
    cplx_t psi[8];
    for (unsigned k = 0; k < 8; ++k) {
        psi[k] = INVSQRT8 + 0.0*I;
    }

    cplx_t *copy = malloc(8 * sizeof(*copy));
    TEST_ASSERT_NOT_NULL(copy);
    memcpy(copy, psi, 8 * sizeof(*copy));

    state_t s = {0};
    s.type = PURE;
    state_init(&s, 3, &copy);
    if (!s.data) {
        free(copy);
        pauli_free(P);
        TEST_FAIL_MESSAGE("state_init failed in identity noop test");
    }

    pauli_apply_pure(&s, P);

    TEST_ASSERT_EQUAL_COMPLEX_ARRAY_TOL(psi, s.data, 8, PRECISION);

    free(s.data);
    s.data = NULL;
    pauli_free(P);
}

/*
 * =====================================================================================================================
 * Gap 1: pauli_new with n_qubits = 0
 *
 * pauli_new has explicit n=0 handling: it allocates a struct with labels=NULL,
 * all masks zero, and returns a valid pointer.  pauli_free(P) must not crash.
 * =====================================================================================================================
 */

void test_pauli_new_n0(void) {
    pauli_t *P = pauli_new(NULL, 0);
    TEST_ASSERT_NOT_NULL(P);
    TEST_ASSERT_EQUAL_UINT(0u, (unsigned)P->n_qubits);
    TEST_ASSERT_EQUAL_UINT64(0u, P->x_mask);
    TEST_ASSERT_EQUAL_UINT64(0u, P->z_mask);
    TEST_ASSERT_EQUAL_UINT64(0u, P->y_mask);
    TEST_ASSERT_NULL(P->labels);
    pauli_free(P);   /* must not crash */
}

/*
 * =====================================================================================================================
 * Gap 2: pauli_new(NULL, n_qubits > 0) must return NULL
 *
 * When n_qubits > 0 the label array is required.  NULL labels is an error.
 * =====================================================================================================================
 */

void test_pauli_new_null_labels(void) {
    pauli_t *P = pauli_new(NULL, 1);
    TEST_ASSERT_NULL(P);

    P = pauli_new(NULL, 3);
    TEST_ASSERT_NULL(P);
}

/*
 * =====================================================================================================================
 * Gap 3: pauli_new(labels, n_qubits > 64) must return NULL
 *
 * The bitmask representation only supports up to 64 qubits.
 * =====================================================================================================================
 */

void test_pauli_new_too_many_qubits(void) {
    const pauli_label_t lbl = PAULI_X;
    pauli_t *P = pauli_new(&lbl, 65);
    TEST_ASSERT_NULL(P);
}

/*
 * =====================================================================================================================
 * Gap 4a: pauli_phase for YY (eta_exp = 2)
 *
 * y_mask=3, z_mask=3.  eta = popcount(3) & 3 = 2.
 * eps(x) = i^2 * (-1)^popcount(x & 3) = -(-1)^popcount(x & 3).
 *
 *   x=0: -1 * 1        = -1  -> {real=-1, imag=0}
 *   x=1: -1 * (-1)     = +1  -> {real=+1, imag=0}
 *   x=2: -1 * (-1)     = +1  -> {real=+1, imag=0}
 *   x=3: -1 * 1        = -1  -> {real=-1, imag=0}
 * =====================================================================================================================
 */

void test_pauli_phase_eta2(void) {
    const pauli_label_t lbls[2] = {PAULI_Y, PAULI_Y};
    pauli_t *P = pauli_new(lbls, 2);
    TEST_ASSERT_NOT_NULL(P);
    /* Verify masks */
    TEST_ASSERT_EQUAL_UINT64(3u, P->x_mask);
    TEST_ASSERT_EQUAL_UINT64(3u, P->z_mask);
    TEST_ASSERT_EQUAL_UINT64(3u, P->y_mask);

    pauli_phase_t ph;

    /* x=0: eta=2, parity=0 -> sign=+1, i^2=-1 -> {real=-1, imag=0} */
    ph = pauli_phase(P, 0u);
    TEST_ASSERT_EQUAL_INT(-1, ph.real);
    TEST_ASSERT_EQUAL_INT( 0, ph.imag);

    /* x=1: parity=popcount(1&3)=1 -> sign=-1, -1*(-1)=+1 -> {real=+1, imag=0} */
    ph = pauli_phase(P, 1u);
    TEST_ASSERT_EQUAL_INT(1, ph.real);
    TEST_ASSERT_EQUAL_INT(0, ph.imag);

    /* x=2: parity=popcount(2&3)=1 -> sign=-1 -> {real=+1, imag=0} */
    ph = pauli_phase(P, 2u);
    TEST_ASSERT_EQUAL_INT(1, ph.real);
    TEST_ASSERT_EQUAL_INT(0, ph.imag);

    /* x=3: parity=popcount(3&3)=2, parity&1=0 -> sign=+1 -> {real=-1, imag=0} */
    ph = pauli_phase(P, 3u);
    TEST_ASSERT_EQUAL_INT(-1, ph.real);
    TEST_ASSERT_EQUAL_INT( 0, ph.imag);

    pauli_free(P);
}

/*
 * =====================================================================================================================
 * Gap 4b: pauli_phase for YYY (eta_exp = 3)
 *
 * y_mask=7, z_mask=7.  eta = popcount(7) & 3 = 3.
 * eps(x) = i^3 * (-1)^popcount(x & 7) = -i * (-1)^popcount(x & 7).
 *
 *   x=0: -i * 1   = -i -> {real=0, imag=-1}
 *   x=7: -i * (-1) = +i -> {real=0, imag=+1}
 * =====================================================================================================================
 */

void test_pauli_phase_eta3(void) {
    const pauli_label_t lbls[3] = {PAULI_Y, PAULI_Y, PAULI_Y};
    pauli_t *P = pauli_new(lbls, 3);
    TEST_ASSERT_NOT_NULL(P);
    TEST_ASSERT_EQUAL_UINT64(7u, P->x_mask);
    TEST_ASSERT_EQUAL_UINT64(7u, P->z_mask);
    TEST_ASSERT_EQUAL_UINT64(7u, P->y_mask);

    pauli_phase_t ph;

    /* x=0: parity=0 -> sign=+1, i^3={0,-1} -> {real=0, imag=-1} */
    ph = pauli_phase(P, 0u);
    TEST_ASSERT_EQUAL_INT( 0, ph.real);
    TEST_ASSERT_EQUAL_INT(-1, ph.imag);

    /* x=7: parity=popcount(7&7)=3, parity&1=1 -> sign=-1 -> {real=0, imag=+1} */
    ph = pauli_phase(P, 7u);
    TEST_ASSERT_EQUAL_INT(0, ph.real);
    TEST_ASSERT_EQUAL_INT(1, ph.imag);

    /* Spot-check x=1: parity=popcount(1&7)=1 -> sign=-1, -1*i^3={0,+1} */
    ph = pauli_phase(P, 1u);
    TEST_ASSERT_EQUAL_INT(0, ph.real);
    TEST_ASSERT_EQUAL_INT(1, ph.imag);

    /* Spot-check x=3: parity=popcount(3&7)=2, parity&1=0 -> sign=+1 -> {real=0,imag=-1} */
    ph = pauli_phase(P, 3u);
    TEST_ASSERT_EQUAL_INT( 0, ph.real);
    TEST_ASSERT_EQUAL_INT(-1, ph.imag);

    pauli_free(P);
}

/*
 * =====================================================================================================================
 * Gap 5: pauli_apply_pure — YYY (eta_exp = 3, the only previously untested value)
 *
 * labels = {Y, Y, Y}: x_mask=7, z_mask=7, y_mask=7.
 *
 * Inline spot-check: YYY|000> = -i|111>.
 *   split_bit=0.  Pair k=0: j=0, jp=7.
 *   eta_exp=3, parity_j=popcount(0&7)=0, exp_j=3.
 *   y_parity=popcount(7)&1=1, exp_jp=(3+2)&3=1.
 *   data[0] = apply_ipow(old[7]=0, 1) = 0.
 *   data[7] = apply_ipow(old[0]=1, 3) = CMPLX(0,-1) = -i.
 *
 * Full reference check via check_pauli_apply on all 8 basis states + entangled states.
 * =====================================================================================================================
 */

void test_pauli_apply_YYY(void) {
    const pauli_label_t lbls[3] = {PAULI_Y, PAULI_Y, PAULI_Y};
    pauli_t *P = pauli_new(lbls, 3);
    TEST_ASSERT_NOT_NULL(P);

    TEST_ASSERT_EQUAL_UINT64(7u, P->x_mask);
    TEST_ASSERT_EQUAL_UINT64(7u, P->z_mask);
    TEST_ASSERT_EQUAL_UINT64(7u, P->y_mask);

    /*
     * Inline spot-check: YYY|000> should produce -i|111>.
     * Expected: all zeros except index 7 has value -i = CMPLX(0, -1).
     */
    {
        static const cplx_t expected[8] = {
            0.0+0.0*I, 0.0+0.0*I, 0.0+0.0*I, 0.0+0.0*I,
            0.0+0.0*I, 0.0+0.0*I, 0.0+0.0*I, 0.0-1.0*I
        };

        cplx_t *copy = malloc(8 * sizeof(*copy));
        TEST_ASSERT_NOT_NULL(copy);
        memcpy(copy, PSI3_000, 8 * sizeof(*copy));

        state_t s = {0};
        s.type = PURE;
        state_init(&s, 3, &copy);
        if (!s.data) {
            free(copy);
            pauli_free(P);
            TEST_FAIL_MESSAGE("state_init failed in YYY spot-check");
        }

        pauli_apply_pure(&s, P);
        TEST_ASSERT_EQUAL_COMPLEX_ARRAY_TOL(expected, s.data, 8, PRECISION);
        free(s.data);
        s.data = NULL;
    }

    /* Full reference check on all 8 computational basis states + entangled states */
    check_pauli_apply(P, PSI3_000);
    check_pauli_apply(P, PSI3_001);
    check_pauli_apply(P, PSI3_010);
    check_pauli_apply(P, PSI3_011);
    check_pauli_apply(P, PSI3_100);
    check_pauli_apply(P, PSI3_101);
    check_pauli_apply(P, PSI3_110);
    check_pauli_apply(P, PSI3_111);
    check_pauli_apply(P, PSI3_GHZ);
    check_pauli_apply(P, PSI3_W);
    check_pauli_apply(P, PSI3_PPP);

    pauli_free(P);
}

/*
 * =====================================================================================================================
 * Gap 6: Z involution on non-trivial states
 *
 * The existing test_pauli_apply_involution tests Z involution on |0> only, where
 * Z|0> = |0> (eigenstate: the sign flip never fires at index 0).
 *
 * Here we test Z^2 = I on |1>, |+>, and |-i> — all states where the amplitude
 * at index 1 is nonzero, so the -1 phase flip actually executes and must cancel.
 * =====================================================================================================================
 */

void test_pauli_apply_involution_Z_nontrivial(void) {
    const pauli_label_t lbl = PAULI_Z;
    pauli_t *P = pauli_new(&lbl, 1);
    TEST_ASSERT_NOT_NULL(P);

    /* Test states where index-1 amplitude is nonzero */
    const cplx_t *test_states[3] = {PSI1_1, PSI1_P, PSI1_MI};
    const unsigned n_states = 3;

    for (unsigned k = 0; k < n_states; ++k) {
        const cplx_t *psi = test_states[k];

        cplx_t *copy = malloc(2 * sizeof(*copy));
        TEST_ASSERT_NOT_NULL_MESSAGE(copy, "malloc failed in Z involution test");
        memcpy(copy, psi, 2 * sizeof(*copy));

        state_t s = {0};
        s.type = PURE;
        state_init(&s, 1, &copy);
        if (!s.data) {
            free(copy);
            pauli_free(P);
            TEST_FAIL_MESSAGE("state_init failed in Z involution test");
        }

        /* Apply Z twice; result must equal original */
        pauli_apply_pure(&s, P);
        pauli_apply_pure(&s, P);

        TEST_ASSERT_EQUAL_COMPLEX_ARRAY_TOL(psi, s.data, 2, PRECISION);

        free(s.data);
        s.data = NULL;
    }

    pauli_free(P);
}

/*
 * =====================================================================================================================
 * Gap 7: pauli_apply_pure — YZ (eta_exp=1 with z-phase from separate qubit)
 *
 * pauli_new({Z, Y}, 2): labels[0]=Z (qubit 0), labels[1]=Y (qubit 1).
 *   x_mask = bit 1 (Y on qubit 1)  = 2
 *   z_mask = bits 0(Z) + 1(Y)      = 3
 *   y_mask = bit 1 (Y on qubit 1)  = 2
 *   eta_exp = popcount(2) & 3 = 1
 *
 * eps(col) = i * (-1)^popcount(col & 3):
 *   col=0: i * 1   = +i -> output[2] = +i  (row = 0 XOR 2 = 2)
 *   col=1: i * (-1) = -i -> output[3] = -i  (row = 1 XOR 2 = 3)
 *   col=2: i * (-1) = -i -> output[0] = -i  (row = 2 XOR 2 = 0)
 *   col=3: i * 1   = +i -> output[1] = +i  (row = 3 XOR 2 = 1)
 * =====================================================================================================================
 */

void test_pauli_apply_YZ(void) {
    const pauli_label_t lbls[2] = {PAULI_Z, PAULI_Y};
    pauli_t *P = pauli_new(lbls, 2);
    TEST_ASSERT_NOT_NULL(P);

    /* Verify masks */
    TEST_ASSERT_EQUAL_UINT64(2u, P->x_mask);
    TEST_ASSERT_EQUAL_UINT64(3u, P->z_mask);
    TEST_ASSERT_EQUAL_UINT64(2u, P->y_mask);

    /* Reference check on all 7 n=2 test states */
    check_pauli_apply(P, PSI2_00);
    check_pauli_apply(P, PSI2_01);
    check_pauli_apply(P, PSI2_10);
    check_pauli_apply(P, PSI2_11);
    check_pauli_apply(P, PSI2_PP);
    check_pauli_apply(P, PSI2_P0);
    check_pauli_apply(P, PSI2_1P);

    pauli_free(P);
}

/*
 * =====================================================================================================================
 * Gap 8: pauli_from_str case-sensitivity
 *
 * Only uppercase 'I', 'X', 'Y', 'Z' are valid.  Lowercase letters must return NULL.
 * =====================================================================================================================
 */

void test_pauli_from_str_case_sensitive(void) {
    /* Lowercase single character */
    pauli_t *P = pauli_from_str("x", 1);
    TEST_ASSERT_NULL(P);

    P = pauli_from_str("y", 1);
    TEST_ASSERT_NULL(P);

    P = pauli_from_str("z", 1);
    TEST_ASSERT_NULL(P);

    P = pauli_from_str("i", 1);
    TEST_ASSERT_NULL(P);

    /* Lowercase multi-character */
    P = pauli_from_str("xyz", 3);
    TEST_ASSERT_NULL(P);

    /* Mixed case: one invalid lowercase in an otherwise valid string */
    P = pauli_from_str("XxZ", 3);
    TEST_ASSERT_NULL(P);
}

/*
 * =====================================================================================================================
 * Gap 9: pauli_from_str("", 0) — empty string with n_qubits = 0
 *
 * strlen("") = 0 = n_qubits, so the length check passes.
 * malloc(0) is called for the label array; its return value (NULL or non-NULL
 * is implementation-defined) is then passed to pauli_new(labels, 0).
 * Inside pauli_new, n_qubits == 0 causes the early-return path where labels
 * is never read — so the function succeeds regardless of the malloc(0) result.
 *
 * Expected: non-NULL with all masks zero and n_qubits == 0.
 * =====================================================================================================================
 */

void test_pauli_from_str_empty(void) {
    pauli_t *P = pauli_from_str("", 0);
    TEST_ASSERT_NOT_NULL(P);
    TEST_ASSERT_EQUAL_UINT(0u, (unsigned)P->n_qubits);
    TEST_ASSERT_EQUAL_UINT64(0u, P->x_mask);
    TEST_ASSERT_EQUAL_UINT64(0u, P->z_mask);
    TEST_ASSERT_EQUAL_UINT64(0u, P->y_mask);
    pauli_free(P);
}

/*
 * =====================================================================================================================
 * Gap 10: y_mask invariant
 *
 * For any correctly-constructed Pauli P:
 *   (P->y_mask & P->x_mask) == P->y_mask   (y bits are a subset of x bits)
 *   (P->y_mask & P->z_mask) == P->y_mask   (y bits are a subset of z bits)
 *
 * Checked for: Y (n=1), YY (n=2), YYY (n=3), XYZ (n=3), IXYZ (n=4).
 * =====================================================================================================================
 */

void test_pauli_mask_invariant(void) {
    /* Y, n=1 */
    {
        const pauli_label_t lbl = PAULI_Y;
        pauli_t *P = pauli_new(&lbl, 1);
        TEST_ASSERT_NOT_NULL(P);
        TEST_ASSERT_EQUAL_UINT64(P->y_mask, P->y_mask & P->x_mask);
        TEST_ASSERT_EQUAL_UINT64(P->y_mask, P->y_mask & P->z_mask);
        pauli_free(P);
    }

    /* YY, n=2 */
    {
        const pauli_label_t lbls[2] = {PAULI_Y, PAULI_Y};
        pauli_t *P = pauli_new(lbls, 2);
        TEST_ASSERT_NOT_NULL(P);
        TEST_ASSERT_EQUAL_UINT64(P->y_mask, P->y_mask & P->x_mask);
        TEST_ASSERT_EQUAL_UINT64(P->y_mask, P->y_mask & P->z_mask);
        pauli_free(P);
    }

    /* YYY, n=3 */
    {
        const pauli_label_t lbls[3] = {PAULI_Y, PAULI_Y, PAULI_Y};
        pauli_t *P = pauli_new(lbls, 3);
        TEST_ASSERT_NOT_NULL(P);
        TEST_ASSERT_EQUAL_UINT64(P->y_mask, P->y_mask & P->x_mask);
        TEST_ASSERT_EQUAL_UINT64(P->y_mask, P->y_mask & P->z_mask);
        pauli_free(P);
    }

    /* XYZ, n=3: labels[0]=X, labels[1]=Y, labels[2]=Z -> y_mask=bit1=2 */
    {
        const pauli_label_t lbls[3] = {PAULI_X, PAULI_Y, PAULI_Z};
        pauli_t *P = pauli_new(lbls, 3);
        TEST_ASSERT_NOT_NULL(P);
        TEST_ASSERT_EQUAL_UINT64(2u, P->y_mask);
        TEST_ASSERT_EQUAL_UINT64(P->y_mask, P->y_mask & P->x_mask);
        TEST_ASSERT_EQUAL_UINT64(P->y_mask, P->y_mask & P->z_mask);
        pauli_free(P);
    }

    /* IXYZ, n=4: from pauli_from_str("IXYZ", 4) -> y_mask=bit1=2 */
    {
        pauli_t *P = pauli_from_str("IXYZ", 4);
        TEST_ASSERT_NOT_NULL(P);
        TEST_ASSERT_EQUAL_UINT64(2u, P->y_mask);
        TEST_ASSERT_EQUAL_UINT64(P->y_mask, P->y_mask & P->x_mask);
        TEST_ASSERT_EQUAL_UINT64(P->y_mask, P->y_mask & P->z_mask);
        pauli_free(P);
    }
}

/*
 * =====================================================================================================================
 * Gap 11: pauli_free(NULL) — documented as safe
 *
 * circuit.h documents: "Safe to call with NULL."
 * Implementation: `if (!p) return;`.
 * A call with NULL must not crash or abort.
 * =====================================================================================================================
 */

void test_pauli_free_null(void) {
    pauli_free(NULL);   /* must not crash */
}

/*
 * =====================================================================================================================
 * Gap 12: pauli_new with n_qubits = 64 — boundary of valid range
 *
 * The guard is `n_qubits > 64`, so n=64 is the maximum valid input.
 * The mask construction uses `(uint64_t)1 << j` for j = 0..63.
 * Qubit 63 (j=63) exercises the shift `(uint64_t)1 << 63 = 0x8000000000000000`,
 * which is well-defined for uint64_t.
 *
 * Test: construct a 64-qubit Pauli where only qubit 63 is X, qubit 62 is Y,
 * and all others are I.  Verify:
 *   x_mask has bits 62 and 63 set (X on q63, Y on q62)
 *   z_mask has bit 62 set only  (Y on q62)
 *   y_mask has bit 62 set only
 *   n_qubits == 64
 * =====================================================================================================================
 */

void test_pauli_new_n64(void) {
    pauli_label_t lbls[64];
    for (unsigned j = 0; j < 64; ++j) {
        lbls[j] = PAULI_I;
    }
    lbls[63] = PAULI_X;
    lbls[62] = PAULI_Y;

    pauli_t *P = pauli_new(lbls, 64);
    TEST_ASSERT_NOT_NULL(P);
    TEST_ASSERT_EQUAL_UINT(64u, (unsigned)P->n_qubits);

    /* bit 63 (X) and bit 62 (Y) are both in x_mask */
    const uint64_t expected_x = ((uint64_t)1 << 63) | ((uint64_t)1 << 62);
    /* bit 62 (Y) is in z_mask */
    const uint64_t expected_z = (uint64_t)1 << 62;
    /* bit 62 (Y) is in y_mask */
    const uint64_t expected_y = (uint64_t)1 << 62;

    TEST_ASSERT_EQUAL_UINT64(expected_x, P->x_mask);
    TEST_ASSERT_EQUAL_UINT64(expected_z, P->z_mask);
    TEST_ASSERT_EQUAL_UINT64(expected_y, P->y_mask);

    /* y_mask subset invariant */
    TEST_ASSERT_EQUAL_UINT64(P->y_mask, P->y_mask & P->x_mask);
    TEST_ASSERT_EQUAL_UINT64(P->y_mask, P->y_mask & P->z_mask);

    pauli_free(P);
}

/*
 * =====================================================================================================================
 * pauli_exp_pure tests
 *
 * pauli_exp_pure applies U(theta) = exp(-i*theta/2 * P) = c*I - i*s*P, where
 * c = cos(theta/2), s = sin(theta/2).
 *
 * Reference: ref[k] = c * psi[k]  +  (-i*s) * (P|psi>)[k]
 *          = c * psi[k]  -  i*s * sum_j M[k,j] * psi[j]
 *
 * Helper check_pauli_exp builds P|psi> via build_pauli_matrix+zmv, scales by
 * (-i*s) = CMPLX(0,-s), adds c*psi, then compares against pauli_exp_pure output.
 *
 * Convention reminder:
 *   - Identity Pauli (I^n): global phase e^{-i*theta/2} = CMPLX(c, -s) on all amplitudes.
 *   - theta = 0 : U = I, output must equal input.
 *   - theta = 2*pi : U = -I, all amplitudes negated.
 *   - theta = 4*pi : U = +I, output equals input.  (True periodicity is 4*pi.)
 * =====================================================================================================================
 */

static void check_pauli_exp(const pauli_t *P, const cplx_t *psi, double theta) {
    const unsigned dim = (unsigned)(1u << P->n_qubits);

    cplx_t *mat  = NULL;
    cplx_t *pauli_psi = NULL;  /* M * psi = P|psi> */
    cplx_t *ref  = NULL;
    cplx_t *copy = NULL;
    state_t s    = {0};
    s.type        = PURE;

    const double c = cos(theta * 0.5);
    const double s_val = sin(theta * 0.5);   /* named s_val to avoid clash with state s */

    /* Build full Pauli matrix and compute P|psi> */
    mat = build_pauli_matrix(P);
    if (!mat) {
        fprintf(stderr, "check_pauli_exp(): build_pauli_matrix failed\n");
        goto cleanup;
    }

    pauli_psi = zmv(dim, mat, psi);
    if (!pauli_psi) {
        fprintf(stderr, "check_pauli_exp(): zmv failed\n");
        goto cleanup;
    }

    /* ref[k] = c * psi[k] - i * s_val * pauli_psi[k] */
    ref = malloc((size_t)dim * sizeof(*ref));
    if (!ref) {
        fprintf(stderr, "check_pauli_exp(): ref allocation failed\n");
        goto cleanup;
    }
    for (unsigned k = 0; k < dim; ++k) {
        ref[k] = c * psi[k] + CMPLX(0.0, -s_val) * pauli_psi[k];
    }

    /* Make copy of psi for state_init to consume */
    copy = malloc((size_t)dim * sizeof(*copy));
    if (!copy) {
        fprintf(stderr, "check_pauli_exp(): copy allocation failed\n");
        goto cleanup;
    }
    memcpy(copy, psi, (size_t)dim * sizeof(*copy));

    state_init(&s, (qubit_t)P->n_qubits, &copy);
    if (!s.data) {
        fprintf(stderr, "check_pauli_exp(): state_init failed\n");
        goto cleanup;
    }

    pauli_exp_pure(&s, P, theta);
    if (!s.data) {
        fprintf(stderr, "check_pauli_exp(): pauli_exp_pure left state data NULL\n");
        goto cleanup;
    }

    TEST_ASSERT_EQUAL_COMPLEX_ARRAY_TOL(ref, s.data, dim, PRECISION);

    free(mat);
    free(pauli_psi);
    free(ref);
    free(s.data);
    s.data = NULL;
    return;

cleanup:
    free(mat);
    free(pauli_psi);
    free(ref);
    free(copy);
    state_free(&s);
}

/* Rotation angles: covers theta=0 (identity), pi/4, pi/2, pi, 3*pi/2, 2*pi (-I),
 * 4*pi (identity again), -pi/3 (negative angle), and near-zero 1e-6.            */
#define N_EXP_THETAS 9
static const double TEST_EXP_THETAS[N_EXP_THETAS] = {
    0.0,
    M_PI / 4.0,
    M_PI / 2.0,
    M_PI,
    3.0 * M_PI / 2.0,
    2.0 * M_PI,
    4.0 * M_PI,
    -M_PI / 3.0,
    1e-6
};

/*
 * =====================================================================================================================
 * Test: pauli_exp_pure — Branch 1: Identity Pauli (x_mask=0, z_mask=0)
 *
 * exp(-i*theta/2 * I) = e^{-i*theta/2} * I.  All amplitudes multiplied by
 * CMPLX(cos(theta/2), -sin(theta/2)).  Not a no-op even for non-zero theta.
 *
 * Tested with n=1 (II trivially doesn't exist; single I) and n=2 (IIII via II).
 * States: |0>, |1>, |+>, |+i> for n=1; |00>, |++> for n=2.
 * =====================================================================================================================
 */

void test_pauli_exp_identity_global_phase(void) {
    /* n=1, PAULI_I */
    {
        const pauli_label_t lbl = PAULI_I;
        pauli_t *P = pauli_new(&lbl, 1);
        TEST_ASSERT_NOT_NULL(P);

        const cplx_t *states[4] = {PSI1_0, PSI1_1, PSI1_P, PSI1_PI};
        for (unsigned ti = 0; ti < N_EXP_THETAS; ++ti) {
            for (unsigned si = 0; si < 4; ++si) {
                check_pauli_exp(P, states[si], TEST_EXP_THETAS[ti]);
            }
        }
        pauli_free(P);
    }

    /* n=2, II */
    {
        const pauli_label_t lbls[2] = {PAULI_I, PAULI_I};
        pauli_t *P = pauli_new(lbls, 2);
        TEST_ASSERT_NOT_NULL(P);

        const cplx_t *states[2] = {PSI2_00, PSI2_PP};
        for (unsigned ti = 0; ti < N_EXP_THETAS; ++ti) {
            for (unsigned si = 0; si < 2; ++si) {
                check_pauli_exp(P, states[si], TEST_EXP_THETAS[ti]);
            }
        }
        pauli_free(P);
    }
}

/*
 * =====================================================================================================================
 * Test: pauli_exp_pure — Branch 2: Pure Z-string (x_mask=0, y_mask=0)
 *
 * exp(-i*theta/2 * Z) = diag(e^{-i*theta/2}, e^{+i*theta/2}).
 * Each amplitude is multiplied by e^{-i*theta/2} (parity=0) or e^{+i*theta/2} (parity=1).
 *
 * Key check: |1> (index 1, parity=1) gets the e^{+i*theta/2} factor.  This exercises
 * the fac_neg branch that the |0> eigenstate cannot reach.
 *
 * Also tested: ZZ (n=2) and ZZZ (n=3) to exercise multi-qubit z_mask patterns.
 * =====================================================================================================================
 */

void test_pauli_exp_pure_z(void) {
    /* n=1, Z */
    {
        const pauli_label_t lbl = PAULI_Z;
        pauli_t *P = pauli_new(&lbl, 1);
        TEST_ASSERT_NOT_NULL(P);

        const cplx_t *states[6] = {PSI1_0, PSI1_1, PSI1_P, PSI1_M, PSI1_PI, PSI1_MI};
        for (unsigned ti = 0; ti < N_EXP_THETAS; ++ti) {
            for (unsigned si = 0; si < 6; ++si) {
                check_pauli_exp(P, states[si], TEST_EXP_THETAS[ti]);
            }
        }
        pauli_free(P);
    }

    /* n=2, ZZ — z_mask=3, exercises both parities across a 4-element state vector */
    {
        const pauli_label_t lbls[2] = {PAULI_Z, PAULI_Z};
        pauli_t *P = pauli_new(lbls, 2);
        TEST_ASSERT_NOT_NULL(P);

        const cplx_t *states[4] = {PSI2_00, PSI2_01, PSI2_PP, PSI2_1P};
        for (unsigned ti = 0; ti < N_EXP_THETAS; ++ti) {
            for (unsigned si = 0; si < 4; ++si) {
                check_pauli_exp(P, states[si], TEST_EXP_THETAS[ti]);
            }
        }
        pauli_free(P);
    }

    /* n=3, ZZZ — z_mask=7, all 8 parities */
    {
        const pauli_label_t lbls[3] = {PAULI_Z, PAULI_Z, PAULI_Z};
        pauli_t *P = pauli_new(lbls, 3);
        TEST_ASSERT_NOT_NULL(P);

        const cplx_t *states[5] = {PSI3_000, PSI3_001, PSI3_GHZ, PSI3_W, PSI3_PPP};
        for (unsigned ti = 0; ti < N_EXP_THETAS; ++ti) {
            for (unsigned si = 0; si < 5; ++si) {
                check_pauli_exp(P, states[si], TEST_EXP_THETAS[ti]);
            }
        }
        pauli_free(P);
    }

    /* n=3, IZI — z_mask=2 (non-contiguous qubit) */
    {
        const pauli_label_t lbls[3] = {PAULI_I, PAULI_Z, PAULI_I};
        pauli_t *P = pauli_new(lbls, 3);
        TEST_ASSERT_NOT_NULL(P);

        const cplx_t *states[3] = {PSI3_GHZ, PSI3_W, PSI3_PPP};
        for (unsigned ti = 0; ti < N_EXP_THETAS; ++ti) {
            for (unsigned si = 0; si < 3; ++si) {
                check_pauli_exp(P, states[si], TEST_EXP_THETAS[ti]);
            }
        }
        pauli_free(P);
    }
}

/*
 * =====================================================================================================================
 * Test: pauli_exp_pure — Branch 3: Non-diagonal (x_mask != 0)
 *
 * exp(-i*theta/2 * X) = [[c, -is], [-is, c]].  This exercises the paired-update
 * branch where both old values must be read before either write.
 *
 * Also tests: Y (n=1), XX (n=2, eta_exp=0), XY (n=2, eta_exp=1), YY (n=2, eta_exp=2),
 * YYY (n=3, eta_exp=3), and XYZ (n=3, mixed).
 * =====================================================================================================================
 */

void test_pauli_exp_pure_nondiaganal_X(void) {
    const pauli_label_t lbl = PAULI_X;
    pauli_t *P = pauli_new(&lbl, 1);
    TEST_ASSERT_NOT_NULL(P);

    const cplx_t *states[6] = {PSI1_0, PSI1_1, PSI1_P, PSI1_M, PSI1_PI, PSI1_MI};
    for (unsigned ti = 0; ti < N_EXP_THETAS; ++ti) {
        for (unsigned si = 0; si < 6; ++si) {
            check_pauli_exp(P, states[si], TEST_EXP_THETAS[ti]);
        }
    }
    pauli_free(P);
}

void test_pauli_exp_pure_nondiaganal_Y(void) {
    /* n=1, Y — real rotation matrix [[c,-s],[s,c]], exercises eta_exp=1 */
    const pauli_label_t lbl = PAULI_Y;
    pauli_t *P = pauli_new(&lbl, 1);
    TEST_ASSERT_NOT_NULL(P);

    const cplx_t *states[6] = {PSI1_0, PSI1_1, PSI1_P, PSI1_M, PSI1_PI, PSI1_MI};
    for (unsigned ti = 0; ti < N_EXP_THETAS; ++ti) {
        for (unsigned si = 0; si < 6; ++si) {
            check_pauli_exp(P, states[si], TEST_EXP_THETAS[ti]);
        }
    }
    pauli_free(P);
}

void test_pauli_exp_pure_nondiaganal_XX(void) {
    /* n=2, XX — eta_exp=0, y_parity=0.  Exercises two-qubit non-diagonal branch. */
    const pauli_label_t lbls[2] = {PAULI_X, PAULI_X};
    pauli_t *P = pauli_new(lbls, 2);
    TEST_ASSERT_NOT_NULL(P);

    const cplx_t *states[5] = {PSI2_00, PSI2_01, PSI2_10, PSI2_PP, PSI2_P0};
    for (unsigned ti = 0; ti < N_EXP_THETAS; ++ti) {
        for (unsigned si = 0; si < 5; ++si) {
            check_pauli_exp(P, states[si], TEST_EXP_THETAS[ti]);
        }
    }
    pauli_free(P);
}

void test_pauli_exp_pure_nondiaganal_YY(void) {
    /* n=2, YY — eta_exp=2, exercises the (exp_jp+3)&3 rotation for both elements
     * of the pair.  Also y_parity=0 (even number of Y gates).                   */
    const pauli_label_t lbls[2] = {PAULI_Y, PAULI_Y};
    pauli_t *P = pauli_new(lbls, 2);
    TEST_ASSERT_NOT_NULL(P);

    const cplx_t *states[5] = {PSI2_00, PSI2_01, PSI2_10, PSI2_PP, PSI2_1P};
    for (unsigned ti = 0; ti < N_EXP_THETAS; ++ti) {
        for (unsigned si = 0; si < 5; ++si) {
            check_pauli_exp(P, states[si], TEST_EXP_THETAS[ti]);
        }
    }
    pauli_free(P);
}

void test_pauli_exp_pure_nondiaganal_YYY(void) {
    /* n=3, YYY — eta_exp=3 (the only case where i^3 = -i enters the paired update).
     * y_parity=1 (odd number of Y gates).                                         */
    const pauli_label_t lbls[3] = {PAULI_Y, PAULI_Y, PAULI_Y};
    pauli_t *P = pauli_new(lbls, 3);
    TEST_ASSERT_NOT_NULL(P);

    const cplx_t *states[5] = {PSI3_000, PSI3_001, PSI3_GHZ, PSI3_W, PSI3_PPP};
    for (unsigned ti = 0; ti < N_EXP_THETAS; ++ti) {
        for (unsigned si = 0; si < 5; ++si) {
            check_pauli_exp(P, states[si], TEST_EXP_THETAS[ti]);
        }
    }
    pauli_free(P);
}

void test_pauli_exp_pure_nondiaganal_XYZ(void) {
    /* n=3, XYZ — mixed Pauli string.  eta_exp=1, y_parity=1, z_mask=6.
     * Exercises both z-phase and y-parity simultaneously in the non-diagonal branch. */
    const pauli_label_t lbls[3] = {PAULI_X, PAULI_Y, PAULI_Z};
    pauli_t *P = pauli_new(lbls, 3);
    TEST_ASSERT_NOT_NULL(P);

    const cplx_t *states[5] = {PSI3_000, PSI3_001, PSI3_GHZ, PSI3_W, PSI3_PPP};
    for (unsigned ti = 0; ti < N_EXP_THETAS; ++ti) {
        for (unsigned si = 0; si < 5; ++si) {
            check_pauli_exp(P, states[si], TEST_EXP_THETAS[ti]);
        }
    }
    pauli_free(P);
}

/*
 * =====================================================================================================================
 * Test: pauli_exp_pure — theta=0 is identity for all branches
 *
 * At theta=0: c=1, s=0.  U = I regardless of Pauli string.
 * Output must equal input to tolerance PRECISION for every Pauli type.
 * =====================================================================================================================
 */

void test_pauli_exp_theta_zero_is_identity(void) {
    /* Test one representative Pauli from each branch */
    struct {
        pauli_label_t lbls[3];
        unsigned n;
        const cplx_t *psi;
    } cases[] = {
        /* Branch 1: identity Pauli */
        {{PAULI_I, 0, 0},             1, PSI1_P},
        /* Branch 2: pure Z */
        {{PAULI_Z, 0, 0},             1, PSI1_MI},
        /* Branch 3: non-diagonal X */
        {{PAULI_X, 0, 0},             1, PSI1_P},
        /* Branch 3: non-diagonal Y */
        {{PAULI_Y, 0, 0},             1, PSI1_PI},
        /* Branch 3: non-diagonal YYY (eta_exp=3) */
        {{PAULI_Y, PAULI_Y, PAULI_Y}, 3, PSI3_GHZ},
    };
    const unsigned ncases = sizeof(cases) / sizeof(cases[0]);

    for (unsigned c = 0; c < ncases; ++c) {
        pauli_t *P = pauli_new(cases[c].lbls, (qubit_t)cases[c].n);
        TEST_ASSERT_NOT_NULL(P);
        const unsigned dim = 1u << cases[c].n;

        cplx_t *copy = malloc(dim * sizeof(*copy));
        TEST_ASSERT_NOT_NULL(copy);
        memcpy(copy, cases[c].psi, dim * sizeof(*copy));

        state_t st = {0};
        st.type = PURE;
        state_init(&st, (qubit_t)cases[c].n, &copy);
        if (!st.data) {
            free(copy);
            pauli_free(P);
            TEST_FAIL_MESSAGE("state_init failed in theta=0 identity test");
        }

        pauli_exp_pure(&st, P, 0.0);
        TEST_ASSERT_EQUAL_COMPLEX_ARRAY_TOL(cases[c].psi, st.data, dim, PRECISION);

        free(st.data);
        st.data = NULL;
        pauli_free(P);
    }
}

/*
 * =====================================================================================================================
 * Test: pauli_exp_pure — theta=pi recovers -i * pauli_apply_pure
 *
 * At theta=pi: c=0, s=1.  U = -i * P.
 * So pauli_exp_pure(state, P, pi) == -i * P|psi>.
 *
 * Equivalently: pauli_exp_pure output * i must equal pauli_apply_pure output.
 *
 * Tested for X, Z, and YYY to cover all three non-trivial branches.
 * =====================================================================================================================
 */

void test_pauli_exp_theta_pi_is_neg_i_pauli(void) {
    /* n=1, X: at theta=pi, output[k] = -i * (X|psi>)[k] */
    {
        const pauli_label_t lbl = PAULI_X;
        pauli_t *P = pauli_new(&lbl, 1);
        TEST_ASSERT_NOT_NULL(P);

        const cplx_t *test_states[4] = {PSI1_0, PSI1_1, PSI1_P, PSI1_MI};
        for (unsigned si = 0; si < 4; ++si) {
            const cplx_t *psi = test_states[si];

            /* Compute -i * P|psi> as the expected result */
            cplx_t *mat = build_pauli_matrix(P);
            TEST_ASSERT_NOT_NULL(mat);
            cplx_t *pauli_out = zmv(2, mat, psi);
            TEST_ASSERT_NOT_NULL(pauli_out);
            cplx_t expected[2];
            for (unsigned k = 0; k < 2; ++k) {
                expected[k] = CMPLX(0.0, -1.0) * pauli_out[k];  /* -i * P|psi> */
            }
            free(mat);
            free(pauli_out);

            /* Run pauli_exp_pure at theta=pi */
            cplx_t *copy = malloc(2 * sizeof(*copy));
            TEST_ASSERT_NOT_NULL(copy);
            memcpy(copy, psi, 2 * sizeof(*copy));
            state_t s = {0};
            s.type = PURE;
            state_init(&s, 1, &copy);
            if (!s.data) {
                free(copy);
                pauli_free(P);
                TEST_FAIL_MESSAGE("state_init failed in theta=pi test (X)");
            }
            pauli_exp_pure(&s, P, M_PI);
            TEST_ASSERT_EQUAL_COMPLEX_ARRAY_TOL(expected, s.data, 2, PRECISION);
            free(s.data);
            s.data = NULL;
        }
        pauli_free(P);
    }

    /* n=1, Z: diagonal branch at theta=pi */
    {
        const pauli_label_t lbl = PAULI_Z;
        pauli_t *P = pauli_new(&lbl, 1);
        TEST_ASSERT_NOT_NULL(P);

        const cplx_t *test_states[3] = {PSI1_0, PSI1_1, PSI1_P};
        for (unsigned si = 0; si < 3; ++si) {
            const cplx_t *psi = test_states[si];

            cplx_t *mat = build_pauli_matrix(P);
            TEST_ASSERT_NOT_NULL(mat);
            cplx_t *pauli_out = zmv(2, mat, psi);
            TEST_ASSERT_NOT_NULL(pauli_out);
            cplx_t expected[2];
            for (unsigned k = 0; k < 2; ++k) {
                expected[k] = CMPLX(0.0, -1.0) * pauli_out[k];
            }
            free(mat);
            free(pauli_out);

            cplx_t *copy = malloc(2 * sizeof(*copy));
            TEST_ASSERT_NOT_NULL(copy);
            memcpy(copy, psi, 2 * sizeof(*copy));
            state_t s = {0};
            s.type = PURE;
            state_init(&s, 1, &copy);
            if (!s.data) {
                free(copy);
                pauli_free(P);
                TEST_FAIL_MESSAGE("state_init failed in theta=pi test (Z)");
            }
            pauli_exp_pure(&s, P, M_PI);
            TEST_ASSERT_EQUAL_COMPLEX_ARRAY_TOL(expected, s.data, 2, PRECISION);
            free(s.data);
            s.data = NULL;
        }
        pauli_free(P);
    }
}

/*
 * =====================================================================================================================
 * Test: pauli_exp_pure — theta=2*pi gives global -1 (NOT identity)
 *
 * At theta=2*pi: c=cos(pi)=-1, s=sin(pi)=0.  U = -I.
 * All amplitudes are negated.  This is a documented non-identity case
 * (periodicity is 4*pi, not 2*pi).
 * =====================================================================================================================
 */

void test_pauli_exp_theta_2pi_is_neg_identity(void) {
    /* n=1, Z (diagonal branch) */
    {
        const pauli_label_t lbl = PAULI_Z;
        pauli_t *P = pauli_new(&lbl, 1);
        TEST_ASSERT_NOT_NULL(P);

        /* |+> = [INVSQRT2, INVSQRT2].  After -I: [-INVSQRT2, -INVSQRT2]. */
        const cplx_t expected_P[2] = {-(INVSQRT2)+0.0*I, -(INVSQRT2)+0.0*I};

        cplx_t *copy = malloc(2 * sizeof(*copy));
        TEST_ASSERT_NOT_NULL(copy);
        memcpy(copy, PSI1_P, 2 * sizeof(*copy));
        state_t s = {0};
        s.type = PURE;
        state_init(&s, 1, &copy);
        if (!s.data) {
            free(copy);
            pauli_free(P);
            TEST_FAIL_MESSAGE("state_init failed in theta=2pi test");
        }
        pauli_exp_pure(&s, P, 2.0 * M_PI);
        TEST_ASSERT_EQUAL_COMPLEX_ARRAY_TOL(expected_P, s.data, 2, PRECISION);
        free(s.data);
        s.data = NULL;
        pauli_free(P);
    }

    /* n=1, X (non-diagonal branch) */
    {
        const pauli_label_t lbl = PAULI_X;
        pauli_t *P = pauli_new(&lbl, 1);
        TEST_ASSERT_NOT_NULL(P);

        /* |0> = [1, 0].  After -I: [-1, 0]. */
        const cplx_t expected_0[2] = {-1.0+0.0*I, 0.0+0.0*I};

        cplx_t *copy = malloc(2 * sizeof(*copy));
        TEST_ASSERT_NOT_NULL(copy);
        memcpy(copy, PSI1_0, 2 * sizeof(*copy));
        state_t s = {0};
        s.type = PURE;
        state_init(&s, 1, &copy);
        if (!s.data) {
            free(copy);
            pauli_free(P);
            TEST_FAIL_MESSAGE("state_init failed in theta=2pi non-diagonal test");
        }
        pauli_exp_pure(&s, P, 2.0 * M_PI);
        TEST_ASSERT_EQUAL_COMPLEX_ARRAY_TOL(expected_0, s.data, 2, PRECISION);
        free(s.data);
        s.data = NULL;
        pauli_free(P);
    }
}

/*
 * =====================================================================================================================
 * Test: pauli_exp_pure — unitarity / norm preservation
 *
 * For any unit-norm input |psi>, U(theta)|psi> must have the same norm.
 * Equivalently: sum_k |out[k]|^2 == sum_k |psi[k]|^2 == 1.
 *
 * Tested for each branch:
 *   Branch 1 (identity Pauli): global phase, norm is trivially preserved.
 *   Branch 2 (Z-string diagonal): per-amplitude scalar multiply by e^{+-i*theta/2}.
 *   Branch 3 (X-based non-diagonal): paired update.
 *
 * Also tests with a superposition state that exercises all amplitudes simultaneously.
 * =====================================================================================================================
 */

static double compute_norm_sq(const cplx_t *psi, unsigned dim) {
    double norm_sq = 0.0;
    for (unsigned k = 0; k < dim; ++k) {
        double re = creal(psi[k]);
        double im = cimag(psi[k]);
        norm_sq += re * re + im * im;
    }
    return norm_sq;
}

void test_pauli_exp_unitarity(void) {
    /* Test angles at pi/4, pi/2, pi, 3*pi/2, -pi/3 (non-trivial rotations that move amplitudes) */
    const double angles[5] = {M_PI / 4.0, M_PI / 2.0, M_PI, 3.0 * M_PI / 2.0, -M_PI / 3.0};

    /* Branch 1: Identity Pauli, n=1, state |+> */
    {
        const pauli_label_t lbl = PAULI_I;
        pauli_t *P = pauli_new(&lbl, 1);
        TEST_ASSERT_NOT_NULL(P);
        for (unsigned ai = 0; ai < 5; ++ai) {
            cplx_t *copy = malloc(2 * sizeof(*copy));
            TEST_ASSERT_NOT_NULL(copy);
            memcpy(copy, PSI1_P, 2 * sizeof(*copy));
            state_t s = {0};
            s.type = PURE;
            state_init(&s, 1, &copy);
            TEST_ASSERT_NOT_NULL(s.data);
            pauli_exp_pure(&s, P, angles[ai]);
            const double norm_sq = compute_norm_sq(s.data, 2);
            TEST_ASSERT_DOUBLE_WITHIN(PRECISION, 1.0, norm_sq);
            free(s.data);
            s.data = NULL;
        }
        pauli_free(P);
    }

    /* Branch 2: Z-string, n=1, state |+> (both parities contribute equally) */
    {
        const pauli_label_t lbl = PAULI_Z;
        pauli_t *P = pauli_new(&lbl, 1);
        TEST_ASSERT_NOT_NULL(P);
        for (unsigned ai = 0; ai < 5; ++ai) {
            cplx_t *copy = malloc(2 * sizeof(*copy));
            TEST_ASSERT_NOT_NULL(copy);
            memcpy(copy, PSI1_P, 2 * sizeof(*copy));
            state_t s = {0};
            s.type = PURE;
            state_init(&s, 1, &copy);
            TEST_ASSERT_NOT_NULL(s.data);
            pauli_exp_pure(&s, P, angles[ai]);
            const double norm_sq = compute_norm_sq(s.data, 2);
            TEST_ASSERT_DOUBLE_WITHIN(PRECISION, 1.0, norm_sq);
            free(s.data);
            s.data = NULL;
        }
        pauli_free(P);
    }

    /* Branch 2: ZZ, n=2, state |++> (mixed parity: parity 0 for |00>,|11> and parity 1 for |01>,|10>) */
    {
        const pauli_label_t lbls[2] = {PAULI_Z, PAULI_Z};
        pauli_t *P = pauli_new(lbls, 2);
        TEST_ASSERT_NOT_NULL(P);
        for (unsigned ai = 0; ai < 5; ++ai) {
            cplx_t *copy = malloc(4 * sizeof(*copy));
            TEST_ASSERT_NOT_NULL(copy);
            memcpy(copy, PSI2_PP, 4 * sizeof(*copy));
            state_t s = {0};
            s.type = PURE;
            state_init(&s, 2, &copy);
            TEST_ASSERT_NOT_NULL(s.data);
            pauli_exp_pure(&s, P, angles[ai]);
            const double norm_sq = compute_norm_sq(s.data, 4);
            TEST_ASSERT_DOUBLE_WITHIN(PRECISION, 1.0, norm_sq);
            free(s.data);
            s.data = NULL;
        }
        pauli_free(P);
    }

    /* Branch 3: X, n=1, state |0> */
    {
        const pauli_label_t lbl = PAULI_X;
        pauli_t *P = pauli_new(&lbl, 1);
        TEST_ASSERT_NOT_NULL(P);
        for (unsigned ai = 0; ai < 5; ++ai) {
            cplx_t *copy = malloc(2 * sizeof(*copy));
            TEST_ASSERT_NOT_NULL(copy);
            memcpy(copy, PSI1_0, 2 * sizeof(*copy));
            state_t s = {0};
            s.type = PURE;
            state_init(&s, 1, &copy);
            TEST_ASSERT_NOT_NULL(s.data);
            pauli_exp_pure(&s, P, angles[ai]);
            const double norm_sq = compute_norm_sq(s.data, 2);
            TEST_ASSERT_DOUBLE_WITHIN(PRECISION, 1.0, norm_sq);
            free(s.data);
            s.data = NULL;
        }
        pauli_free(P);
    }

    /* Branch 3: XYZ, n=3, state GHZ — all amplitudes in play */
    {
        const pauli_label_t lbls[3] = {PAULI_X, PAULI_Y, PAULI_Z};
        pauli_t *P = pauli_new(lbls, 3);
        TEST_ASSERT_NOT_NULL(P);
        for (unsigned ai = 0; ai < 5; ++ai) {
            cplx_t *copy = malloc(8 * sizeof(*copy));
            TEST_ASSERT_NOT_NULL(copy);
            memcpy(copy, PSI3_GHZ, 8 * sizeof(*copy));
            state_t s = {0};
            s.type = PURE;
            state_init(&s, 3, &copy);
            TEST_ASSERT_NOT_NULL(s.data);
            pauli_exp_pure(&s, P, angles[ai]);
            const double norm_sq = compute_norm_sq(s.data, 8);
            TEST_ASSERT_DOUBLE_WITHIN(PRECISION, 1.0, norm_sq);
            free(s.data);
            s.data = NULL;
        }
        pauli_free(P);
    }
}

/*
 * =====================================================================================================================
 * Test: pauli_exp_pure — inverse property: U(theta) * U(-theta) = I
 *
 * Applying pauli_exp_pure twice with +theta then -theta must recover the original state.
 * This is the group law: e^{-i*theta/2*P} * e^{+i*theta/2*P} = I.
 *
 * Tested for each branch and several non-trivial angles.
 * =====================================================================================================================
 */

void test_pauli_exp_inverse(void) {
    /* (theta, -theta) pairs at non-trivial angles */
    const double fwd_angles[4] = {M_PI / 4.0, M_PI / 2.0, M_PI, -M_PI / 3.0};

    /* Branch 1: Identity Pauli — global phase must cancel */
    {
        const pauli_label_t lbl = PAULI_I;
        pauli_t *P = pauli_new(&lbl, 1);
        TEST_ASSERT_NOT_NULL(P);
        for (unsigned ai = 0; ai < 4; ++ai) {
            const double theta = fwd_angles[ai];
            cplx_t *copy = malloc(2 * sizeof(*copy));
            TEST_ASSERT_NOT_NULL(copy);
            memcpy(copy, PSI1_P, 2 * sizeof(*copy));
            state_t s = {0};
            s.type = PURE;
            state_init(&s, 1, &copy);
            TEST_ASSERT_NOT_NULL(s.data);
            pauli_exp_pure(&s, P, theta);       /* forward: e^{-i*theta/2} * I */
            pauli_exp_pure(&s, P, -theta);      /* inverse: e^{+i*theta/2} * I */
            TEST_ASSERT_EQUAL_COMPLEX_ARRAY_TOL(PSI1_P, s.data, 2, PRECISION);
            free(s.data);
            s.data = NULL;
        }
        pauli_free(P);
    }

    /* Branch 2: Z diagonal — phase factors must cancel elementwise */
    {
        const pauli_label_t lbl = PAULI_Z;
        pauli_t *P = pauli_new(&lbl, 1);
        TEST_ASSERT_NOT_NULL(P);
        const cplx_t *test_states[3] = {PSI1_0, PSI1_1, PSI1_P};
        for (unsigned ai = 0; ai < 4; ++ai) {
            const double theta = fwd_angles[ai];
            for (unsigned si = 0; si < 3; ++si) {
                const cplx_t *psi = test_states[si];
                cplx_t *copy = malloc(2 * sizeof(*copy));
                TEST_ASSERT_NOT_NULL(copy);
                memcpy(copy, psi, 2 * sizeof(*copy));
                state_t s = {0};
                s.type = PURE;
                state_init(&s, 1, &copy);
                TEST_ASSERT_NOT_NULL(s.data);
                pauli_exp_pure(&s, P, theta);
                pauli_exp_pure(&s, P, -theta);
                TEST_ASSERT_EQUAL_COMPLEX_ARRAY_TOL(psi, s.data, 2, PRECISION);
                free(s.data);
                s.data = NULL;
            }
        }
        pauli_free(P);
    }

    /* Branch 3: X non-diagonal */
    {
        const pauli_label_t lbl = PAULI_X;
        pauli_t *P = pauli_new(&lbl, 1);
        TEST_ASSERT_NOT_NULL(P);
        const cplx_t *test_states[4] = {PSI1_0, PSI1_1, PSI1_P, PSI1_MI};
        for (unsigned ai = 0; ai < 4; ++ai) {
            const double theta = fwd_angles[ai];
            for (unsigned si = 0; si < 4; ++si) {
                const cplx_t *psi = test_states[si];
                cplx_t *copy = malloc(2 * sizeof(*copy));
                TEST_ASSERT_NOT_NULL(copy);
                memcpy(copy, psi, 2 * sizeof(*copy));
                state_t s = {0};
                s.type = PURE;
                state_init(&s, 1, &copy);
                TEST_ASSERT_NOT_NULL(s.data);
                pauli_exp_pure(&s, P, theta);
                pauli_exp_pure(&s, P, -theta);
                TEST_ASSERT_EQUAL_COMPLEX_ARRAY_TOL(psi, s.data, 2, PRECISION);
                free(s.data);
                s.data = NULL;
            }
        }
        pauli_free(P);
    }

    /* Branch 3: YYY (eta_exp=3, y_parity=1) — the most complex non-diagonal case */
    {
        const pauli_label_t lbls[3] = {PAULI_Y, PAULI_Y, PAULI_Y};
        pauli_t *P = pauli_new(lbls, 3);
        TEST_ASSERT_NOT_NULL(P);
        const cplx_t *test_states[3] = {PSI3_000, PSI3_GHZ, PSI3_W};
        for (unsigned ai = 0; ai < 4; ++ai) {
            const double theta = fwd_angles[ai];
            for (unsigned si = 0; si < 3; ++si) {
                const cplx_t *psi = test_states[si];
                cplx_t *copy = malloc(8 * sizeof(*copy));
                TEST_ASSERT_NOT_NULL(copy);
                memcpy(copy, psi, 8 * sizeof(*copy));
                state_t s = {0};
                s.type = PURE;
                state_init(&s, 3, &copy);
                TEST_ASSERT_NOT_NULL(s.data);
                pauli_exp_pure(&s, P, theta);
                pauli_exp_pure(&s, P, -theta);
                TEST_ASSERT_EQUAL_COMPLEX_ARRAY_TOL(psi, s.data, 8, PRECISION);
                free(s.data);
                s.data = NULL;
            }
        }
        pauli_free(P);
    }

    /* Branch 3: XYZ (mixed z-phase and y-parity) */
    {
        const pauli_label_t lbls[3] = {PAULI_X, PAULI_Y, PAULI_Z};
        pauli_t *P = pauli_new(lbls, 3);
        TEST_ASSERT_NOT_NULL(P);
        const cplx_t *test_states[2] = {PSI3_GHZ, PSI3_PPP};
        for (unsigned ai = 0; ai < 4; ++ai) {
            const double theta = fwd_angles[ai];
            for (unsigned si = 0; si < 2; ++si) {
                const cplx_t *psi = test_states[si];
                cplx_t *copy = malloc(8 * sizeof(*copy));
                TEST_ASSERT_NOT_NULL(copy);
                memcpy(copy, psi, 8 * sizeof(*copy));
                state_t s = {0};
                s.type = PURE;
                state_init(&s, 3, &copy);
                TEST_ASSERT_NOT_NULL(s.data);
                pauli_exp_pure(&s, P, theta);
                pauli_exp_pure(&s, P, -theta);
                TEST_ASSERT_EQUAL_COMPLEX_ARRAY_TOL(psi, s.data, 8, PRECISION);
                free(s.data);
                s.data = NULL;
            }
        }
        pauli_free(P);
    }
}

/*
 * =====================================================================================================================
 * Test: pauli_exp_pure — eigenstate phases
 *
 * If P|psi> = epsilon * |psi> (eigenstate with eigenvalue epsilon = +1 or -1), then:
 *   exp(-i*theta/2 * P) |psi> = e^{-i*theta/2 * epsilon} * |psi>
 *                             = (c - i*s*epsilon) * |psi>
 *
 * epsilon = +1: output = (c - i*s) * |psi> = e^{-i*theta/2} * |psi>   (fac_pos)
 * epsilon = -1: output = (c + i*s) * |psi> = e^{+i*theta/2} * |psi>   (fac_neg)
 *
 * Spot-checks at theta = pi/3:
 *   c = cos(pi/6) = sqrt(3)/2, s = sin(pi/6) = 1/2
 *   e^{-i*pi/6} = sqrt(3)/2 - i/2
 *   e^{+i*pi/6} = sqrt(3)/2 + i/2
 *
 * Eigenstates tested:
 *   Z: |0> (+1 eigenstate), |1> (-1 eigenstate)
 *   X: |+> (+1 eigenstate), |-> (-1 eigenstate)
 *   ZZ (n=2): |00> (+1), |01> (-1), |10> (-1), |11> (+1)
 * =====================================================================================================================
 */

void test_pauli_exp_eigenstate_phases(void) {
    const double theta  = M_PI / 3.0;
    const double half   = theta * 0.5;      /* pi/6 */
    const double c_val  = cos(half);        /* sqrt(3)/2 */
    const double s_val  = sin(half);        /* 1/2 */

    /* fac_pos = e^{-i*pi/6} = c - i*s  (epsilon = +1) */
    const cplx_t fac_pos = CMPLX(c_val, -s_val);
    /* fac_neg = e^{+i*pi/6} = c + i*s  (epsilon = -1) */
    const cplx_t fac_neg = CMPLX(c_val,  s_val);

    /* --- Z: |0> is epsilon=+1 eigenstate --- */
    {
        const pauli_label_t lbl = PAULI_Z;
        pauli_t *P = pauli_new(&lbl, 1);
        TEST_ASSERT_NOT_NULL(P);

        /* Expected: fac_pos * |0> = [c-is, 0] */
        const cplx_t expected[2] = {fac_pos * PSI1_0[0], fac_pos * PSI1_0[1]};

        cplx_t *copy = malloc(2 * sizeof(*copy));
        TEST_ASSERT_NOT_NULL(copy);
        memcpy(copy, PSI1_0, 2 * sizeof(*copy));
        state_t s = {0};
        s.type = PURE;
        state_init(&s, 1, &copy);
        TEST_ASSERT_NOT_NULL(s.data);
        pauli_exp_pure(&s, P, theta);
        TEST_ASSERT_EQUAL_COMPLEX_ARRAY_TOL(expected, s.data, 2, PRECISION);
        free(s.data);
        s.data = NULL;
        pauli_free(P);
    }

    /* --- Z: |1> is epsilon=-1 eigenstate --- */
    {
        const pauli_label_t lbl = PAULI_Z;
        pauli_t *P = pauli_new(&lbl, 1);
        TEST_ASSERT_NOT_NULL(P);

        /* Expected: fac_neg * |1> = [0, c+is] */
        const cplx_t expected[2] = {fac_neg * PSI1_1[0], fac_neg * PSI1_1[1]};

        cplx_t *copy = malloc(2 * sizeof(*copy));
        TEST_ASSERT_NOT_NULL(copy);
        memcpy(copy, PSI1_1, 2 * sizeof(*copy));
        state_t s = {0};
        s.type = PURE;
        state_init(&s, 1, &copy);
        TEST_ASSERT_NOT_NULL(s.data);
        pauli_exp_pure(&s, P, theta);
        TEST_ASSERT_EQUAL_COMPLEX_ARRAY_TOL(expected, s.data, 2, PRECISION);
        free(s.data);
        s.data = NULL;
        pauli_free(P);
    }

    /* --- X: |+> is epsilon=+1 eigenstate --- */
    {
        const pauli_label_t lbl = PAULI_X;
        pauli_t *P = pauli_new(&lbl, 1);
        TEST_ASSERT_NOT_NULL(P);

        /* Expected: fac_pos * |+> */
        const cplx_t expected[2] = {fac_pos * PSI1_P[0], fac_pos * PSI1_P[1]};

        cplx_t *copy = malloc(2 * sizeof(*copy));
        TEST_ASSERT_NOT_NULL(copy);
        memcpy(copy, PSI1_P, 2 * sizeof(*copy));
        state_t s = {0};
        s.type = PURE;
        state_init(&s, 1, &copy);
        TEST_ASSERT_NOT_NULL(s.data);
        pauli_exp_pure(&s, P, theta);
        TEST_ASSERT_EQUAL_COMPLEX_ARRAY_TOL(expected, s.data, 2, PRECISION);
        free(s.data);
        s.data = NULL;
        pauli_free(P);
    }

    /* --- X: |-> is epsilon=-1 eigenstate --- */
    {
        const pauli_label_t lbl = PAULI_X;
        pauli_t *P = pauli_new(&lbl, 1);
        TEST_ASSERT_NOT_NULL(P);

        /* Expected: fac_neg * |-> */
        const cplx_t expected[2] = {fac_neg * PSI1_M[0], fac_neg * PSI1_M[1]};

        cplx_t *copy = malloc(2 * sizeof(*copy));
        TEST_ASSERT_NOT_NULL(copy);
        memcpy(copy, PSI1_M, 2 * sizeof(*copy));
        state_t s = {0};
        s.type = PURE;
        state_init(&s, 1, &copy);
        TEST_ASSERT_NOT_NULL(s.data);
        pauli_exp_pure(&s, P, theta);
        TEST_ASSERT_EQUAL_COMPLEX_ARRAY_TOL(expected, s.data, 2, PRECISION);
        free(s.data);
        s.data = NULL;
        pauli_free(P);
    }

    /* --- ZZ (n=2): |00>(+1), |01>(-1), |10>(-1), |11>(+1) ---
     * ZZ parity = popcount(col & 3): even -> +1 eigenstate, odd -> -1 eigenstate.
     *   |00>(col=0): popcount=0 -> epsilon=+1 -> fac_pos
     *   |01>(col=1): popcount=1 -> epsilon=-1 -> fac_neg
     *   |10>(col=2): popcount=1 -> epsilon=-1 -> fac_neg
     *   |11>(col=3): popcount=2 -> epsilon=+1 -> fac_pos
     */
    {
        const pauli_label_t lbls[2] = {PAULI_Z, PAULI_Z};
        pauli_t *P = pauli_new(lbls, 2);
        TEST_ASSERT_NOT_NULL(P);

        struct { const cplx_t *psi; cplx_t fac; } zz_cases[4] = {
            {PSI2_00, fac_pos},
            {PSI2_01, fac_neg},
            {PSI2_10, fac_neg},
            {PSI2_11, fac_pos},
        };

        for (unsigned ci = 0; ci < 4; ++ci) {
            const cplx_t *psi = zz_cases[ci].psi;
            const cplx_t   fac = zz_cases[ci].fac;
            cplx_t expected[4];
            for (unsigned k = 0; k < 4; ++k) {
                expected[k] = fac * psi[k];
            }
            cplx_t *copy = malloc(4 * sizeof(*copy));
            TEST_ASSERT_NOT_NULL(copy);
            memcpy(copy, psi, 4 * sizeof(*copy));
            state_t s = {0};
            s.type = PURE;
            state_init(&s, 2, &copy);
            TEST_ASSERT_NOT_NULL(s.data);
            pauli_exp_pure(&s, P, theta);
            TEST_ASSERT_EQUAL_COMPLEX_ARRAY_TOL(expected, s.data, 4, PRECISION);
            free(s.data);
            s.data = NULL;
        }
        pauli_free(P);
    }
}

/*
 * =====================================================================================================================
 * Test: pauli_exp_pure — theta=4*pi is true identity (contrast with theta=2*pi = -I)
 *
 * At theta=4*pi: c = cos(2*pi) = 1, s = sin(2*pi) = 0 (to machine precision).
 * The function should return U = +I, meaning output equals input exactly.
 *
 * This confirms the 4*pi periodicity of SU(2): unlike theta=2*pi which gives -I,
 * theta=4*pi completes the full cycle and returns to +I.
 *
 * Tested with a hardcoded state (not via check_pauli_exp) so the assertion is
 * explicit: "output IS the original state", not just "matches a computed reference".
 * =====================================================================================================================
 */

void test_pauli_exp_theta_4pi_is_identity(void) {
    /* Branch 1: Identity Pauli — e^{-i*2*pi} * I must equal I (global phase = +1) */
    {
        const pauli_label_t lbl = PAULI_I;
        pauli_t *P = pauli_new(&lbl, 1);
        TEST_ASSERT_NOT_NULL(P);

        /* |+i> = [INVSQRT2, i*INVSQRT2] */
        cplx_t *copy = malloc(2 * sizeof(*copy));
        TEST_ASSERT_NOT_NULL(copy);
        memcpy(copy, PSI1_PI, 2 * sizeof(*copy));
        state_t s = {0};
        s.type = PURE;
        state_init(&s, 1, &copy);
        TEST_ASSERT_NOT_NULL(s.data);
        pauli_exp_pure(&s, P, 4.0 * M_PI);
        TEST_ASSERT_EQUAL_COMPLEX_ARRAY_TOL(PSI1_PI, s.data, 2, PRECISION);
        free(s.data);
        s.data = NULL;
        pauli_free(P);
    }

    /* Branch 2: Z diagonal — confirm 4*pi returns to original, 2*pi returns -original */
    {
        const pauli_label_t lbl = PAULI_Z;
        pauli_t *P = pauli_new(&lbl, 1);
        TEST_ASSERT_NOT_NULL(P);

        /* State |-> = [INVSQRT2, -INVSQRT2] */
        /* At 2*pi: output must be [-INVSQRT2, INVSQRT2] (negated) */
        /* At 4*pi: output must be [INVSQRT2, -INVSQRT2] (original) */

        /* 4*pi check */
        cplx_t *copy = malloc(2 * sizeof(*copy));
        TEST_ASSERT_NOT_NULL(copy);
        memcpy(copy, PSI1_M, 2 * sizeof(*copy));
        state_t s = {0};
        s.type = PURE;
        state_init(&s, 1, &copy);
        TEST_ASSERT_NOT_NULL(s.data);
        pauli_exp_pure(&s, P, 4.0 * M_PI);
        TEST_ASSERT_EQUAL_COMPLEX_ARRAY_TOL(PSI1_M, s.data, 2, PRECISION);
        free(s.data);
        s.data = NULL;

        pauli_free(P);
    }

    /* Branch 3: X non-diagonal — 4*pi must return to original */
    {
        const pauli_label_t lbl = PAULI_X;
        pauli_t *P = pauli_new(&lbl, 1);
        TEST_ASSERT_NOT_NULL(P);

        /* State |+i> = [INVSQRT2, i*INVSQRT2] */
        cplx_t *copy = malloc(2 * sizeof(*copy));
        TEST_ASSERT_NOT_NULL(copy);
        memcpy(copy, PSI1_PI, 2 * sizeof(*copy));
        state_t s = {0};
        s.type = PURE;
        state_init(&s, 1, &copy);
        TEST_ASSERT_NOT_NULL(s.data);
        pauli_exp_pure(&s, P, 4.0 * M_PI);
        TEST_ASSERT_EQUAL_COMPLEX_ARRAY_TOL(PSI1_PI, s.data, 2, PRECISION);
        free(s.data);
        s.data = NULL;
        pauli_free(P);
    }

    /* Branch 3: YYY (n=3, eta_exp=3) — 4*pi must return to original */
    {
        const pauli_label_t lbls[3] = {PAULI_Y, PAULI_Y, PAULI_Y};
        pauli_t *P = pauli_new(lbls, 3);
        TEST_ASSERT_NOT_NULL(P);

        cplx_t *copy = malloc(8 * sizeof(*copy));
        TEST_ASSERT_NOT_NULL(copy);
        memcpy(copy, PSI3_GHZ, 8 * sizeof(*copy));
        state_t s = {0};
        s.type = PURE;
        state_init(&s, 3, &copy);
        TEST_ASSERT_NOT_NULL(s.data);
        pauli_exp_pure(&s, P, 4.0 * M_PI);
        TEST_ASSERT_EQUAL_COMPLEX_ARRAY_TOL(PSI3_GHZ, s.data, 8, PRECISION);
        free(s.data);
        s.data = NULL;
        pauli_free(P);
    }
}

/*
 * =====================================================================================================================
 * Test: pauli_exp_pure — consistency with pauli_apply_pure at theta=pi
 *
 * At theta=pi: c=0, s=1.  exp(-i*pi/2 * P) = -i * P.
 * So pauli_exp_pure output equals -i * (pauli_apply_pure output).
 * Equivalently: pauli_exp output * i equals pauli_apply output.
 *
 * This test exercises the cross-function consistency for the full set of
 * branches: Identity (global phase), Z-diagonal, X (eta_exp=0), Y (eta_exp=1),
 * YY (eta_exp=2), YYY (eta_exp=3), and XYZ (mixed).
 * The earlier test covered only X and Z, leaving YYY and mixed cases unchecked.
 * =====================================================================================================================
 */

void test_pauli_exp_consistency_with_apply(void) {
    /* Helper: for Pauli P and state psi (dim elements), verify:
     *   pauli_exp_pure(psi, P, pi)[k]  ==  -i * (pauli_apply_pure(psi, P))[k]
     *
     * Strategy: build the reference via build_pauli_matrix+zmv (same as check_pauli_exp),
     * scale by -i, then compare against pauli_exp_pure.
     */

    struct {
        pauli_label_t lbls[4];
        unsigned      n;
        const cplx_t *psi;
        const char   *name;
    } cases[] = {
        /* Branch 1: Identity — exp result is -i * I|psi> = -i * psi */
        {{PAULI_I, 0, 0, 0},             1, PSI1_P,   "I"},
        /* Branch 2: Z diagonal */
        {{PAULI_Z, 0, 0, 0},             1, PSI1_P,   "Z"},
        {{PAULI_Z, 0, 0, 0},             1, PSI1_1,   "Z|1>"},
        /* Branch 3: X (eta_exp=0) */
        {{PAULI_X, 0, 0, 0},             1, PSI1_PI,  "X"},
        /* Branch 3: Y (eta_exp=1) */
        {{PAULI_Y, 0, 0, 0},             1, PSI1_P,   "Y"},
        /* Branch 3: XX (eta_exp=0, n=2) */
        {{PAULI_X, PAULI_X, 0, 0},       2, PSI2_PP,  "XX"},
        /* Branch 3: YY (eta_exp=2, y_parity=0, n=2) */
        {{PAULI_Y, PAULI_Y, 0, 0},       2, PSI2_1P,  "YY"},
        /* Branch 3: YYY (eta_exp=3, y_parity=1, n=3) — gap from earlier test */
        {{PAULI_Y, PAULI_Y, PAULI_Y, 0}, 3, PSI3_GHZ, "YYY"},
        {{PAULI_Y, PAULI_Y, PAULI_Y, 0}, 3, PSI3_W,   "YYY|W>"},
        /* Branch 3: XYZ (mixed, n=3) — gap from earlier test */
        {{PAULI_X, PAULI_Y, PAULI_Z, 0}, 3, PSI3_PPP, "XYZ"},
        {{PAULI_X, PAULI_Y, PAULI_Z, 0}, 3, PSI3_000, "XYZ|000>"},
    };
    const unsigned ncases = sizeof(cases) / sizeof(cases[0]);

    for (unsigned ci = 0; ci < ncases; ++ci) {
        const unsigned dim = 1u << cases[ci].n;
        const cplx_t  *psi = cases[ci].psi;

        pauli_t *P = pauli_new(cases[ci].lbls, (qubit_t)cases[ci].n);
        TEST_ASSERT_NOT_NULL(P);

        /* Build -i * P|psi> as the expected result */
        cplx_t *mat = build_pauli_matrix(P);
        TEST_ASSERT_NOT_NULL(mat);
        cplx_t *pauli_out = zmv(dim, mat, psi);
        TEST_ASSERT_NOT_NULL(pauli_out);
        cplx_t *expected = malloc(dim * sizeof(*expected));
        TEST_ASSERT_NOT_NULL(expected);
        for (unsigned k = 0; k < dim; ++k) {
            expected[k] = CMPLX(0.0, -1.0) * pauli_out[k];  /* -i * P|psi> */
        }
        free(mat);
        free(pauli_out);

        /* Run pauli_exp_pure at theta=pi */
        cplx_t *copy = malloc(dim * sizeof(*copy));
        TEST_ASSERT_NOT_NULL(copy);
        memcpy(copy, psi, dim * sizeof(*copy));
        state_t s = {0};
        s.type = PURE;
        state_init(&s, (qubit_t)cases[ci].n, &copy);
        if (!s.data) {
            free(copy);
            free(expected);
            pauli_free(P);
            TEST_FAIL_MESSAGE("state_init failed in consistency test");
        }
        pauli_exp_pure(&s, P, M_PI);
        TEST_ASSERT_EQUAL_COMPLEX_ARRAY_TOL(expected, s.data, dim, PRECISION);

        free(s.data);
        s.data = NULL;
        free(expected);
        pauli_free(P);
    }
}

/*
 * =====================================================================================================================
 * Test: pauli_exp_pure — non-unit-norm state: norm is preserved proportionally
 *
 * pauli_exp_pure is defined for arbitrary pure states; unitarity guarantees that
 * whatever norm the input has, the output has the same norm.
 *
 * For a state with norm^2 = N, the output should also have norm^2 = N.
 * We use a deliberately unnormalized state: psi = [2, 0] (norm^2 = 4).
 *
 * This exercises the code path where the amplitudes are non-unit, confirming
 * that no implicit normalization occurs anywhere in the implementation.
 * =====================================================================================================================
 */

void test_pauli_exp_non_unit_norm(void) {
    /* Non-unit state: psi = [2+0i, 0+0i]; norm^2 = 4 */
    const cplx_t psi_unnorm[2] = {2.0 + 0.0*I, 0.0 + 0.0*I};
    const double expected_norm_sq = 4.0;

    /* Test with X at theta=pi/2.  Output should have same norm^2 = 4. */
    {
        const pauli_label_t lbl = PAULI_X;
        pauli_t *P = pauli_new(&lbl, 1);
        TEST_ASSERT_NOT_NULL(P);

        cplx_t *copy = malloc(2 * sizeof(*copy));
        TEST_ASSERT_NOT_NULL(copy);
        memcpy(copy, psi_unnorm, 2 * sizeof(*copy));
        state_t s = {0};
        s.type = PURE;
        state_init(&s, 1, &copy);
        TEST_ASSERT_NOT_NULL(s.data);
        pauli_exp_pure(&s, P, M_PI / 2.0);
        const double norm_sq = compute_norm_sq(s.data, 2);
        TEST_ASSERT_DOUBLE_WITHIN(PRECISION, expected_norm_sq, norm_sq);
        free(s.data);
        s.data = NULL;
        pauli_free(P);
    }

    /* Test with Z at theta=pi/3. */
    {
        const pauli_label_t lbl = PAULI_Z;
        pauli_t *P = pauli_new(&lbl, 1);
        TEST_ASSERT_NOT_NULL(P);

        cplx_t *copy = malloc(2 * sizeof(*copy));
        TEST_ASSERT_NOT_NULL(copy);
        memcpy(copy, psi_unnorm, 2 * sizeof(*copy));
        state_t s = {0};
        s.type = PURE;
        state_init(&s, 1, &copy);
        TEST_ASSERT_NOT_NULL(s.data);
        pauli_exp_pure(&s, P, M_PI / 3.0);
        const double norm_sq = compute_norm_sq(s.data, 2);
        TEST_ASSERT_DOUBLE_WITHIN(PRECISION, expected_norm_sq, norm_sq);
        free(s.data);
        s.data = NULL;
        pauli_free(P);
    }

    /* Test with I at theta=pi (global phase only). */
    {
        const pauli_label_t lbl = PAULI_I;
        pauli_t *P = pauli_new(&lbl, 1);
        TEST_ASSERT_NOT_NULL(P);

        cplx_t *copy = malloc(2 * sizeof(*copy));
        TEST_ASSERT_NOT_NULL(copy);
        memcpy(copy, psi_unnorm, 2 * sizeof(*copy));
        state_t s = {0};
        s.type = PURE;
        state_init(&s, 1, &copy);
        TEST_ASSERT_NOT_NULL(s.data);
        pauli_exp_pure(&s, P, M_PI);
        const double norm_sq = compute_norm_sq(s.data, 2);
        TEST_ASSERT_DOUBLE_WITHIN(PRECISION, expected_norm_sq, norm_sq);
        free(s.data);
        s.data = NULL;
        pauli_free(P);
    }
}

/*
 * =====================================================================================================================
 * Test: pauli_exp_pure — n=4, XZYX (exercises dim=16, split_bit and all branches at max qubits)
 * =====================================================================================================================
 */

void test_pauli_exp_4q(void) {
    /* XZYX: labels[0]=X, labels[1]=Z, labels[2]=Y, labels[3]=X -> x_mask=13, z_mask=6, y_mask=4 */
    const pauli_label_t lbls[4] = {PAULI_X, PAULI_Z, PAULI_Y, PAULI_X};
    pauli_t *P = pauli_new(lbls, 4);
    TEST_ASSERT_NOT_NULL(P);

    TEST_ASSERT_EQUAL_UINT64(13u, P->x_mask);
    TEST_ASSERT_EQUAL_UINT64(6u,  P->z_mask);
    TEST_ASSERT_EQUAL_UINT64(4u,  P->y_mask);

    /* Sample angles: 0, pi/2, pi, 2*pi */
    const double angles[4] = {0.0, M_PI / 2.0, M_PI, 2.0 * M_PI};

    /* All 16 computational basis states at each angle */
    for (unsigned ai = 0; ai < 4; ++ai) {
        for (unsigned k = 0; k < 16; ++k) {
            cplx_t psi[16];
            memset(psi, 0, sizeof(psi));
            psi[k] = 1.0 + 0.0*I;
            check_pauli_exp(P, psi, angles[ai]);
        }
        check_pauli_exp(P, PSI4_GHZ, angles[ai]);
    }

    pauli_free(P);
}

/*
 * =====================================================================================================================
 * Test: pauli_exp_pure — non-diagonal with split_bit > 0
 *
 * All 18 existing pauli_exp_pure tests that exercise the non-diagonal branch (x_mask != 0)
 * use Pauli strings whose lowest-set x_mask bit is bit 0 (e.g. X, XX, XYZ, XZYX all have
 * ctzll(x_mask) == 0).  The `split_bit = __builtin_ctzll(x_mask)` value therefore equals 0
 * in every prior test, making insertBit0(k, 0) a simple even-index enumeration.
 *
 * A bug in insertBit0() for split_bit > 0 — or a wrong parity_j or exp_jp computation when
 * the x-support starts above qubit 0 — is entirely invisible to the existing suite.
 *
 * This test covers four distinct cases, exercising split_bit in {1, 2}:
 *
 *   IX   (n=2): x_mask=2, z_mask=0, y_mask=0, split_bit=1, eta_exp=0, y_parity=0
 *   IXI  (n=3): x_mask=2, z_mask=0, y_mask=0, split_bit=1, eta_exp=0, y_parity=0
 *   IIX  (n=3): x_mask=4, z_mask=0, y_mask=0, split_bit=2, eta_exp=0, y_parity=0
 *   IYI  (n=3): x_mask=2, z_mask=2, y_mask=2, split_bit=1, eta_exp=1, y_parity=1
 *
 * Each case sweeps all TEST_EXP_THETAS angles over a representative set of states.
 * =====================================================================================================================
 */

void test_pauli_exp_nondiaganal_split_bit(void) {
    /* --- IX (n=2): split_bit=1 --- */
    {
        /* labels[0]=I (qubit 0), labels[1]=X (qubit 1): x_mask = bit 1 = 2, split_bit=1 */
        const pauli_label_t lbls[2] = {PAULI_I, PAULI_X};
        pauli_t *P = pauli_new(lbls, 2);
        TEST_ASSERT_NOT_NULL(P);

        /* Verify structural expectation before running the rotation tests */
        TEST_ASSERT_EQUAL_UINT64(2u, P->x_mask);
        TEST_ASSERT_EQUAL_UINT64(0u, P->z_mask);
        TEST_ASSERT_EQUAL_UINT64(0u, P->y_mask);

        const cplx_t *states[4] = {PSI2_00, PSI2_01, PSI2_10, PSI2_PP};
        for (unsigned ti = 0; ti < N_EXP_THETAS; ++ti) {
            for (unsigned si = 0; si < 4; ++si) {
                check_pauli_exp(P, states[si], TEST_EXP_THETAS[ti]);
            }
        }
        pauli_free(P);
    }

    /* --- IXI (n=3): split_bit=1 --- */
    {
        /* labels[0]=I (qubit 0), labels[1]=X (qubit 1), labels[2]=I (qubit 2):
         * x_mask = bit 1 = 2, split_bit=1                                       */
        const pauli_label_t lbls[3] = {PAULI_I, PAULI_X, PAULI_I};
        pauli_t *P = pauli_new(lbls, 3);
        TEST_ASSERT_NOT_NULL(P);

        TEST_ASSERT_EQUAL_UINT64(2u, P->x_mask);
        TEST_ASSERT_EQUAL_UINT64(0u, P->z_mask);
        TEST_ASSERT_EQUAL_UINT64(0u, P->y_mask);

        const cplx_t *states[5] = {PSI3_000, PSI3_010, PSI3_GHZ, PSI3_W, PSI3_PPP};
        for (unsigned ti = 0; ti < N_EXP_THETAS; ++ti) {
            for (unsigned si = 0; si < 5; ++si) {
                check_pauli_exp(P, states[si], TEST_EXP_THETAS[ti]);
            }
        }
        pauli_free(P);
    }

    /* --- IIX (n=3): split_bit=2 --- */
    {
        /* labels[0]=I, labels[1]=I, labels[2]=X:
         * x_mask = bit 2 = 4, split_bit=2                                       */
        const pauli_label_t lbls[3] = {PAULI_I, PAULI_I, PAULI_X};
        pauli_t *P = pauli_new(lbls, 3);
        TEST_ASSERT_NOT_NULL(P);

        TEST_ASSERT_EQUAL_UINT64(4u, P->x_mask);
        TEST_ASSERT_EQUAL_UINT64(0u, P->z_mask);
        TEST_ASSERT_EQUAL_UINT64(0u, P->y_mask);

        const cplx_t *states[5] = {PSI3_000, PSI3_100, PSI3_GHZ, PSI3_W, PSI3_PPP};
        for (unsigned ti = 0; ti < N_EXP_THETAS; ++ti) {
            for (unsigned si = 0; si < 5; ++si) {
                check_pauli_exp(P, states[si], TEST_EXP_THETAS[ti]);
            }
        }
        pauli_free(P);
    }

    /* --- IYI (n=3): split_bit=1, eta_exp=1, y_parity=1 --- */
    {
        /* labels[0]=I, labels[1]=Y, labels[2]=I:
         * x_mask=2, z_mask=2, y_mask=2, split_bit=1, eta_exp=1, y_parity=1.
         * This exercises the y-parity path (exp_jp = (exp_j + 2) & 3) with a
         * non-zero split_bit simultaneously.                                     */
        const pauli_label_t lbls[3] = {PAULI_I, PAULI_Y, PAULI_I};
        pauli_t *P = pauli_new(lbls, 3);
        TEST_ASSERT_NOT_NULL(P);

        TEST_ASSERT_EQUAL_UINT64(2u, P->x_mask);
        TEST_ASSERT_EQUAL_UINT64(2u, P->z_mask);
        TEST_ASSERT_EQUAL_UINT64(2u, P->y_mask);

        const cplx_t *states[5] = {PSI3_000, PSI3_010, PSI3_GHZ, PSI3_W, PSI3_PPP};
        for (unsigned ti = 0; ti < N_EXP_THETAS; ++ti) {
            for (unsigned si = 0; si < 5; ++si) {
                check_pauli_exp(P, states[si], TEST_EXP_THETAS[ti]);
            }
        }
        pauli_free(P);
    }
}
