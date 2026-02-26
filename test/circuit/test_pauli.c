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
