// test_pauli.c - Unit tests for pauli_apply_pure and related Pauli string functions.
//
// Reference oracle: check_pauli_apply and check_pauli_exp each receive a big-endian
// Pauli string and build two independent paths:
//   - Kronecker path:  string → mat_pauli_str (gatemat.c, no masks) → zmv → ref
//   - Mask path:       string → pauli_from_str → pauli_apply_pure / pauli_exp_pure
// Both paths follow goto-cleanup for leak-free failure handling.

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
#include "gatemat.h"
#include "test_pure_states.h"

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/*
 * =====================================================================================================================
 * Helper: apply a Pauli string (given as big-endian char string) to a copy of psi
 * and compare against the Kronecker-product reference.
 *
 * Two completely independent paths share only the string input:
 *   - Reference:   pauli_str → mat_pauli_str (Kronecker products) → zmv → ref
 *   - Under test:  pauli_str → pauli_from_str (bitmasks) → pauli_apply_pure → result
 * =====================================================================================================================
 */

static void check_pauli_apply(const char *pauli_str, unsigned n_qubits, const cplx_t *psi) {
    const unsigned dim = (unsigned)(1u << n_qubits);

    cplx_t *mat  = NULL;
    cplx_t *ref  = NULL;
    cplx_t *copy = NULL;
    pauli_t *P   = NULL;
    state_t st   = {0};
    st.type       = PURE;

    /* Reference path: string → Kronecker matrix → M·psi */
    mat = mat_pauli_str(pauli_str, n_qubits);
    if (!mat) {
        fprintf(stderr, "check_pauli_apply(): mat_pauli_str failed\n");
        goto cleanup;
    }

    ref = zmv(dim, mat, psi);
    if (!ref) {
        fprintf(stderr, "check_pauli_apply(): zmv failed\n");
        goto cleanup;
    }

    /* Mask path: string → pauli_t → pauli_apply_pure */
    P = pauli_from_str(pauli_str, (qubit_t)n_qubits);

    copy = malloc((size_t)dim * sizeof(*copy));
    if (!copy) {
        fprintf(stderr, "check_pauli_apply(): copy allocation failed\n");
        goto cleanup;
    }
    memcpy(copy, psi, (size_t)dim * sizeof(*copy));

    state_init(&st, (qubit_t)n_qubits, &copy);
    if (!st.data) {
        fprintf(stderr, "check_pauli_apply(): state_init failed\n");
        goto cleanup;
    }

    pauli_apply_pure(&st, P);
    if (!st.data) {
        fprintf(stderr, "check_pauli_apply(): pauli_apply_pure left state data NULL\n");
        goto cleanup;
    }

    TEST_ASSERT_EQUAL_COMPLEX_ARRAY_TOL(ref, st.data, dim, PRECISION);

    free(mat);
    free(ref);
    pauli_free(P);
    free(st.data);
    st.data = NULL;
    return;

cleanup:
    free(mat);
    free(ref);
    free(copy);
    pauli_free(P);
    state_free(&st);
}

/* Test state vectors are generated per-test via test_mk_states_pure(n, &nvecs)
 * and freed via test_rm_states_pure(n, states).  No global state arrays needed. */

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
 * Test: pauli_to_str — serialise masks back to big-endian character string
 * =====================================================================================================================
 */

void test_pauli_to_str(void) {
    /* n=1 single-qubit cases */
    {
        const pauli_label_t lbl = PAULI_I;
        pauli_t *P = pauli_new(&lbl, 1);
        TEST_ASSERT_NOT_NULL(P);
        char *s = pauli_to_str(P);
        TEST_ASSERT_NOT_NULL(s);
        TEST_ASSERT_EQUAL_STRING("I", s);
        free(s);
        pauli_free(P);
    }
    {
        const pauli_label_t lbl = PAULI_X;
        pauli_t *P = pauli_new(&lbl, 1);
        TEST_ASSERT_NOT_NULL(P);
        char *s = pauli_to_str(P);
        TEST_ASSERT_NOT_NULL(s);
        TEST_ASSERT_EQUAL_STRING("X", s);
        free(s);
        pauli_free(P);
    }
    {
        const pauli_label_t lbl = PAULI_Y;
        pauli_t *P = pauli_new(&lbl, 1);
        TEST_ASSERT_NOT_NULL(P);
        char *s = pauli_to_str(P);
        TEST_ASSERT_NOT_NULL(s);
        TEST_ASSERT_EQUAL_STRING("Y", s);
        free(s);
        pauli_free(P);
    }
    {
        const pauli_label_t lbl = PAULI_Z;
        pauli_t *P = pauli_new(&lbl, 1);
        TEST_ASSERT_NOT_NULL(P);
        char *s = pauli_to_str(P);
        TEST_ASSERT_NOT_NULL(s);
        TEST_ASSERT_EQUAL_STRING("Z", s);
        free(s);
        pauli_free(P);
    }
    /* n=4, round-trip through pauli_from_str */
    {
        pauli_t *P = pauli_from_str("IXYZ", 4);
        TEST_ASSERT_NOT_NULL(P);
        char *s = pauli_to_str(P);
        TEST_ASSERT_NOT_NULL(s);
        TEST_ASSERT_EQUAL_STRING("IXYZ", s);
        free(s);
        pauli_free(P);
    }
    /* n=3, constructed via pauli_new (labels in little-endian order) */
    {
        const pauli_label_t lbls[3] = {PAULI_X, PAULI_Y, PAULI_Z};
        pauli_t *P = pauli_new(lbls, 3);
        TEST_ASSERT_NOT_NULL(P);
        char *s = pauli_to_str(P);
        TEST_ASSERT_NOT_NULL(s);
        /* labels[0]=X (q0), labels[1]=Y (q1), labels[2]=Z (q2)
         * big-endian: s[0]=q2='Z', s[1]=q1='Y', s[2]=q0='X' */
        TEST_ASSERT_EQUAL_STRING("ZYX", s);
        free(s);
        pauli_free(P);
    }
    /* n=0: empty string */
    {
        pauli_t *P = pauli_new(NULL, 0);
        TEST_ASSERT_NOT_NULL(P);
        char *s = pauli_to_str(P);
        TEST_ASSERT_NOT_NULL(s);
        TEST_ASSERT_EQUAL_STRING("", s);
        free(s);
        pauli_free(P);
    }
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
 * test_pauli_apply_all — all Pauli string application tests
 *
 * Each row of APPLY_CASES is one Pauli string (big-endian) and its qubit count.
 * For every row the full set of test states from test_mk_states_pure is iterated
 * and verified against the Kronecker-product reference in check_pauli_apply.
 * =====================================================================================================================
 */

static const struct { const char *str; unsigned n; } APPLY_CASES[] = {
    /* n=1 */
    {"I",    1}, {"X",    1}, {"Y",    1}, {"Z",    1},
    /* n=2 — big-endian strings; e.g. "ZX" = Z on q1, X on q0 */
    {"XX",   2}, {"ZZ",   2}, {"ZX",   2}, {"YY",   2},
    {"IZ",   2}, {"XI",   2}, {"YI",   2}, {"YZ",   2},
    /* n=3 */
    {"ZYX",  3}, {"ZZZ",  3}, {"IXI",  3}, {"XII",  3},
    {"IZI",  3}, {"ZIZ",  3}, {"YYY",  3}, {"III",  3},
    /* n=4 */
    {"XYZX", 4}, {"IIXI", 4}, {"IXII", 4}, {"XIII", 4},
};

void test_pauli_apply_all(void) {
    const unsigned ncases = sizeof(APPLY_CASES) / sizeof(APPLY_CASES[0]);
    for (unsigned c = 0; c < ncases; ++c) {
        const unsigned n = APPLY_CASES[c].n;
        unsigned nvecs; cplx_t **s = test_mk_states_pure(n, &nvecs);
        for (unsigned i = 0; i < nvecs; ++i)
            check_pauli_apply(APPLY_CASES[c].str, n, s[i]);
        test_rm_states_pure(n, s);
    }
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
    pauli_t *P_str = pauli_from_str("ZYX", 3);

    TEST_ASSERT_EQUAL_UINT64(P_new->x_mask, P_str->x_mask);
    TEST_ASSERT_EQUAL_UINT64(P_new->z_mask, P_str->z_mask);
    TEST_ASSERT_EQUAL_UINT64(P_new->y_mask, P_str->y_mask);

    unsigned nvecs;
    cplx_t **states = test_mk_states_pure(3, &nvecs);

    for (unsigned si = 0; si < nvecs; ++si) {
        cplx_t *cn = malloc(8 * sizeof(*cn));
        TEST_ASSERT_NOT_NULL(cn);
        memcpy(cn, states[si], 8 * sizeof(*cn));
        state_t sn = {0}; sn.type = PURE;
        state_init(&sn, 3, &cn);
        if (!sn.data) { free(cn); pauli_free(P_new); pauli_free(P_str); test_rm_states_pure(3, states); TEST_FAIL_MESSAGE("state_init failed (P_new)"); }
        pauli_apply_pure(&sn, P_new);

        cplx_t *cs = malloc(8 * sizeof(*cs));
        TEST_ASSERT_NOT_NULL(cs);
        memcpy(cs, states[si], 8 * sizeof(*cs));
        state_t ss = {0}; ss.type = PURE;
        state_init(&ss, 3, &cs);
        if (!ss.data) { free(cs); free(sn.data); pauli_free(P_new); pauli_free(P_str); test_rm_states_pure(3, states); TEST_FAIL_MESSAGE("state_init failed (P_str)"); }
        pauli_apply_pure(&ss, P_str);

        TEST_ASSERT_EQUAL_COMPLEX_ARRAY_TOL(sn.data, ss.data, 8, PRECISION);
        free(sn.data); sn.data = NULL;
        free(ss.data); ss.data = NULL;
    }

    test_rm_states_pure(3, states);
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
    struct { const char *str; unsigned n; } cases[] = {
        {"X", 1}, {"Y", 1}, {"Z", 1}, {"XX", 2}, {"ZYX", 3},
    };
    const unsigned ncases = sizeof(cases) / sizeof(cases[0]);

    for (unsigned c = 0; c < ncases; ++c) {
        const unsigned n = cases[c].n;
        const unsigned dim = 1u << n;
        unsigned nvecs;
        cplx_t **states = test_mk_states_pure(n, &nvecs);
        pauli_t *P = pauli_from_str(cases[c].str, (qubit_t)n);

        for (unsigned si = 0; si < nvecs; ++si) {
            cplx_t *copy = malloc(dim * sizeof(*copy));
            TEST_ASSERT_NOT_NULL(copy);
            memcpy(copy, states[si], dim * sizeof(*copy));
            state_t s = {0}; s.type = PURE;
            state_init(&s, (qubit_t)n, &copy);
            if (!s.data) { free(copy); pauli_free(P); test_rm_states_pure(n, states); TEST_FAIL_MESSAGE("state_init failed in involution test"); }
            pauli_apply_pure(&s, P);
            pauli_apply_pure(&s, P);
            TEST_ASSERT_EQUAL_COMPLEX_ARRAY_TOL(states[si], s.data, dim, PRECISION);
            free(s.data); s.data = NULL;
        }

        pauli_free(P);
        test_rm_states_pure(n, states);
    }
}

/*
 * =====================================================================================================================
 * Test: III on n=3 uniform superposition is identity
 *
 * I^(x)3 applied to |+++> must return |+++> unchanged.
 * =====================================================================================================================
 */


/*
 * =====================================================================================================================
 * Gap 1: pauli_new with n_qubits = 0
 *
 * pauli_new with n=0 allocates a struct with all masks zero and returns a valid
 * pointer.  pauli_free(P) must not crash.
 * =====================================================================================================================
 */

void test_pauli_new_n0(void) {
    pauli_t *P = pauli_new(NULL, 0);
    TEST_ASSERT_NOT_NULL(P);
    TEST_ASSERT_EQUAL_UINT(0u, (unsigned)P->n_qubits);
    TEST_ASSERT_EQUAL_UINT64(0u, P->x_mask);
    TEST_ASSERT_EQUAL_UINT64(0u, P->z_mask);
    TEST_ASSERT_EQUAL_UINT64(0u, P->y_mask);
    pauli_free(P);   /* must not crash */
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
    /* Z²=I verified on all n=1 states (subsumes the original 3-state check) */
    unsigned nvecs;
    cplx_t **states = test_mk_states_pure(1, &nvecs);
    pauli_t *P = pauli_from_str("Z", 1);

    for (unsigned si = 0; si < nvecs; ++si) {
        cplx_t *copy = malloc(2 * sizeof(*copy));
        TEST_ASSERT_NOT_NULL_MESSAGE(copy, "malloc failed in Z involution test");
        memcpy(copy, states[si], 2 * sizeof(*copy));
        state_t s = {0}; s.type = PURE;
        state_init(&s, 1, &copy);
        if (!s.data) { free(copy); pauli_free(P); test_rm_states_pure(1, states); TEST_FAIL_MESSAGE("state_init failed in Z involution test"); }
        pauli_apply_pure(&s, P);
        pauli_apply_pure(&s, P);
        TEST_ASSERT_EQUAL_COMPLEX_ARRAY_TOL(states[si], s.data, 2, PRECISION);
        free(s.data); s.data = NULL;
    }

    pauli_free(P);
    test_rm_states_pure(1, states);
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
 * Helper check_pauli_exp takes the Pauli string and builds both paths independently:
 *   - Reference:   pauli_str → mat_pauli_str (Kronecker) → P|psi> → ref
 *   - Under test:  pauli_str → pauli_from_str (bitmasks)  → pauli_exp_pure → result
 *
 * Convention reminder:
 *   - Identity Pauli (I^n): global phase e^{-i*theta/2} = CMPLX(c, -s) on all amplitudes.
 *   - theta = 0 : U = I, output must equal input.
 *   - theta = 2*pi : U = -I, all amplitudes negated.
 *   - theta = 4*pi : U = +I, output equals input.  (True periodicity is 4*pi.)
 * =====================================================================================================================
 */

static void check_pauli_exp(const char *pauli_str, unsigned n_qubits,
                             const cplx_t *psi, double theta) {
    const unsigned dim = (unsigned)(1u << n_qubits);

    cplx_t *mat       = NULL;
    cplx_t *pauli_psi = NULL;
    cplx_t *ref       = NULL;
    cplx_t *copy      = NULL;
    pauli_t *P        = NULL;
    state_t st        = {0};
    st.type            = PURE;

    const double c  = cos(theta * 0.5);
    const double sv = sin(theta * 0.5);

    /* Reference path: string → Kronecker matrix → P|psi> → ref */
    mat = mat_pauli_str(pauli_str, n_qubits);
    if (!mat) {
        fprintf(stderr, "check_pauli_exp(): mat_pauli_str failed\n");
        goto cleanup;
    }

    pauli_psi = zmv(dim, mat, psi);
    if (!pauli_psi) {
        fprintf(stderr, "check_pauli_exp(): zmv failed\n");
        goto cleanup;
    }

    ref = malloc((size_t)dim * sizeof(*ref));
    if (!ref) {
        fprintf(stderr, "check_pauli_exp(): ref allocation failed\n");
        goto cleanup;
    }
    for (unsigned k = 0; k < dim; ++k) {
        ref[k] = c * psi[k] + CMPLX(0.0, -sv) * pauli_psi[k];
    }

    /* Mask path: string → pauli_t → pauli_exp_pure */
    P = pauli_from_str(pauli_str, (qubit_t)n_qubits);

    copy = malloc((size_t)dim * sizeof(*copy));
    if (!copy) {
        fprintf(stderr, "check_pauli_exp(): copy allocation failed\n");
        goto cleanup;
    }
    memcpy(copy, psi, (size_t)dim * sizeof(*copy));

    state_init(&st, (qubit_t)n_qubits, &copy);
    if (!st.data) {
        fprintf(stderr, "check_pauli_exp(): state_init failed\n");
        goto cleanup;
    }

    pauli_exp_pure(&st, P, theta);
    if (!st.data) {
        fprintf(stderr, "check_pauli_exp(): pauli_exp_pure left state data NULL\n");
        goto cleanup;
    }

    TEST_ASSERT_EQUAL_COMPLEX_ARRAY_TOL(ref, st.data, dim, PRECISION);

    free(mat);
    free(pauli_psi);
    free(ref);
    pauli_free(P);
    free(st.data);
    st.data = NULL;
    return;

cleanup:
    free(mat);
    free(pauli_psi);
    free(ref);
    free(copy);
    pauli_free(P);
    state_free(&st);
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
    /* n=1, I — test all n=1 states */
    {
        unsigned nvecs; cplx_t **states = test_mk_states_pure(1, &nvecs);
        for (unsigned ti = 0; ti < N_EXP_THETAS; ++ti)
            for (unsigned si = 0; si < nvecs; ++si)
                check_pauli_exp("I", 1, states[si], TEST_EXP_THETAS[ti]);
        test_rm_states_pure(1, states);
    }
    /* n=2, II — test all n=2 states */
    {
        unsigned nvecs; cplx_t **states = test_mk_states_pure(2, &nvecs);
        for (unsigned ti = 0; ti < N_EXP_THETAS; ++ti)
            for (unsigned si = 0; si < nvecs; ++si)
                check_pauli_exp("II", 2, states[si], TEST_EXP_THETAS[ti]);
        test_rm_states_pure(2, states);
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
    const char *strs[] = {"Z", "ZZ", "ZZZ", "IZI"};
    const unsigned ns[] = {1, 2, 3, 3};
    for (unsigned ci = 0; ci < 4; ++ci) {
        unsigned nvecs; cplx_t **states = test_mk_states_pure(ns[ci], &nvecs);
        for (unsigned ti = 0; ti < N_EXP_THETAS; ++ti)
            for (unsigned si = 0; si < nvecs; ++si)
                check_pauli_exp(strs[ci], ns[ci], states[si], TEST_EXP_THETAS[ti]);
        test_rm_states_pure(ns[ci], states);
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
    unsigned nvecs; cplx_t **s = test_mk_states_pure(1, &nvecs);
    for (unsigned ti = 0; ti < N_EXP_THETAS; ++ti)
        for (unsigned si = 0; si < nvecs; ++si)
            check_pauli_exp("X", 1, s[si], TEST_EXP_THETAS[ti]);
    test_rm_states_pure(1, s);
}

void test_pauli_exp_pure_nondiaganal_Y(void) {
    unsigned nvecs; cplx_t **s = test_mk_states_pure(1, &nvecs);
    for (unsigned ti = 0; ti < N_EXP_THETAS; ++ti)
        for (unsigned si = 0; si < nvecs; ++si)
            check_pauli_exp("Y", 1, s[si], TEST_EXP_THETAS[ti]);
    test_rm_states_pure(1, s);
}

void test_pauli_exp_pure_nondiaganal_XX(void) {
    unsigned nvecs; cplx_t **s = test_mk_states_pure(2, &nvecs);
    for (unsigned ti = 0; ti < N_EXP_THETAS; ++ti)
        for (unsigned si = 0; si < nvecs; ++si)
            check_pauli_exp("XX", 2, s[si], TEST_EXP_THETAS[ti]);
    test_rm_states_pure(2, s);
}

void test_pauli_exp_pure_nondiaganal_YY(void) {
    unsigned nvecs; cplx_t **s = test_mk_states_pure(2, &nvecs);
    for (unsigned ti = 0; ti < N_EXP_THETAS; ++ti)
        for (unsigned si = 0; si < nvecs; ++si)
            check_pauli_exp("YY", 2, s[si], TEST_EXP_THETAS[ti]);
    test_rm_states_pure(2, s);
}

void test_pauli_exp_pure_nondiaganal_YYY(void) {
    unsigned nvecs; cplx_t **s = test_mk_states_pure(3, &nvecs);
    for (unsigned ti = 0; ti < N_EXP_THETAS; ++ti)
        for (unsigned si = 0; si < nvecs; ++si)
            check_pauli_exp("YYY", 3, s[si], TEST_EXP_THETAS[ti]);
    test_rm_states_pure(3, s);
}

void test_pauli_exp_pure_nondiaganal_XYZ(void) {
    /* labels[0]=X(q0),Y(q1),Z(q2) → big-endian "ZYX" */
    unsigned nvecs; cplx_t **s = test_mk_states_pure(3, &nvecs);
    for (unsigned ti = 0; ti < N_EXP_THETAS; ++ti)
        for (unsigned si = 0; si < nvecs; ++si)
            check_pauli_exp("ZYX", 3, s[si], TEST_EXP_THETAS[ti]);
    test_rm_states_pure(3, s);
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
    /* U(0) = I for all Paulis and all states */
    const char *strs[] = {"I", "Z", "X", "Y", "YYY"};
    const unsigned ns[] = {1, 1, 1, 1, 3};
    for (unsigned c = 0; c < 5; ++c) {
        const unsigned n = ns[c]; const unsigned dim = 1u << n;
        unsigned nvecs; cplx_t **states = test_mk_states_pure(n, &nvecs);
        pauli_t *P = pauli_from_str(strs[c], (qubit_t)n);
        for (unsigned si = 0; si < nvecs; ++si) {
            cplx_t *copy = malloc(dim * sizeof(*copy));
            TEST_ASSERT_NOT_NULL(copy);
            memcpy(copy, states[si], dim * sizeof(*copy));
            state_t st = {0}; st.type = PURE;
            state_init(&st, (qubit_t)n, &copy);
            if (!st.data) { free(copy); pauli_free(P); test_rm_states_pure(n, states); TEST_FAIL_MESSAGE("state_init failed"); }
            pauli_exp_pure(&st, P, 0.0);
            TEST_ASSERT_EQUAL_COMPLEX_ARRAY_TOL(states[si], st.data, dim, PRECISION);
            free(st.data); st.data = NULL;
        }
        pauli_free(P); test_rm_states_pure(n, states);
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
    /* U(pi)|psi> = -i*P|psi> for all states; tested for X (non-diagonal) and Z (diagonal) */
    const char *strs[] = {"X", "Z"};
    for (unsigned c = 0; c < 2; ++c) {
        unsigned nvecs; cplx_t **states = test_mk_states_pure(1, &nvecs);
        pauli_t *P = pauli_from_str(strs[c], 1);
        for (unsigned si = 0; si < nvecs; ++si) {
            const cplx_t *psi = states[si];
            /* Reference: -i * P|psi> */
            cplx_t *rc = malloc(2 * sizeof(*rc)); TEST_ASSERT_NOT_NULL(rc);
            memcpy(rc, psi, 2 * sizeof(*rc));
            state_t sr = {0}; sr.type = PURE; state_init(&sr, 1, &rc);
            if (!sr.data) { free(rc); pauli_free(P); test_rm_states_pure(1, states); TEST_FAIL_MESSAGE("ref init failed"); }
            pauli_apply_pure(&sr, P);
            cplx_t expected[2];
            for (unsigned k = 0; k < 2; ++k) expected[k] = CMPLX(0.0, -1.0) * sr.data[k];
            free(sr.data); sr.data = NULL;
            /* Under test */
            cplx_t *copy = malloc(2 * sizeof(*copy)); TEST_ASSERT_NOT_NULL(copy);
            memcpy(copy, psi, 2 * sizeof(*copy));
            state_t s = {0}; s.type = PURE; state_init(&s, 1, &copy);
            if (!s.data) { free(copy); pauli_free(P); test_rm_states_pure(1, states); TEST_FAIL_MESSAGE("test init failed"); }
            pauli_exp_pure(&s, P, M_PI);
            TEST_ASSERT_EQUAL_COMPLEX_ARRAY_TOL(expected, s.data, 2, PRECISION);
            free(s.data); s.data = NULL;
        }
        pauli_free(P); test_rm_states_pure(1, states);
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
    /* U(2pi) = -I for all Paulis and states: output[k] == -psi[k] */
    const char *strs[] = {"I", "Z", "X", "YYY"};
    const unsigned ns[] = {1, 1, 1, 3};
    for (unsigned c = 0; c < 4; ++c) {
        const unsigned n = ns[c]; const unsigned dim = 1u << n;
        unsigned nvecs; cplx_t **states = test_mk_states_pure(n, &nvecs);
        pauli_t *P = pauli_from_str(strs[c], (qubit_t)n);
        for (unsigned si = 0; si < nvecs; ++si) {
            cplx_t *expected = malloc(dim * sizeof(*expected)); TEST_ASSERT_NOT_NULL(expected);
            for (unsigned k = 0; k < dim; ++k) expected[k] = -states[si][k];
            cplx_t *copy = malloc(dim * sizeof(*copy)); TEST_ASSERT_NOT_NULL(copy);
            memcpy(copy, states[si], dim * sizeof(*copy));
            state_t s = {0}; s.type = PURE; state_init(&s, (qubit_t)n, &copy);
            if (!s.data) { free(copy); free(expected); pauli_free(P); test_rm_states_pure(n, states); TEST_FAIL_MESSAGE("state_init failed"); }
            pauli_exp_pure(&s, P, 2.0 * M_PI);
            TEST_ASSERT_EQUAL_COMPLEX_ARRAY_TOL(expected, s.data, dim, PRECISION);
            free(expected); free(s.data); s.data = NULL;
        }
        pauli_free(P); test_rm_states_pure(n, states);
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
    const double angles[5] = {M_PI/4.0, M_PI/2.0, M_PI, 3.0*M_PI/2.0, -M_PI/3.0};
    const char *strs[] = {"I", "Z", "ZZ", "X", "ZYX"};
    const unsigned ns[] = {1, 1, 2, 1, 3};
    for (unsigned c = 0; c < 5; ++c) {
        const unsigned n = ns[c]; const unsigned dim = 1u << n;
        unsigned nvecs; cplx_t **states = test_mk_states_pure(n, &nvecs);
        pauli_t *P = pauli_from_str(strs[c], (qubit_t)n);
        for (unsigned ai = 0; ai < 5; ++ai) {
            for (unsigned si = 0; si < nvecs; ++si) {
                cplx_t *copy = malloc(dim * sizeof(*copy)); TEST_ASSERT_NOT_NULL(copy);
                memcpy(copy, states[si], dim * sizeof(*copy));
                state_t s = {0}; s.type = PURE; state_init(&s, (qubit_t)n, &copy);
                TEST_ASSERT_NOT_NULL(s.data);
                pauli_exp_pure(&s, P, angles[ai]);
                TEST_ASSERT_DOUBLE_WITHIN(PRECISION, 1.0, compute_norm_sq(s.data, dim));
                free(s.data); s.data = NULL;
            }
        }
        pauli_free(P); test_rm_states_pure(n, states);
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
    /* U(theta)*U(-theta) = I for all Paulis and all states */
    const double fwd_angles[4] = {M_PI/4.0, M_PI/2.0, M_PI, -M_PI/3.0};
    const char *strs[] = {"I", "Z", "X", "YYY", "ZYX"};
    const unsigned ns[] = {1, 1, 1, 3, 3};
    for (unsigned c = 0; c < 5; ++c) {
        const unsigned n = ns[c]; const unsigned dim = 1u << n;
        unsigned nvecs; cplx_t **states = test_mk_states_pure(n, &nvecs);
        pauli_t *P = pauli_from_str(strs[c], (qubit_t)n);
        for (unsigned ai = 0; ai < 4; ++ai) {
            const double theta = fwd_angles[ai];
            for (unsigned si = 0; si < nvecs; ++si) {
                cplx_t *copy = malloc(dim * sizeof(*copy)); TEST_ASSERT_NOT_NULL(copy);
                memcpy(copy, states[si], dim * sizeof(*copy));
                state_t s = {0}; s.type = PURE; state_init(&s, (qubit_t)n, &copy);
                TEST_ASSERT_NOT_NULL(s.data);
                pauli_exp_pure(&s, P, theta);
                pauli_exp_pure(&s, P, -theta);
                TEST_ASSERT_EQUAL_COMPLEX_ARRAY_TOL(states[si], s.data, dim, PRECISION);
                free(s.data); s.data = NULL;
            }
        }
        pauli_free(P); test_rm_states_pure(n, states);
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

    /* Eigenstates constructed inline — no global PSI arrays needed */
    const cplx_t e0[2] = {1.0+0.0*I, 0.0+0.0*I};             /* |0> */
    const cplx_t e1[2] = {0.0+0.0*I, 1.0+0.0*I};             /* |1> */
    const cplx_t ep[2] = {INVSQRT2+0.0*I,  INVSQRT2+0.0*I};  /* |+> */
    const cplx_t em[2] = {INVSQRT2+0.0*I, -INVSQRT2+0.0*I};  /* |-> */
    const cplx_t e00[4] = {1,0,0,0}; const cplx_t e01[4] = {0,1,0,0};
    const cplx_t e10[4] = {0,0,1,0}; const cplx_t e11[4] = {0,0,0,1};

    /* Helper macro: apply exp(P, theta) to psi, check == fac*psi */
#define CHECK_EIGEN(P, psi, dim, fac) do { \
    cplx_t expected_[dim]; \
    for (unsigned k_ = 0; k_ < (dim); ++k_) expected_[k_] = (fac) * (psi)[k_]; \
    cplx_t *copy_ = malloc((dim) * sizeof(*copy_)); TEST_ASSERT_NOT_NULL(copy_); \
    memcpy(copy_, (psi), (dim) * sizeof(*copy_)); \
    state_t s_ = {0}; s_.type = PURE; state_init(&s_, (qubit_t)__builtin_ctz(dim), &copy_); \
    TEST_ASSERT_NOT_NULL(s_.data); \
    pauli_exp_pure(&s_, (P), theta); \
    TEST_ASSERT_EQUAL_COMPLEX_ARRAY_TOL(expected_, s_.data, (dim), PRECISION); \
    free(s_.data); s_.data = NULL; \
} while(0)

    { pauli_t *P = pauli_from_str("Z", 1);
      CHECK_EIGEN(P, e0, 2, fac_pos); CHECK_EIGEN(P, e1, 2, fac_neg); pauli_free(P); }
    { pauli_t *P = pauli_from_str("X", 1);
      CHECK_EIGEN(P, ep, 2, fac_pos); CHECK_EIGEN(P, em, 2, fac_neg); pauli_free(P); }
    { pauli_t *P = pauli_from_str("ZZ", 2);
      CHECK_EIGEN(P, e00, 4, fac_pos); CHECK_EIGEN(P, e01, 4, fac_neg);
      CHECK_EIGEN(P, e10, 4, fac_neg); CHECK_EIGEN(P, e11, 4, fac_pos); pauli_free(P); }
#undef CHECK_EIGEN
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
    /* U(4pi) = +I for all Paulis and states (contrast: 2pi gives -I) */
    const char *strs[] = {"I", "Z", "X", "YYY"};
    const unsigned ns[] = {1, 1, 1, 3};
    for (unsigned c = 0; c < 4; ++c) {
        const unsigned n = ns[c]; const unsigned dim = 1u << n;
        unsigned nvecs; cplx_t **states = test_mk_states_pure(n, &nvecs);
        pauli_t *P = pauli_from_str(strs[c], (qubit_t)n);
        for (unsigned si = 0; si < nvecs; ++si) {
            cplx_t *copy = malloc(dim * sizeof(*copy)); TEST_ASSERT_NOT_NULL(copy);
            memcpy(copy, states[si], dim * sizeof(*copy));
            state_t s = {0}; s.type = PURE; state_init(&s, (qubit_t)n, &copy);
            TEST_ASSERT_NOT_NULL(s.data);
            pauli_exp_pure(&s, P, 4.0 * M_PI);
            TEST_ASSERT_EQUAL_COMPLEX_ARRAY_TOL(states[si], s.data, dim, PRECISION);
            free(s.data); s.data = NULL;
        }
        pauli_free(P); test_rm_states_pure(n, states);
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
     * Reference: apply P via pauli_apply_pure, scale by -i, compare against pauli_exp_pure.
     */

    /* exp(pi) = -i*P: test all states for each Pauli covering all branches */
    const char *strs[] = {"I", "Z", "X", "Y", "XX", "YY", "YYY", "ZYX"};
    const unsigned ns[] = {1, 1, 1, 1, 2, 2, 3, 3};
    const unsigned ncases = sizeof(strs) / sizeof(strs[0]);

    for (unsigned ci = 0; ci < ncases; ++ci) {
        const unsigned n = ns[ci]; const unsigned dim = 1u << n;
        unsigned nvecs; cplx_t **states = test_mk_states_pure(n, &nvecs);
        pauli_t *P = pauli_from_str(strs[ci], (qubit_t)n);

        for (unsigned si = 0; si < nvecs; ++si) {
            const cplx_t *psi = states[si];
            /* Reference: -i * P|psi> */
            cplx_t *rc = malloc(dim * sizeof(*rc)); TEST_ASSERT_NOT_NULL(rc);
            memcpy(rc, psi, dim * sizeof(*rc));
            state_t sr = {0}; sr.type = PURE; state_init(&sr, (qubit_t)n, &rc);
            if (!sr.data) { free(rc); pauli_free(P); test_rm_states_pure(n, states); TEST_FAIL_MESSAGE("ref init failed"); }
            pauli_apply_pure(&sr, P);
            cplx_t *expected = malloc(dim * sizeof(*expected)); TEST_ASSERT_NOT_NULL(expected);
            for (unsigned k = 0; k < dim; ++k) expected[k] = CMPLX(0.0, -1.0) * sr.data[k];
            free(sr.data); sr.data = NULL;
            /* Under test */
            cplx_t *copy = malloc(dim * sizeof(*copy)); TEST_ASSERT_NOT_NULL(copy);
            memcpy(copy, psi, dim * sizeof(*copy));
            state_t s = {0}; s.type = PURE; state_init(&s, (qubit_t)n, &copy);
            if (!s.data) { free(copy); free(expected); pauli_free(P); test_rm_states_pure(n, states); TEST_FAIL_MESSAGE("test init failed"); }
            pauli_exp_pure(&s, P, M_PI);
            TEST_ASSERT_EQUAL_COMPLEX_ARRAY_TOL(expected, s.data, dim, PRECISION);
            free(expected); free(s.data); s.data = NULL;
        }
        pauli_free(P); test_rm_states_pure(n, states);
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
    /* labels[0]=X(q0),Z(q1),Y(q2),X(q3) → big-endian "XYZX"; x_mask=13, y_mask=4 */
    const double angles[4] = {0.0, M_PI / 2.0, M_PI, 2.0 * M_PI};
    unsigned nvecs; cplx_t **s = test_mk_states_pure(4, &nvecs);
    for (unsigned ai = 0; ai < 4; ++ai)
        for (unsigned si = 0; si < nvecs; ++si)
            check_pauli_exp("XYZX", 4, s[si], angles[ai]);
    test_rm_states_pure(4, s);
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
    const char *strs[] = {"XI", "IXI", "XII", "IYI"};
    const unsigned ns[] = {2, 3, 3, 3};
    for (unsigned ci = 0; ci < 4; ++ci) {
        unsigned nvecs; cplx_t **s = test_mk_states_pure(ns[ci], &nvecs);
        for (unsigned ti = 0; ti < N_EXP_THETAS; ++ti)
            for (unsigned si = 0; si < nvecs; ++si)
                check_pauli_exp(strs[ci], ns[ci], s[si], TEST_EXP_THETAS[ti]);
        test_rm_states_pure(ns[ci], s);
    }
}
