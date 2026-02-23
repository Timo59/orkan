// test_gate_tiled_2q_mat.c  - Test harness for two_from_mat on tiled mixed states.
//                             The reference is given by similarity transformation UρU†
//                             computed on the full density matrix.
//
// Branch coverage targets (with LOG_TILE_DIM=3):
//   Branch A (hi < 3)  : n=2 (q1=0,q2=1); n=3 (q1=0,q2=2)
//   Branch B (lo<3<=hi): n=4 (q1=0,q2=3), (q1=2,q2=3)
//   Branch C (lo >= 3) : n=5 (q1=3,q2=4); n=6 (q1=3,q2=5), (q1=4,q2=5)
// All cases are covered automatically by the full (q1,q2) sweep over all nqubits in [2,6].

/*
 * =====================================================================================================================
 * includes
 * =====================================================================================================================
 */

#ifndef TEST_GATE_H
#include "test_gate.h"
#endif

#ifndef GATEMAT_H
#include "gatemat.h"
#endif

#ifndef LINALG_H
#include "linalg.h"
#endif

#ifndef TEST_MIXED_STATES_H
#include "test_mixed_states.h"
#endif

#ifndef UTILS_H
#include "utils.h"
#endif

#ifndef TEST_MIXED_UTILS_H
#include "test_mixed_utils.h"
#endif

#include <stdlib.h>

/* Forward declaration: two_from_mat is an internal tiled-gate function
 * not yet declared in gate_tiled.h (not on the test include path). */
void two_from_mat(state_t *state, qubit_t q1, qubit_t q2, const cplx_t *mat);

/*
 * =====================================================================================================================
 * Additional two-qubit gate matrices (column-major, hi⊗lo convention)
 * Column-major, hi⊗lo convention: mat[row + col*4] where row = 2*hi_bit + lo_bit.
 * Defined non-static so test_gate.c can reference them via the extern declarations
 * in test_gate.h.
 * =====================================================================================================================
 */

/* iSWAP: |00>->|00>, |01>->i|10>, |10>->i|01>, |11>->|11>
 * Column-major, hi⊗lo convention: mat[row + col*4] where row = 2*hi_bit + lo_bit.
 * col 0: [1,0,0,0]^T  col 1: [0,0,i,0]^T  col 2: [0,i,0,0]^T  col 3: [0,0,0,1]^T */
const cplx_t ISWAPMAT[16] = {
    1.0+0.0*I, 0.0+0.0*I, 0.0+0.0*I, 0.0+0.0*I,  /* col 0: |00>->{1}|00>  */
    0.0+0.0*I, 0.0+0.0*I, 0.0+1.0*I, 0.0+0.0*I,  /* col 1: |01>->{i}|10>  */
    0.0+0.0*I, 0.0+1.0*I, 0.0+0.0*I, 0.0+0.0*I,  /* col 2: |10>->{i}|01>  */
    0.0+0.0*I, 0.0+0.0*I, 0.0+0.0*I, 1.0+0.0*I   /* col 3: |11>->{1}|11>  */
};

/* H⊗H: Hadamard on hi tensor Hadamard on lo.
 * Column-major, hi⊗lo convention: mat[row + col*4] where row = 2*hi_bit + lo_bit.
 * H⊗H = (1/2)*[[1,1,1,1],[1,-1,1,-1],[1,1,-1,-1],[1,-1,-1,1]] */
const cplx_t HHMAT[16] = {
     0.5+0.0*I,  0.5+0.0*I,  0.5+0.0*I,  0.5+0.0*I,  /* col 0 */
     0.5+0.0*I, -0.5+0.0*I,  0.5+0.0*I, -0.5+0.0*I,  /* col 1 */
     0.5+0.0*I,  0.5+0.0*I, -0.5+0.0*I, -0.5+0.0*I,  /* col 2 */
     0.5+0.0*I, -0.5+0.0*I, -0.5+0.0*I,  0.5+0.0*I   /* col 3 */
};

/*
 * =====================================================================================================================
 * Test: two_from_mat (single application, tiled)
 *
 * For every nqubits in [2, MAXQUBITS_TILED] and every ordered pair (q1,q2) with q1 != q2:
 *   - Compute hi=max(q1,q2), lo=min(q1,q2)
 *   - Build the N-qubit reference matrix via mat_two_qubit_gate(nqubits, mat4x4, hi, lo)
 *     (hi as q1, lo as q2: matches two_from_mat's hi⊗lo convention)
 *   - Compute reference: ref_full = U * rho * U†  via zumu
 *   - Build a tiled state from rho_full, apply two_from_mat(&state, q1, q2, mat4x4)
 *   - Compare element-by-element with assert_tiled_equals_full
 *
 * Both orderings (q1<q2) and (q1>q2) are visited, so the internal sort in
 * two_from_mat is exercised in both directions for every gate.
 * =====================================================================================================================
 */

void testTwoFromMat(const cplx_t *mat4x4) {
    state_t test_state = {0};
    test_state.type = MIXED_TILED;
    cplx_t **test_rhos = NULL, *gateMat = NULL;
    cplx_t *rho_full = NULL, *ref_full = NULL;
    unsigned nmixed = 0;

    for (unsigned nqubits = 2; nqubits <= MAXQUBITS_TILED; ++nqubits) {
        const unsigned dim = POW2(nqubits, dim_t);

        for (unsigned q1 = 0; q1 < nqubits; ++q1) {
            for (unsigned q2 = 0; q2 < nqubits; ++q2) {
                if (q1 == q2) continue;

                const unsigned hi = (q1 > q2) ? q1 : q2;
                const unsigned lo = (q1 < q2) ? q1 : q2;

                if (!((test_rhos = test_mk_states_mixed(nqubits, &nmixed)))) {
                    fprintf(stderr, "testTwoFromMat(): test_rhos initialization failed\n");
                    goto cleanup;
                }

                /* Build N-qubit gate matrix in hi⊗lo convention */
                if (!((gateMat = mat_two_qubit_gate(nqubits, mat4x4, hi, lo)))) goto cleanup;

                for (unsigned i = 0; i < nmixed; ++i) {
                    rho_full = density_unpack(dim, test_rhos[i]);
                    if (!rho_full) {
                        fprintf(stderr, "testTwoFromMat(): density_unpack failed\n");
                        goto cleanup;
                    }

                    /* Reference: U * rho * U† */
                    ref_full = zumu(dim, gateMat, rho_full);
                    if (!ref_full) {
                        fprintf(stderr, "testTwoFromMat(): zumu failed\n");
                        goto cleanup;
                    }

                    tiled_state_from_full(&test_state, nqubits, rho_full);
                    if (!test_state.data) {
                        fprintf(stderr, "testTwoFromMat(): tiled state creation failed\n");
                        goto cleanup;
                    }

                    two_from_mat(&test_state, q1, q2, mat4x4);

                    assert_tiled_equals_full(&test_state, ref_full, dim);

                    free(rho_full);
                    rho_full = NULL;
                    free(ref_full);
                    ref_full = NULL;

                    state_free(&test_state);
                }

                free(gateMat);
                gateMat = NULL;

                for (unsigned j = 0; j < nmixed; ++j) {
                    free(test_rhos[j]);
                }
                free(test_rhos);
                test_rhos = NULL;
            }
        }
    }

    return;

    cleanup:
        if (test_rhos) {
            for (unsigned i = 0; i < nmixed; ++i) {
                free(test_rhos[i]);
            }
        }
        free(test_rhos);
        test_rhos = NULL;
        free(gateMat);
        gateMat = NULL;
        free(rho_full);
        rho_full = NULL;
        free(ref_full);
        ref_full = NULL;
        state_free(&test_state);
}

/*
 * =====================================================================================================================
 * Test: two_from_mat (double sequential application of the same gate, tiled)
 *
 * Applies mat4x4 twice via two_from_mat on the same qubit pair (q1,q2) and compares
 * against U²ρ(U†)² = U(UρU†)U†.  This catches stale off-diagonal mirror bugs that a
 * single application can miss: if the first call writes an incorrect conjugate into a
 * mirror slot, the second call reads that stale value from the lower triangle's mirror
 * and propagates the error; the final assert_tiled_equals_full then detects it.
 *
 * Reference is computed by calling zumu twice with the same N-qubit gate matrix.
 * =====================================================================================================================
 */

void testTwoFromMatDouble(const cplx_t *mat4x4) {
    state_t test_state = {0};
    test_state.type = MIXED_TILED;
    cplx_t **test_rhos = NULL, *gateMat = NULL;
    cplx_t *rho_full = NULL, *ref1 = NULL, *ref2 = NULL;
    unsigned nmixed = 0;

    for (unsigned nqubits = 2; nqubits <= MAXQUBITS_TILED; ++nqubits) {
        const unsigned dim = POW2(nqubits, dim_t);

        for (unsigned q1 = 0; q1 < nqubits; ++q1) {
            for (unsigned q2 = 0; q2 < nqubits; ++q2) {
                if (q1 == q2) continue;

                const unsigned hi = (q1 > q2) ? q1 : q2;
                const unsigned lo = (q1 < q2) ? q1 : q2;

                if (!((test_rhos = test_mk_states_mixed(nqubits, &nmixed)))) {
                    fprintf(stderr, "testTwoFromMatDouble(): test_rhos initialization failed\n");
                    goto cleanup_double;
                }

                /* Build N-qubit gate matrix in hi⊗lo convention (same gate applied twice) */
                if (!((gateMat = mat_two_qubit_gate(nqubits, mat4x4, hi, lo)))) goto cleanup_double;

                for (unsigned i = 0; i < nmixed; ++i) {
                    rho_full = density_unpack(dim, test_rhos[i]);
                    if (!rho_full) {
                        fprintf(stderr, "testTwoFromMatDouble(): density_unpack failed\n");
                        goto cleanup_double;
                    }

                    /* Reference: U (U rho U†) U† */
                    ref1 = zumu(dim, gateMat, rho_full);
                    if (!ref1) {
                        fprintf(stderr, "testTwoFromMatDouble(): zumu (first application) failed\n");
                        goto cleanup_double;
                    }

                    ref2 = zumu(dim, gateMat, ref1);
                    if (!ref2) {
                        fprintf(stderr, "testTwoFromMatDouble(): zumu (second application) failed\n");
                        goto cleanup_double;
                    }

                    free(ref1);
                    ref1 = NULL;

                    tiled_state_from_full(&test_state, nqubits, rho_full);
                    if (!test_state.data) {
                        fprintf(stderr, "testTwoFromMatDouble(): tiled state creation failed\n");
                        goto cleanup_double;
                    }

                    two_from_mat(&test_state, q1, q2, mat4x4);
                    two_from_mat(&test_state, q1, q2, mat4x4);

                    assert_tiled_equals_full(&test_state, ref2, dim);

                    free(rho_full);
                    rho_full = NULL;
                    free(ref2);
                    ref2 = NULL;

                    state_free(&test_state);
                }

                free(gateMat);
                gateMat = NULL;

                for (unsigned j = 0; j < nmixed; ++j) {
                    free(test_rhos[j]);
                }
                free(test_rhos);
                test_rhos = NULL;
            }
        }
    }

    return;

    cleanup_double:
        if (test_rhos) {
            for (unsigned i = 0; i < nmixed; ++i) {
                free(test_rhos[i]);
            }
        }
        free(test_rhos);
        test_rhos = NULL;
        free(gateMat);
        gateMat = NULL;
        free(rho_full);
        rho_full = NULL;
        free(ref1);
        ref1 = NULL;
        free(ref2);
        ref2 = NULL;
        state_free(&test_state);
}

/*
 * =====================================================================================================================
 * Test: two_from_mat argument-order symmetry (tiled)
 *
 * Verifies that calling two_from_mat(state, q1, q2, mat) and two_from_mat(state, q2, q1, mat)
 * on identical initial states produces the same result.  This tests that the internal sort
 * (hi = max(q1,q2), lo = min(q1,q2)) is correctly applied regardless of argument order.
 *
 * For each unordered pair {q1, q2} (loop q1 < q2 to visit each pair once):
 *   - Both orderings are compared against the same zumu reference
 *   - state_a: initial rho_full, apply two_from_mat(q1, q2, mat), compare with ref_full
 *   - state_b: same initial rho_full, apply two_from_mat(q2, q1, mat), compare with ref_full
 *   - If the sort is correct, both comparisons pass; if not, at least one fails.
 * =====================================================================================================================
 */

void testTwoFromMatSymmetry(const cplx_t *mat4x4) {
    state_t state_a = {0}, state_b = {0};
    state_a.type = MIXED_TILED;
    state_b.type = MIXED_TILED;
    cplx_t **test_rhos = NULL, *gateMat = NULL;
    cplx_t *rho_full = NULL, *ref_full = NULL;
    unsigned nmixed = 0;

    for (unsigned nqubits = 2; nqubits <= MAXQUBITS_TILED; ++nqubits) {
        const unsigned dim = POW2(nqubits, dim_t);

        /* q1 < q2: visit each unordered pair exactly once */
        for (unsigned q1 = 0; q1 < nqubits; ++q1) {
            for (unsigned q2 = q1 + 1; q2 < nqubits; ++q2) {
                /* hi=q2, lo=q1 since q1 < q2 */
                if (!((test_rhos = test_mk_states_mixed(nqubits, &nmixed)))) {
                    fprintf(stderr, "testTwoFromMatSymmetry(): test_rhos initialization failed\n");
                    goto cleanup_sym;
                }

                /* Build N-qubit reference matrix in hi⊗lo convention */
                if (!((gateMat = mat_two_qubit_gate(nqubits, mat4x4, q2, q1)))) goto cleanup_sym;

                for (unsigned i = 0; i < nmixed; ++i) {
                    rho_full = density_unpack(dim, test_rhos[i]);
                    if (!rho_full) {
                        fprintf(stderr, "testTwoFromMatSymmetry(): density_unpack failed\n");
                        goto cleanup_sym;
                    }

                    /* Reference: U * rho * U† */
                    ref_full = zumu(dim, gateMat, rho_full);
                    if (!ref_full) {
                        fprintf(stderr, "testTwoFromMatSymmetry(): zumu failed\n");
                        goto cleanup_sym;
                    }

                    /* state_a: apply two_from_mat(q1, q2) */
                    tiled_state_from_full(&state_a, nqubits, rho_full);
                    if (!state_a.data) {
                        fprintf(stderr, "testTwoFromMatSymmetry(): tiled state_a creation failed\n");
                        goto cleanup_sym;
                    }

                    two_from_mat(&state_a, q1, q2, mat4x4);
                    assert_tiled_equals_full(&state_a, ref_full, dim);

                    /* state_b: same initial state, apply two_from_mat(q2, q1) */
                    tiled_state_from_full(&state_b, nqubits, rho_full);
                    if (!state_b.data) {
                        fprintf(stderr, "testTwoFromMatSymmetry(): tiled state_b creation failed\n");
                        goto cleanup_sym;
                    }

                    two_from_mat(&state_b, q2, q1, mat4x4);
                    assert_tiled_equals_full(&state_b, ref_full, dim);

                    free(rho_full);
                    rho_full = NULL;
                    free(ref_full);
                    ref_full = NULL;

                    state_free(&state_a);
                    state_free(&state_b);
                }

                free(gateMat);
                gateMat = NULL;

                for (unsigned j = 0; j < nmixed; ++j) {
                    free(test_rhos[j]);
                }
                free(test_rhos);
                test_rhos = NULL;
            }
        }
    }

    return;

    cleanup_sym:
        if (test_rhos) {
            for (unsigned i = 0; i < nmixed; ++i) {
                free(test_rhos[i]);
            }
        }
        free(test_rhos);
        test_rhos = NULL;
        free(gateMat);
        gateMat = NULL;
        free(rho_full);
        rho_full = NULL;
        free(ref_full);
        ref_full = NULL;
        state_free(&state_a);
        state_free(&state_b);
}
