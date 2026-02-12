// test_gate_tiled.c  - Test harness for quantum gates on tiled mixed states. The reference is given by
//                      similarity transformation UρU† computed on the full density matrix.

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

#include <math.h>
#include <stdlib.h>

/*
 * Maximum qubits for tiled tests. Set to 6 (dim=64) to exercise cross-tile
 * code paths, which require dim > TILE_DIM (32 with LOG_TILE_DIM=5).
 * For 6-qubit mixed states the full density matrix is 64x64 = 4096 elements,
 * and the reference computation (U rho U^dag) runs in O(dim^3) which is
 * fast enough for testing.
 */
#define MAXQUBITS_TILED 6

/*
 * =====================================================================================================================
 * Tiled storage utilities
 * =====================================================================================================================
 */

/*
 * @brief   Unpack a packed lower-triangular density matrix to full matrix
 */
static cplx_t* density_unpack(const unsigned n, const cplx_t *packed) {
    cplx_t *full = NULL;

    if (!((full = malloc(n * n * sizeof(*full))))) {
        fprintf(stderr, "density_unpack(): full allocation failed\n");
        return full;
    }

    unsigned packed_idx = 0;
    for (unsigned col = 0; col < n; ++col) {
        for (unsigned row = col; row < n; ++row) {
            full[row + col * n] = packed[packed_idx];
            if (row != col) {
                full[col + row * n] = conj(packed[packed_idx]);
            }
            ++packed_idx;
        }
    }

    return full;
}

/*
 * @brief   Create a MIXED_TILED state from a full density matrix
 *
 * Allocates a tiled state and populates it element by element using state_set().
 */
static void tiled_state_from_full(state_t *state, unsigned nqubits, const cplx_t *full) {
    const unsigned dim = POW2(nqubits, dim_t);

    state->type = MIXED_TILED;
    state_init(state, nqubits, NULL);  // Allocate zero-initialized tiled storage

    if (!state->data) return;

    // Populate via state_set -- only need lower triangle (Hermitian)
    for (unsigned col = 0; col < dim; ++col) {
        for (unsigned row = col; row < dim; ++row) {
            state_set(state, row, col, full[row + col * dim]);
        }
    }
}

/*
 * @brief   Compare a tiled state against a full reference matrix element by element
 */
static void assert_tiled_equals_full(const state_t *state, const cplx_t *ref_full, unsigned dim) {
    for (unsigned col = 0; col < dim; ++col) {
        for (unsigned row = 0; row < dim; ++row) {
            cplx_t got = state_get(state, row, col);
            cplx_t exp = ref_full[row + col * dim];
            double re_diff = creal(got) - creal(exp);
            double im_diff = cimag(got) - cimag(exp);
            if (re_diff > PRECISION || re_diff < -PRECISION ||
                im_diff > PRECISION || im_diff < -PRECISION) {
                char msg[256];
                snprintf(msg, sizeof(msg),
                    "Tiled mismatch at (%u,%u): got (%.15e, %.15e) expected (%.15e, %.15e)",
                    row, col, creal(got), cimag(got), creal(exp), cimag(exp));
                TEST_FAIL_MESSAGE(msg);
            }
        }
    }
}

extern cplx_t* mat_two_qubit_gate(unsigned nqubits, const cplx_t *gate, unsigned q1, unsigned q2);

/*
 * =====================================================================================================================
 * Test: Single qubit gates (tiled)
 * =====================================================================================================================
 */

void testSingleQubitGateTiled(const single_qubit_gate gate, const cplx_t *mat) {
    state_t test_state = {0};
    test_state.type = MIXED_TILED;
    cplx_t **test_rhos = NULL, *gateMat = NULL;
    cplx_t *rho_full = NULL, *ref_full = NULL;
    unsigned nmixed = 0;

    for (unsigned nqubits = 1; nqubits <= MAXQUBITS_TILED; ++nqubits) {
        const unsigned dim = POW2(nqubits, dim_t);

        for (unsigned pos = 0; pos < nqubits; ++pos) {
            if (!((test_rhos = test_mk_states_mixed(nqubits, &nmixed)))) {
                fprintf(stderr, "testSingleQubitGateTiled(): test_rhos initialization failed\n");
                goto cleanup;
            }

            if (!((gateMat = mat_single_qubit_gate(nqubits, mat, pos)))) goto cleanup;

            for (unsigned i = 0; i < nmixed; ++i) {
                // Unpack packed density matrix to full
                rho_full = density_unpack(dim, test_rhos[i]);
                if (!rho_full) {
                    fprintf(stderr, "testSingleQubitGateTiled(): density_unpack failed\n");
                    goto cleanup;
                }

                // Compute reference: ref_full = U * rho * U†
                ref_full = zumu(dim, gateMat, rho_full);
                if (!ref_full) {
                    fprintf(stderr, "testSingleQubitGateTiled(): zumu failed\n");
                    goto cleanup;
                }

                // Create tiled state from full density matrix
                tiled_state_from_full(&test_state, nqubits, rho_full);
                if (!test_state.data) {
                    fprintf(stderr, "testSingleQubitGateTiled(): tiled state creation failed\n");
                    goto cleanup;
                }

                // Apply gate
                gate(&test_state, pos);
                if (!test_state.data) {
                    fprintf(stderr, "testSingleQubitGateTiled(): gate application failed\n");
                    goto cleanup;
                }

                // Compare element by element
                assert_tiled_equals_full(&test_state, ref_full, dim);

                free(rho_full);
                rho_full = NULL;
                free(ref_full);
                ref_full = NULL;

                state_free(&test_state);
            }

            free(gateMat);
            gateMat = NULL;

            // Free test density matrix array (packed format data)
            for (unsigned j = 0; j < nmixed; ++j) {
                free(test_rhos[j]);
            }
            free(test_rhos);
            test_rhos = NULL;
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
 * Test: Rotation gates (tiled)
 * =====================================================================================================================
 */

void testRotationGateTiled(const rotation_gate gate, void (*mat_fn)(double theta, cplx_t *mat)) {
    state_t test_state = {0};
    test_state.type = MIXED_TILED;
    cplx_t **test_rhos = NULL, *gateMat = NULL;
    cplx_t *rho_full = NULL, *ref_full = NULL;
    cplx_t mat2x2[4];
    unsigned nmixed = 0;

    for (unsigned t = 0; t < NUM_THETAS; ++t) {
        double theta = TEST_THETAS[t];

        mat_fn(theta, mat2x2);

        for (unsigned nqubits = 1; nqubits <= MAXQUBITS_TILED; ++nqubits) {
            const unsigned dim = POW2(nqubits, dim_t);

            for (unsigned pos = 0; pos < nqubits; ++pos) {
                if (!((test_rhos = test_mk_states_mixed(nqubits, &nmixed)))) {
                    fprintf(stderr, "testRotationGateTiled(): test_rhos initialization failed\n");
                    goto cleanup;
                }

                if (!((gateMat = mat_single_qubit_gate(nqubits, mat2x2, pos)))) goto cleanup;

                for (unsigned i = 0; i < nmixed; ++i) {
                    rho_full = density_unpack(dim, test_rhos[i]);
                    if (!rho_full) {
                        fprintf(stderr, "testRotationGateTiled(): density_unpack failed\n");
                        goto cleanup;
                    }

                    ref_full = zumu(dim, gateMat, rho_full);
                    if (!ref_full) {
                        fprintf(stderr, "testRotationGateTiled(): zumu failed\n");
                        goto cleanup;
                    }

                    tiled_state_from_full(&test_state, nqubits, rho_full);
                    if (!test_state.data) {
                        fprintf(stderr, "testRotationGateTiled(): tiled state creation failed\n");
                        goto cleanup;
                    }

                    gate(&test_state, pos, theta);
                    if (!test_state.data) {
                        fprintf(stderr, "testRotationGateTiled(): gate application failed\n");
                        goto cleanup;
                    }

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
 * Test: Two-qubit gates (tiled)
 * =====================================================================================================================
 */

void testTwoQubitGateTiled(const two_qubit_gate gate, const cplx_t *mat) {
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

                if (!((test_rhos = test_mk_states_mixed(nqubits, &nmixed)))) {
                    fprintf(stderr, "testTwoQubitGateTiled(): test_rhos initialization failed\n");
                    goto cleanup;
                }

                if (!((gateMat = mat_two_qubit_gate(nqubits, mat, q1, q2)))) goto cleanup;

                for (unsigned i = 0; i < nmixed; ++i) {
                    rho_full = density_unpack(dim, test_rhos[i]);
                    if (!rho_full) {
                        fprintf(stderr, "testTwoQubitGateTiled(): density_unpack failed\n");
                        goto cleanup;
                    }

                    ref_full = zumu(dim, gateMat, rho_full);
                    if (!ref_full) {
                        fprintf(stderr, "testTwoQubitGateTiled(): zumu failed\n");
                        goto cleanup;
                    }

                    tiled_state_from_full(&test_state, nqubits, rho_full);
                    if (!test_state.data) {
                        fprintf(stderr, "testTwoQubitGateTiled(): tiled state creation failed\n");
                        goto cleanup;
                    }

                    gate(&test_state, q1, q2);
                    if (!test_state.data) {
                        fprintf(stderr, "testTwoQubitGateTiled(): gate application failed\n");
                        goto cleanup;
                    }

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

