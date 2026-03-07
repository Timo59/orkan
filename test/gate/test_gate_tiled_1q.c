// test_gate_tiled_1q.c  - Test harness for single-qubit quantum gates on tiled mixed states. The reference is given by
//                         similarity transformation UρU† computed on the full density matrix.

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

#include <math.h>
#include <stdlib.h>

/* Forward declaration: single_from_mat is an internal tiled-gate function
 * declared in gate_tiled.h (not on the test include path). */
void single_from_mat(state_t *state, qubit_t target, const cplx_t *mat);

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
        const unsigned dim = POW2(nqubits, unsigned);

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
            const unsigned dim = POW2(nqubits, unsigned);

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
 * Test: single_from_mat (single application, tiled)
 * =====================================================================================================================
 */

void testSingleQubitMatGateTiled(const cplx_t *mat) {
    state_t test_state = {0};
    test_state.type = MIXED_TILED;
    cplx_t **test_rhos = NULL, *gateMat = NULL;
    cplx_t *rho_full = NULL, *ref_full = NULL;
    unsigned nmixed = 0;

    for (unsigned nqubits = 1; nqubits <= MAXQUBITS_TILED; ++nqubits) {
        const unsigned dim = POW2(nqubits, unsigned);

        for (unsigned pos = 0; pos < nqubits; ++pos) {
            if (!((test_rhos = test_mk_states_mixed(nqubits, &nmixed)))) {
                fprintf(stderr, "testSingleQubitMatGateTiled(): test_rhos initialization failed\n");
                goto cleanup_mat;
            }

            if (!((gateMat = mat_single_qubit_gate(nqubits, mat, pos)))) goto cleanup_mat;

            for (unsigned i = 0; i < nmixed; ++i) {
                rho_full = density_unpack(dim, test_rhos[i]);
                if (!rho_full) {
                    fprintf(stderr, "testSingleQubitMatGateTiled(): density_unpack failed\n");
                    goto cleanup_mat;
                }

                ref_full = zumu(dim, gateMat, rho_full);
                if (!ref_full) {
                    fprintf(stderr, "testSingleQubitMatGateTiled(): zumu failed\n");
                    goto cleanup_mat;
                }

                tiled_state_from_full(&test_state, nqubits, rho_full);
                if (!test_state.data) {
                    fprintf(stderr, "testSingleQubitMatGateTiled(): tiled state creation failed\n");
                    goto cleanup_mat;
                }

                single_from_mat(&test_state, pos, mat);

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

    return;

    cleanup_mat:
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
 * Test: single_from_mat (double sequential application, tiled)
 *
 * Applies mat1 then mat2 and compares against U2(U1 rho U1†)U2†.
 * Catches stale diagonal-tile mirrors that assert_tiled_equals_full cannot
 * detect via state_tiled_get (which reads from the lower triangle and conjugates).
 * =====================================================================================================================
 */

void testSingleQubitMatGateTiledDouble(const cplx_t *mat1, const cplx_t *mat2) {
    state_t test_state = {0};
    test_state.type = MIXED_TILED;
    cplx_t **test_rhos = NULL;
    cplx_t *gateMat1 = NULL, *gateMat2 = NULL;
    cplx_t *rho_full = NULL, *ref1 = NULL, *ref2 = NULL;
    unsigned nmixed = 0;

    for (unsigned nqubits = 1; nqubits <= MAXQUBITS_TILED; ++nqubits) {
        const unsigned dim = POW2(nqubits, unsigned);

        for (unsigned pos = 0; pos < nqubits; ++pos) {
            if (!((test_rhos = test_mk_states_mixed(nqubits, &nmixed)))) {
                fprintf(stderr, "testSingleQubitMatGateTiledDouble(): test_rhos initialization failed\n");
                goto cleanup_double;
            }

            if (!((gateMat1 = mat_single_qubit_gate(nqubits, mat1, pos)))) goto cleanup_double;
            if (!((gateMat2 = mat_single_qubit_gate(nqubits, mat2, pos)))) goto cleanup_double;

            for (unsigned i = 0; i < nmixed; ++i) {
                rho_full = density_unpack(dim, test_rhos[i]);
                if (!rho_full) {
                    fprintf(stderr, "testSingleQubitMatGateTiledDouble(): density_unpack failed\n");
                    goto cleanup_double;
                }

                /* Reference: U2 (U1 rho U1†) U2† */
                ref1 = zumu(dim, gateMat1, rho_full);
                if (!ref1) {
                    fprintf(stderr, "testSingleQubitMatGateTiledDouble(): zumu (mat1) failed\n");
                    goto cleanup_double;
                }

                ref2 = zumu(dim, gateMat2, ref1);
                if (!ref2) {
                    fprintf(stderr, "testSingleQubitMatGateTiledDouble(): zumu (mat2) failed\n");
                    goto cleanup_double;
                }

                free(ref1);
                ref1 = NULL;

                tiled_state_from_full(&test_state, nqubits, rho_full);
                if (!test_state.data) {
                    fprintf(stderr, "testSingleQubitMatGateTiledDouble(): tiled state creation failed\n");
                    goto cleanup_double;
                }

                single_from_mat(&test_state, pos, mat1);
                single_from_mat(&test_state, pos, mat2);

                assert_tiled_equals_full(&test_state, ref2, dim);

                free(rho_full);
                rho_full = NULL;
                free(ref2);
                ref2 = NULL;

                state_free(&test_state);
            }

            free(gateMat1);
            gateMat1 = NULL;
            free(gateMat2);
            gateMat2 = NULL;

            for (unsigned j = 0; j < nmixed; ++j) {
                free(test_rhos[j]);
            }
            free(test_rhos);
            test_rhos = NULL;
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
        free(gateMat1);
        gateMat1 = NULL;
        free(gateMat2);
        gateMat2 = NULL;
        free(rho_full);
        rho_full = NULL;
        free(ref1);
        ref1 = NULL;
        free(ref2);
        ref2 = NULL;
        state_free(&test_state);
}
