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
