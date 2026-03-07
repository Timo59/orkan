// test_gate_tiled_3q.c  - Test harness for three-qubit quantum gates on tiled mixed states. The reference is given by
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

#include <stdlib.h>

/*
 * =====================================================================================================================
 * Test: Three-qubit gates (tiled)
 * =====================================================================================================================
 */

void testThreeQubitGateTiled(const three_qubit_gate gate, const cplx_t *mat) {
    state_t test_state = {0};
    test_state.type = MIXED_TILED;
    cplx_t **test_rhos = NULL, *gateMat = NULL;
    cplx_t *rho_full = NULL, *ref_full = NULL;
    unsigned nmixed = 0;

    for (unsigned nqubits = 3; nqubits <= MAXQUBITS_TILED; ++nqubits) {
        const unsigned dim = POW2(nqubits, unsigned);

        for (unsigned q1 = 0; q1 < nqubits; ++q1) {
            for (unsigned q2 = 0; q2 < nqubits; ++q2) {
                if (q1 == q2) continue;

                for (unsigned q3 = 0; q3 < nqubits; ++q3) {
                    if (q1 == q3 || q2 == q3) continue;

                    if (!((test_rhos = test_mk_states_mixed(nqubits, &nmixed)))) {
                        fprintf(stderr, "testThreeQubitGateTiled(): test_rhos initialization failed\n");
                        goto cleanup;
                    }

                    if (!((gateMat = mat_three_qubit_gate(nqubits, mat, q1, q2, q3)))) goto cleanup;

                    for (unsigned i = 0; i < nmixed; ++i) {
                        rho_full = density_unpack(dim, test_rhos[i]);
                        if (!rho_full) {
                            fprintf(stderr, "testThreeQubitGateTiled(): density_unpack failed\n");
                            goto cleanup;
                        }

                        ref_full = zumu(dim, gateMat, rho_full);
                        if (!ref_full) {
                            fprintf(stderr, "testThreeQubitGateTiled(): zumu failed\n");
                            goto cleanup;
                        }

                        tiled_state_from_full(&test_state, nqubits, rho_full);
                        if (!test_state.data) {
                            fprintf(stderr, "testThreeQubitGateTiled(): tiled state creation failed\n");
                            goto cleanup;
                        }

                        gate(&test_state, q1, q2, q3);
                        if (!test_state.data) {
                            fprintf(stderr, "testThreeQubitGateTiled(): gate application failed\n");
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
