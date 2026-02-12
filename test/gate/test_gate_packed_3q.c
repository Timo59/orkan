// test_gate_packed_3q.c  - Test harness for three-qubit quantum gates on packed mixed states. The reference is given by
//                         matrix-matrix multiplication of the density matrix with the matrix representation of the quantum gate.

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
 * Test: Three-qubit gates
 * =====================================================================================================================
 */

void testThreeQubitGateMixed(const three_qubit_gate gate, const cplx_t *mat) {
    state_t test_state = {0};
    test_state.type = MIXED_PACKED;
    cplx_t **test_rhos = NULL, *gateMat = NULL;
    cplx_t *rho_full = NULL, *ref_full = NULL;
    unsigned nmixed = 0;

    // Iterate the number of qubits (minimum 3 for three-qubit gates)
    for (unsigned nqubits = 3; nqubits <= MAXQUBITS; ++nqubits) {
        const unsigned dim = POW2(nqubits, dim_t);

        // Iterate all (q1, q2, q3) triples where all three are distinct
        for (unsigned q1 = 0; q1 < nqubits; ++q1) {
            for (unsigned q2 = 0; q2 < nqubits; ++q2) {
                if (q1 == q2) continue;

                for (unsigned q3 = 0; q3 < nqubits; ++q3) {
                    if (q1 == q3 || q2 == q3) continue;

                    // Generate test density matrices
                    if (!((test_rhos = test_mk_states_mixed(nqubits, &nmixed)))) {
                        fprintf(stderr, "testThreeQubitGateMixed(): test_rhos initialization failed\n");
                        goto cleanup;
                    }

                    // Initialize matrix representation of the gate acting on q1, q2, and q3
                    if (!((gateMat = mat_three_qubit_gate(nqubits, mat, q1, q2, q3)))) goto cleanup;

                    // Iterate the test density matrices
                    for (unsigned i = 0; i < nmixed; ++i) {
                        // Unpack density matrix from packed storage to full matrix
                        rho_full = density_unpack(dim, test_rhos[i]);
                        if (!rho_full) {
                            fprintf(stderr, "testThreeQubitGateMixed(): density_unpack failed\n");
                            goto cleanup;
                        }

                        // Compute reference: ref_full = U * rho * U†
                        ref_full = zumu(dim, gateMat, rho_full);
                        if (!ref_full) {
                            fprintf(stderr, "testThreeQubitGateMixed(): zumu failed\n");
                            goto cleanup;
                        }

                        // Pack reference result to packed storage
                        cplx_t *ref = density_pack(dim, ref_full);
                        if (!ref) {
                            fprintf(stderr, "testThreeQubitGateMixed(): density_pack failed\n");
                            goto cleanup;
                        }

                        // Initialize test state with the density matrix
                        state_init(&test_state, nqubits, &test_rhos[i]);
                        if (!test_state.data) {
                            fprintf(stderr, "testThreeQubitGateMixed(): test state data initialization failed\n");
                            goto cleanup;
                        }

                        // Apply test function
                        gate(&test_state, q1, q2, q3);
                        if (!test_state.data) {
                            fprintf(stderr, "testThreeQubitGateMixed(): gate application failed\n");
                            goto cleanup;
                        }

                        // Compare density matrices in packed format
                        const unsigned packed_len = dim * (dim + 1) / 2;
                        TEST_ASSERT_EQUAL_COMPLEX_ARRAY_TOL(ref, test_state.data, packed_len, PRECISION);

                        // Free temporary matrices
                        free(ref);
                        ref = NULL;
                        free(rho_full);
                        rho_full = NULL;
                        free(ref_full);
                        ref_full = NULL;

                        // Free test state data
                        free(test_state.data);
                        test_state.data = NULL;
                    }

                    // Free matrix representation
                    free(gateMat);
                    gateMat = NULL;

                    // Free test density matrix array
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
                test_rhos[i] = NULL;
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
