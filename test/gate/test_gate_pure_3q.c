// test_gate_pure_3q.c  - Test harness for three-qubit quantum gates on pure states. The reference is given by
//                        matrix-vector multiplication of the state vector with the matrix representation of the gate.

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

#ifndef TEST_PURE_STATES_H
#include "test_pure_states.h"
#endif

#ifndef UTILS_H
#include "utils.h"
#endif

#include <stdlib.h>

/*
 * =====================================================================================================================
 * Test: Three-qubit gates
 * =====================================================================================================================
 */

void testThreeQubitGate(const three_qubit_gate gate, const cplx_t *mat) {
    state_t test_state = {0};
    test_state.type = PURE;
    cplx_t **test_vecs = NULL, *gateMat = NULL;
    unsigned nvecs = 0;

    // Iterate the number of qubits (minimum 3 for three-qubit gates)
    for (unsigned nqubits = 3; nqubits <= MAXQUBITS; ++nqubits) {
        const unsigned dim = POW2(nqubits, unsigned);

        // Iterate all (q1, q2, q3) triples where all three are distinct
        for (unsigned q1 = 0; q1 < nqubits; ++q1) {
            for (unsigned q2 = 0; q2 < nqubits; ++q2) {
                if (q1 == q2) continue;

                for (unsigned q3 = 0; q3 < nqubits; ++q3) {
                    if (q1 == q3 || q2 == q3) continue;

                    // Generate test state vectors
                    if (!((test_vecs = test_mk_states_pure(nqubits, &nvecs)))) {
                        fprintf(stderr, "testThreeQubitGate(): test_vecs initialization failed\n");
                        goto cleanup;
                    }

                    // Initialize matrix representation of the gate acting on q1, q2, and q3
                    if (!((gateMat = mat_three_qubit_gate(nqubits, mat, q1, q2, q3)))) goto cleanup;

                    // Iterate the test state vectors
                    for (unsigned i = 0; i < nvecs; ++i) {
                        // Reference state vector by matrix-vector multiplication
                        cplx_t* ref = zmv(dim, gateMat, test_vecs[i]);
                        if (!ref) {
                            fprintf(stderr, "testThreeQubitGate(): matrix-vector multiplication failed\n");
                            goto cleanup;
                        }

                        // Initialize test state with the state vector
                        state_init(&test_state, nqubits, &test_vecs[i]);
                        if (!test_state.data) {
                            fprintf(stderr, "testThreeQubitGate(): test state data initialization failed\n");
                            goto cleanup;
                        }

                        // Apply test function
                        gate(&test_state, q1, q2, q3);
                        if (!test_state.data) {
                            fprintf(stderr, "testThreeQubitGate(): gate application failed\n");
                            goto cleanup;
                        }

                        // Compare state vectors
                        TEST_ASSERT_EQUAL_COMPLEX_ARRAY_TOL(ref, test_state.data, dim, PRECISION);

                        // Free the reference result
                        free(ref);
                        ref = NULL;

                        // Free test state vector
                        free(test_state.data);
                        test_state.data = NULL;
                    }

                    // Free matrix representation
                    free(gateMat);
                    gateMat = NULL;

                    // Free test vector array
                    free(test_vecs);
                    test_vecs = NULL;
                }
            }
        }
    }

    return;

    cleanup:
        if (test_vecs) {
            for (unsigned i = 0; i < nvecs; ++i) {
                free(test_vecs[i]);
                test_vecs[i] = NULL;
            }
        }
        free(test_vecs);
        test_vecs = NULL;
        free(gateMat);
        gateMat = NULL;
        state_free(&test_state);
}
