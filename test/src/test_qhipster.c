// test_qhipster.c  - Test harness for quantum gates on pure states. The reference is given by matrix-vector
//                    multiplication of the state vector with the matrix representation of the quantum gate.

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
 * Test: Single qubit gates
 * =====================================================================================================================
 */

void testSingleQubitGate(const single_qubit_gate gate, const cplx_t *mat) {
    state_t test_state = {0}; // Pure test state
    test_state.type = PURE;
    cplx_t **test_vecs = NULL, *gateMat = NULL;  // Test state vectors, matrix representation of the quantum gate
    unsigned nvecs = 0;    // Number of test state vectors

    // Iterate the number of qubits
    for (unsigned nqubits = 1; nqubits <= MAXQUBITS; ++nqubits) {
        const unsigned dim = POW2(nqubits, dim_t);

        // Iterate target qubits
        for (unsigned pos = 0; pos < nqubits; ++pos) {
            // Generate test state vectors
            if (!((test_vecs = test_mk_states_pure(nqubits, &nvecs)))) {
                fprintf(stderr, "testSingleQubitGate(): test_vecs initialization failed\n");
                goto cleanup;
            }

            // Initialize matrix representation of the gate acting on the target qubit
            if (!((gateMat = mat_single_qubit_gate(nqubits, mat, pos)))) goto cleanup;

            // Iterate the test state vectors
            for (unsigned i = 0; i < nvecs; ++i) {
                // Reference state vector by matrix-vector multiplication
                cplx_t* ref = zmv(dim, gateMat, test_vecs[i]);
                if (!ref) {
                    fprintf(stderr, "testSingleQubitGate(): matrix-vector multiplication failed\n");
                    goto cleanup;
                }

                // Initialize test state with the state vector (previous state vector gets freed inside state_init())
                state_init(&test_state, nqubits, &test_vecs[i]);
                if (!test_state.data) {
                    fprintf(stderr, "testSingleQubitGate(): test state data initialization failed\n");
                    goto cleanup;
                }

                // Apply test function
                gate(&test_state, pos);
                if (!test_state.data) {
                    fprintf(stderr, "testSingleQubitGate(): gate application failed\n");
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

    return;

    cleanup:
        // Free test state vectors
        if (test_vecs) {
            for (unsigned i = 0; i < nvecs; ++i) {
                free(test_vecs[i]);
                test_vecs[i] = NULL;
            }
        }
        free(test_vecs);
        test_vecs = NULL;

        // Free gate matrix representation
        free(gateMat);
        gateMat = NULL;

        // Free test state
        state_free(&test_state);
}
