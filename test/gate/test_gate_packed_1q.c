// test_gate_packed_1q.c  - Test harness for single-qubit quantum gates on packed mixed states. The reference is given by
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

#include <math.h>
#include <stdlib.h>

/*
 * =====================================================================================================================
 * Test: Single qubit gates
 * =====================================================================================================================
 */

void testSingleQubitGateMixed(const single_qubit_gate gate, const cplx_t *mat) {
    state_t test_state = {0}; // Mixed test state
    test_state.type = MIXED_PACKED;
    cplx_t **test_rhos = NULL, *gateMat = NULL;  // Test density matrices, matrix representation of the quantum gate
    cplx_t *rho_full = NULL, *ref_full = NULL;  // Temporary matrices for reference computation
    unsigned nmixed = 0;    // Number of test density matrices

    // Iterate the number of qubits (start from 1 to cover 1-qubit mixed states)
    for (unsigned nqubits = 1; nqubits <= MAXQUBITS; ++nqubits) {
        const unsigned dim = POW2(nqubits, dim_t);

        // Iterate target qubits
        for (unsigned pos = 0; pos < nqubits; ++pos) {
            // Generate test density matrices
            if (!((test_rhos = test_mk_states_mixed(nqubits, &nmixed)))) {
                fprintf(stderr, "testSingleQubitGateMixed(): test_rhos initialization failed\n");
                goto cleanup;
            }

            // Initialize matrix representation of the gate acting on the target qubit
            if (!((gateMat = mat_single_qubit_gate(nqubits, mat, pos)))) goto cleanup;

            // Iterate the test density matrices
            for (unsigned i = 0; i < nmixed; ++i) {
                // Unpack density matrix from packed storage to full matrix
                rho_full = density_unpack(dim, test_rhos[i]);
                if (!rho_full) {
                    fprintf(stderr, "testSingleQubitGateMixed(): density_unpack failed\n");
                    goto cleanup;
                }

                // Compute reference: ref_full = U * rho * U†
                ref_full = zumu(dim, gateMat, rho_full);
                if (!ref_full) {
                    fprintf(stderr, "testSingleQubitGateMixed(): zumu failed\n");
                    goto cleanup;
                }

                // Pack reference result to packed storage
                cplx_t *ref = density_pack(dim, ref_full);
                if (!ref) {
                    fprintf(stderr, "testSingleQubitGateMixed(): density_pack failed\n");
                    goto cleanup;
                }

                // Initialize test state with the density matrix (previous density matrix gets freed inside state_init())
                state_init(&test_state, nqubits, &test_rhos[i]);
                if (!test_state.data) {
                    fprintf(stderr, "testSingleQubitGateMixed(): test state data initialization failed\n");
                    goto cleanup;
                }

                // Apply test function
                gate(&test_state, pos);
                if (!test_state.data) {
                    fprintf(stderr, "testSingleQubitGateMixed(): gate application failed\n");
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

    return;

    cleanup:
        // Free test density matrices
        if (test_rhos) {
            for (unsigned i = 0; i < nmixed; ++i) {
                free(test_rhos[i]);
                test_rhos[i] = NULL;
            }
        }
        free(test_rhos);
        test_rhos = NULL;

        // Free gate matrix representation
        free(gateMat);
        gateMat = NULL;

        // Free temporary matrices
        free(rho_full);
        rho_full = NULL;
        free(ref_full);
        ref_full = NULL;

        // Free test state
        state_free(&test_state);
}

/*
 * =====================================================================================================================
 * Test: Rotation gates (mixed states)
 * =====================================================================================================================
 */

void testRotationGateMixed(const rotation_gate gate, void (*mat_fn)(double theta, cplx_t *mat)) {
    state_t test_state = {0};
    test_state.type = MIXED_PACKED;
    cplx_t **test_rhos = NULL, *gateMat = NULL;
    cplx_t *rho_full = NULL, *ref_full = NULL;
    cplx_t mat2x2[4];
    unsigned nmixed = 0;

    // Iterate test angles
    for (unsigned t = 0; t < NUM_THETAS; ++t) {
        double theta = TEST_THETAS[t];

        // Build 2x2 matrix for this angle
        mat_fn(theta, mat2x2);

        // Iterate the number of qubits (start from 1 to cover 1-qubit mixed states)
        for (unsigned nqubits = 1; nqubits <= MAXQUBITS; ++nqubits) {
            const unsigned dim = POW2(nqubits, dim_t);

            // Iterate target qubits
            for (unsigned pos = 0; pos < nqubits; ++pos) {
                // Generate test density matrices
                if (!((test_rhos = test_mk_states_mixed(nqubits, &nmixed)))) {
                    fprintf(stderr, "testRotationGateMixed(): test_rhos initialization failed\n");
                    goto cleanup;
                }

                // Initialize matrix representation of the gate acting on the target qubit
                if (!((gateMat = mat_single_qubit_gate(nqubits, mat2x2, pos)))) goto cleanup;

                // Iterate the test density matrices
                for (unsigned i = 0; i < nmixed; ++i) {
                    // Unpack density matrix from packed storage to full matrix
                    rho_full = density_unpack(dim, test_rhos[i]);
                    if (!rho_full) {
                        fprintf(stderr, "testRotationGateMixed(): density_unpack failed\n");
                        goto cleanup;
                    }

                    // Compute reference: ref_full = U * rho * U†
                    ref_full = zumu(dim, gateMat, rho_full);
                    if (!ref_full) {
                        fprintf(stderr, "testRotationGateMixed(): zumu failed\n");
                        goto cleanup;
                    }

                    // Pack reference result to packed storage
                    cplx_t *ref = density_pack(dim, ref_full);
                    if (!ref) {
                        fprintf(stderr, "testRotationGateMixed(): density_pack failed\n");
                        goto cleanup;
                    }

                    // Initialize test state with the density matrix
                    state_init(&test_state, nqubits, &test_rhos[i]);
                    if (!test_state.data) {
                        fprintf(stderr, "testRotationGateMixed(): test state data initialization failed\n");
                        goto cleanup;
                    }

                    // Apply test function
                    gate(&test_state, pos, theta);
                    if (!test_state.data) {
                        fprintf(stderr, "testRotationGateMixed(): gate application failed\n");
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
