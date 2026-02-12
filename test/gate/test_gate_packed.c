// test_gate_packed.c  - Test harness for quantum gates on packed mixed states. The reference is given by
//                       matrix-matrix multiplication of the density matrix with the matrix representation of the quantum gate.

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
 * =====================================================================================================================
 * Packed storage utilities
 * =====================================================================================================================
 */

static cplx_t* density_unpack(const unsigned n, const cplx_t *packed) {
    cplx_t *full = NULL;

    // Allocate memory for full n×n matrix
    if (!((full = malloc(n * n * sizeof(*full))))) {
        fprintf(stderr, "density_unpack(): full allocation failed\n");
        return full;
    }

    // Unpack lower triangle from packed storage (column-major)
    // Packed format: column 0 (n elements), column 1 (n-1 elements), ..., column n-1 (1 element)
    unsigned packed_idx = 0;
    for (unsigned col = 0; col < n; ++col) {
        for (unsigned row = col; row < n; ++row) {
            full[row + col * n] = packed[packed_idx];  // Lower triangle
            if (row != col) {
                // Hermitian property: H[col,row] = conj(H[row,col])
                full[col + row * n] = conj(packed[packed_idx]);  // Upper triangle
            }
            ++packed_idx;
        }
    }

    return full;
}


static cplx_t* density_pack(const unsigned n, const cplx_t *full) {
    cplx_t *packed = NULL;

    // Allocate memory for packed storage
    const unsigned packed_len = n * (n + 1) / 2;
    if (!((packed = malloc(packed_len * sizeof(*packed))))) {
        fprintf(stderr, "density_pack(): packed allocation failed\n");
        return packed;
    }

    // Pack lower triangle into packed storage (column-major)
    unsigned packed_idx = 0;
    for (unsigned col = 0; col < n; ++col) {
        for (unsigned row = col; row < n; ++row) {
            packed[packed_idx] = full[row + col * n];
            ++packed_idx;
        }
    }

    return packed;
}

/*
 * @brief   Returns the matrix representing the application of a three-qubit gate
 */
extern cplx_t* mat_three_qubit_gate(unsigned nqubits, const cplx_t *gate, unsigned q1, unsigned q2, unsigned q3);

/*
 * @brief   Returns the matrix representing the application of a two-qubit gate
 */
extern cplx_t* mat_two_qubit_gate(unsigned nqubits, const cplx_t *gate, unsigned q1, unsigned q2);

/*
 * =====================================================================================================================
 * Test: Single qubit gates
 * =====================================================================================================================
 */

void testTwoQubitGateMixed(const two_qubit_gate gate, const cplx_t *mat) {
    state_t test_state = {0};
    test_state.type = MIXED_PACKED;
    cplx_t **test_rhos = NULL, *gateMat = NULL;
    cplx_t *rho_full = NULL, *ref_full = NULL;
    unsigned nmixed = 0;

    // Iterate the number of qubits (minimum 2 for two-qubit gates)
    for (unsigned nqubits = 2; nqubits <= MAXQUBITS; ++nqubits) {
        const unsigned dim = POW2(nqubits, dim_t);

        // Iterate all (q1, q2) pairs where q1 != q2
        for (unsigned q1 = 0; q1 < nqubits; ++q1) {
            for (unsigned q2 = 0; q2 < nqubits; ++q2) {
                if (q1 == q2) continue;

                // Generate test density matrices
                if (!((test_rhos = test_mk_states_mixed(nqubits, &nmixed)))) {
                    fprintf(stderr, "testTwoQubitGateMixed(): test_rhos initialization failed\n");
                    goto cleanup;
                }

                // Initialize matrix representation of the gate acting on q1 and q2
                if (!((gateMat = mat_two_qubit_gate(nqubits, mat, q1, q2)))) goto cleanup;

                // Iterate the test density matrices
                for (unsigned i = 0; i < nmixed; ++i) {
                    // Unpack density matrix from packed storage to full matrix
                    rho_full = density_unpack(dim, test_rhos[i]);
                    if (!rho_full) {
                        fprintf(stderr, "testTwoQubitGateMixed(): density_unpack failed\n");
                        goto cleanup;
                    }

                    // Compute reference: ref_full = U * rho * U†
                    ref_full = zumu(dim, gateMat, rho_full);
                    if (!ref_full) {
                        fprintf(stderr, "testTwoQubitGateMixed(): zumu failed\n");
                        goto cleanup;
                    }

                    // Pack reference result to packed storage
                    cplx_t *ref = density_pack(dim, ref_full);
                    if (!ref) {
                        fprintf(stderr, "testTwoQubitGateMixed(): density_pack failed\n");
                        goto cleanup;
                    }

                    // Initialize test state with the density matrix
                    state_init(&test_state, nqubits, &test_rhos[i]);
                    if (!test_state.data) {
                        fprintf(stderr, "testTwoQubitGateMixed(): test state data initialization failed\n");
                        goto cleanup;
                    }

                    // Apply test function
                    gate(&test_state, q1, q2);
                    if (!test_state.data) {
                        fprintf(stderr, "testTwoQubitGateMixed(): gate application failed\n");
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