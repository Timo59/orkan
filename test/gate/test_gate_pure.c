// test_gate_pure.c  - Test harness for quantum gates on pure states. The reference is given by matrix-vector
//                     multiplication of the state vector with the matrix representation of the quantum gate.

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

#include <math.h>
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

/*
 * =====================================================================================================================
 * Rotation gate matrix builders
 * =====================================================================================================================
 */

void mat_rx(double theta, cplx_t *mat) {
    double c = cos(theta / 2.0);
    double s = sin(theta / 2.0);
    // [[cos(θ/2), -i·sin(θ/2)], [-i·sin(θ/2), cos(θ/2)]]
    // Column-major: [0,0], [1,0], [0,1], [1,1]
    mat[0] = c;
    mat[1] = -s * I;
    mat[2] = -s * I;
    mat[3] = c;
}

void mat_ry(double theta, cplx_t *mat) {
    double c = cos(theta / 2.0);
    double s = sin(theta / 2.0);
    // [[cos(θ/2), -sin(θ/2)], [sin(θ/2), cos(θ/2)]]
    // Column-major: [0,0], [1,0], [0,1], [1,1]
    mat[0] = c;
    mat[1] = s;
    mat[2] = -s;
    mat[3] = c;
}

void mat_rz(double theta, cplx_t *mat) {
    double c = cos(theta / 2.0);
    double s = sin(theta / 2.0);
    // [[e^(-iθ/2), 0], [0, e^(iθ/2)]]
    // Column-major: [0,0], [1,0], [0,1], [1,1]
    mat[0] = c - s * I;
    mat[1] = 0;
    mat[2] = 0;
    mat[3] = c + s * I;
}

void mat_p(double theta, cplx_t *mat) {
    // [[1, 0], [0, e^(iθ)]]
    // Column-major: [0,0], [1,0], [0,1], [1,1]
    mat[0] = 1.0;
    mat[1] = 0.0;
    mat[2] = 0.0;
    mat[3] = cos(theta) + sin(theta) * I;
}

/*
 * =====================================================================================================================
 * Test: Rotation gates
 * =====================================================================================================================
 */

/*
 * =====================================================================================================================
 * Test: Two-qubit gates
 * =====================================================================================================================
 */

void testTwoQubitGate(const two_qubit_gate gate, const cplx_t *mat) {
    state_t test_state = {0};
    test_state.type = PURE;
    cplx_t **test_vecs = NULL, *gateMat = NULL;
    unsigned nvecs = 0;

    // Iterate the number of qubits (minimum 2 for two-qubit gates)
    for (unsigned nqubits = 2; nqubits <= MAXQUBITS; ++nqubits) {
        const unsigned dim = POW2(nqubits, dim_t);

        // Iterate all (q1, q2) pairs where q1 != q2
        for (unsigned q1 = 0; q1 < nqubits; ++q1) {
            for (unsigned q2 = 0; q2 < nqubits; ++q2) {
                if (q1 == q2) continue;

                // Generate test state vectors
                if (!((test_vecs = test_mk_states_pure(nqubits, &nvecs)))) {
                    fprintf(stderr, "testTwoQubitGate(): test_vecs initialization failed\n");
                    goto cleanup;
                }

                // Initialize matrix representation of the gate acting on q1 and q2
                if (!((gateMat = mat_two_qubit_gate(nqubits, mat, q1, q2)))) goto cleanup;

                // Iterate the test state vectors
                for (unsigned i = 0; i < nvecs; ++i) {
                    // Reference state vector by matrix-vector multiplication
                    cplx_t* ref = zmv(dim, gateMat, test_vecs[i]);
                    if (!ref) {
                        fprintf(stderr, "testTwoQubitGate(): matrix-vector multiplication failed\n");
                        goto cleanup;
                    }

                    // Initialize test state with the state vector
                    state_init(&test_state, nqubits, &test_vecs[i]);
                    if (!test_state.data) {
                        fprintf(stderr, "testTwoQubitGate(): test state data initialization failed\n");
                        goto cleanup;
                    }

                    // Apply test function
                    gate(&test_state, q1, q2);
                    if (!test_state.data) {
                        fprintf(stderr, "testTwoQubitGate(): gate application failed\n");
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

// Test angles: representative values covering identity, small, π/2, π, and negative
static const double TEST_THETAS[] = {0.0, M_PI/4, M_PI/2, M_PI, 3*M_PI/2, -M_PI/3};
static const unsigned NUM_THETAS = sizeof(TEST_THETAS) / sizeof(TEST_THETAS[0]);

void testRotationGate(const rotation_gate gate, void (*mat_fn)(double theta, cplx_t *mat)) {
    state_t test_state = {0};
    test_state.type = PURE;
    cplx_t **test_vecs = NULL, *gateMat = NULL;
    cplx_t mat2x2[4];
    unsigned nvecs = 0;

    // Iterate test angles
    for (unsigned t = 0; t < NUM_THETAS; ++t) {
        double theta = TEST_THETAS[t];

        // Build 2x2 matrix for this angle
        mat_fn(theta, mat2x2);

        // Iterate the number of qubits
        for (unsigned nqubits = 1; nqubits <= MAXQUBITS; ++nqubits) {
            const unsigned dim = POW2(nqubits, dim_t);

            // Iterate target qubits
            for (unsigned pos = 0; pos < nqubits; ++pos) {
                // Generate test state vectors
                if (!((test_vecs = test_mk_states_pure(nqubits, &nvecs)))) {
                    fprintf(stderr, "testRotationGate(): test_vecs initialization failed\n");
                    goto cleanup;
                }

                // Initialize matrix representation of the gate acting on the target qubit
                if (!((gateMat = mat_single_qubit_gate(nqubits, mat2x2, pos)))) goto cleanup;

                // Iterate the test state vectors
                for (unsigned i = 0; i < nvecs; ++i) {
                    // Reference state vector by matrix-vector multiplication
                    cplx_t* ref = zmv(dim, gateMat, test_vecs[i]);
                    if (!ref) {
                        fprintf(stderr, "testRotationGate(): matrix-vector multiplication failed\n");
                        goto cleanup;
                    }

                    // Initialize test state with the state vector
                    state_init(&test_state, nqubits, &test_vecs[i]);
                    if (!test_state.data) {
                        fprintf(stderr, "testRotationGate(): test state data initialization failed\n");
                        goto cleanup;
                    }

                    // Apply test function
                    gate(&test_state, pos, theta);
                    if (!test_state.data) {
                        fprintf(stderr, "testRotationGate(): gate application failed\n");
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
