/*
 * =====================================================================================================================
 *                                                      includes
 * =====================================================================================================================
 */
#ifndef GATELIB_H
#include "qhipster.h"
#endif

#ifndef HELPFUNC_H
#include "utils.h"
#endif

#ifndef GATEMAT_H
#include "gateMat.h"
#endif

#ifndef UNITY_FRAMEWORK_H
#include "unity.h"
#endif

/*
 * =====================================================================================================================
 *                                                      macros
 * =====================================================================================================================
 */
#define PRECISION               1e-8
#define MAXQUBITS               6

/*
 * =====================================================================================================================
 *                                                  setUp and tearDown
 * =====================================================================================================================
 */

void setUp(void) {}

void tearDown(void) {}

/*
 * =====================================================================================================================
 *                                              test applyX
 * =====================================================================================================================
 */
/*
 * This function checks whether the execution of the Pauli-X gate by applyX coincides with matrix multiplication
 */
void testApplyX(void) {
    for (qubit_t qubits = 1; qubits <= MAXQUBITS; ++qubits) {       // For each number of qubits up to MAXQUBITS:
        dim_t dim = POW2(qubits, dim_t);                            // Hilbert space dimension
        state_t testState;                                          // Initialize a pure quantum state with specified
        stateInitEmpty(&testState, qubits);                         // number of qubits
        cplx_t** vecs = generateTestVectors(qubits);                // State vectors for matrix multiplication

        for (qubit_t pos = 0; pos < qubits; ++pos) {                    // For each target qubit from the right:
            cplx_t* gateMat = singleQubitGateMat(XMAT, qubits, pos);    // generate the matrix representation from the
                                                                        // Kronecker product
            for (dim_t i = 0; i < dim + 1; ++i) {                       // For each state vector:
                stateInitVector(&testState, vecs[i]);                   // Copy state vector to testState's vector,
                applyX(&testState, pos);                                // apply test function
                cplx_t* ref = cmatVecMul(gateMat, vecs[i], dim);    // Matrix multiplication as reference
                TEST_ASSERT_TRUE(cvectorAlmostEqual(ref, testState.vec, dim, PRECISION));
                free(ref);
            }
            free(gateMat);
        }
        stateFreeVector(&testState);
        freeTestVectors(vecs, qubits);
    }
}

/*
 * =====================================================================================================================
 *                                              test applyY
 * =====================================================================================================================
 */

/*
 * This function checks whether the execution of the Pauli-Y gate by applyY coincides with matrix multiplication
 */
void testApplyY(void) {
    for (qubit_t qubits = 1; qubits <= MAXQUBITS; ++qubits) {       // For each number of qubits up to MAXQUBITS:
        dim_t dim = POW2(qubits, dim_t);                            // Hilbert space dimension
        state_t testState;                                          // Initialize a pure quantum state with specified
        stateInitEmpty(&testState, qubits);                         // number of qubits
        cplx_t** vecs = generateTestVectors(qubits);                // State vectors for matrix multiplication

        for (qubit_t pos = 0; pos < qubits; ++pos) {                    // For each target qubit from the right:
            cplx_t* gateMat = singleQubitGateMat(YMAT, qubits, pos);    // generate the matrix representation from the
                                                                        // Kronecker product
            for (dim_t i = 0; i < dim + 1; ++i) {                       // For each state vector:
                stateInitVector(&testState, vecs[i]);                   // Copy state vector to testState's vector,
                applyY(&testState, pos);                                // apply test function
                cplx_t* ref = cmatVecMul(gateMat, vecs[i], dim);    // Matrix multiplication as reference
                TEST_ASSERT_TRUE(cvectorAlmostEqual(ref, testState.vec, dim, PRECISION));
                free(ref);
            }
            free(gateMat);
        }
        stateFreeVector(&testState);
        freeTestVectors(vecs, qubits);
    }
}

/*
 * =====================================================================================================================
 *                                              test applyZ
 * =====================================================================================================================
 */
/*
 * This function checks whether the execution of the Pauli-Z gate by applyZ coincides with matrix multiplication
 */
void testApplyZ(void) {
    for (qubit_t qubits = 1; qubits <= MAXQUBITS; ++qubits) {       // For each number of qubits up to MAXQUBITS:
        dim_t dim = POW2(qubits, dim_t);                            // Hilbert space dimension
        state_t testState;                                          // Initialize a pure quantum state with specified
        stateInitEmpty(&testState, qubits);                         // number of qubits
        cplx_t** vecs = generateTestVectors(qubits);                // State vectors for matrix multiplication

        for (qubit_t pos = 0; pos < qubits; ++pos) {                    // For each target qubit from the right:
            cplx_t* gateMat = singleQubitGateMat(ZMAT, qubits, pos);    // generate the matrix representation from the
            // Kronecker product
            for (dim_t i = 0; i < dim + 1; ++i) {                       // For each state vector:
                stateInitVector(&testState, vecs[i]);                   // Copy state vector to testState's vector,
                applyZ(&testState, pos);                                // apply test function
                cplx_t* ref = cmatVecMul(gateMat, vecs[i], dim);    // Matrix multiplication as reference
                TEST_ASSERT_TRUE(cvectorAlmostEqual(ref, testState.vec, dim, PRECISION));
                free(ref);
            }
            free(gateMat);
        }
        stateFreeVector(&testState);
        freeTestVectors(vecs, qubits);
    }
}

/*
 * =====================================================================================================================
 *                                              test applyS(dagger)
 * =====================================================================================================================
 */

/*
 * This function checks whether the execution of the S-gate by applyS and its inverse coincides with matrix
 * multiplication
 */
void testApplyS(void) {
    for (qubit_t qubits = 1; qubits <= MAXQUBITS; ++qubits) {       // For each number of qubits up to MAXQUBITS:
        dim_t dim = POW2(qubits, dim_t);                            // Hilbert space dimension
        state_t testState;                                          // Initialize a pure quantum state with specified
        stateInitEmpty(&testState, qubits);                         // number of qubits
        cplx_t** vecs = generateTestVectors(qubits);                // State vectors for matrix multiplication

        for (qubit_t pos = 0; pos < qubits; ++pos) {                    // For each target qubit from the right:
            cplx_t* gateMat = singleQubitGateMat(SMAT, qubits, pos);    // generate the matrix representation from the
            // Kronecker product
            for (dim_t i = 0; i < dim + 1; ++i) {                       // For each state vector:
                stateInitVector(&testState, vecs[i]);                   // Copy state vector to testState's vector,
                applyS(&testState, pos);                                // apply test function
                cplx_t* ref = cmatVecMul(gateMat, vecs[i], dim);    // Matrix multiplication as reference

                TEST_ASSERT_TRUE(cvectorAlmostEqual(ref, testState.vec, dim, PRECISION));
                free(ref);

                applySdagger(&testState, pos);                          // Apply the inverse of S
                TEST_ASSERT_TRUE(cvectorAlmostEqual(vecs[i], testState.vec, dim, PRECISION));
            }
            free(gateMat);
        }
        stateFreeVector(&testState);
        freeTestVectors(vecs, qubits);
    }
}

/*
 * =====================================================================================================================
 *                                              test applyH
 * =====================================================================================================================
 */
/*
 * This function checks whether the execution of the Pauli-Z gate by applyZ coincides with matrix multiplication
*/
void testApplyH(void) {
    for (qubit_t qubits = 1; qubits <= MAXQUBITS; ++qubits) {       // For each number of qubits up to MAXQUBITS:
        dim_t dim = POW2(qubits, dim_t);                            // Hilbert space dimension
        state_t testState;                                          // Initialize a pure quantum state with specified
        stateInitEmpty(&testState, qubits);                         // number of qubits
        cplx_t** vecs = generateTestVectors(qubits);                // State vectors for matrix multiplication

        for (qubit_t pos = 0; pos < qubits; ++pos) {                    // For each target qubit from the right:
            cplx_t* gateMat = singleQubitGateMat(HMAT, qubits, pos);    // generate the matrix representation from the
                                                                        // Kronecker product
            for (dim_t i = 0; i < dim + 1; ++i) {                       // For each state vector:
                stateInitVector(&testState, vecs[i]);                   // Copy state vector to testState's vector,
                applyH(&testState, pos);                                // apply test function
                cplx_t* ref = cmatVecMul(gateMat, vecs[i], dim);    // Matrix multiplication as reference
                TEST_ASSERT_TRUE(cvectorAlmostEqual(ref, testState.vec, dim, PRECISION));
                free(ref);
            }
            free(gateMat);
        }
        stateFreeVector(&testState);
        freeTestVectors(vecs, qubits);
    }
}

/*
 * =====================================================================================================================
 *                                              test applyT(dagger)
 * =====================================================================================================================
 */

/*
 * This function checks whether the execution of the T-gate by applyT and its inverse coincides with matrix
 * multiplication
 *
 * Input:
 *      -
 *      Output:
 *      There is no return value. Only checks whether the state vector after applyS coincides with the state vector
 *      after matrix multiplication for all target qubits and whether the subsequent application of T^dagger leads to
 *      the starting vector.
 */
void testApplyT(void) {
    for (qubit_t qubits = 1; qubits <= MAXQUBITS; ++qubits) {       // For each number of qubits up to MAXQUBITS:
        dim_t dim = POW2(qubits, dim_t);                            // Hilbert space dimension
        state_t testState;                                          // Initialize a pure quantum state with specified
        stateInitEmpty(&testState, qubits);                         // number of qubits
        cplx_t** vecs = generateTestVectors(qubits);                // State vectors for matrix multiplication

        for (qubit_t pos = 0; pos < qubits; ++pos) {                    // For each target qubit from the right:
            cplx_t* gateMat = singleQubitGateMat(TMAT, qubits, pos);    // generate the matrix representation from the
                                                                        // Kronecker product
            for (dim_t i = 0; i < dim + 1; ++i) {                       // For each state vector:
                stateInitVector(&testState, vecs[i]);                   // Copy state vector to testState's vector,
                applyT(&testState, pos);                                // apply test function
                cplx_t* ref = cmatVecMul(gateMat, vecs[i], dim);    // Matrix multiplication as reference

                TEST_ASSERT_TRUE(cvectorAlmostEqual(ref, testState.vec, dim, PRECISION));
                free(ref);

                applyTdagger(&testState, pos);                          // Apply the inverse of S
                TEST_ASSERT_TRUE(cvectorAlmostEqual(vecs[i], testState.vec, dim, PRECISION));
            }
            free(gateMat);
        }
        stateFreeVector(&testState);
        freeTestVectors(vecs, qubits);
    }
}

/*
 * =====================================================================================================================
 *                                              test applyCX
 * =====================================================================================================================
 */
/*
 * This function checks whether the execution of the Pauli-X gate on one qubit, controlled by another, coincides with
 * matrix multiplication
 */
void testApplyCX(void) {
    for (qubit_t qubits = 2; qubits <= MAXQUBITS; ++qubits) {
        dim_t dim = POW2(qubits, dim_t);                            // Hilbert space dimension
        state_t testState;                                          // Initialize a pure quantum state with specified
        stateInitEmpty(&testState, qubits);                         // number of qubits
        cplx_t** vecs = generateTestVectors(qubits);                // State vectors for matrix multiplication

        for (qubit_t control = 0; control < qubits; ++control) {        // For each qubit as a control qubit
            for (qubit_t target = 0; target < qubits; ++target) {       // For each qubit as a target qubit
                if (control == target) {                                // Skip if control and target coincide
                    continue;
                }
                cplx_t* gateMat = controlledGateMat(XMAT, qubits, control, target);     // Controlled operation's matrix
                                                                                        // representation
                for (dim_t i = 0; i < dim + 1; ++i) {                   // For each state vector
                    stateInitVector(&testState, vecs[i]);                   // Copy state vector to testState's vector,
                    applyCX(&testState, control, target);                   // apply test function
                    cplx_t* reference = cmatVecMul(gateMat, vecs[i], dim);  // Matrix multiplication as reference
                    TEST_ASSERT_TRUE(cvectorAlmostEqual(reference, testState.vec, dim, PRECISION));
                    free(reference);
                }

                free(gateMat);
            }
        }
        freeTestVectors(vecs, qubits);
        stateFreeVector(&testState);
    }
}

/*
 * =====================================================================================================================
 *                                              test applyCY
 * =====================================================================================================================
 */
/*
 * This function checks whether the execution of the Pauli-Y gate on one qubit, controlled by another, coincides with
 * matrix multiplication.
 */
void testApplyCY(void) {
    for (qubit_t qubits = 2; qubits <= MAXQUBITS; ++qubits) {
        dim_t dim = POW2(qubits, dim_t);                            // Hilbert space dimension
        state_t testState;                                          // Initialize a pure quantum state with specified
        stateInitEmpty(&testState, qubits);                         // number of qubits
        cplx_t** vecs = generateTestVectors(qubits);                // State vectors for matrix multiplication

        for (qubit_t control = 0; control < qubits; ++control) {        // For each qubit as a control qubit
            for (qubit_t target = 0; target < qubits; ++target) {       // For each qubit as a target qubit
                if (control == target) {                                // Skip if control and target coincide
                    continue;
                }
                cplx_t* gateMat = controlledGateMat(YMAT, qubits, control, target);     // Controlled operation's matrix
                                                                                        // representation
                for (dim_t i = 0; i < dim + 1; ++i) {                   // For each state vector
                    stateInitVector(&testState, vecs[i]);                   // Copy state vector to testState's vector,
                    applyCY(&testState, control, target);                   // apply test function
                    cplx_t* reference = cmatVecMul(gateMat, vecs[i], dim);  // Matrix multiplication as reference
                    TEST_ASSERT_TRUE(cvectorAlmostEqual(reference, testState.vec, dim, PRECISION));
                    free(reference);
                }

                free(gateMat);
            }
        }
        freeTestVectors(vecs, qubits);
        stateFreeVector(&testState);
    }
}

/*
 * =====================================================================================================================
 *                                              test applyCZ
 * =====================================================================================================================
 */
/*
 * This function checks whether the execution of the Pauli-Z gate on one qubit, controlled by another, coincides with
 * matrix multiplication.
 */
void testApplyCZ(void) {
    for (qubit_t qubits = 2; qubits <= MAXQUBITS; ++qubits) {
        dim_t dim = POW2(qubits, dim_t);                            // Hilbert space dimension
        state_t testState;                                          // Initialize a pure quantum state with specified
        stateInitEmpty(&testState, qubits);                         // number of qubits
        cplx_t** vecs = generateTestVectors(qubits);                // State vectors for matrix multiplication

        for (qubit_t control = 0; control < qubits; ++control) {        // For each qubit as a control qubit
            for (qubit_t target = 0; target < qubits; ++target) {       // For each qubit as a target qubit
                if (control == target) {                                // Skip if control and target coincide
                    continue;
                }
                cplx_t* gateMat = controlledGateMat(ZMAT, qubits, control, target);     // Controlled operation's matrix
                                                                                        // representation
                for (dim_t i = 0; i < dim + 1; ++i) {                   // For each state vector
                    stateInitVector(&testState, vecs[i]);                   // Copy state vector to testState's vector,
                    applyCZ(&testState, control, target);                   // apply test function
                    cplx_t* reference = cmatVecMul(gateMat, vecs[i], dim);  // Matrix multiplication as reference
                    TEST_ASSERT_TRUE(cvectorAlmostEqual(reference, testState.vec, dim, PRECISION));
                    free(reference);
                }

                free(gateMat);
            }
        }
        freeTestVectors(vecs, qubits);
        stateFreeVector(&testState);
    }
}

/*
 * =====================================================================================================================
 *                                              test applyCS
 * =====================================================================================================================
 */
/*
 * This function checks whether the execution of the S-gate on one qubit, controlled by another and its inverse,
 * coincides with matrix multiplication.
 */
void testApplyCS(void) {
    for (qubit_t qubits = 2; qubits <= MAXQUBITS; ++qubits) {
        dim_t dim = POW2(qubits, dim_t);                            // Hilbert space dimension
        state_t testState;                                          // Initialize a pure quantum state with specified
        stateInitEmpty(&testState, qubits);                         // number of qubits
        cplx_t** vecs = generateTestVectors(qubits);                // State vectors for matrix multiplication

        for (qubit_t control = 0; control < qubits; ++control) {        // For each qubit as a control qubit
            for (qubit_t target = 0; target < qubits; ++target) {       // For each qubit as a target qubit
                if (control == target) {                                // Skip if control and target coincide
                    continue;
                }
                cplx_t* gateMat = controlledGateMat(SMAT, qubits, control, target);     // Controlled operation's matrix
                                                                                        // representation
                for (dim_t i = 0; i < dim + 1; ++i) {                   // For each state vector
                    stateInitVector(&testState, vecs[i]);                   // Copy state vector to testState's vector,
                    applyCS(&testState, control, target);                   // apply test function
                    cplx_t* reference = cmatVecMul(gateMat, vecs[i], dim);  // Matrix multiplication as reference
                    TEST_ASSERT_TRUE(cvectorAlmostEqual(reference, testState.vec, dim, PRECISION));
                    free(reference);
                }

                free(gateMat);
            }
        }
        freeTestVectors(vecs, qubits);
        stateFreeVector(&testState);
    }
}

/*
 * =====================================================================================================================
 *                                              test applyCH
 * =====================================================================================================================
 */
/*
 * This function checks whether the execution of the Pauli-X gate on one qubit, controlled by another, coincides with
 * matrix multiplication
 */
void testApplyCH(void) {
    for (qubit_t qubits = 2; qubits <= MAXQUBITS; ++qubits) {
        dim_t dim = POW2(qubits, dim_t);                            // Hilbert space dimension
        state_t testState;                                          // Initialize a pure quantum state with specified
        stateInitEmpty(&testState, qubits);                         // number of qubits
        cplx_t** vecs = generateTestVectors(qubits);                // State vectors for matrix multiplication

        for (qubit_t control = 0; control < qubits; ++control) {        // For each qubit as a control qubit
            for (qubit_t target = 0; target < qubits; ++target) {       // For each qubit as a target qubit
                if (control == target) {                                // Skip if control and target coincide
                    continue;
                }
                cplx_t* gateMat = controlledGateMat(HMAT, qubits, control, target);     // Controlled operation's matrix
                                                                                        // representation
                for (dim_t i = 0; i < dim + 1; ++i) {                   // For each state vector
                    stateInitVector(&testState, vecs[i]);                   // Copy state vector to testState's vector,
                    applyCH(&testState, control, target);                   // apply test function
                    cplx_t* reference = cmatVecMul(gateMat, vecs[i], dim);  // Matrix multiplication as reference
                    TEST_ASSERT_TRUE(cvectorAlmostEqual(reference, testState.vec, dim, PRECISION));
                    free(reference);
                }

                free(gateMat);
            }
        }
        freeTestVectors(vecs, qubits);
        stateFreeVector(&testState);
    }
}

/*
 * =====================================================================================================================
 *                                              test applyCS
 * =====================================================================================================================
 */

/*
 * This function checks whether the execution of the S-gate on one qubit, controlled by another and its inverse,
 * coincides with matrix multiplication
 */

void testApplyCT(void) {
    for (qubit_t qubits = 2; qubits <= MAXQUBITS; ++qubits) {
        dim_t dim = POW2(qubits, dim_t);                            // Hilbert space dimension
        state_t testState;                                          // Initialize a pure quantum state with specified
        stateInitEmpty(&testState, qubits);                         // number of qubits
        cplx_t** vecs = generateTestVectors(qubits);                // State vectors for matrix multiplication

        for (qubit_t control = 0; control < qubits; ++control) {        // For each qubit as a control qubit
            for (qubit_t target = 0; target < qubits; ++target) {       // For each qubit as a target qubit
                if (control == target) {                                // Skip if control and target coincide
                    continue;
                }
                cplx_t* gateMat = controlledGateMat(TMAT, qubits, control, target);     // Controlled operation's matrix
                                                                                        // representation
                for (dim_t i = 0; i < dim + 1; ++i) {                   // For each state vector
                    stateInitVector(&testState, vecs[i]);                   // Copy state vector to testState's vector,
                    applyCT(&testState, control, target);                   // apply test function
                    cplx_t* reference = cmatVecMul(gateMat, vecs[i], dim);  // Matrix multiplication as reference
                    TEST_ASSERT_TRUE(cvectorAlmostEqual(reference, testState.vec, dim, PRECISION));
                    free(reference);
                }

                free(gateMat);
            }
        }
        freeTestVectors(vecs, qubits);
        stateFreeVector(&testState);
    }
}

/*
 * =====================================================================================================================
 *                                              test applyRX
 * =====================================================================================================================
 */
/*
 * This function checks whether the execution of the X rotation by applyRX and its inverse coincide with matrix
 * multiplication
 */
void testApplyRX(void) {
    for (qubit_t qubits = 1; qubits <= MAXQUBITS; ++qubits) {
        dim_t dim = POW2(qubits, dim_t);                            // Hilbert space dimension
        state_t testState;                                          // Initialize a pure quantum state with specified
        stateInitEmpty(&testState, qubits);                         // number of qubits
        cplx_t** vecs = generateTestVectors(qubits);                // State vectors for matrix multiplication
        double angle = qubits * 0.25;

        for (qubit_t pos = 0; pos < qubits; ++pos) {                    // For each target qubit from the right
            cplx_t* gateMat = RXGateMat(qubits, pos, angle);            // generate the matrix representation from the
                                                                        // Kronecker product
            for (dim_t i = 0; i < dim + 1; ++i) {                           // For each state vector
                stateInitVector(&testState, vecs[i]);                       // Copy state vector to testState's vector,
                applyRX(&testState, pos, angle);                            // apply test function
                cplx_t* ref = cmatVecMul(gateMat, vecs[i], dim);    // Matrix multiplication as reference
                TEST_ASSERT_TRUE(cvectorAlmostEqual(ref, testState.vec, dim, PRECISION));
                free(ref);

                applyRXdagger(&testState, pos, angle);                    // Apply the inverse of RX
                TEST_ASSERT_TRUE(cvectorAlmostEqual(vecs[i], testState.vec, dim, PRECISION));
            }
            free(gateMat);
        }
        freeTestVectors(vecs, qubits);
        stateFreeVector(&testState);
    }
}

/*
 * =====================================================================================================================
 *                                              test applyRY
 * =====================================================================================================================
 */
/*
 * This function checks whether the execution of the Y rotation by applyRY and its inverse coincide with matrix
 * multiplication
 */
void testApplyRY(void) {
    for (qubit_t qubits = 1; qubits <= MAXQUBITS; ++qubits) {
        dim_t dim = POW2(qubits, dim_t);                            // Hilbert space dimension
        state_t testState;                                          // Initialize a pure quantum state with specified
        stateInitEmpty(&testState, qubits);                         // number of qubits
        cplx_t** vecs = generateTestVectors(qubits);                // State vectors for matrix multiplication
        double angle = qubits * 0.25;

        for (qubit_t pos = 0; pos < qubits; ++pos) {                    // For each target qubit from the right
            cplx_t* gateMat = RYGateMat(qubits, pos, angle);            // generate the matrix representation from the
            // Kronecker product
            for (dim_t i = 0; i < dim + 1; ++i) {                           // For each state vector
                stateInitVector(&testState, vecs[i]);                       // Copy state vector to testState's vector,
                applyRY(&testState, pos, angle);                            // apply test function
                cplx_t* ref = cmatVecMul(gateMat, vecs[i], dim);    // Matrix multiplication as reference
                TEST_ASSERT_TRUE(cvectorAlmostEqual(ref, testState.vec, dim, PRECISION));
                free(ref);

                applyRYdagger(&testState, pos, angle);                    // Apply the inverse of RX
                TEST_ASSERT_TRUE(cvectorAlmostEqual(vecs[i], testState.vec, dim, PRECISION));
            }
            free(gateMat);
        }
        freeTestVectors(vecs, qubits);
        stateFreeVector(&testState);
    }
}

/*
 * =====================================================================================================================
 *                                              test applyRZ
 * =====================================================================================================================
 */
/*
 * This function checks whether the execution of the X rotation by applyRZ and its inverse coincide with matrix
 * multiplication
 */
void testApplyRZ(void) {
    for (qubit_t qubits = 1; qubits <= MAXQUBITS; ++qubits) {
        dim_t dim = POW2(qubits, dim_t);                            // Hilbert space dimension
        state_t testState;                                          // Initialize a pure quantum state with specified
        stateInitEmpty(&testState, qubits);                         // number of qubits
        cplx_t** vecs = generateTestVectors(qubits);                // State vectors for matrix multiplication
        double angle = qubits * 0.25;

        for (qubit_t pos = 0; pos < qubits; ++pos) {                    // For each target qubit from the right
            cplx_t* gateMat = RZGateMat(qubits, pos, angle);            // generate the matrix representation from the
                                                                            // Kronecker product
            for (dim_t i = 0; i < dim + 1; ++i) {                           // For each state vector
                stateInitVector(&testState, vecs[i]);                       // Copy state vector to testState's vector,
                applyRZ(&testState, pos, angle);                            // apply test function
                cplx_t* ref = cmatVecMul(gateMat, vecs[i], dim);    // Matrix multiplication as reference
                TEST_ASSERT_TRUE(cvectorAlmostEqual(ref, testState.vec, dim, PRECISION));
                free(ref);

                applyRZdagger(&testState, pos, angle);                    // Apply the inverse of RX
                TEST_ASSERT_TRUE(cvectorAlmostEqual(vecs[i], testState.vec, dim, PRECISION));
            }
            free(gateMat);
        }
        freeTestVectors(vecs, qubits);
        stateFreeVector(&testState);
    }
}
/*
 * =====================================================================================================================
 *                                                      test SWAP
 * =====================================================================================================================
 */
/*
 * This function checks whether the application of a SWAP gate by applySWAP coincides with the matrix multiplication
 */
void testApplySWAP(void) {
    for (qubit_t qubits = 2; qubits <= MAXQUBITS; ++qubits) {
        dim_t dim = POW2(qubits, dim_t);                            // Hilbert space dimension
        state_t testState;                                          // Initialize a pure quantum state with specified
        stateInitEmpty(&testState, qubits);                         // number of qubits
        cplx_t** vecs = generateTestVectors(qubits);                // State vectors for matrix multiplication

        for (qubit_t right = 0; right < qubits - 1 ; ++right) {         // For each constituent in the swap
                                                                        // up to second to last qubit
            for (qubit_t left = right + 1; left < qubits; ++left) {     // For each qubit at least one qubit away
                cplx_t* gateMat = swapGateMat(qubits, right, left);

                for (dim_t i = 0; i < dim + 1; ++i) {
                    stateInitVector(&testState, vecs[i]);
                    applySWAP(&testState, right, left);
                    cplx_t* ref = cmatVecMul(gateMat, vecs[i], dim);
                    TEST_ASSERT_TRUE(cvectorAlmostEqual(ref, testState.vec, dim, PRECISION));
                    free(ref);
                }
                free(gateMat);
            }
        }
        freeTestVectors(vecs, qubits);
        stateFreeVector(&testState);
    }
}

/*
 * =====================================================================================================================
 *                                              test Toffoli
 * =====================================================================================================================
 */
/*
 * This function checks whether the application of a Toffoli-gate by applyToffoli coincides with the matrix
 * multiplication
 */
void testApplyToffoli(void) {
    for (qubit_t qubits = 3; qubits <= MAXQUBITS; ++qubits) {
        dim_t dim = POW2(qubits, dim_t);
        state_t testState;                                          // Initialize a pure quantum state with specified
        stateInitEmpty(&testState, qubits);                         // number of qubits
        cplx_t** vecs = generateTestVectors(qubits);

        for (qubit_t control1 = 0; control1 < qubits; ++control1) {
            for (qubit_t control2 = 0; control2 < qubits; ++control2) {
                for(qubit_t target = 0; target < qubits; ++target) {
                    /*
                    When only one of the three positions coincide continue with the next triple
                    without executing the test.
                    */
                    if (control1 == target || control2 == target || control1 == control2) {
                        continue;
                    }
                    cplx_t* gateMat = doublyControlledGateMat(XMAT, qubits, control1, control2, target);
                    /*
                    Starting at the first test and reference vector, respectively, check whether the
                    application of Toffoli by matrix multiplication and qhipster coincide.
                    */
                    for (dim_t i = 0; i < dim + 1; ++i) {
                        stateInitVector(&testState, vecs[i]);
                        cplx_t* ref = cmatVecMul(gateMat, vecs[i], dim);
                        applyToffoli(&testState, control1, control2, target);

                        TEST_ASSERT_TRUE(cvectorAlmostEqual(ref, testState.vec, dim, PRECISION));
                        free(ref);
                    }
                    free(gateMat);
                }
            }
        }
        freeTestVectors(vecs, qubits);
        stateFreeVector(&testState);
    }
}

/*
 * =====================================================================================================================
 *                                                      main
 * =====================================================================================================================
 */

int main(void) {
    UNITY_BEGIN();
    RUN_TEST(testApplyX);
    RUN_TEST(testApplyY);
    RUN_TEST(testApplyZ);
    RUN_TEST(testApplyS);
    RUN_TEST(testApplyH);
    RUN_TEST(testApplyT);
    // RUN_TEST(testApplyP);
    RUN_TEST(testApplyCX);
    RUN_TEST(testApplyCY);
    RUN_TEST(testApplyCZ);
    RUN_TEST(testApplyCS);
    RUN_TEST(testApplyCH);
    RUN_TEST(testApplyCT);
    // RUN_TEST(testApplyCP);
    RUN_TEST(testApplyRX);
    RUN_TEST(testApplyRY);
    RUN_TEST(testApplyRZ);
    RUN_TEST(testApplySWAP);
    RUN_TEST(testApplyToffoli);
    return UNITY_END();
}
