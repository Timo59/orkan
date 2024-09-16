/*
 * =====================================================================================================================
 *                                                      includes
 * =====================================================================================================================
 */

#include "helpfunc.h"
#include "gatelib.h"
#include "cplxutil.h"
#include "unity.h"

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
 * Input:
 *      -
 * Output:
 *      There is no return value. Only checks whether the state vector after applyX coincides with the state vector
 *      after matrix multiplication for all target qubits.
*/
void testApplyX(void) {
    state_t testState;

    for (qubit_t qubits = 1; qubits <= MAXQUBITS; ++qubits) {

        dim_t dim = POW2(qubits, dim_t);                                    // Hilbert space dimension
        cplx_t** refVectors = generateTestVectors(qubits);                  // State vectors for matrix multiplication

        for (qubit_t pos = 0; pos < qubits; ++pos) {                        // For each target qubit from the right
            cplx_t** testVectors = generateTestVectors(qubits);             // generate state vector to test on
            cplx_t* gateMat = singleQubitGateMat(XMAT, qubits, pos);   // generate the matrix representation from
                                                                            // the Kronecker product

            for (dim_t i = 0; i < dim + 1; ++i) {                                       // For each state vector
                stateInitVector(&testState, testVectors[i], qubits);       // attach vector to test state
                applyX(&testState, pos);                                    // apply test function
                cplx_t* reference = cmatVecMul(gateMat, refVectors[i], dim);    // Matrix multiplication as
                                                                                        // reference
                TEST_ASSERT_TRUE(cvectorAlmostEqual(reference, testState.vec, dim, PRECISION));
                free(reference);
            }
            freeTestVectors(testVectors, qubits);
            free(gateMat);
        }
        freeTestVectors(refVectors, qubits);
    }
}

void testApplyXblas(void) {
    state_t testState;

    for (qubit_t qubits = 1; qubits <= MAXQUBITS; ++qubits) {

        dim_t dim = POW2(qubits, dim_t);                                    // Hilbert space dimension
        cplx_t** refVectors = generateTestVectors(qubits);                  // State vectors for matrix multiplication

        for (qubit_t pos = 0; pos < qubits; ++pos) {                        // For each target qubit from the right
            cplx_t** testVectors = generateTestVectors(qubits);             // generate state vector to test on
            cplx_t* gateMat = singleQubitGateMat(XMAT, qubits, pos);   // generate the matrix representation from
            // the Kronecker product

            for (dim_t i = 0; i < dim + 1; ++i) {                                       // For each state vector
                stateInitVector(&testState, testVectors[i], qubits);       // attach vector to test state
                applyXblas(&testState, pos);                                // apply test function
                cplx_t* reference = cmatVecMul(gateMat, refVectors[i], dim);    // Matrix multiplication as
                // reference
                TEST_ASSERT_TRUE(cvectorAlmostEqual(reference, testState.vec, dim, PRECISION));
                free(reference);
            }
            freeTestVectors(testVectors, qubits);
            free(gateMat);
        }
        freeTestVectors(refVectors, qubits);
    }
}

void testApplyXomp(void) {
    state_t testState;

    for (qubit_t qubits = 1; qubits <= MAXQUBITS; ++qubits) {

        dim_t dim = POW2(qubits, dim_t);                                    // Hilbert space dimension
        cplx_t** refVectors = generateTestVectors(qubits);                  // State vectors for matrix multiplication

        for (qubit_t pos = 0; pos < qubits; ++pos) {                        // For each target qubit from the right
            cplx_t** testVectors = generateTestVectors(qubits);             // generate state vector to test on
            cplx_t* gateMat = singleQubitGateMat(XMAT, qubits, pos);   // generate the matrix representation from
            // the Kronecker product

            for (dim_t i = 0; i < dim + 1; ++i) {                                       // For each state vector
                stateInitVector(&testState, testVectors[i], qubits);       // attach vector to test state
                applyXomp(&testState, pos);                                    // apply test function
                cplx_t* reference = cmatVecMul(gateMat, refVectors[i], dim);    // Matrix multiplication as
                // reference
                TEST_ASSERT_TRUE(cvectorAlmostEqual(reference, testState.vec, dim, PRECISION));
                free(reference);
            }
            freeTestVectors(testVectors, qubits);
            free(gateMat);
        }
        freeTestVectors(refVectors, qubits);
    }
}

void testApplyXomp_blas(void) {
    state_t testState;

    for (qubit_t qubits = 1; qubits <= MAXQUBITS; ++qubits) {

        dim_t dim = POW2(qubits, dim_t);                                    // Hilbert space dimension
        cplx_t** refVectors = generateTestVectors(qubits);                  // State vectors for matrix multiplication

        for (qubit_t pos = 0; pos < qubits; ++pos) {                        // For each target qubit from the right
            cplx_t** testVectors = generateTestVectors(qubits);             // generate state vector to test on
            cplx_t* gateMat = singleQubitGateMat(XMAT, qubits, pos);   // generate the matrix representation from
            // the Kronecker product

            for (dim_t i = 0; i < dim + 1; ++i) {                                       // For each state vector
                stateInitVector(&testState, testVectors[i], qubits);       // attach vector to test state
                applyXomp_blas(&testState, pos);                                    // apply test function
                cplx_t* reference = cmatVecMul(gateMat, refVectors[i], dim);    // Matrix multiplication as
                // reference
                TEST_ASSERT_TRUE(cvectorAlmostEqual(reference, testState.vec, dim, PRECISION));
                free(reference);
            }
            freeTestVectors(testVectors, qubits);
            free(gateMat);
        }
        freeTestVectors(refVectors, qubits);
    }
}

void testApplyXgcd(void) {
    state_t testState;

    for (qubit_t qubits = 1; qubits <= MAXQUBITS; ++qubits) {

        dim_t dim = POW2(qubits, dim_t);                                    // Hilbert space dimension
        cplx_t** refVectors = generateTestVectors(qubits);                  // State vectors for matrix multiplication

        for (qubit_t pos = 0; pos < qubits; ++pos) {                        // For each target qubit from the right
            cplx_t** testVectors = generateTestVectors(qubits);             // generate state vector to test on
            cplx_t* gateMat = singleQubitGateMat(XMAT, qubits, pos);   // generate the matrix representation from
            // the Kronecker product

            for (dim_t i = 0; i < dim + 1; ++i) {                                       // For each state vector
                stateInitVector(&testState, testVectors[i], qubits);       // attach vector to test state
                applyXgcd(&testState, pos);                                    // apply test function
                cplx_t* reference = cmatVecMul(gateMat, refVectors[i], dim);    // Matrix multiplication as
                // reference
                TEST_ASSERT_TRUE(cvectorAlmostEqual(reference, testState.vec, dim, PRECISION));
                free(reference);
            }
            freeTestVectors(testVectors, qubits);
            free(gateMat);
        }
        freeTestVectors(refVectors, qubits);
    }
}

/*
 * =====================================================================================================================
 *                                              test applyY
 * =====================================================================================================================
 */

/*
 * This function checks whether the execution of the Pauli-Y gate by applyY coincides with matrix multiplication
 * Input:
 *      -
 * Output:
 *      There is no return value. Only checks whether the state vector after applyY coincides with the state vector
 *      after matrix multiplication for all target qubits.
*/
void testApplyY(void) {
    state_t testState;

    for (qubit_t qubits = 1; qubits <= MAXQUBITS; ++qubits) {

        dim_t dim = POW2(qubits, dim_t);                                    // Hilbert space dimension
        cplx_t** refVectors = generateTestVectors(qubits);                  // State vectors for matrix multiplication

        for (qubit_t pos = 0; pos < qubits; ++pos) {                        // For each target qubit from the right
            cplx_t** testVectors = generateTestVectors(qubits);             // generate state vector to test on
            cplx_t* gateMat = singleQubitGateMat(YMAT, qubits, pos);   // generate the matrix representation from
            // the Kronecker product

            for (dim_t i = 0; i < dim + 1; ++i) {                                       // For each state vector
                stateInitVector(&testState, testVectors[i], qubits);       // attach vector to test state
                applyY(&testState, pos);                                    // apply test function
                cplx_t* reference = cmatVecMul(gateMat, refVectors[i], dim);   // Matrix multiplication as
                                                                                        // reference
                TEST_ASSERT_TRUE(cvectorAlmostEqual(reference, testState.vec, dim, PRECISION));
                free(reference);
            }
            freeTestVectors(testVectors, qubits);
            free(gateMat);
        }
        freeTestVectors(refVectors, qubits);
    }
}

void testApplyYblas(void) {
    state_t testState;

    for (qubit_t qubits = 1; qubits <= MAXQUBITS; ++qubits) {

        dim_t dim = POW2(qubits, dim_t);                                    // Hilbert space dimension
        cplx_t** refVectors = generateTestVectors(qubits);                  // State vectors for matrix multiplication

        for (qubit_t pos = 0; pos < qubits; ++pos) {                        // For each target qubit from the right
            cplx_t** testVectors = generateTestVectors(qubits);             // generate state vector to test on
            cplx_t* gateMat = singleQubitGateMat(YMAT, qubits, pos);   // generate the matrix representation from
            // the Kronecker product

            for (dim_t i = 0; i < dim + 1; ++i) {                                       // For each state vector
                stateInitVector(&testState, testVectors[i], qubits);       // attach vector to test state
                applyYblas(&testState, pos);                                    // apply test function
                cplx_t* reference = cmatVecMul(gateMat, refVectors[i], dim);   // Matrix multiplication as
                // reference
                TEST_ASSERT_TRUE(cvectorAlmostEqual(reference, testState.vec, dim, PRECISION));
                free(reference);
            }
            freeTestVectors(testVectors, qubits);
            free(gateMat);
        }
        freeTestVectors(refVectors, qubits);
    }
}

void testApplyYomp(void) {
    state_t testState;

    for (qubit_t qubits = 1; qubits <= MAXQUBITS; ++qubits) {

        dim_t dim = POW2(qubits, dim_t);                                    // Hilbert space dimension
        cplx_t** refVectors = generateTestVectors(qubits);                  // State vectors for matrix multiplication

        for (qubit_t pos = 0; pos < qubits; ++pos) {                        // For each target qubit from the right
            cplx_t** testVectors = generateTestVectors(qubits);             // generate state vector to test on
            cplx_t* gateMat = singleQubitGateMat(YMAT, qubits, pos);   // generate the matrix representation from
            // the Kronecker product

            for (dim_t i = 0; i < dim + 1; ++i) {                                       // For each state vector
                stateInitVector(&testState, testVectors[i], qubits);       // attach vector to test state
                applyYomp(&testState, pos);                                    // apply test function
                cplx_t* reference = cmatVecMul(gateMat, refVectors[i], dim);   // Matrix multiplication as
                // reference
                TEST_ASSERT_TRUE(cvectorAlmostEqual(reference, testState.vec, dim, PRECISION));
                free(reference);
            }
            freeTestVectors(testVectors, qubits);
            free(gateMat);
        }
        freeTestVectors(refVectors, qubits);
    }
}

void testApplyYomp_blas(void) {
    state_t testState;

    for (qubit_t qubits = 1; qubits <= MAXQUBITS; ++qubits) {

        dim_t dim = POW2(qubits, dim_t);                                    // Hilbert space dimension
        cplx_t** refVectors = generateTestVectors(qubits);                  // State vectors for matrix multiplication

        for (qubit_t pos = 0; pos < qubits; ++pos) {                        // For each target qubit from the right
            cplx_t** testVectors = generateTestVectors(qubits);             // generate state vector to test on
            cplx_t* gateMat = singleQubitGateMat(YMAT, qubits, pos);   // generate the matrix representation from
            // the Kronecker product

            for (dim_t i = 0; i < dim + 1; ++i) {                                       // For each state vector
                stateInitVector(&testState, testVectors[i], qubits);       // attach vector to test state
                applyYomp_blas(&testState, pos);                                    // apply test function
                cplx_t* reference = cmatVecMul(gateMat, refVectors[i], dim);   // Matrix multiplication as
                // reference
                TEST_ASSERT_TRUE(cvectorAlmostEqual(reference, testState.vec, dim, PRECISION));
                free(reference);
            }
            freeTestVectors(testVectors, qubits);
            free(gateMat);
        }
        freeTestVectors(refVectors, qubits);
    }
}

void testApplyYgcd(void) {
    state_t testState;

    for (qubit_t qubits = 1; qubits <= MAXQUBITS; ++qubits) {

        dim_t dim = POW2(qubits, dim_t);                                    // Hilbert space dimension
        cplx_t** refVectors = generateTestVectors(qubits);                  // State vectors for matrix multiplication

        for (qubit_t pos = 0; pos < qubits; ++pos) {                        // For each target qubit from the right
            cplx_t** testVectors = generateTestVectors(qubits);             // generate state vector to test on
            cplx_t* gateMat = singleQubitGateMat(YMAT, qubits, pos);   // generate the matrix representation from
            // the Kronecker product

            for (dim_t i = 0; i < dim + 1; ++i) {                                       // For each state vector
                stateInitVector(&testState, testVectors[i], qubits);       // attach vector to test state
                applyYgcd(&testState, pos);                                    // apply test function
                cplx_t* reference = cmatVecMul(gateMat, refVectors[i], dim);   // Matrix multiplication as
                // reference
                TEST_ASSERT_TRUE(cvectorAlmostEqual(reference, testState.vec, dim, PRECISION));
                free(reference);
            }
            freeTestVectors(testVectors, qubits);
            free(gateMat);
        }
        freeTestVectors(refVectors, qubits);
    }
}

/*
 * =====================================================================================================================
 *                                              test applyZ
 * =====================================================================================================================
 */

/*
 * This function checks whether the execution of the Pauli-Z gate by applyZ coincides with matrix multiplication
 * Input:
 *      -
 * Output:
 *      There is no return value. Only checks whether the state vector after applyZ coincides with the state vector
 *      after matrix multiplication for all target qubits.
*/
void testApplyZ(void) {
    state_t testState;

    for (qubit_t qubits = 1; qubits <= MAXQUBITS; ++qubits) {

        dim_t dim = POW2(qubits, dim_t);                                    // Hilbert space dimension
        cplx_t** refVectors = generateTestVectors(qubits);                  // State vectors for matrix multiplication

        for (qubit_t pos = 0; pos < qubits; ++pos) {                        // For each target qubit from the right
            cplx_t** testVectors = generateTestVectors(qubits);             // generate state vector to test on
            cplx_t* gateMat = singleQubitGateMat(ZMAT, qubits, pos);   // generate the matrix representation from
                                                                            // the Kronecker product

            for (dim_t i = 0; i < dim + 1; ++i) {                                       // For each state vector
                stateInitVector(&testState, testVectors[i], qubits);        // attach vector to test state
                applyZ(&testState, pos);                                    // apply test function
                cplx_t* reference = cmatVecMul(gateMat, refVectors[i], dim);    // Matrix multiplication as
                                                                                        // reference
                TEST_ASSERT_TRUE(cvectorAlmostEqual(reference, testState.vec, dim, PRECISION));
                free(reference);
            }
            freeTestVectors(testVectors, qubits);
            free(gateMat);
        }
        freeTestVectors(refVectors, qubits);
    }
}

void testApplyZblas(void) {
    state_t testState;

    for (qubit_t qubits = 1; qubits <= MAXQUBITS; ++qubits) {

        dim_t dim = POW2(qubits, dim_t);                                    // Hilbert space dimension
        cplx_t** refVectors = generateTestVectors(qubits);                  // State vectors for matrix multiplication

        for (qubit_t pos = 0; pos < qubits; ++pos) {                        // For each target qubit from the right
            cplx_t** testVectors = generateTestVectors(qubits);             // generate state vector to test on
            cplx_t* gateMat = singleQubitGateMat(ZMAT, qubits, pos);   // generate the matrix representation from
            // the Kronecker product

            for (dim_t i = 0; i < dim + 1; ++i) {                                       // For each state vector
                stateInitVector(&testState, testVectors[i], qubits);        // attach vector to test state
                applyZblas(&testState, pos);                                    // apply test function
                cplx_t* reference = cmatVecMul(gateMat, refVectors[i], dim);    // Matrix multiplication as
                // reference
                TEST_ASSERT_TRUE(cvectorAlmostEqual(reference, testState.vec, dim, PRECISION));
                free(reference);
            }
            freeTestVectors(testVectors, qubits);
            free(gateMat);
        }
        freeTestVectors(refVectors, qubits);
    }
}

void testApplyZomp(void) {
    state_t testState;

    for (qubit_t qubits = 1; qubits <= MAXQUBITS; ++qubits) {

        dim_t dim = POW2(qubits, dim_t);                                    // Hilbert space dimension
        cplx_t** refVectors = generateTestVectors(qubits);                  // State vectors for matrix multiplication

        for (qubit_t pos = 0; pos < qubits; ++pos) {                        // For each target qubit from the right
            cplx_t** testVectors = generateTestVectors(qubits);             // generate state vector to test on
            cplx_t* gateMat = singleQubitGateMat(ZMAT, qubits, pos);   // generate the matrix representation from
            // the Kronecker product

            for (dim_t i = 0; i < dim + 1; ++i) {                                       // For each state vector
                stateInitVector(&testState, testVectors[i], qubits);        // attach vector to test state
                applyZomp(&testState, pos);                                    // apply test function
                cplx_t* reference = cmatVecMul(gateMat, refVectors[i], dim);    // Matrix multiplication as
                // reference
                TEST_ASSERT_TRUE(cvectorAlmostEqual(reference, testState.vec, dim, PRECISION));
                free(reference);
            }
            freeTestVectors(testVectors, qubits);
            free(gateMat);
        }
        freeTestVectors(refVectors, qubits);
    }
}

void testApplyZomp_blas(void) {
    state_t testState;

    for (qubit_t qubits = 1; qubits <= MAXQUBITS; ++qubits) {

        dim_t dim = POW2(qubits, dim_t);                                    // Hilbert space dimension
        cplx_t** refVectors = generateTestVectors(qubits);                  // State vectors for matrix multiplication

        for (qubit_t pos = 0; pos < qubits; ++pos) {                        // For each target qubit from the right
            cplx_t** testVectors = generateTestVectors(qubits);             // generate state vector to test on
            cplx_t* gateMat = singleQubitGateMat(ZMAT, qubits, pos);   // generate the matrix representation from
            // the Kronecker product

            for (dim_t i = 0; i < dim + 1; ++i) {                                       // For each state vector
                stateInitVector(&testState, testVectors[i], qubits);        // attach vector to test state
                applyZomp_blas(&testState, pos);                                    // apply test function
                cplx_t* reference = cmatVecMul(gateMat, refVectors[i], dim);    // Matrix multiplication as
                // reference
                TEST_ASSERT_TRUE(cvectorAlmostEqual(reference, testState.vec, dim, PRECISION));
                free(reference);
            }
            freeTestVectors(testVectors, qubits);
            free(gateMat);
        }
        freeTestVectors(refVectors, qubits);
    }
}

void testApplyZgcd(void) {
    state_t testState;

    for (qubit_t qubits = 1; qubits <= MAXQUBITS; ++qubits) {

        dim_t dim = POW2(qubits, dim_t);                                    // Hilbert space dimension
        cplx_t** refVectors = generateTestVectors(qubits);                  // State vectors for matrix multiplication

        for (qubit_t pos = 0; pos < qubits; ++pos) {                        // For each target qubit from the right
            cplx_t** testVectors = generateTestVectors(qubits);             // generate state vector to test on
            cplx_t* gateMat = singleQubitGateMat(ZMAT, qubits, pos);   // generate the matrix representation from
            // the Kronecker product

            for (dim_t i = 0; i < dim + 1; ++i) {                                       // For each state vector
                stateInitVector(&testState, testVectors[i], qubits);        // attach vector to test state
                applyZgcd(&testState, pos);                                    // apply test function
                cplx_t* reference = cmatVecMul(gateMat, refVectors[i], dim);    // Matrix multiplication as
                // reference
                TEST_ASSERT_TRUE(cvectorAlmostEqual(reference, testState.vec, dim, PRECISION));
                free(reference);
            }
            freeTestVectors(testVectors, qubits);
            free(gateMat);
        }
        freeTestVectors(refVectors, qubits);
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
 *
 * Input:
 *      -
 *      Output:
 *      There is no return value. Only checks whether the state vector after applyS coincides with the state vector
 *      after matrix multiplication for all target qubits and whether the subsequent application of S^dagger leads to
 *      the starting vector.
*/
void testApplyS(void) {
    state_t testState;

    for (qubit_t qubits = 1; qubits <= MAXQUBITS; ++qubits) {

        dim_t dim = POW2(qubits, dim_t);                                    // Hilbert space dimension
        cplx_t** refVectors = generateTestVectors(qubits);                  // State vectors for matrix multiplication

        for (qubit_t pos = 0; pos < qubits; ++pos) {                        // For each target qubit from the right
            cplx_t** testVectors = generateTestVectors(qubits);             // generate state vector to test on
            cplx_t* gateMat = singleQubitGateMat(SMAT, qubits, pos);   // generate the matrix representation from
                                                                            // the Kronecker product

            for (dim_t i = 0; i < dim + 1; ++i) {                                       // For each state vector
                stateInitVector(&testState, testVectors[i], qubits);        // attach vector to test state
                applyS(&testState, pos);                                    // apply test function
                cplx_t* reference = cmatVecMul(gateMat, refVectors[i], dim);    // Matrix multiplication as
                                                                                        // reference
                TEST_ASSERT_TRUE(cvectorAlmostEqual(reference, testState.vec, dim, PRECISION));
                free(reference);

                applySdagger(&testState, pos);                  // Apply the inverse of S
                TEST_ASSERT_TRUE(cvectorAlmostEqual(refVectors[i], testState.vec, dim, PRECISION));
            }
            freeTestVectors(testVectors, qubits);
            free(gateMat);
        }
        freeTestVectors(refVectors, qubits);
    }
}

/*
 * =====================================================================================================================
 *                                              test applyH
 * =====================================================================================================================
 */

/*
 * This function checks whether the execution of the Pauli-Z gate by applyZ coincides with matrix multiplication
 * Input:
 *      -
 * Output:
 *      There is no return value. Only checks whether the state vector after applyZ coincides with the state vector
 *      after matrix multiplication for all target qubits.
*/
void testApplyH(void) {
    state_t testState;

    for (qubit_t qubits = 1; qubits <= MAXQUBITS; ++qubits) {

        dim_t dim = POW2(qubits, dim_t);                                    // Hilbert space dimension
        cplx_t** refVectors = generateTestVectors(qubits);                  // State vectors for matrix multiplication

        for (qubit_t pos = 0; pos < qubits; ++pos) {                        // For each target qubit from the right
            cplx_t** testVectors = generateTestVectors(qubits);             // generate state vector to test on
            cplx_t* gateMat = singleQubitGateMat(HMAT, qubits, pos);   // generate the matrix representation from
            // the Kronecker product

            for (dim_t i = 0; i < dim + 1; ++i) {                                       // For each state vector
                stateInitVector(&testState, testVectors[i], qubits);        // attach vector to test state
                applyH(&testState, pos);                                    // apply test function
                cplx_t* reference = cmatVecMul(gateMat, refVectors[i], dim);    // Matrix multiplication as
                // reference
                TEST_ASSERT_TRUE(cvectorAlmostEqual(reference, testState.vec, dim, PRECISION));
                free(reference);
            }
            freeTestVectors(testVectors, qubits);
            free(gateMat);
        }
        freeTestVectors(refVectors, qubits);
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
    state_t testState;

    for (qubit_t qubits = 1; qubits <= MAXQUBITS; ++qubits) {

        dim_t dim = POW2(qubits, dim_t);                                    // Hilbert space dimension
        cplx_t** refVectors = generateTestVectors(qubits);                  // State vectors for matrix multiplication

        for (qubit_t pos = 0; pos < qubits; ++pos) {                        // For each target qubit from the right
            cplx_t** testVectors = generateTestVectors(qubits);             // generate state vector to test on
            cplx_t* gateMat = singleQubitGateMat(TMAT, qubits, pos);   // generate the matrix representation from
                                                                            // the Kronecker product

            for (dim_t i = 0; i < dim + 1; ++i) {                                       // For each state vector
                stateInitVector(&testState, testVectors[i], qubits);        // attach vector to test state
                applyT(&testState, pos);                                    // apply test function
                cplx_t* reference = cmatVecMul(gateMat, refVectors[i], dim);    // Matrix multiplication as
                                                                                        // reference
                TEST_ASSERT_TRUE(cvectorAlmostEqual(reference, testState.vec, dim, PRECISION));
                free(reference);

                applyTdagger(&testState, pos);                  // Apply the inverse of S
                TEST_ASSERT_TRUE(cvectorAlmostEqual(refVectors[i], testState.vec, dim, PRECISION));
            }
            freeTestVectors(testVectors, qubits);
            free(gateMat);
        }
        freeTestVectors(refVectors, qubits);
    }
}

/*
 * =====================================================================================================================
 *                                              test applyCX
 * =====================================================================================================================
 */

/*
 * This function checks whether the execution of the Pauli-X gate on one qubit, controlled by another, coincides with
 * matrix multiplication.
 * Input:
 *      -
 * Output:
 *      There is no return value. Only checks whether the state vector after applyCX coincides with the state vector
 *      after matrix multiplication for all target qubits.
*/

void testApplyCX(void) {
    state_t testState;

    for (qubit_t qubits = 2; qubits <= MAXQUBITS; ++qubits) {
        dim_t dim = POW2(qubits, dim_t);                                    // Hilbert space dimension
        cplx_t** refVectors = generateTestVectors(qubits);                  // State vectors for matrix multiplication

        for (qubit_t control = 0; control < qubits; ++control) {                // For each qubit as a control qubit
            for (qubit_t target = 0; target < qubits; ++target) {               // For each qubit as a target qubit
                if (control == target) {                                        // qubits coincide
                    continue;
                }
                cplx_t** testVectors = generateTestVectors(qubits);             // Do nothing if control and target
                cplx_t* gateMat = controlledGateMat(XMAT, qubits, control, target); // Controlled operation's matrix
                                                                                        // representation
                for (dim_t i = 0; i < dim + 1; ++i) {                                   // For each state vector
                    stateInitVector(&testState, testVectors[i], qubits);    // attach vector to test state
                    applyCX(&testState, control, target);                           // apply test function
                    cplx_t* reference = cmatVecMul(gateMat, refVectors[i], dim);    // Matrix multiplication as
                                                                                        // reference
                    TEST_ASSERT_TRUE(cvectorAlmostEqual(reference, testState.vec, dim, PRECISION));
                    free(reference);
                }
                freeTestVectors(testVectors, qubits);
                free(gateMat);
            }
        }
        freeTestVectors(refVectors, qubits);
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
 * Input:
 *      -
 * Output:
 *      There is no return value. Only checks whether the state vector after applyCY coincides with the state vector
 *      after matrix multiplication for all target qubits.
*/

void testApplyCY(void) {
    state_t testState;

    for (qubit_t qubits = 2; qubits <= MAXQUBITS; ++qubits) {
        dim_t dim = POW2(qubits, dim_t);                                    // Hilbert space dimension
        cplx_t** refVectors = generateTestVectors(qubits);                  // State vectors for matrix multiplication

        for (qubit_t control = 0; control < qubits; ++control) {                // For each qubit as a control qubit
            for (qubit_t target = 0; target < qubits; ++target) {               // For each qubit as a target qubit
                if (control == target) {                                        // Do nothing if control and target
                    continue;                                                   // qubits coincide
                }
                cplx_t** testVectors = generateTestVectors(qubits);
                cplx_t* gateMat = controlledGateMat(YMAT, qubits, control, target); // Controlled operation's matrix
                                                                                        // representation
                for (dim_t i = 0; i < dim + 1; ++i) {                                   // For each state vector
                    stateInitVector(&testState, testVectors[i], qubits);    // attach vector to test state
                    applyCY(&testState, control, target);                           // apply test function
                    cplx_t* reference = cmatVecMul(gateMat, refVectors[i], dim);    // Matrix multiplication as
                                                                                        // reference
                    TEST_ASSERT_TRUE(cvectorAlmostEqual(reference, testState.vec, dim, PRECISION));
                    free(reference);
                }
                freeTestVectors(testVectors, qubits);
                free(gateMat);
            }
        }
        freeTestVectors(refVectors, qubits);
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
 * Input:
 *      -
 * Output:
 *      There is no return value. Only checks whether the state vector after applyCZ coincides with the state vector
 *      after matrix multiplication for all target qubits.
*/

void testApplyCZ(void) {
    state_t testState;

    for (qubit_t qubits = 2; qubits <= MAXQUBITS; ++qubits) {
        dim_t dim = POW2(qubits, dim_t);                                    // Hilbert space dimension
        cplx_t** refVectors = generateTestVectors(qubits);                  // State vectors for matrix multiplication

        for (qubit_t control = 0; control < qubits; ++control) {                // For each qubit as a control qubit
            for (qubit_t target = 0; target < qubits; ++target) {               // For each qubit as a target qubit
                if (control == target) {                                        // Do nothing if control and target
                    continue;                                                   // qubits coincide
                }
                cplx_t** testVectors = generateTestVectors(qubits);
                cplx_t* gateMat = controlledGateMat(ZMAT, qubits, control, target); // Controlled operation's matrix
                                                                                        // representation
                for (dim_t i = 0; i < dim + 1; ++i) {                                   // For each state vector
                    stateInitVector(&testState, testVectors[i], qubits);    // attach vector to test state
                    applyCZ(&testState, control, target);                           // apply test function
                    cplx_t* reference = cmatVecMul(gateMat, refVectors[i], dim);    // Matrix multiplication as
                    // reference
                    TEST_ASSERT_TRUE(cvectorAlmostEqual(reference, testState.vec, dim, PRECISION));
                    free(reference);
                }
                freeTestVectors(testVectors, qubits);
                free(gateMat);
            }
        }
        freeTestVectors(refVectors, qubits);
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
 * Input:
 *      -
 * Output:
 *      There is no return value. Only checks whether the state vector after applyCS coincides with the state vector
 *      after matrix multiplication for all target qubits and whether the subsequent application of the inverse returns
 *      the starting state vector.
*/

void testApplyCS(void) {
    state_t testState;

    for (qubit_t qubits = 2; qubits <= MAXQUBITS; ++qubits) {
        dim_t dim = POW2(qubits, dim_t);                                    // Hilbert space dimension
        cplx_t** refVectors = generateTestVectors(qubits);                  // State vectors for matrix multiplication

        for (qubit_t control = 0; control < qubits; ++control) {                // For each qubit as a control qubit
            for (qubit_t target = 0; target < qubits; ++target) {               // For each qubit as a target qubit
                if (control == target) {                                        // Do nothing if control and target
                    continue;                                                   // qubits coincide
                }
                cplx_t** testVectors = generateTestVectors(qubits);
                cplx_t* gateMat = controlledGateMat(SMAT, qubits, control, target); // Controlled operation's matrix
                                                                                        // representation
                for (dim_t i = 0; i < dim + 1; ++i) {                                   // For each state vector
                    stateInitVector(&testState, testVectors[i], qubits);    // attach vector to test state
                    applyCS(&testState, control, target);                           // apply test function
                    cplx_t* reference = cmatVecMul(gateMat, refVectors[i], dim);    // Matrix multiplication as
                                                                                        // reference
                    TEST_ASSERT_TRUE(cvectorAlmostEqual(reference, testState.vec, dim, PRECISION));
                    free(reference);

                    applyCSdagger(&testState, control, target);                     // Apply the inverse test fucntion
                    TEST_ASSERT_TRUE(cvectorAlmostEqual(refVectors[i], testState.vec, dim, PRECISION));
                }
                freeTestVectors(testVectors, qubits);
                free(gateMat);
            }
        }
        freeTestVectors(refVectors, qubits);
    }
}

/*
 * =====================================================================================================================
 *                                              test applyCH
 * =====================================================================================================================
 */

/*
 * This function checks whether the execution of the Pauli-X gate on one qubit, controlled by another, coincides with
 * matrix multiplication.
 * Input:
 *      -
 * Output:
 *      There is no return value. Only checks whether the state vector after applyCX coincides with the state vector
 *      after matrix multiplication for all target qubits.
*/

void testApplyCH(void) {
    state_t testState;

    for (qubit_t qubits = 2; qubits <= MAXQUBITS; ++qubits) {
        dim_t dim = POW2(qubits, dim_t);                                    // Hilbert space dimension
        cplx_t** refVectors = generateTestVectors(qubits);                  // State vectors for matrix multiplication

        for (qubit_t control = 0; control < qubits; ++control) {                // For each qubit as a control qubit
            for (qubit_t target = 0; target < qubits; ++target) {               // For each qubit as a target qubit
                if (control == target) {                                        // Do nothing if control and target
                    continue;                                                   // qubits coincide
                }
                cplx_t** testVectors = generateTestVectors(qubits);
                cplx_t* gateMat = controlledGateMat(HMAT, qubits, control, target); // Controlled operation's matrix
                                                                                        // representation
                for (dim_t i = 0; i < dim + 1; ++i) {                                   // For each state vector
                    stateInitVector(&testState, testVectors[i], qubits);    // attach vector to test state
                    applyCH(&testState, control, target);                           // apply test function
                    cplx_t* reference = cmatVecMul(gateMat, refVectors[i], dim);    // Matrix multiplication as
                    // reference
                    TEST_ASSERT_TRUE(cvectorAlmostEqual(reference, testState.vec, dim, PRECISION));
                    free(reference);
                }
                freeTestVectors(testVectors, qubits);
                free(gateMat);
            }
        }
        freeTestVectors(refVectors, qubits);
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
 * Input:
 *      -
 * Output:
 *      There is no return value. Only checks whether the state vector after applyCS coincides with the state vector
 *      after matrix multiplication for all target qubits and whether the subsequent application of the inverse returns
 *      the starting state vector.
*/

void testApplyCT(void) {
    state_t testState;

    for (qubit_t qubits = 2; qubits <= MAXQUBITS; ++qubits) {
        dim_t dim = POW2(qubits, dim_t);                                    // Hilbert space dimension
        cplx_t** refVectors = generateTestVectors(qubits);                  // State vectors for matrix multiplication

        for (qubit_t control = 0; control < qubits; ++control) {                // For each qubit as a control qubit
            for (qubit_t target = 0; target < qubits; ++target) {               // For each qubit as a target qubit
                if (control == target) {                                        // Do nothing if control and target
                    continue;                                                   // qubits coincide
                }
                cplx_t** testVectors = generateTestVectors(qubits);
                cplx_t* gateMat = controlledGateMat(TMAT, qubits, control, target); // Controlled operation's matrix
                // representation
                for (dim_t i = 0; i < dim + 1; ++i) {                                   // For each state vector
                    stateInitVector(&testState, testVectors[i], qubits);    // attach vector to test state
                    applyCT(&testState, control, target);                           // apply test function
                    cplx_t* reference = cmatVecMul(gateMat, refVectors[i], dim);    // Matrix multiplication as
                    // reference
                    TEST_ASSERT_TRUE(cvectorAlmostEqual(reference, testState.vec, dim, PRECISION));
                    free(reference);

                    applyCTdagger(&testState, control, target);                     // Apply the inverse test fucntion
                    TEST_ASSERT_TRUE(cvectorAlmostEqual(refVectors[i], testState.vec, dim, PRECISION));
                }
                freeTestVectors(testVectors, qubits);
                free(gateMat);
            }
        }
        freeTestVectors(refVectors, qubits);
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
 *
 * Input:
 *      -
 *      Output:
 *      There is no return value. Only checks whether the state vector after applyRX coincides with the state vector
 *      after matrix multiplication for all target qubits and whether the subsequent application of Rx^dagger leads to
 *      the starting vector.
*/
void testApplyRX(void) {
    state_t testState;

    for (qubit_t qubits = 1; qubits <= MAXQUBITS; ++qubits) {

        dim_t dim = POW2(qubits, dim_t);                                    // Hilbert space dimension
        cplx_t** refVectors = generateTestVectors(qubits);                  // State vectors for matrix multiplication
        double angle = qubits * 0.25;

        for (qubit_t pos = 0; pos < qubits; ++pos) {                        // For each target qubit from the right
            cplx_t** testVectors = generateTestVectors(qubits);             // generate state vector to test on
            cplx_t* gateMat = RXGateMat(qubits, pos, angle);                // generate the matrix representation from
            // the Kronecker product

            for (dim_t i = 0; i < dim + 1; ++i) {                                       // For each state vector
                stateInitVector(&testState, testVectors[i], qubits);        // attach vector to test state
                applyRX(&testState, pos, angle);                            // apply test function
                cplx_t* reference = cmatVecMul(gateMat, refVectors[i], dim);    // Matrix multiplication as
                // reference
                TEST_ASSERT_TRUE(cvectorAlmostEqual(reference, testState.vec, dim, PRECISION));
                free(reference);

                applyRXdagger(&testState, pos, angle);                    // Apply the inverse of RX
                TEST_ASSERT_TRUE(cvectorAlmostEqual(refVectors[i], testState.vec, dim, PRECISION));
            }
            freeTestVectors(testVectors, qubits);
            free(gateMat);
        }
        freeTestVectors(refVectors, qubits);
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
 *
 * Input:
 *      -
 *      Output:
 *      There is no return value. Only checks whether the state vector after applyRY coincides with the state vector
 *      after matrix multiplication for all target qubits and whether the subsequent application of Ry^dagger leads to
 *      the starting vector.
*/
void testApplyRY(void) {
    state_t testState;

    for (qubit_t qubits = 1; qubits <= MAXQUBITS; ++qubits) {

        dim_t dim = POW2(qubits, dim_t);                                    // Hilbert space dimension
        cplx_t** refVectors = generateTestVectors(qubits);                  // State vectors for matrix multiplication
        double angle = qubits * 0.25;

        for (qubit_t pos = 0; pos < qubits; ++pos) {                        // For each target qubit from the right
            cplx_t** testVectors = generateTestVectors(qubits);             // generate state vector to test on
            cplx_t* gateMat = RYGateMat(qubits, pos, angle);                // generate the matrix representation from
            // the Kronecker product

            for (dim_t i = 0; i < dim + 1; ++i) {                                       // For each state vector
                stateInitVector(&testState, testVectors[i], qubits);        // attach vector to test state
                applyRY(&testState, pos, angle);                            // apply test function
                cplx_t* reference = cmatVecMul(gateMat, refVectors[i], dim);    // Matrix multiplication as
                // reference
                TEST_ASSERT_TRUE(cvectorAlmostEqual(reference, testState.vec, dim, PRECISION));
                free(reference);

                applyRYdagger(&testState, pos, angle);                    // Apply the inverse of RX
                TEST_ASSERT_TRUE(cvectorAlmostEqual(refVectors[i], testState.vec, dim, PRECISION));
            }
            freeTestVectors(testVectors, qubits);
            free(gateMat);
        }
        freeTestVectors(refVectors, qubits);
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
 *
 * Input:
 *      -
 *      Output:
 *      There is no return value. Only checks whether the state vector after applyRZ coincides with the state vector
 *      after matrix multiplication for all target qubits and whether the subsequent application of Rz^dagger leads to
 *      the starting vector.
*/
void testApplyRZ(void) {
    state_t testState;

    for (qubit_t qubits = 1; qubits <= MAXQUBITS; ++qubits) {

        dim_t dim = POW2(qubits, dim_t);                                    // Hilbert space dimension
        cplx_t** refVectors = generateTestVectors(qubits);                  // State vectors for matrix multiplication
        double angle = qubits * 0.25;

        for (qubit_t pos = 0; pos < qubits; ++pos) {                        // For each target qubit from the right
            cplx_t** testVectors = generateTestVectors(qubits);             // generate state vector to test on
            cplx_t* gateMat = RZGateMat(qubits, pos, angle);                // generate the matrix representation from
                                                                            // the Kronecker product

            for (dim_t i = 0; i < dim + 1; ++i) {                                       // For each state vector
                stateInitVector(&testState, testVectors[i], qubits);        // attach vector to test state
                applyRZ(&testState, pos, angle);                            // apply test function
                cplx_t* reference = cmatVecMul(gateMat, refVectors[i], dim);    // Matrix multiplication as
                // reference
                TEST_ASSERT_TRUE(cvectorAlmostEqual(reference, testState.vec, dim, PRECISION));
                free(reference);

                applyRZdagger(&testState, pos, angle);                    // Apply the inverse of RX
                TEST_ASSERT_TRUE(cvectorAlmostEqual(refVectors[i], testState.vec, dim, PRECISION));
            }
            freeTestVectors(testVectors, qubits);
            free(gateMat);
        }
        freeTestVectors(refVectors, qubits);
    }
}

/*
 * =====================================================================================================================
 *                                                      test SWAP
 * =====================================================================================================================
 */

/*
 * This function checks whether the application of a SWAP gate by applySWAP coincides with the matrix multiplication
 * Input:
 *      -
 * Output:
 *      There is no return value. Only checks whether the state vector after applySWAP coincides with the state vector
 *      from matrix multiplication to some precision.
*/
void testApplySWAP(void) {
    state_t testState;

    for (qubit_t qubits = 2; qubits <= MAXQUBITS; ++qubits) {

        dim_t dim = POW2(qubits, dim_t);                        // Hilbert space dimension
        cplx_t** refVectors = generateTestVectors(qubits);      // State vectors for matrix multiplication

        for (qubit_t right_qubit = 0; right_qubit < qubits - 1 ; ++right_qubit) {  // For each constituent in the swap
                                                                                    // up to second to last qubit

            for (qubit_t left_qubit = right_qubit + 1; left_qubit < qubits; ++left_qubit) {    // For each qubit at
                                                                                                // least one qubit away
                cplx_t** testVectors = generateTestVectors(qubits);
                cplx_t* gateMat = swapGateMat(qubits, right_qubit, left_qubit);
                for (dim_t i = 0; i < dim + 1; ++i) {
                    stateInitVector(&testState, testVectors[i], qubits);
                    applySWAP(&testState, right_qubit, left_qubit);
                    cplx_t* reference = cmatVecMul(gateMat, refVectors[i], dim);
                    TEST_ASSERT_TRUE(cvectorAlmostEqual(reference, testState.vec, dim, PRECISION));

                    free(reference);
                }
                freeTestVectors(testVectors, qubits);
                free(gateMat);
            }
        }
        freeTestVectors(refVectors, qubits);
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
 * Input:
 *      -
 * Output:
 *      There is no return value. Only checks whether the state vector after applyToffoli coincides with the state
 *      vector from matrix multiplication to some precision.
*/

void testApplyToffoli(void) {
    state_t testState;

    for (qubit_t qubits = 3; qubits <= MAXQUBITS; ++qubits) {
        dim_t dim = POW2(qubits, dim_t);
        cplx_t** refVectors = generateTestVectors(qubits);

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
                    cplx_t** testVectors = generateTestVectors(qubits);
                    cplx_t* gateMat = doublyControlledGateMat(XMAT, qubits, control1, control2, target);
                    /*
                    Starting at the first test and reference vector, respectively, check whether the
                    application of Toffoli by matrix multiplication and qhipster coincide.
                    */
                    for (dim_t i = 0; i < dim + 1; ++i) {
                        stateInitVector(&testState, testVectors[i], qubits);
                        cplx_t* reference = cmatVecMul(gateMat, testVectors[i], dim);
                        applyToffoli(&testState, control1, control2, target);

                        TEST_ASSERT_TRUE(cvectorAlmostEqual(reference, testState.vec, dim, PRECISION));
                        free(reference);
                    }
                    freeTestVectors(testVectors, qubits);
                    free(gateMat);
                }
            }
        }
        freeTestVectors(refVectors, qubits);
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
    RUN_TEST(testApplyXblas);
    RUN_TEST(testApplyXomp);
    RUN_TEST(testApplyXomp_blas);
    RUN_TEST(testApplyXgcd);
    RUN_TEST(testApplyY);
    RUN_TEST(testApplyYblas);
    RUN_TEST(testApplyYomp);
    RUN_TEST(testApplyYomp_blas);
    RUN_TEST(testApplyYgcd);
    RUN_TEST(testApplyZ);
    RUN_TEST(testApplyZblas);
    RUN_TEST(testApplyZomp);
    RUN_TEST(testApplyZomp_blas);
    RUN_TEST(testApplyZgcd);
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
