/*
 * =================================================================================================
 *                                              includes
 * =================================================================================================
 */

#include "helpfunc.h"
#include "gatelib.h"
#include "cplxutil.h"
#include "unity.h"

/*
 * =================================================================================================
 *                                              macros
 * =================================================================================================
 */

#define PRECISION               1e-8
#define MAXQUBITS               6

/*
 * =================================================================================================
 *                                              setUp and tearDown
 * =================================================================================================
 */

void setUp(void) {}

void tearDown(void) {}

/*
 * =================================================================================================
 *                                              test applyX
 * =================================================================================================
 */

/*
This function compares the execution of the Pauli-X gate of the tested method with the matrix
multiplication method.

Input:
    -
Output:
    There is no return value. Only checks whether the tested state vector and the reference state
    vector are equal to some precision.
*/
void testApplyX(void) {
    dim_t dim;
    cplx_t** testVectors;       // double complex pointer to an array of pointers defining the set
                                // of reference vectors

    cplx_t** stateVectors;      // double complex pointer to an array of pointers defining the set
                                // of test vectors

    state_t testState;          // custom struct holding the double complex test state vector, the
                                // number of qubits and the Hilbert's space dimension

    cplx_t* gateMatrix;         // double complex pointer to the reference matrix emerging from 
                                // Kronecker products of identity matrices and the Pauli-X

    cplx_t* result;             // double complex pointer holding the result of the reference vector
                                // from matrix vector multiplication

    /*
    Starting with one qubit, the instances are created up to the number MAXQUBITS defined in the
    macros.
    */
    for (qubit_t qubits = 1; qubits <= MAXQUBITS; ++qubits) {

        dim = POW2(qubits, dim_t);                                          // unsigned integer
                                                                            // holding the Hilbert
                                                                            // space dimension
                                                                            // according to the
                                                                            // number of qubits                    
        
        testVectors = generateTestVectors(qubits);                          // generates the set of
                                                                            // reference vectors and
                                                                            // attaches them to 
                                                                            // their pointer
        /*
        Starting with the rightmost qubit the Pauli-X gate is applied to each qubit for all test
        vectors.
        */
        for (qubit_t position = 0; position < qubits; ++position) {
            stateVectors = generateTestVectors(qubits);                     // generates the set of
                                                                            // test vectors and
                                                                            // attaches them to
                                                                            // their pointer
    
            gateMatrix = singleQubitGateMat(XMAT, qubits, position);        // initiates the matrix
                                                                            // corresponding to the
                                                                            // Pauli-X gate on the  
                                                                            // current qubit
            /*
            Starting at the first reference and test vector, respectively, goes through all test
            vectors checking if they coincide step by step.
            */
            for (dim_t i = 0; i < dim + 1; ++i) {
                stateInitVector(&testState, stateVectors[i], qubits);       // initialises the test
                                                                            // state with the test
                                                                            // vector corresponding
                                                                            // to the current i

                result = cmatVecMul(gateMatrix, testVectors[i], dim);       // multiplies the 
                                                                            // reference vector
                                                                            // corresponding to the
                                                                            // current i with the 
                                                                            // gate matrix from 
                                                                            // Kronecker product

                applyXgcd(&testState, position);                 // applies the Pauli-X
                                                                            // gate according to the
                                                                            // tested method

                TEST_ASSERT_TRUE(cvectorAlmostEqual(result, testState.vector, dim, PRECISION));
                                                                            // evaluate if the state
                                                                            // vector obtained by
                                                                            // matrix multiplication
                                                                            // coincides with the
                                                                            // one obtained by the
                                                                            // tested method

                free(result);                                               // free the memory
                                                                            // allocated to the
                                                                            // results pointer for
                                                                            // the next reference
                                                                            // vector
            }
            freeTestVectors(stateVectors, qubits);                          // free all test vectors
                                                                            // altered by the tested
                                                                            // method for the next
                                                                            // position of the
                                                                            // Pauli-X gate                          
        }
        freeTestVectors(testVectors, qubits);                               // free all reference
                                                                            // vectors used for
                                                                            // matrix vector 
                                                                            // multiplication for 
                                                                            // a higher number of 
                                                                            // qubits
    }
}

/*
 * =================================================================================================
 *                                              test applyY
 * =================================================================================================
 */

void testApplyY(void) {
    dim_t dim;
    cplx_t** testVectors;
    cplx_t** stateVectors;
    state_t testState;
    cplx_t* gateMatrix;
    cplx_t* result;
    for (qubit_t qubits = 1; qubits <= MAXQUBITS; ++qubits) {
        dim = POW2(qubits, dim_t);
        testVectors = generateTestVectors(qubits);
        for (qubit_t position = 0; position < qubits; ++position) {
            stateVectors = generateTestVectors(qubits);
            gateMatrix = singleQubitGateMat(YMAT, qubits, position);
            for (dim_t i = 0; i < dim + 1; ++i) {
                stateInitVector(&testState, stateVectors[i], qubits);
                result = cmatVecMul(gateMatrix, testVectors[i], dim);
                applyY(&testState, position);
                TEST_ASSERT_TRUE(cvectorAlmostEqual(result, testState.vector, dim, PRECISION));
                free(result);
            }
            freeTestVectors(stateVectors, qubits);
        }
        freeTestVectors(testVectors, qubits);
    }
}

/*
 * =================================================================================================
 *                                              test applyZ
 * =================================================================================================
 */

void testApplyZ(void) {
    dim_t dim;
    cplx_t** testVectors;
    cplx_t** stateVectors;
    state_t testState;
    cplx_t* gateMatrix;
    cplx_t* result;
    for (qubit_t qubits = 1; qubits <= MAXQUBITS; ++qubits) {
        dim = POW2(qubits, dim_t);
        testVectors = generateTestVectors(qubits);
        for (qubit_t position = 0; position < qubits; ++position) {
            stateVectors = generateTestVectors(qubits);
            gateMatrix = singleQubitGateMat(ZMAT, qubits, position);
            for (dim_t i = 0; i < dim + 1; ++i) {
                stateInitVector(&testState, stateVectors[i], qubits);
                result = cmatVecMul(gateMatrix, testVectors[i], dim);
                applyZ(&testState, position);
                TEST_ASSERT_TRUE(cvectorAlmostEqual(result, testState.vector, dim, PRECISION));
                free(result);
            }
            freeTestVectors(stateVectors, qubits);
        }
        freeTestVectors(testVectors, qubits);
    }
}

/*
 * =================================================================================================
 *                                              test applyS(dagger)
 * =================================================================================================
 */

void testApplyS(void) {
    dim_t dim;
    cplx_t** testVectors;
    cplx_t** stateVectors;
    state_t testState;
    cplx_t* gateMatrix;
    cplx_t* result;
    for (qubit_t qubits = 1; qubits <= MAXQUBITS; ++qubits) {
        dim = POW2(qubits, dim_t);
        testVectors = generateTestVectors(qubits);
        for (qubit_t position = 0; position < qubits; ++position) {
            stateVectors = generateTestVectors(qubits);
            gateMatrix = singleQubitGateMat(SMAT, qubits, position);
            for (dim_t i = 0; i < dim + 1; ++i) {
                stateInitVector(&testState, stateVectors[i], qubits);
                result = cmatVecMul(gateMatrix, testVectors[i], dim);
                applyS(&testState, position);
                TEST_ASSERT_TRUE(cvectorAlmostEqual(result, testState.vector, dim, PRECISION));
                applySdagger(&testState, position);
                TEST_ASSERT_TRUE(cvectorAlmostEqual(testVectors[i], testState.vector, dim, \
                                                    PRECISION));
                free(result);
            }
            freeTestVectors(stateVectors, qubits);
        }
        freeTestVectors(testVectors, qubits);
    }
}

/*
 * =================================================================================================
 *                                              test applyH
 * =================================================================================================
 */

void testApplyH(void) {
    dim_t dim;
    cplx_t** testVectors;
    cplx_t** stateVectors;
    state_t testState;
    cplx_t* gateMatrix;
    cplx_t* result;
    for (qubit_t qubits = 1; qubits <= MAXQUBITS; ++qubits) {
        dim = POW2(qubits, dim_t);
        testVectors = generateTestVectors(qubits);
        for (qubit_t position = 0; position < qubits; ++position) {
            stateVectors = generateTestVectors(qubits);
            gateMatrix = singleQubitGateMat(HMAT, qubits, position);
            for (dim_t i = 0; i < dim + 1; ++i) {
               stateInitVector(&testState, stateVectors[i], qubits);
                result = cmatVecMul(gateMatrix, testVectors[i], dim);
                applyH(&testState, position);
                TEST_ASSERT_TRUE(cvectorAlmostEqual(result, testState.vector, dim, PRECISION));
                free(result);
            }
            freeTestVectors(stateVectors, qubits);
        }
        freeTestVectors(testVectors, qubits);
    }
}

/*
 * =================================================================================================
 *                                              test applyT(dagger)
 * =================================================================================================
 */

void testApplyT(void) {
    dim_t dim;
    cplx_t** testVectors;
    cplx_t** stateVectors;
    state_t testState;
    cplx_t* gateMatrix;
    cplx_t* result;
    for (qubit_t qubits = 1; qubits <= MAXQUBITS; ++qubits) {
        dim = POW2(qubits, dim_t);
        testVectors = generateTestVectors(qubits);
        for (qubit_t position = 0; position < qubits; ++position) {
            stateVectors = generateTestVectors(qubits);
            gateMatrix = singleQubitGateMat(TMAT, qubits, position);
            for (dim_t i = 0; i < dim + 1; ++i) {
                stateInitVector(&testState, stateVectors[i], qubits);
                result = cmatVecMul(gateMatrix, testVectors[i], dim);
                applyT(&testState, position);
                TEST_ASSERT_TRUE(cvectorAlmostEqual(result, testState.vector, dim, PRECISION));
                applyTdagger(&testState, position);
                TEST_ASSERT_TRUE(cvectorAlmostEqual(testVectors[i], testState.vector, dim, \
                                                    PRECISION));
                free(result);
            }
            freeTestVectors(stateVectors, qubits);
        }
        freeTestVectors(testVectors, qubits);
    }
}

/*
 * =================================================================================================
 *                                              test applyCX
 * =================================================================================================
 */

void testApplyCX(void) {
    dim_t dim;
    cplx_t** testVectors;
    cplx_t** stateVectors;
    state_t testState;
    cplx_t* gateMatrix;
    cplx_t* result;
    for (qubit_t qubits = 2; qubits <= MAXQUBITS; ++qubits) {
        dim = POW2(qubits, dim_t);
        testVectors = generateTestVectors(qubits);
        for (qubit_t control = 0; control < qubits; ++control) {
            for (qubit_t target = 0; target < qubits; ++target) {
                stateVectors = generateTestVectors(qubits);
                if (control == target) {
                    continue;
                }
                gateMatrix = controlledGateMat(XMAT, qubits, control, target);
                for (dim_t i = 0; i < dim + 1; ++i) {
                    stateInitVector(&testState, stateVectors[i], qubits);
                    result = cmatVecMul(gateMatrix, testVectors[i], dim);
                    applyCX(&testState, control, target);
                    TEST_ASSERT_TRUE(cvectorAlmostEqual(result, testState.vector, dim, PRECISION));
                    free(result);
                }
                freeTestVectors(stateVectors, qubits);
            }
        }
        freeTestVectors(testVectors, qubits);
    }
}

/*
 * =================================================================================================
 *                                              test applyCY
 * =================================================================================================
 */

void testApplyCY(void) {
    dim_t dim;
    cplx_t** testVectors;
    cplx_t** stateVectors;
    state_t testState;
    cplx_t* gateMatrix;
    cplx_t* result;
    for (qubit_t qubits = 2; qubits <= MAXQUBITS; ++qubits) {
        dim = POW2(qubits, dim_t);
        testVectors = generateTestVectors(qubits);
        for (qubit_t control = 0; control < qubits; ++control) {
            for (qubit_t target = 0; target < qubits; ++target) {
                stateVectors = generateTestVectors(qubits);
                if (control == target) {
                    continue;
                }
                gateMatrix = controlledGateMat(YMAT, qubits, control, target);
                for (dim_t i = 0; i < dim + 1; ++i) {
                    stateInitVector(&testState, stateVectors[i], qubits);
                    result = cmatVecMul(gateMatrix, testVectors[i], dim);
                    applyCY(&testState, control, target);
                    TEST_ASSERT_TRUE(cvectorAlmostEqual(result, testState.vector, dim, PRECISION));
                    free(result);
                }
                freeTestVectors(stateVectors, qubits);
            }
        }
        freeTestVectors(testVectors, qubits);
    }
}

/*
 * =================================================================================================
 *                                              test applyCZ
 * =================================================================================================
 */

void testApplyCZ(void) {
    dim_t dim;
    cplx_t** testVectors;
    cplx_t** stateVectors;
    state_t testState;
    cplx_t* gateMatrix;
    cplx_t* result;
    for (qubit_t qubits = 2; qubits <= MAXQUBITS; ++qubits) {
        dim = POW2(qubits, dim_t);
        testVectors = generateTestVectors(qubits);
        for (qubit_t control = 0; control < qubits; ++control) {
            for (qubit_t target = 0; target < qubits; ++target) {
                stateVectors = generateTestVectors(qubits);
                if (control == target) {
                    continue;
                }
                gateMatrix = controlledGateMat(ZMAT, qubits, control, target);
                for (dim_t i = 0; i < dim + 1; ++i) {
                    stateInitVector(&testState, stateVectors[i], qubits);
                    result = cmatVecMul(gateMatrix, testVectors[i], dim);
                    applyCZ(&testState, control, target);
                    TEST_ASSERT_TRUE(cvectorAlmostEqual(result, testState.vector, dim, PRECISION));
                    free(result);
                }
                freeTestVectors(stateVectors, qubits);
            }
        }
        freeTestVectors(testVectors, qubits);
    }
}

/*
 * =================================================================================================
 *                                              test applyCS(dagger)
 * =================================================================================================
 */

void testApplyCS(void) {
    dim_t dim;
    cplx_t** testVectors;
    cplx_t** stateVectors;
    state_t testState;
    cplx_t* gateMatrix;
    cplx_t* result;
    for (qubit_t qubits = 2; qubits <= MAXQUBITS; ++qubits) {
        dim = POW2(qubits, dim_t);
        testVectors = generateTestVectors(qubits);
        for (qubit_t control = 0; control < qubits; ++control) {
            for (qubit_t target = 0; target < qubits; ++target) {
                stateVectors = generateTestVectors(qubits);
                if (control == target) {
                    continue;
                }
                gateMatrix = controlledGateMat(SMAT, qubits, control, target);
                for (dim_t i = 0; i < dim + 1; ++i) {
                    stateInitVector(&testState, stateVectors[i], qubits);
                    result = cmatVecMul(gateMatrix, testVectors[i], dim);
                    applyCS(&testState, control, target);
                    TEST_ASSERT_TRUE(cvectorAlmostEqual(result, testState.vector, dim, PRECISION));
                    applyCSdagger(&testState, control, target);
                    TEST_ASSERT_TRUE(cvectorAlmostEqual(testVectors[i], testState.vector, dim, \
                                                        PRECISION));
                    free(result);
                }
                freeTestVectors(stateVectors, qubits);
            }
        }
        freeTestVectors(testVectors, qubits);
    }
}

/*
 * =================================================================================================
 *                                              test applyCH
 * =================================================================================================
 */

void testApplyCH(void) {
    dim_t dim;
    cplx_t** testVectors;
    cplx_t** stateVectors;
    state_t testState;
    cplx_t* gateMatrix;
    cplx_t* result;
    for (qubit_t qubits = 2; qubits <= MAXQUBITS; ++qubits) {
        dim = POW2(qubits, dim_t);
        testVectors = generateTestVectors(qubits);
        for (qubit_t control = 0; control < qubits; ++control) {
            for (qubit_t target = 0; target < qubits; ++target) {
                stateVectors = generateTestVectors(qubits);
                if (control == target) {
                    continue;
                }
                gateMatrix = controlledGateMat(HMAT, qubits, control, target);
                for (dim_t i = 0; i < dim + 1; ++i) {
                    stateInitVector(&testState, stateVectors[i], qubits);
                    result = cmatVecMul(gateMatrix, testVectors[i], dim);
                    applyCH(&testState, control, target);
                    TEST_ASSERT_TRUE(cvectorAlmostEqual(result, testState.vector, dim, PRECISION));
                    free(result);
                }
                freeTestVectors(stateVectors, qubits);
            }
        }
        freeTestVectors(testVectors, qubits);
    }
}

/*
 * =================================================================================================
 *                                              test applyCT(dagger)
 * =================================================================================================
 */

void testApplyCT(void) {
    dim_t dim;
    cplx_t** testVectors;
    cplx_t** stateVectors;
    state_t testState;
    cplx_t* gateMatrix;
    cplx_t* result;
    for (qubit_t qubits = 2; qubits <= MAXQUBITS; ++qubits) {
        dim = POW2(qubits, dim_t);
        testVectors = generateTestVectors(qubits);
        for (qubit_t control = 0; control < qubits; ++control) {
            for (qubit_t target = 0; target < qubits; ++target) {
                stateVectors = generateTestVectors(qubits);
                if (control == target) {
                    continue;
                }
                gateMatrix = controlledGateMat(TMAT, qubits, control, target);
                for (dim_t i = 0; i < dim + 1; ++i) {
                    stateInitVector(&testState, stateVectors[i], qubits);
                    result = cmatVecMul(gateMatrix, testVectors[i], dim);
                    applyCT(&testState, control, target);
                    TEST_ASSERT_TRUE(cvectorAlmostEqual(result, testState.vector, dim, PRECISION));
                    applyCTdagger(&testState, control, target);
                    TEST_ASSERT_TRUE(cvectorAlmostEqual(testVectors[i], testState.vector, dim, \
                                                        PRECISION));
                    free(result);
                }
                freeTestVectors(stateVectors, qubits);
            }
        }
        freeTestVectors(testVectors, qubits);
    }
}

/*
 * =================================================================================================
 *                                              test applyRX
 * =================================================================================================
 */

void testApplyRX(void) {
    dim_t dim;
    cplx_t** testVectors;
    cplx_t** stateVectors;
    state_t testState;
    cplx_t* gateMatrix;
    cplx_t* result;
    double angle;
    for (qubit_t qubits = 1; qubits <= MAXQUBITS; ++qubits) {
        angle = qubits * 0.25;
        dim = POW2(qubits, dim_t);
        testVectors = generateTestVectors(qubits);
        for (qubit_t position = 0; position < qubits; ++position) {
            stateVectors = generateTestVectors(qubits);
            gateMatrix = RXGateMat(qubits, position, angle);
            for (dim_t i = 0; i < dim + 1; ++i) {
                stateInitVector(&testState, stateVectors[i], qubits);
                result = cmatVecMul(gateMatrix, testVectors[i], dim);
                applyRX(&testState, position, angle);
                TEST_ASSERT_TRUE(cvectorAlmostEqual(result, testState.vector, dim, PRECISION));
                free(result);
                applyRXdagger(&testState, position, angle);
                TEST_ASSERT_TRUE(cvectorAlmostEqual(testVectors[i], testState.vector, dim, \
                                                    PRECISION));
            }
            freeTestVectors(stateVectors, qubits);
        }
        freeTestVectors(testVectors, qubits);
    }
}

/*
 * =================================================================================================
 *                                              test applyRY
 * =================================================================================================
 */

void testApplyRY(void) {
    dim_t dim;
    cplx_t** testVectors;
    cplx_t** stateVectors;
    state_t testState;
    cplx_t* gateMatrix;
    cplx_t* result;
    double angle;
    for (qubit_t qubits = 1; qubits <= MAXQUBITS; ++qubits) {
        angle = qubits * 0.25;
        dim = POW2(qubits, dim_t);
        testVectors = generateTestVectors(qubits);
        for (qubit_t position = 0; position < qubits; ++position) {
            stateVectors = generateTestVectors(qubits);
            gateMatrix = RYGateMat(qubits, position, angle);
            for (dim_t i = 0; i < dim + 1; ++i) {
                stateInitVector(&testState, stateVectors[i], qubits);
                result = cmatVecMul(gateMatrix, testVectors[i], dim);
                applyRY(&testState, position, angle);
                TEST_ASSERT_TRUE(cvectorAlmostEqual(result, testState.vector, dim, PRECISION));
                free(result);
                applyRYdagger(&testState, position, angle);
                TEST_ASSERT_TRUE(cvectorAlmostEqual(testVectors[i], testState.vector, dim, \
                                                    PRECISION));
            }
            freeTestVectors(stateVectors, qubits);
        }
        freeTestVectors(testVectors, qubits);
    }
}

/*
 * =================================================================================================
 *                                              test applyRZ
 * =================================================================================================
 */

void testApplyRZ(void) {
    dim_t dim;
    cplx_t** testVectors;
    cplx_t** stateVectors;
    state_t testState;
    cplx_t* gateMatrix;
    cplx_t* result;
    double angle;
    for (qubit_t qubits = 1; qubits <= MAXQUBITS; ++qubits) {
        angle = qubits * 0.25;
        dim = POW2(qubits, dim_t);
        testVectors = generateTestVectors(qubits);
        for (qubit_t position = 0; position < qubits; ++position) {
            stateVectors = generateTestVectors(qubits);
            gateMatrix = RZGateMat(qubits, position, angle);
            for (dim_t i = 0; i < dim + 1; ++i) {
                stateInitVector(&testState, stateVectors[i], qubits);
                result = cmatVecMul(gateMatrix, testVectors[i], dim);
                applyRZ(&testState, position, angle);
                TEST_ASSERT_TRUE(cvectorAlmostEqual(result, testState.vector, dim, PRECISION));
                free(result);
                applyRZdagger(&testState, position, angle);
                TEST_ASSERT_TRUE(cvectorAlmostEqual(testVectors[i], testState.vector, dim, \
                                                    PRECISION));
            }
            freeTestVectors(stateVectors, qubits);
        }
        freeTestVectors(testVectors, qubits);
    }
}

/*
 * =================================================================================================
 *                                              test SWAP
 * =================================================================================================
 */

/*
This function compares the execution of the SWAP gate of the tested method to the matrix
multiplication outcome.

Input:
    -
Output:
    There is no return value. Only checks whether the tested state vector and the reference state
    vector are equal to some precision.
*/
void testApplySWAP(void) {
    dim_t dim;

    cplx_t** testVectors;       // double complex pointer to an array of pointers defining the set
                                // of reference vectors

    cplx_t** stateVectors;      // double complex pointer to an array of pointers defining the set
                                // of test vectors

    state_t testState;          // custom struct holding the double complex test state vector, the
                                // number of qubits and the Hilbert's space dimension

    cplx_t* gateMatrix;         // double complex pointer to the reference matrix emerging from 
                                // Kronecker products of identity matrices and the Pauli-X

    cplx_t* result;             // double complex pointer holding the result of the reference vector
                                // from matrix vector multiplication

    /*
    Starting with two qubit (SWAP needs at least two qubits), the instances are created up to the
    number MAXQUBITS defined in the macros.
    */
    for (squbit_t qubits = 2; qubits <= MAXQUBITS; ++qubits) {

        dim = POW2(qubits, dim_t);                      // unsigned integer holding the Hilbert
                                                        // space dimension according to the number
                                                        // of qubits                    
        
        testVectors = generateTestVectors(qubits);      // generates the set of reference vectors
                                                        // and attaches them to their pointer
        /*
        Starting with the two rightmost qubits, the swap is executed on each possible pair of
        qubits. The integer for the smaller index can only go up to the second to last qubit,
        because there is no other qubit left to it to be swapped with.
        */
        for (squbit_t right_qubit = 0; right_qubit < qubits - 1 ; ++right_qubit) {
            /*
            Starting at the position left to the right qubit, the larger index goes up to the last
            possible qubit.
            */
            for (squbit_t left_qubit = right_qubit + 1; left_qubit < qubits; ++left_qubit) {
                stateVectors = generateTestVectors(qubits);                 // generates the set of
                                                                            // test vectors and
                                                                            // attaches them to
                                                                            // their pointer
    
                gateMatrix = swapGateMat(qubits, right_qubit, left_qubit);  // initiates the matrix
                                                                            // corresponding to the
                                                                            // Pauli-X gate on the  
                                                                            // current qubit
                
                /*
                Starting at the first reference and test vector, respectively, goes through all test
                vectors checking if they coincide step by step.
                */
                for (dim_t i = 0; i < dim + 1; ++i) {
                    stateInitVector(&testState, stateVectors[i], qubits);   // initialises the test
                                                                            // state with the test
                                                                            // vector corresponding
                                                                            // to the current i

                    result = cmatVecMul(gateMatrix, testVectors[i], dim);   // multiplies the 
                                                                            // reference vector
                                                                            // corresponding to the
                                                                            // current i with the 
                                                                            // gate matrix from 
                                                                            // Kronecker product

                    applySWAP(&testState, right_qubit, left_qubit);          // applies the SWAP
                                                                            // gate according to the
                                                                            // tested method

                    TEST_ASSERT_TRUE(cvectorAlmostEqual(result, testState.vector, dim, PRECISION));
                                                                            // evaluate if the state
                                                                            // vector obtained by
                                                                            // matrix multiplication
                                                                            // coincides with the
                                                                            // one obtained by the
                                                                            // tested method

                    free(result);                                           // free the memory
                                                                            // allocated to the
                                                                            // results pointer for
                                                                            // the next reference
                                                                            // vector
                }
                freeTestVectors(stateVectors, qubits);                      // free all test vectors
                                                                            // altered by the tested
                                                                            // method for the next
                                                                            // pair to be swapped
            }                         
        }
        freeTestVectors(testVectors, qubits);                               // free all reference
                                                                            // vectors used for
                                                                            // matrix vector 
                                                                            // multiplication for 
                                                                            // a higher number of 
                                                                            // qubits
    }
}

/*
 * =================================================================================================
 *                                              test Toffoli
 * =================================================================================================
 */

void testApplyToffoli(void) {
    dim_t dim;

    cplx_t** testVectors;       // double complex pointer to an array of pointers defining the set
                                // of reference vectors

    cplx_t** stateVectors;      // double complex pointer to an array of pointers defining the set
                                // of test vectors

    state_t testState;          // custom struct holding the double complex test state vector, the
                                // number of qubits and the Hilbert's space dimension

    cplx_t* gateMatrix;         // double complex pointer to the reference matrix emerging from 
                                // Kronecker products of identity matrices and the Pauli-X

    cplx_t* result;             // double complex pointer holding the result of the reference vector
                                // from matrix vector multiplication

    /*
    Starting with two qubit (SWAP needs at least two qubits), the instances are created up to the
    number MAXQUBITS defined in the macros.
    */
    for (qubit_t qubits = 3; qubits <= MAXQUBITS; ++qubits) {
        dim = POW2(qubits, dim_t);                      // unsigned integer holding the Hilbert
                                                        // space dimension according to the number
                                                        // of qubits                    
        
        testVectors = generateTestVectors(qubits);      // generates the set of reference vectors
                                                        // and attaches them to their pointer
        /*
        Starting with the three rightmost qubits, the Toffoli is executed on each possible triple of
        qubits.
        */
        for (qubit_t control1 = 0; control1 < qubits; ++control1) {
            for (qubit_t control2 = 0; control2 < qubits; ++control2) {
                for(qubit_t target = 0; target < qubits; ++target) {
                    stateVectors = generateTestVectors(qubits);
                    /*
                    When only one of the three positions coincide continue with the next triple
                    without executing the test.
                    */
                    if (control1 == target || control2 == target || control1 == control2) {
                        continue;
                    }
                    gateMatrix = doublyControlledGateMat(XMAT, qubits, control1, control2, target);
                    /*
                    Starting at the first test and reference vector, respectively, check whether the
                    application of Toffoli by matrix multiplication and qhipster coincide.
                    */
                    for (dim_t i = 0; i < dim + 1; ++i) {
                        stateInitVector(&testState, stateVectors[i], qubits);
                        result = cmatVecMul(gateMatrix, testVectors[i], dim);
                        applyToffoli(&testState, control1, control2, target);

                        TEST_ASSERT_TRUE(cvectorAlmostEqual(result, testState.vector, dim, \
                                         PRECISION));
                        free(result);
                    }
                }
                freeTestVectors(stateVectors, qubits);
            }
        }
        freeTestVectors(testVectors, qubits);
    }
}

/*
 * =================================================================================================
 *                                              main
 * =================================================================================================
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
