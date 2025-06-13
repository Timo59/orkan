//
// Created by Timo Ziegler on 13.03.25.
//
#ifndef TESTEXACT_H
#include "testExact.h"
#endif

/*
 * =====================================================================================================================
 *                                                  testMeanObs
 * =====================================================================================================================
 */
void testMeanDiagObs(void) {
    for (qubit_t qubits = 1; qubits <= MAXQUBITS; qubits++) {
        const dim_t dim = POW2(qubits, dim_t);

        // Define test state and initialize test state vectors
        stateInitEmpty(&testState, qubits);
        vecs = generateTestVectors(qubits);

        for (uint8_t i = 0; i < dim + 1; i++) {
            stateInitVector(&testState, vecs[i]);
            const double test = meanObs(&testState, diagObs);

            double ref = 0;
            for (dim_t j = 0; j < dim; j++) {
                ref += pow(cabs(vecs[i][j]), 2) * diagObs[j];
            }

            TEST_ASSERT_TRUE(ralmostEqual(ref, test, PRECISION));
        }
        freeTestVectors(vecs, qubits);
        vecs = NULL;
        stateFreeVector(&testState);
        testState.vec = NULL;
    }
}

void testMeanObsHerm(void) {
    if (atexit(cleanup)) {
        fprintf(stderr, "Failed to register cleanup function\n");
        exit(EXIT_FAILURE);
    }

    for (qubit_t qubits = 2; qubits <= MAXQUBITS; qubits++) {
        const dim_t dim = POW2(qubits, dim_t);

        // Define test state and initialize test state vectors
        stateInitEmpty(&testState, qubits);
        vecs = generateTestVectors(qubits);

        // Initialize the hermitian operator struct from the composite quantum gates in testExact.h
        testHerm.len = 5;
        testHerm.comp = cg + 5 * (qubits - 2);
        testHerm.weight = dcoeff;

        // Initialize the matrix representation of the hermitian operator
        testObsMat[0] = calloc(dim * dim, sizeof(cplx_t));
        if (testObsMat[0] == NULL) {
            fprintf(stderr, "testMeanObsHerm(): testObsMat[0] allocation failed\n");
            exit(EXIT_FAILURE);
        }
        for (uint8_t j = 0; j < 5; ++j) {
            cplx_t* gateMat = cgMat[j](qubits);
            cscalarMatMulInPlace(dcoeff[j], gateMat, dim);
            cmatAddInPlace(testObsMat[0], gateMat, dim);
            free(gateMat);
        }

        for (uint8_t i = 0; i < dim + 1; i++) {
            stateInitVector(&testState, vecs[i]);
            const double test = meanObs(&testState, &testHerm);

            cplx_t* tmp = malloc(dim * sizeof(cplx_t));
            if (tmp == NULL) {
                exit(EXIT_FAILURE);
            }
            for (dim_t k = 0; k < dim; k++) {
                tmp[k] = vecs[i][k];
            }

            cmatVecMulInPlace(testObsMat[0], tmp, dim);
            const double ref = creal(cInner(tmp, vecs[i], dim));
            free(tmp);

            TEST_ASSERT_TRUE(ralmostEqual(ref, test, PRECISION));
        }
        free(testObsMat[0]);
        testObsMat[0] = NULL;
        freeTestVectors(vecs, qubits);
        vecs = NULL;
        stateFreeVector(&testState);
        testState.vec = NULL;
    }
}

/*
 * =====================================================================================================================
 *                                                  testMeanOp
 * =====================================================================================================================
 */

void testMeanCG(void) {
    if (atexit(cleanup)) {                                          // Register cleanup function
        fprintf(stderr, "Failed to register cleanup function\n");
        exit(EXIT_FAILURE);
    }

    for (qubit_t qubits = 2; qubits <= MAXQUBITS; qubits++) {
        const dim_t dim = POW2(qubits, dim_t);

        // Define test state and initialize test state vectors
        stateInitEmpty(&testState, qubits);
        vecs = generateTestVectors(qubits);

        // Initialize the array of hermitian composite quantum gates serving as observables
        const applyCG *testOp = cg + 5 * (qubits - 2);

        // Iterate the observables and start by defining their matrix representation
        for (uint8_t j = 0; j < 4; ++j) {
            testOpMat = cgMat[j](qubits);

            for (uint8_t i = 0; i < dim + 1; i++) {
                stateInitVector(&testState, vecs[i]);
                const double test = meanCG(&testState, testOp[j]);

                cplx_t* tmp = malloc(dim * sizeof(cplx_t));
                if (tmp == NULL) {
                    exit(EXIT_FAILURE);
                }
                for (dim_t k = 0; k < dim; k++) {
                    tmp[k] = vecs[i][k];
                }

                cmatVecMulInPlace(testOpMat, tmp, dim);
                const cplx_t ref = cInner(tmp, vecs[i], dim);
                free(tmp);

                TEST_ASSERT_TRUE(calmostEqual(ref, test, PRECISION));
            }
            free(testOpMat);
            testOpMat = NULL;
        }
        freeTestVectors(vecs, qubits);
        vecs = NULL;
        stateFreeVector(&testState);
        testState.vec = NULL;
    }
}

/*
 * =====================================================================================================================
 *                                                 testGradientPQC
 * =====================================================================================================================
 */
static double* finDiffMeth(const cplx_t* x, unsigned char d, const cplx_t* H, const matPCG U[], unsigned char parc,
    double par[]);

/*
 * Tests gradPQC with diagonal operator as observable and a PQC comprising evolutions of composite quantum gates
 */
void testGradPQCDiag(void) {
    if (atexit(cleanup)) {
        fprintf(stderr, "Failed to register cleanup function\n");
        exit(EXIT_FAILURE);
    }

    for (qubit_t qubits = 2; qubits < MAXQUBITS; ++qubits) {
        const dim_t dim = POW2(qubits, dim_t);

        // Define test state and initialize test state vectors
        stateInitEmpty(&testState, qubits);
        vecs = generateTestVectors(qubits);

        // Define the test observable as the diagonal entries of a diagonal hermitian operator
        double* testObs = diagObs;

        // Define the matrix representation of the test observable
        testObsMat[0] = diagMat(qubits);

        // Define test PQC and the generators of the unitary transformations
        const applyPCG* testPQC = pqc + 5 * (qubits - 2);
        const applyCG* testGen = cg + 5 * (qubits - 2);

        // Define matrix representations of the test PQC
        const matPCG* testPQCMat = pqcMat;

        /* TESTING */
        for (dim_t i = 0; i < dim + 1; ++i) {
            stateInitVector(&testState, vecs[i]);

            // Ascertain the exact gradient by tested method
            testVec = gradPQC(&testState, testObs, 5, testPQC, randPar, testGen);

            // Ascertain the reference gradient by finite difference method
            refVec = finDiffMeth(vecs[i], qubits, testObsMat[0], testPQCMat, 5, randPar);

            // Check if the entries of both gradients are at most 1e-4 apart
            TEST_ASSERT_TRUE(rvectorAlmostEqual(refVec, testVec, 5, APPROXPRECISION));

            free(testVec);
            testVec = NULL;
            free(refVec);
            refVec = NULL;
        }

        free(testObsMat[0]);
        testObsMat[0] = NULL;
        freeTestVectors(vecs, qubits);
        vecs = NULL;
        stateFreeVector(&testState);
        testState.vec = NULL;
    }
}

/*
 * Tests gradPQC with general hermitian operator as observable and a PQC comprising evolutions of composite quantum
 * gates
 */
void testGradPQCHerm(void) {
    // Register cleanup function
    if (atexit(cleanup)) {
        fprintf(stderr, "Failed to register cleanup function\n");
        exit(EXIT_FAILURE);
    }

    for (qubit_t qubits = 2; qubits < MAXQUBITS; ++qubits) {
        const dim_t dim = POW2(qubits, dim_t);

        // Define test state and initialize test state vectors
        stateInitEmpty(&testState, qubits);
        vecs = generateTestVectors(qubits);

        // Define the test observable as a general hermitian operator
        herm_t testObs;
        testObs.len = 5;                                // xi, yi, zi, swapi, diag for i = qubits defined in testExact.h
        testObs.comp = cg + testObs.len * (qubits - 2);
        testObs.weight = dcoeff;

        // Define the matrix representation of the test observable
        testObsMat[0] = calloc(dim * dim, sizeof(cplx_t));
        if (testObsMat[0] == NULL) {
            fprintf(stderr, "testGradPQCHerm(): testObsMat[0] allocation failed\n");
            exit(EXIT_FAILURE);
        }
        for (uint8_t j = 0; j < 5; ++j) {
            cplx_t* gateMat = cgMat[j](qubits);
            cscalarMatMulInPlace(dcoeff[j], gateMat, dim);
            cmatAddInPlace(testObsMat[0], gateMat, dim);
            free(gateMat);
        }

        // Define test PQC and the generators of the unitary transformations
        const applyPCG* testPQC = pqc + 5 * (qubits - 2);
        const applyCG* testGen = cg + 5 * (qubits - 2);

        // Define matrix representations of the test PQC
        const matPCG* testPQCMat = pqcMat;

        /* TESTING */
        for (dim_t i = 0; i < dim + 1; ++i) {
            stateInitVector(&testState, vecs[i]);

            // Ascertain the exact gradient by tested method
            testVec = gradPQC(&testState, &testObs, 5, testPQC, randPar, testGen);

            // Ascertain the reference gradient by finite difference method
            refVec = finDiffMeth(vecs[i], qubits, testObsMat[0], testPQCMat, 5, randPar);

            // Check if the entries of both gradients are at most 1e-4 apart
            TEST_ASSERT_TRUE(rvectorAlmostEqual(refVec, testVec, 5, APPROXPRECISION));

            free(testVec);
            testVec = NULL;
            free(refVec);
            refVec = NULL;
        }

        free(testObsMat[0]);
        testObsMat[0] = NULL;
        freeTestVectors(vecs, qubits);
        vecs = NULL;
        stateFreeVector(&testState);
        testState.vec = NULL;
    }
}

/*
 * Tests gradPQC with diagonal operator as observable and a PQC comprising a composition of parametrized quantum gates
 */
void testGradPQCDiag2(void) {
    if (atexit(cleanup)) {
        fprintf(stderr, "Failed to register cleanup function\n");
        exit(EXIT_FAILURE);
    }

    for (qubit_t qubits = 2; qubits < MAXQUBITS; ++qubits) {
        const dim_t dim = POW2(qubits, dim_t);

        // Define test state and initialize test state vectors
        stateInitEmpty(&testState, qubits);
        vecs = generateTestVectors(qubits);

        // Define the test observable as the diagonal entries of a diagonal hermitian operator
        double* testObs = diagObs;

        // Define the matrix representation of the test observable
        testObsMat[0] = diagMat(qubits);

        // Define test PQC and the generators of the unitary transformations
        const applyPCG* testPQC = rpqc + 4 * (qubits - 2);
        const applyCG* testGen = gen + 4 * (qubits - 2);

        // Define matrix representations of the test PQC
        const matPCG* testPQCMat = rpqcMat;

        /* TESTING */
        for (dim_t i = 0; i < dim + 1; ++i) {
            stateInitVector(&testState, vecs[i]);

            // Ascertain the exact gradient by tested method
            testVec = gradPQC(&testState, testObs, 4, testPQC, randPar, testGen);

            // Ascertain the reference gradient by finite difference method
            refVec = finDiffMeth(vecs[i], qubits, testObsMat[0], testPQCMat, 4, randPar);

            // Check if the entries of both gradients are at most 1e-4 apart
            TEST_ASSERT_TRUE(rvectorAlmostEqual(refVec, testVec, 1, APPROXPRECISION));

            free(testVec);
            testVec = NULL;
            free(refVec);
            refVec = NULL;
        }

        free(testObsMat[0]);
        testObsMat[0] = NULL;
        freeTestVectors(vecs, qubits);
        vecs = NULL;
        stateFreeVector(&testState);
        testState.vec = NULL;
    }
}

/*
 * Tests gradPQC with general hermitian operator  as observable and a PQC comprising a composition of parametrized
 * quantum gates
 */
void testGradPQCHerm2(void) {
    if (atexit(cleanup)) {
        fprintf(stderr, "Failed to register cleanup function\n");
        exit(EXIT_FAILURE);
    }

    for (qubit_t qubits = 2; qubits < MAXQUBITS; ++qubits) {
        const dim_t dim = POW2(qubits, dim_t);

        // Define test state and initialize test state vectors
        stateInitEmpty(&testState, qubits);
        vecs = generateTestVectors(qubits);

        // Define the test observable as a general hermitian operator
        herm_t testObs;
        testObs.len = 5;                                // xi, yi, zi, swapi, diag for i = qubits defined in testExact.h
        testObs.comp = cg + testObs.len * (qubits - 2);
        testObs.weight = dcoeff;

        // Define the matrix representation of the test observable
        testObsMat[0] = calloc(dim * dim, sizeof(cplx_t));
        if (testObsMat[0] == NULL) {
            fprintf(stderr, "testGradPQCHerm(): testObsMat[0] allocation failed\n");
            exit(EXIT_FAILURE);
        }
        for (uint8_t j = 0; j < 5; ++j) {
            cplx_t* gateMat = cgMat[j](qubits);
            cscalarMatMulInPlace(dcoeff[j], gateMat, dim);
            cmatAddInPlace(testObsMat[0], gateMat, dim);
            free(gateMat);
        }

        // Define test PQC and the generators of the unitary transformations
        const applyPCG* testPQC = rpqc + 4 * (qubits - 2);
        const applyCG* testGen = gen + 4 * (qubits - 2);

        // Define matrix representations of the test PQC
        const matPCG* testPQCMat = rpqcMat;

        /* TESTING */
        for (dim_t i = 0; i < dim + 1; ++i) {
            stateInitVector(&testState, vecs[i]);

            // Ascertain the exact gradient by tested method
            testVec = gradPQC(&testState, &testObs, 4, testPQC, randPar, testGen);

            // Ascertain the reference gradient by finite difference method
            refVec = finDiffMeth(vecs[i], qubits, testObsMat[0], testPQCMat, 4, randPar);

            // Check if the entries of both gradients are at most 1e-4 apart
            TEST_ASSERT_TRUE(rvectorAlmostEqual(refVec, testVec, 1, APPROXPRECISION));

            free(testVec);
            testVec = NULL;
            free(refVec);
            refVec = NULL;
        }

        free(testObsMat[0]);
        testObsMat[0] = NULL;
        freeTestVectors(vecs, qubits);
        vecs = NULL;
        stateFreeVector(&testState);
        testState.vec = NULL;
    }
}

/*
 * =====================================================================================================================
 *                                                 testMMseq
 * =====================================================================================================================
 */
void testMMseq(void) {
    if (atexit(cleanup)) {
        fprintf(stderr, "Failed to register cleanup function\n");
        exit(EXIT_FAILURE);
    }

    for (qubit_t qubits = 2; qubits < MAXQUBITS; ++qubits) {
        const dim_t dim = POW2(qubits, dim_t);

        // Define test state
        stateInitEmpty(&testState, qubits);

        // Define the test observables, after initialising the test hermitian operator
        testHerm.len = 5;
        testHerm.comp = cg + 5 * (qubits - 2);
        testHerm.weight = dcoeff;
        const applyCG* testObs = obs;

        // Define the matrix representation of the test observables
        testObsMat[0] = diagMat(qubits);
        testObsMat[1] = calloc(dim * dim, sizeof(cplx_t));
        if (testObsMat[1] == NULL) {
            fprintf(stderr, "testMMseq(): testObsMat[1] allocation failed\n");
            exit(EXIT_FAILURE);
        }
        for (uint8_t j = 0; j < 5; ++j) {
            cplx_t* gateMat = cgMat[j](qubits);
            cscalarMatMulInPlace(dcoeff[j], gateMat, dim);
            cmatAddInPlace(testObsMat[1], gateMat, dim);
            free(gateMat);
        }

        // Define the sequence of LCUs
        lccg_t testLCU[3];
        for (unsigned int i = 0; i < 3; ++i) {
            testLCU[i].len = 5;
            testLCU[i].comp = unitary + 15 * (qubits - 2) + i * 5;
            testLCU[i].weight = zcoeff + i * 5;
        }

        // Initialize the matrix representations of LCU components
        for (uint8_t i = 0; i < 15; ++i) {
            testLCUMat[i] = unitaryMat[i](qubits);
        }

        // Iterate the channels of the sequence, each time calculating the moment matrix of the specified LCU
        for (uint8_t j = 0; j < 3; ++j) {
            // Initialize the test state vectors
            vecs = generateTestVectors(qubits);

            for (dim_t i = 0; i < dim + 1; ++i) {
                stateInitVector(&testState, vecs[i]);

                // Define the variable to hold the moment matrices
                for (uint8_t k = 0; k < 2; ++k) {
                    if ((testMomMat[k] = malloc(5 * 5 * sizeof (cplx_t))) == NULL) {
                        fprintf(stderr, "testMMseq(): testMomMat[%d] allocation failed\n", k);
                        exit(EXIT_FAILURE);
                    }
                }
                mmseq(&testState, 2, testObs, 3, testLCU, j, testMomMat);

                for (uint8_t k = 0; k < 2; ++k) {
                    if ((testMat[k] = malloc(5 * 5 * sizeof (cplx_t))) == NULL) {
                        fprintf(stderr, "testMMseq(): testMat[%d] allocation failed\n", k);
                        exit(EXIT_FAILURE);
                    }

                    if ((refMat[k] = malloc(5 * 5 * sizeof (cplx_t))) == NULL) {
                        fprintf(stderr, "testMMseq(): refMat[%d] allocation failed\n", k);
                        exit(EXIT_FAILURE);
                    }
                }

                // Bring testMomMat from column major format to row major
                for (uint8_t k = 0; k < 5; ++k) {
                    for (uint8_t mat = 0; mat < 2; ++mat) {
                        testMat[mat][k + k * 5] = testMomMat[mat][k + k * 5];
                        for (uint8_t l = 0; l < 5; ++l) {
                            testMat[mat][l + 5 * k] = testMomMat[mat][k + 5 * l];
                        }
                    }
                }
                // Free testMomMat as it is not needed anymore
                for (uint8_t k = 0; k < 2; ++k) {
                    free(testMomMat[k]);
                    testMomMat[k] = NULL;
                }

                // Apply the matrix representations of LCU channels prior to the j-th channel to the reference vector
                for (uint8_t k = 0; k < j; ++k) {
                    cplx_t* sum = calloc(dim, sizeof (cplx_t));
                    if (sum == NULL) {
                        fprintf(stderr, "testMMseq(): %d-th sum allocation prior to link failed\n", k);
                        exit(EXIT_FAILURE);
                    }
                    for (uint8_t l = 0; l < 5; ++l) {
                        cplx_t* tmp = cmatVecMul(testLCUMat[k * 5 + l], vecs[i], dim);
                        cscalarVecMulInPlace(zcoeff[k * 5 + l], tmp, dim);
                        cvecAddInPlace(sum, tmp, dim);
                        free(tmp);
                    }
                    free(vecs[i]);
                    vecs[i] = sum;
                }

                // Calculate the vectors spanning the reference matrix; i.e., the LCU components applied to the input
                for (uint8_t k = 0; k < 5; ++k) {
                    if ((vecsCopy[k] = malloc(dim * sizeof (cplx_t))) == NULL) {
                        fprintf(stderr, "testMMseq(): vecsCopy[%d] allocation failed\n", k);
                        exit(EXIT_FAILURE);
                    }

                    for (dim_t l = 0; l < dim; ++l) {
                        vecsCopy[k][l] = vecs[i][l];
                    }
                    printf("refVec[%d] =", k);
                    vectorPrint(vecsCopy[k], dim);

                    // Apply respective search unitary
                    cmatVecMulInPlace(testLCUMat[j * 5 + k], vecsCopy[k], dim);
                    printf("testLCUMat[%d] =\n", j * 5 + k);
                    matrixPrint(testLCUMat[j * 5 + k], dim);

                    // Apply the remaining LCUs to vecsCopy[k]
                    for (uint8_t l = j + 1; l < 3; ++l) {
                        cplx_t* sum = calloc(dim, sizeof (cplx_t));
                        if (sum == NULL) {
                            fprintf(stderr, "testMMseq: %d-th sum allocation after link failed\n", l);
                            exit(EXIT_FAILURE);
                        }

                        for (uint8_t m = 0; m < 5; ++m) {
                            cplx_t* tmp = cmatVecMul(testLCUMat[j * 5 + m], vecsCopy[k], dim);
                            cscalarVecMulInPlace(zcoeff[l * 5 + m], tmp, dim);
                            cvecAddInPlace(sum, tmp, dim);
                            free(tmp);
                        }
                        free(vecsCopy[k]);
                        vecsCopy[k] = sum;
                    } // for LCU channels
                    printf("refVec[%d] =", k);
                    vectorPrint(vecsCopy[k], dim);
                } // for reference vector

                // Calculate the matrix elements of the reference moment matrix
                for (uint8_t k = 0; k < 5; ++k) {
                    for (uint8_t m = 0; m < 2; ++m) {
                        cplx_t* tmp = cmatVecMul(testObsMat[m], vecsCopy[k], dim);

                        // Iterate columns of the reference moment matrix
                        for (uint8_t l = 0; l < 5; ++l) {
                            refMat[m][l + k * 5] = cInner(tmp, vecsCopy[l], dim);
                        } // for column of moment matrix

                        free(tmp);
                    } // for moment matrix
                } // for row of reference moment matrix

                for (uint8_t k = 0; k < 2; ++k) {
                    printf("test :\n");
                    cmatrixPrint(testMat[k], dim);
                    printf("ref :\n");
                    cmatrixPrint(refMat[k], dim);
                    TEST_ASSERT_TRUE(cvectorAlmostEqual(refMat[k], testMat[k], 25, APPROXPRECISION));
                }

                for (uint8_t k = 0; k < 2; ++k) {
                    free(refMat[k]);
                    refMat[k] = NULL;
                    free(testMat[k]);
                    testMat[k] = NULL;
                }
                for (uint8_t k = 0; k < 5; ++k) {
                    free(vecsCopy[k]);
                    vecsCopy[k] = NULL;
                }
            } // for test vector
            freeTestVectors(vecs, qubits);
            vecs = NULL;
        } // for link
        for (unsigned char i = 0; i < 15; ++i) {
            free(testLCUMat[i]);
            testLCUMat[i] = NULL;
        }
        for (unsigned char i = 0; i < 2; ++i) {
            free(testObsMat[i]);
            testObsMat[i] = NULL;
        }
        stateFreeVector(&testState);
        testState.vec = NULL;
    } // for qubits
}

/*
 * =====================================================================================================================
 *                                                      main
 * =====================================================================================================================
 */
int main(void) {
    UNITY_BEGIN();
//    RUN_TEST(testMeanDiagObs);
//    RUN_TEST(testMeanObsHerm);
//    RUN_TEST(testMeanCG);
//    RUN_TEST(testGradPQCDiag);
//    RUN_TEST(testGradPQCHerm);
//    RUN_TEST(testGradPQCDiag2);
//    RUN_TEST(testGradPQCHerm2);
    RUN_TEST(testMMseq);
    return UNITY_END();
}

/*
 * =====================================================================================================================
 *                                                 Helper functions
 * =====================================================================================================================
 */
/*
 * @brief   Returns the gradient of the quadratic form associated with the hermitian matrix using finite difference
 * method; i.e.,
 *                                                  <Ux|H|Ux>
 *
 * @param[in,out] x : Vector transformed with the unitary matrix
 * @param[in] d     : Logarithm of the vector space dimension; i.e., dim = 2**d
 * @param[in,out] H : Hermitian matrix
 * @param[in,out] U : Sequence of parametrized unitary transformation
 * @param[in] par   : Parameter setting
 *
 * @returns The approximate gradient of the quadratic form of H with respect to x after performing a sequence of
 *          parametrized unitary transformations at the parameter setting
 */
static double* finDiffMeth(const cplx_t* x, const unsigned char d, const cplx_t* H, const matPCG U[],
                           const unsigned char parc, double par[]) {
    const unsigned long dim = POW2(d, unsigned long);

    double* out = malloc(parc * sizeof(double));
    if (out == NULL) {
        fprintf(stderr, "finDiffMeth(): out allocation failed\n");
        exit(EXIT_FAILURE);
    }

    // Copy the input vector to one temporary to conduct matrix multiplications to
    cplx_t *bra, *ket;
    if ((bra = malloc(dim * sizeof(cplx_t))) == NULL) {
        fprintf(stderr, "finDiffMeth(): bra allocation failed\n");
        free(out);
        exit(EXIT_FAILURE);
    }
    if ((ket = malloc(dim * sizeof(cplx_t))) == NULL) {
        fprintf(stderr, "finDiffMeth(): ket allocation failed\n");
        free(out);
        free(bra);
        exit(EXIT_FAILURE);
    }
    for (unsigned long i = 0; i < dim; ++i) {
        ket[i] = x[i];
    }

    // Initialize the unitary transformation matrices for the input parameters, conduct matrix vector multiplication
    // with ket and free the matrices afterward
    cplx_t* uMat;
    for (unsigned char i = 0; i < parc; ++i) {
        uMat = U[i](d, par[i]);
        cmatVecMulInPlace(uMat, ket, dim);
        free(uMat);
    }

    // Copy ket to bra, conduct matrix vector multiplication of H with ket and calculate the quadratic form by inner
    // product
    for (unsigned long i = 0; i < dim; ++i) {
        bra[i] = ket[i];
    }
    cmatVecMulInPlace(H, ket, dim);
    double current;
    if (cimag(cInner(bra, ket, dim)) < 1e-9) {
        current = creal(cInner(bra, ket, dim));
    } else {
        fprintf(stderr, "finDiffMeth(): Quadratic form is off being real\n");
        free(out);
        free(bra);
        free(ket);
        exit(EXIT_FAILURE);
    }

    // For each parameter entry shift the value by EPSILON, compute the quadratic form at that point and add the
    // difference with the current value to the return pointer
    for (unsigned char i = 0; i < parc; ++i) {
        par[i] += EPSILON;

        // Copy input vector to ket
        for (unsigned long j = 0; j < dim; ++j) {
            ket[j] = x[j];
        }

        // Conduct matrix vector multiplication of the shifted unitary transformations and ket
        for (unsigned char j = 0; j < parc; ++j) {
            uMat = U[j](d, par[j]);
            cmatVecMulInPlace(uMat, ket, dim);
            free(uMat);
        }

        // Copy ket to bra, conduct matrix vector multiplication of H with ket and calculate the quadratic
        for (unsigned long j = 0; j < dim; ++j) {
            bra[j] = ket[j];
        }
        cmatVecMulInPlace(H, ket, dim);
        double shifted;
        if (cimag(cInner(bra, ket, dim)) < 1e-9) {
            shifted = creal(cInner(bra, ket, dim));
        } else {
            fprintf(stderr, "finDiffMeth(): Quadratic form is off being real\n");
            free(out);
            free(bra);
            free(ket);
            exit(EXIT_FAILURE);
        }

        // Add the difference quotient to the return pointer and reset the parameter entry
        out[i] = (1. / EPSILON) * (shifted - current);
        par[i] -= EPSILON;
    }

    free(bra);
    free(ket);
    return out;
}
