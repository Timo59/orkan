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
        applyCG *testObs = cg + 5 * (qubits - 2);

        // Iterate the observables and start by defining their matrix representation
        for (uint8_t j = 0; j < 4; ++j) {
            testObsMat = cgMat[j](qubits);

            for (uint8_t i = 0; i < dim + 1; i++) {
                stateInitVector(&testState, vecs[i]);
                const double test = meanObs(&testState, testObs[j]);

                cplx_t* tmp = malloc(dim * sizeof(cplx_t));
                if (tmp == NULL) {
                    exit(EXIT_FAILURE);
                }
                for (dim_t k = 0; k < dim; k++) {
                    tmp[k] = vecs[i][k];
                }

                cmatVecMulInPlace(testObsMat, tmp, dim);
                const double ref = creal(cInner(tmp, vecs[i], dim));
                free(tmp);

                TEST_ASSERT_TRUE(ralmostEqual(ref, test, PRECISION));
            }
            free(testObsMat);
            testObsMat = NULL;
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
        testObsMat = calloc(dim * dim, sizeof(cplx_t));
        if (testObsMat == NULL) {
            fprintf(stderr, "testMeanObsHerm(): testObsMat allocation failed\n");
            exit(EXIT_FAILURE);
        }
        for (uint8_t j = 0; j < 5; ++j) {
            cplx_t* gateMat = cgMat[j](qubits);
            cscalarMatMulInPlace(dcoeff[j], gateMat, dim);
            cmatAddInPlace(testObsMat, gateMat, dim);
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

            cmatVecMulInPlace(testObsMat, tmp, dim);
            const double ref = creal(cInner(tmp, vecs[i], dim));
            free(tmp);

            TEST_ASSERT_TRUE(ralmostEqual(ref, test, PRECISION));
        }
        free(testObsMat);
        testObsMat = NULL;
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
        testObsMat = diagMat(qubits);

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
            refVec = finDiffMeth(vecs[i], qubits, testObsMat, testPQCMat, 5, randPar);

            // Check if the entries of both gradients are at most 1e-4 apart
            TEST_ASSERT_TRUE(rvectorAlmostEqual(refVec, testVec, 5, APPROXPRECISION));

            free(testVec);
            testVec = NULL;
            free(refVec);
            refVec = NULL;
        }

        free(testObsMat);
        testObsMat = NULL;
        freeTestVectors(vecs, qubits);
        vecs = NULL;
        stateFreeVector(&testState);
        testState.vec = NULL;
    }
}

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
        testObsMat = calloc(dim * dim, sizeof(cplx_t));
        if (testObsMat == NULL) {
            fprintf(stderr, "testGradPQCHerm(): testObsMat allocation failed\n");
            exit(EXIT_FAILURE);
        }
        for (uint8_t j = 0; j < 5; ++j) {
            cplx_t* gateMat = cgMat[j](qubits);
            cscalarMatMulInPlace(dcoeff[j], gateMat, dim);
            cmatAddInPlace(testObsMat, gateMat, dim);
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
            refVec = finDiffMeth(vecs[i], qubits, testObsMat, testPQCMat, 5, randPar);

            // Check if the entries of both gradients are at most 1e-4 apart
            TEST_ASSERT_TRUE(rvectorAlmostEqual(refVec, testVec, 5, APPROXPRECISION));

            free(testVec);
            testVec = NULL;
            free(refVec);
            refVec = NULL;
        }

        free(testObsMat);
        testObsMat = NULL;
        freeTestVectors(vecs, qubits);
        vecs = NULL;
        stateFreeVector(&testState);
        testState.vec = NULL;
    }
}

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
        testObsMat = diagMat(qubits);

        // Define test PQC and the generators of the unitary transformations
        const applyPCG* testPQC = rpqc + 4 * (qubits - 2);
        const applyCG* testGen = gen + 4 * (qubits - 2);

        // Define matrix representations of the test PQC
        const matPCG* testPQCMat = rpqcMat;

        /* TESTING */
        for (dim_t i = dim; i < dim + 1; ++i) {
            stateInitVector(&testState, vecs[i]);

            // Ascertain the exact gradient by tested method
            testVec = gradPQC(&testState, testObs, 1, testPQC, randPar, testGen);

            // Ascertain the reference gradient by finite difference method
            refVec = finDiffMeth(vecs[i], qubits, testObsMat, testPQCMat, 1, randPar);

            printf("qubits = %d\n", qubits);
            printf("test = ");
            vectorPrint(testVec, 1);
            printf("ref = ");
            vectorPrint(refVec, 1);

            // Check if the entries of both gradients are at most 1e-4 apart
            TEST_ASSERT_TRUE(rvectorAlmostEqual(refVec, testVec, 1, APPROXPRECISION));

            free(testVec);
            testVec = NULL;
            free(refVec);
            refVec = NULL;
        }

        free(testObsMat);
        testObsMat = NULL;
        freeTestVectors(vecs, qubits);
        vecs = NULL;
        stateFreeVector(&testState);
        testState.vec = NULL;
    }
}

// /*
//  * =====================================================================================================================
//  *                                                 testMMseq
//  * =====================================================================================================================
//  */
// state_t testState = {                                               // State for the method to be tested
//     .qubits = 0,
//     .dim = 0,
//     .vec = NULL
// };
// cplx_t *observableMat[3] = {NULL, NULL, NULL};                      // Matrix representations of observables
// cplx_t** uMat[3] = {NULL, NULL, NULL};                              // Matrix representations of the channels' unitaries
// cplx_t** vecs = NULL;                                               // Test state vectors;
// cplx_t* momMat[3] = {NULL, NULL, NULL};                             // Moment matrices as the result of mmseq
// cplx_t *test[3] = {NULL, NULL, NULL};                               // Moment matrices in dense row major form
// cplx_t *ref[3] = {NULL, NULL, NULL};                                // Reference from matrix multiplication
// cplx_t *refVec[5] = {NULL, NULL, NULL, NULL, NULL};                 // Reference state vecs to determine moment matrix
//
// void cleanupTestMMseq(void) {
//     if (testState.vec != NULL) {
//         stateFreeVector(&testState);
//     }
//     for (uint8_t i = 0; i < 3; ++i) {
//         if (observableMat[i] != NULL) {
//             free(observableMat[i]);
//         }
//         if (uMat[i] != NULL) {
//             for (uint8_t j = 0; j < 5; ++j) {
//                 if (uMat[i][j] != NULL) {
//                     free(uMat[i][j]);
//                 }
//             }
//             free(uMat[i]);
//         }
//         if (momMat[i] != NULL) {
//             free(momMat[i]);
//         }
//         if (test[i] != NULL) {
//             free(test[i]);
//         }
//         if (ref[i] != NULL) {
//             free(ref[i]);
//         }
//     }
//     if (vecs != NULL) {
//         freeTestVectors(vecs, testState.qubits);
//     }
//     for (uint8_t i = 0; i < 5; ++i) {
//         if (refVec[i] != NULL) {
//             free(refVec[i]);
//         }
//     }
//
// }
//
// void testMMseq(void) {
//     if (atexit(cleanupTestMMseq)) {                                 // Register cleanup function
//         fprintf(stderr, "Failed to register cleanup function\n");
//         exit(EXIT_FAILURE);
//     }
//
//     for (qubit_t qubits = 2; qubits < MAXQUBITS; ++qubits) {        // Iterate number of qubits
//         const dim_t dim = POW2(qubits, dim_t);
//         stateInitEmpty(&testState, qubits);                         // Define the quantum state
//
//         const applyQB* observable = obs + 3 * (qubits - 2);         // diagObs, all-Y and rotational SWAP blocks
//
//         for (uint8_t i = 0; i < 3; ++i) {                           // Matrix representation of observables
//             observableMat[i] = obsMat[i](qubits);
//         }
//
//         const applyQB* u = channel + 15 * (qubits - 2);             // Concatenation of search unitaries
//
//         for (uint8_t i = 0; i < 3; ++i) {
//             if ((uMat[i] = malloc(5 * sizeof (cgMat))) == NULL) {
//                 fprintf(stderr, "testMMseq: uMat[%d] allocation failed\n", i);
//                 exit(EXIT_FAILURE);
//             }
//             for (uint8_t j = 0; j < 5; ++j) {                       // Initilialize the matrix representations of the
//                 uMat[i][j] = channelMat[i][j](qubits);              // channels' unitaries
//             }
//         }
//
//         for (uint8_t j = 0; j < 3; ++j) {                           // Iterate the channels of the sequence
//             vecs = generateTestVectors(qubits);
//             for (dim_t i = 0; i < dim + 1; ++i) {                   // Iterate the test vectors
//                 stateInitVector(&testState, vecs[i]);
//
//                 for (uint8_t k = 0; k < 3; ++k) {
//                     if ((momMat[k] = malloc(5 * 5 * sizeof (cplx_t))) == NULL) {
//                         fprintf(stderr, "testMMseq: momMat[%d] allocation failed\n", k);
//                         exit(EXIT_FAILURE);
//                     }
//                 }
//                 mmseq(&testState, 3, observable, 3, zcoeff, 5, u, j, momMat);
//
//                 for (uint8_t k = 0; k < 3; ++k) {
//                     if ((test[k] = malloc(5 * 5 * sizeof (cplx_t))) == NULL) {
//                         fprintf(stderr, "testMMseq: test[%d] allocation failed\n", k);
//                         exit(EXIT_FAILURE);
//                     }
//
//                     if ((ref[k] = malloc(5 * 5 * sizeof (cplx_t))) == NULL) {
//                         fprintf(stderr, "testMMseq: ref[%d] allocation failed\n", k);
//                         exit(EXIT_FAILURE);
//                     }
//                 }
//                 for (uint8_t k = 0; k < 5; ++k) {                   // Iterate the columns of the column major form
//                     for (uint8_t mat = 0; mat < 3; ++mat) {         // Iterate the moment matrices
//                         test[mat][k + k * 5] = momMat[mat][k + k * 5];  // Copy the diagonal elements
//                         for (uint8_t l = 0; l < 5; ++l) {               // Iterate rows of the lower triangle
//                             test[mat][l + 5 * k] = momMat[mat][k + 5 * l];
//                         }
//                     }
//                 }
//                 for (uint8_t k = 0; k < 3; ++k) {
//                     free(momMat[k]);
//                     momMat[k] = NULL;
//                 }
//
//                 for (uint8_t k = 0; k < j; ++k) {                   // Apply the matrix representation of LCU channels
//                     cplx_t* sum = calloc(dim, sizeof (cplx_t));     // prior to the j-th channel to the reference vector
//                     if (sum == NULL) {
//                         fprintf(stderr, "testMMseq: %d-th sum allocation prior to link failed\n", k);
//                         exit(EXIT_FAILURE);
//                     }
//                     for (uint8_t l = 0; l < 5; ++l) {
//                         cplx_t* tmp = cmatVecMul(uMat[k][l], vecs[i], dim);
//                         cscalarVecMulInPlace(zcoeff[k * 5 + l], tmp, dim);
//                         cvecAddInPlace(sum, tmp, dim);
//                         free(tmp);
//                     }
//                     free(vecs[i]);
//                     vecs[i] = sum;
//                 }
//
//                 for (uint8_t k = 0; k < 5; ++k) {                   // Calculate the reference vectors spanning the
//                     if ((refVec[k] = malloc(dim * sizeof (cplx_t))) == NULL) {  // moment matrix
//                         fprintf(stderr, "testMMseq: refVec[%d] allocation failed\n", k);
//                         exit(EXIT_FAILURE);
//                     }
//
//                     for (dim_t l = 0; l < dim; ++l) {               // Copy current reference state vector
//                         refVec[k][l] = vecs[i][l];
//                     }
//
//                     cmatVecMulInPlace(uMat[j][k], refVec[k], dim);  // Apply the according search unitary
//
//                     for (uint8_t l = j + 1; l < 3; ++l) {           // Apply the remaining LCU channels to refVec[k]
//                         cplx_t* sum = calloc(dim, sizeof (cplx_t));
//                         if (sum == NULL) {
//                             fprintf(stderr, "testMMseq: %d-th sum allocation after link failed\n", l);
//                             exit(EXIT_FAILURE);
//                         }
//
//                         for (uint8_t m = 0; m < 5; ++m) {
//                             cplx_t* tmp = cmatVecMul(uMat[l][m], refVec[k], dim);
//                             cscalarVecMulInPlace(zcoeff[l * 5 + m], tmp, dim);
//                             cvecAddInPlace(sum, tmp, dim);
//                             free(tmp);
//                         }
//                         free(refVec[k]);
//                         refVec[k] = sum;
//                     } // for LCU channels
//                 } // for reference vector
//
//
//                 for (uint8_t k = 0; k < 5; ++k) {                   // Iterate the rows of the reference moment matrices
//                     for (uint8_t m = 0; m < 3; ++m) {               // Iterate the moment matrices
//                         cplx_t* tmp = cmatVecMul(observableMat[m], refVec[k], dim);
//
//                         for (uint8_t l = 0; l < 5; ++l) {               // Iterate columns of the reference moment matrix
//                             ref[m][l + k * 5] = cInner(tmp, refVec[l], dim);
//                         } // for column of moment matrix
//
//                         free(tmp);
//                     } // for moment matrix
//                 } // for row of reference moment matrix
//
//                 for (uint8_t k = 0; k < 3; ++k) {
//                     TEST_ASSERT_TRUE(cvectorAlmostEqual(ref[k], test[k], 25, APPROXPRECISION));
//                 }
//
//                 for (uint8_t k = 0; k < 3; ++k) {
//                     free(ref[k]);
//                     ref[k] = NULL;
//                     free(test[k]);
//                     test[k] = NULL;
//                 }
//                 for (uint8_t k = 0; k < 5; ++k) {
//                     free(refVec[k]);
//                     refVec[k] = NULL;
//                 }
//             } // for test vector
//             freeTestVectors(vecs, qubits);
//             vecs = NULL;
//         } // for link
//         for (uint8_t k = 0; k < 3; ++k) {
//             for (uint8_t l = 0; l < 5; ++l) {
//                 free(uMat[k][l]);
//                 uMat[k][l] = NULL;
//             }
//             free(uMat[k]);
//             uMat[k] = NULL;
//             free(observableMat[k]);
//             observableMat[k] = NULL;
//         }
//         stateFreeVector(&testState);
//         testState.vec = NULL;
//     } // for qubits
// }

/*
 * =====================================================================================================================
 *                                                      main
 * =====================================================================================================================
 */
int main(void) {
    UNITY_BEGIN();
//    RUN_TEST(testMeanDiagObs);
//    RUN_TEST(testMeanCG);
//    RUN_TEST(testMeanObsHerm);
//    RUN_TEST(testGradPQCDiag);
//    RUN_TEST(testGradPQCHerm);
    RUN_TEST(testGradPQCDiag2);
    // RUN_TEST(testMMseq);
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
