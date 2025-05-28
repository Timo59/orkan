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
    state_t testState;
    for (qubit_t qubits = 1; qubits <= MAXQUBITS; qubits++) {
        const dim_t dim = POW2(qubits, dim_t);
        stateInitEmpty(&testState, qubits);
        cplx_t** vecs = generateTestVectors(qubits);

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
        stateFreeVector(&testState);
    }
}

void testMeanObsQb(void) {
    state_t testState;
    for (qubit_t qubits = 2; qubits <= MAXQUBITS; qubits++) {
        const dim_t dim = POW2(qubits, dim_t);
        stateInitEmpty(&testState, qubits);
        cplx_t** vecs = generateTestVectors(qubits);

        applyQB *observables = qb + 5 * (qubits - 2);

        for (uint8_t j = 0; j < 4; ++j) {
            cplx_t* observableMat = qbMat[5 * (qubits - 2) + j](qubits);
            for (uint8_t i = 0; i < dim + 1; i++) {
                stateInitVector(&testState, vecs[i]);
                const double test = meanObs(&testState, observables[j]);

                cplx_t* tmp = malloc(dim * sizeof(cplx_t));
                if (tmp == NULL) {
                    free(observableMat);
                    stateFreeVector(&testState);
                    freeTestVectors(vecs, qubits);
                    exit(EXIT_FAILURE);
                }
                for (dim_t k = 0; k < dim; k++) {
                    tmp[k] = vecs[i][k];
                }

                cmatVecMulInPlace(observableMat, tmp, dim);
                double ref = creal(cInner(tmp, vecs[i], dim));

                TEST_ASSERT_TRUE(ralmostEqual(ref, test, PRECISION));
            }
            free(observableMat);
        }
        freeTestVectors(vecs, qubits);
        stateFreeVector(&testState);
    }
}

/*
 * =====================================================================================================================
 *                                                 testGradientPQC
 * =====================================================================================================================
 */

void testGradPQC1(void) {
    state_t testState;
    for (qubit_t qubits = 2; qubits < MAXQUBITS; ++qubits) {
        const dim_t dim = POW2(qubits, dim_t);
        cplx_t** vecs = generateTestVectors(qubits);

        stateInitEmpty(&testState, qubits);
        const applyQB* blocks = qb + 5 * (qubits - 2);
        const applyPQB* circuit = pqc + 5 * (qubits - 2);

        cplx_t* pqbMat[5];                                          // Matrix representations of parametrized blocks

        /* TESTING */
        for (dim_t i = 0; i < dim + 1; ++i) {
            stateInitVector(&testState, vecs[i]);
            double* test = malloc(5 * sizeof (double));
            gradPQC(&testState, 5, circuit, randPar, blocks, diagObs, test);

            double ref[5];
            cplx_t* tmp = malloc(dim * sizeof (cplx_t));
            if (!tmp) {
                fprintf(stderr, "testGradPQC1: tmp allocation failed\n");
                return;
            }
            for (dim_t j = 0; j < dim; ++j) {                       // Copy current test vector to tmp
                tmp[j] = vecs[i][j];
            }

            for (uint8_t j = 0; j < 5; ++j) {                       // Calculate the evolution matrices with the current
                pqbMat[j] = evoMat[j](qubits, randPar[j]);          // parameters
            }

            for (int8_t j = 0; j < 5; ++j) {                        // Evolve the temporary reference vector by matrix
                cmatVecMulInPlace(pqbMat[j], tmp, dim);             // multiplication
            }
            double mean = 0;                                        // Mean value of the observable at the current
            for (dim_t j = 0; j < dim; ++j) {                       // parameters
                mean += pow(cabs(tmp[j]), 2) * diagObs[j];
            }


            double ref2[5];
            cplx_t* tmpBra = malloc(dim * sizeof (cplx_t));
            if (!tmpBra) {
                fprintf(stderr, "testGradPQC1: tmpBra allocation failed\n");
                return;
            }
            for (dim_t j = 0; j < dim; ++j) {
                tmpBra[j] = vecs[i][j];
            }
            for (uint8_t j = 0; j < 5; ++j) {
                cmatVecMulInPlace(pqbMat[j], tmpBra, dim);
            }

            for (uint8_t j = 0; j < 5; ++j) {
                for (dim_t k = 0; k < dim; ++k) {
                    tmp[k] = vecs[i][k];
                }

                for (uint8_t k = 0; k <= j; ++k) {
                    cmatVecMulInPlace(pqbMat[k], tmp, dim);
                }

                cplx_t* blockMat = qbMat[j](qubits);
                cmatVecMulInPlace(blockMat, tmp, dim);
                free(blockMat);

                for (uint8_t k = j + 1; k < 5; ++k) {
                    cmatVecMulInPlace(pqbMat[k], tmp, dim);
                }

                ref2[j] = 2 * cimag(cInner(tmpBra, tmp, dim));
            }
            free(tmpBra);

            for (uint8_t j = 0; j < 5; ++j) {
                free(pqbMat[j]);
            }


            for (depth_t j = 0; j < 5; ++j) {                       // For all directions in parameter space
                randPar[j] += EPSILON;                              // Slightly shift the j-th parameter
                for (dim_t k = 0; k < dim; ++k) {                   // Copy current test vector to tmp
                    tmp[k] = vecs[i][k];
                }

                for (uint8_t k = 0; k < 5; ++k) {                   // Calculate the evolution matrices with the shifted
                    pqbMat[k] = evoMat[k](qubits, randPar[k]);      // parameters
                }

                for (int8_t k = 0; k < 5; ++k) {                    // Evolve the temporary reference vector by matrix
                    cmatVecMulInPlace(pqbMat[k], tmp, dim);         // multiplication
                }
                double val = 0;                                     // Mean value of the observable at the shifted
                for (dim_t k = 0; k < dim; ++k) {                   // parameters
                    val += pow(cabs(tmp[k]), 2) * diagObs[k];
                }

                ref[j] = 1/EPSILON * (val - mean);
                randPar[j] -= EPSILON;

                for (uint8_t k = 0; k < 5; ++k) {
                    free(pqbMat[k]);
                }

            }

            TEST_ASSERT_TRUE(rvectorAlmostEqual(ref, test, 5, APPROXPRECISION));

            free(tmp);
            free(test);
        }
        stateFreeVector(&testState);
        freeTestVectors(vecs, qubits);
    }
}

//void testGradPQC2(void) {
//    state_t testState;
//    for (qubit_t qubits = 2; qubits < MAXQUBITS; ++qubits) {
//        const dim_t dim = POW2(qubits, dim_t);
//        cplx_t **vecs = generateTestVectors(qubits);
//
//        stateInitEmpty(&testState, qubits);
//        const applyQB* blocks = rqbs + 5 * (qubits - 2);
//        const applyPQB* circuit = rpqc + 5 * (qubits - 2);
//
//        cplx_t* pqbMat[5];                                          // Matrix representations of parametrized blocks
//        const uint8_t swapc = qubits / 2;                           // Number of adjacent swaps
//
//
//    }
//}

/*
 * =====================================================================================================================
 *                                                 testMMseq
 * =====================================================================================================================
 */
state_t testState = {                                               // State for the method to be tested
    .qubits = 0,
    .dim = 0,
    .vec = NULL
};
cplx_t *observableMat[3] = {NULL, NULL, NULL};                      // Matrix representations of observables
cplx_t** uMat[3] = {NULL, NULL, NULL};                              // Matrix representations of the channels' unitaries
cplx_t** vecs = NULL;                                               // Test state vectors;
cplx_t* momMat[3] = {NULL, NULL, NULL};                             // Moment matrices as the result of mmseq
cplx_t *test[3] = {NULL, NULL, NULL};                               // Moment matrices in dense row major form
cplx_t *ref[3] = {NULL, NULL, NULL};                                // Reference from matrix multiplication
cplx_t *refVec[5] = {NULL, NULL, NULL, NULL, NULL};                 // Reference state vecs to determine moment matrix

void cleanupTestMMseq(void) {
    if (testState.vec != NULL) {
        stateFreeVector(&testState);
    }
    for (uint8_t i = 0; i < 3; ++i) {
        if (observableMat[i] != NULL) {
            free(observableMat[i]);
        }
        if (uMat[i] != NULL) {
            for (uint8_t j = 0; j < 5; ++j) {
                if (uMat[i][j] != NULL) {
                    free(uMat[i][j]);
                }
            }
            free(uMat[i]);
        }
        if (momMat[i] != NULL) {
            free(momMat[i]);
        }
        if (test[i] != NULL) {
            free(test[i]);
        }
        if (ref[i] != NULL) {
            free(ref[i]);
        }
    }
    if (vecs != NULL) {
        freeTestVectors(vecs, testState.qubits);
    }
    for (uint8_t i = 0; i < 5; ++i) {
        if (refVec[i] != NULL) {
            free(refVec[i]);
        }
    }

}

void testMMseq(void) {
    if (atexit(cleanupTestMMseq)) {                                 // Register cleanup function
        fprintf(stderr, "Failed to register cleanup function\n");
        exit(EXIT_FAILURE);
    }

    for (qubit_t qubits = 2; qubits < MAXQUBITS; ++qubits) {        // Iterate number of qubits
        const dim_t dim = POW2(qubits, dim_t);
        stateInitEmpty(&testState, qubits);                         // Define the quantum state

        const applyQB* observable = obs + 3 * (qubits - 2);         // diagObs, all-Y and rotational SWAP blocks

        for (uint8_t i = 0; i < 3; ++i) {                           // Matrix representation of observables
            observableMat[i] = obsMat[i](qubits);
        }

        const applyQB* u = channel + 15 * (qubits - 2);             // Concatenation of search unitaries

        for (uint8_t i = 0; i < 3; ++i) {
            if ((uMat[i] = malloc(5 * sizeof (qbMat))) == NULL) {
                fprintf(stderr, "testMMseq: uMat[%d] allocation failed\n", i);
                exit(EXIT_FAILURE);
            }
            for (uint8_t j = 0; j < 5; ++j) {                       // Initilialize the matrix representations of the
                uMat[i][j] = channelMat[i][j](qubits);              // channels' unitaries
            }
        }

        for (uint8_t j = 0; j < 3; ++j) {                           // Iterate the channels of the sequence
            vecs = generateTestVectors(qubits);
            for (dim_t i = 0; i < dim + 1; ++i) {                   // Iterate the test vectors
                stateInitVector(&testState, vecs[i]);

                for (uint8_t k = 0; k < 3; ++k) {
                    if ((momMat[k] = malloc(5 * 5 * sizeof (cplx_t))) == NULL) {
                        fprintf(stderr, "testMMseq: momMat[%d] allocation failed\n", k);
                        exit(EXIT_FAILURE);
                    }
                }
                mmseq(&testState, 3, observable, 3, zcoeff, 5, u, j, momMat);

                for (uint8_t k = 0; k < 3; ++k) {
                    if ((test[k] = malloc(5 * 5 * sizeof (cplx_t))) == NULL) {
                        fprintf(stderr, "testMMseq: test[%d] allocation failed\n", k);
                        exit(EXIT_FAILURE);
                    }

                    if ((ref[k] = malloc(5 * 5 * sizeof (cplx_t))) == NULL) {
                        fprintf(stderr, "testMMseq: ref[%d] allocation failed\n", k);
                        exit(EXIT_FAILURE);
                    }
                }
                for (uint8_t k = 0; k < 5; ++k) {                   // Iterate the columns of the column major form
                    for (uint8_t mat = 0; mat < 3; ++mat) {         // Iterate the moment matrices
                        test[mat][k + k * 5] = momMat[mat][k + k * 5];  // Copy the diagonal elements
                        for (uint8_t l = 0; l < 5; ++l) {               // Iterate rows of the lower triangle
                            test[mat][l + 5 * k] = momMat[mat][k + 5 * l];
                        }
                    }
                }
                for (uint8_t k = 0; k < 3; ++k) {
                    free(momMat[k]);
                    momMat[k] = NULL;
                }

                for (uint8_t k = 0; k < j; ++k) {                   // Apply the matrix representation of LCU channels
                    cplx_t* sum = calloc(dim, sizeof (cplx_t));     // prior to the j-th channel to the reference vector
                    if (sum == NULL) {
                        fprintf(stderr, "testMMseq: %d-th sum allocation prior to link failed\n", k);
                        exit(EXIT_FAILURE);
                    }
                    for (uint8_t l = 0; l < 5; ++l) {
                        cplx_t* tmp = cmatVecMul(uMat[k][l], vecs[i], dim);
                        cscalarVecMulInPlace(zcoeff[k * 5 + l], tmp, dim);
                        cvecAddInPlace(sum, tmp, dim);
                        free(tmp);
                    }
                    free(vecs[i]);
                    vecs[i] = sum;
                }

                for (uint8_t k = 0; k < 5; ++k) {                   // Calculate the reference vectors spanning the
                    if ((refVec[k] = malloc(dim * sizeof (cplx_t))) == NULL) {  // moment matrix
                        fprintf(stderr, "testMMseq: refVec[%d] allocation failed\n", k);
                        exit(EXIT_FAILURE);
                    }

                    for (dim_t l = 0; l < dim; ++l) {               // Copy current reference state vector
                        refVec[k][l] = vecs[i][l];
                    }

                    cmatVecMulInPlace(uMat[j][k], refVec[k], dim);  // Apply the according search unitary

                    for (uint8_t l = j + 1; l < 3; ++l) {           // Apply the remaining LCU channels to refVec[k]
                        cplx_t* sum = calloc(dim, sizeof (cplx_t));
                        if (sum == NULL) {
                            fprintf(stderr, "testMMseq: %d-th sum allocation after link failed\n", l);
                            exit(EXIT_FAILURE);
                        }

                        for (uint8_t m = 0; m < 5; ++m) {
                            cplx_t* tmp = cmatVecMul(uMat[l][m], refVec[k], dim);
                            cscalarVecMulInPlace(zcoeff[l * 5 + m], tmp, dim);
                            cvecAddInPlace(sum, tmp, dim);
                            free(tmp);
                        }
                        free(refVec[k]);
                        refVec[k] = sum;
                    } // for LCU channels
                } // for reference vector


                for (uint8_t k = 0; k < 5; ++k) {                   // Iterate the rows of the reference moment matrices
                    for (uint8_t m = 0; m < 3; ++m) {               // Iterate the moment matrices
                        cplx_t* tmp = cmatVecMul(observableMat[m], refVec[k], dim);

                        for (uint8_t l = 0; l < 5; ++l) {               // Iterate columns of the reference moment matrix
                            ref[m][l + k * 5] = cInner(tmp, refVec[l], dim);
                        } // for column of moment matrix

                        free(tmp);
                    } // for moment matrix
                } // for row of reference moment matrix

                for (uint8_t k = 0; k < 3; ++k) {
                    TEST_ASSERT_TRUE(cvectorAlmostEqual(ref[k], test[k], 25, APPROXPRECISION));
                }

                for (uint8_t k = 0; k < 3; ++k) {
                    free(ref[k]);
                    ref[k] = NULL;
                    free(test[k]);
                    test[k] = NULL;
                }
                for (uint8_t k = 0; k < 5; ++k) {
                    free(refVec[k]);
                    refVec[k] = NULL;
                }
            } // for test vector
            freeTestVectors(vecs, qubits);
            vecs = NULL;
        } // for link
        for (uint8_t k = 0; k < 3; ++k) {
            for (uint8_t l = 0; l < 5; ++l) {
                free(uMat[k][l]);
                uMat[k][l] = NULL;
            }
            free(uMat[k]);
            uMat[k] = NULL;
            free(observableMat[k]);
            observableMat[k] = NULL;
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
    RUN_TEST(testMeanDiagObs);
    RUN_TEST(testMeanObsQb);
    RUN_TEST(testGradPQC1);
    RUN_TEST(testMMseq);
    return UNITY_END();
}