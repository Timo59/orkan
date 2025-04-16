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
void testMeanObs(void) {
    state_t testState;
    for (qubit_t qubits = 1; qubits <= MAXQUBITS; qubits++) {
        const dim_t dim = POW2(qubits, dim_t);
        stateInitEmpty(&testState, qubits);
        cplx_t** vecs = generateTestVectors(qubits);

        for (uint8_t i = 0; i < dim + 1; i++) {
            stateInitVector(&testState, vecs[i]);
            double test = meanObs(&testState, diagObs);

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
cplx_t *refBra = NULL, *refKet = NULL;                              // Reference bra and ket to determine moment matrix

void cleanupTestMMseq(void) {
    if (testState.vec != NULL) {
        printf("Freeing testState...\n");
        stateFreeVector(&testState);
    }
    for (uint8_t i = 0; i < 3; ++i) {
        if (observableMat[i] != NULL) {
            printf("Freeing observableMat[%d]...\n", i);
            free(observableMat[i]);
        }
        if (uMat[i] != NULL) {
            for (uint8_t j = 0; j < 5; ++j) {
                if (uMat[i][j] != NULL) {
                    printf("Freeing uMat[%d][%d]...\n", i, j);
                    free(uMat[i][j]);
                }
            }
            printf("Freeing uMat[%d]...\n", i);
            free(uMat[i]);
        }
        if (momMat[i] != NULL) {
            printf("Freeing momMat[%d]...\n", i);
            free(momMat[i]);
        }
        if (test[i] != NULL) {
            printf("Freeing test[%d]...\n", i);
            free(test[i]);
        }
        if (ref[i] != NULL) {
            printf("Freeing ref[%d]...\n", i);
            free(ref[i]);
        }
    }
    if (vecs != NULL) {
        printf("Freeing vecs...\n");
        freeTestVectors(vecs, testState.qubits);
    }
    if (refBra != NULL) {
        printf("Freeing refBra...\n");
        free(refBra);
    }
    if (refKet != NULL) {
        printf("Freeing refKet...\n");
        free(refKet);
    }
}

void testMMseq(void) {
    if (atexit(cleanupTestMMseq)) {                                 // Register cleanup function
        fprintf(stderr, "Failed to register cleanup function\n");
        exit(EXIT_FAILURE);
    }

    for (qubit_t qubits = 2; qubits < MAXQUBITS; ++qubits) {
        const dim_t dim = POW2(qubits, dim_t);
        stateInitEmpty(&testState, qubits);

        const applyQB* observable = obs + 3 * (qubits - 2);         // diagObs, all-Y and rotational SWAP blocks

        for (uint8_t i = 0; i < 3; ++i) {
            observableMat[i] = obsMat[i](qubits);
        }

        const applyQB* u = channel + 15 * (qubits - 2);             // Concatenation of search unitaries

        for (uint8_t i = 0; i < 3; ++i) {
            if ((uMat[i] = malloc(5 * sizeof (cplx_t*))) == NULL) {
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
                    if ((momMat[k] = malloc(5 * 6 / 2 * sizeof (cplx_t))) == NULL) {
                        fprintf(stderr, "testMMseq: momMat[%d] allocation failed\n", k);
                        exit(EXIT_FAILURE);
                    }
                }
                mmseq(&testState, 3, observable, 3, coeff, 5, u, j, momMat);

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
                for (uint8_t k = 0; k < 5; ++k) {                   // Iterate the columns of packed column major form
                    const uint16_t offset = 5 * k - k * (k - 1) / 2;   // Elements prior to the k-th column
                    for (uint8_t l = 0; l < 5 - k; ++l) {               // Iterate rows of the lower triangle
                        for (uint8_t mat = 0; mat < 3; ++mat) {
                            test[mat][l + k + 5 * k] = conj(momMat[mat][l + offset]);
                            test[mat][k + 5 * (l + k)] = momMat[mat][l + offset];
                        }
                    }
                }
                for (uint8_t k = 0; k < 3; ++k) {
                    free(momMat[k]);
                    momMat[k] = NULL;
                }

                for (uint8_t k = 0; k < j; ++k) {                   // Matrix representation of LCU up to channel prior
                    cplx_t* sum = calloc(dim, sizeof (cplx_t));
                    if (sum == NULL) {
                        fprintf(stderr, "testMMseq: %d-th sum allocation prior to link failed\n", k);
                        exit(EXIT_FAILURE);
                    }
                    for (uint8_t l = 0; l < 5; ++l) {
                        cplx_t* tmp = cmatVecMul(uMat[k][l], vecs[i], dim);
                        cscalarVecMulInPlace(coeff[k * 5 + l], tmp, dim);
                        cvecAddInPlace(sum, tmp, dim);
                        free(tmp);
                    }
                    free(vecs[i]);
                    vecs[i] = sum;
                }

                if ((refBra = malloc(dim * sizeof (cplx_t))) == NULL) {
                    fprintf(stderr, "testMMseq: refBra/refKet allocation failed\n");
                    exit(EXIT_FAILURE);
                }
                if ((refKet = malloc(dim * sizeof (cplx_t))) == NULL) {
                    fprintf(stderr, "testMMseq: refBra/refKet allocation failed\n");
                    exit(EXIT_FAILURE);
                }

                for (uint8_t k = 0; k < 5; ++k) {                   // Iterate the rows of the reference moment matrices
                    for (uint8_t l = 0; l < 5; ++l) {               // Iterate columns of the reference moment matrix
                        for (uint8_t m = 0; m < 3; ++m) {           // Iterate the reference moment matrices themselves

                            for (dim_t n = 0; n < dim; ++n) {       // Copy current state vector to bra and ket
                                refBra[n] = vecs[i][n];
                                refKet[n] = vecs[i][n];
                            }

                            for (uint8_t n = j + 1; n < 3; ++n) {   // Apply the remaining LCU channels to refBra
                                cplx_t* sum = calloc(dim, sizeof (cplx_t));
                                if (sum == NULL) {
                                    fprintf(stderr, "testMMseq: %d-th sum allocation after link failed\n", n);
                                    exit(EXIT_FAILURE);
                                }
                                for (uint8_t s = 0; s < 5; ++s) {
                                    cplx_t* tmp = cmatVecMul(uMat[n][s], refBra, dim);
                                    cscalarVecMulInPlace(coeff[n * 5 + s], tmp, dim);
                                    cvecAddInPlace(sum, tmp, dim);
                                    free(tmp);
                                }
                                free(refBra);
                                refBra = sum;
                            }

                            if (m == 0 && l == 0) {
                                cmatVecMulInPlace(uMat[j][k], refBra, dim);
                                printf("Reference after (%d, %d)-th unitary: ", j, k);
                                vectorPrint(refBra, dim);
                            }

                            cmatVecMulInPlace(observableMat[m], refBra, dim);

                            cmatVecMulInPlace(uMat[j][l], refKet, dim);
                            for (uint8_t n = j + 1; n < 3; ++n) {   // Apply the remaining LCU channels to refBra
                                cplx_t* sum = calloc(dim, sizeof (cplx_t));
                                if (sum == NULL) {
                                    fprintf(stderr, "testMMseq: %d-th sum allocation after link failed\n", n);
                                    exit(EXIT_FAILURE);
                                }
                                for (uint8_t s = 0; s < 5; ++s) {
                                    cplx_t* tmp = cmatVecMul(uMat[n][s], refKet, dim);
                                    cscalarVecMulInPlace(coeff[n * 5 + s], tmp, dim);
                                    cvecAddInPlace(sum, tmp, dim);
                                    free(tmp);
                                }
                                free(refKet);
                                refKet = sum;
                            }
                            ref[m][l + k * 5] = cInner(refBra, refKet, dim);
                        } // for moment matrix
                    } // for column of moment matrix
                } // for row of reference moment matrix

                printf("Qubits = %d, vector = %ld\n", qubits, i);
                for (uint8_t k = 0; k < 3; ++k) {
                    printf("ref[%d] =\n", k);
                    cmatrixPrint(ref[k], 5);
                    printf("test[%d] =\n", k);
                    cmatrixPrint(test[k], 5);
                    TEST_ASSERT_TRUE(cvectorAlmostEqual(ref[k], test[k], 25, PRECISION));
                }

                for (uint8_t k = 0; k < 3; ++k) {
                    free(ref[k]);
                    ref[k] = NULL;
                    free(test[k]);
                    test[k] = NULL;
                }
                free(refBra);
                refBra = NULL;
                free(refKet);
                refKet = NULL;
            } // for test vector
            freeTestVectors(vecs, qubits);
            vecs = NULL;
        } // for link
        for (uint8_t k = 0; k < 3; ++k) {
            for (uint8_t l = 0; l < 5; ++l) {
                free(uMat[k][l]);
                uMat[k][l] = NULL;
            }
            free(observableMat[k]);
            observableMat[k] = NULL;
        }
        stateFreeVector(&testState);
        testState.vec = NULL;
    } // for qubits
    exit(EXIT_SUCCESS);
}
/*
 * =====================================================================================================================
 *                                                      main
 * =====================================================================================================================
 */
int main(void) {
    UNITY_BEGIN();
    // RUN_TEST(testMeanObs);
    // RUN_TEST(testGradPQC1);
    RUN_TEST(testMMseq);
    return UNITY_END();
}