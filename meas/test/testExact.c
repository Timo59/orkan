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

            for (uint8_t j = 0; j < 5; ++j) {                       // Calculate the evolution matrices at the current
                cplx_t* blockMat = qbMat[j](qubits);                // parameters
                pqbMat[j] = zexpm(blockMat, randPar[j], dim);
                free(blockMat);
            }

            for (int8_t j = 4; j >= 0; --j) {                       // Evolve the temporary reference vector by matrix
                cmatVecMulInPlace(pqbMat[j], tmp, dim);             // multiplication
            }
            double mean = 0;                                        // Mean value of the observable at the current
            for (dim_t j = 0; j < dim; ++j) {                       // parameters
                mean += pow(cabs(tmp[j]), 2) * diagObs[j];
            }

            for (uint8_t j = 0; j < 5; ++j) {
                free(pqbMat[j]);
            }

            for (depth_t j = 0; j < 5; ++j) {                       // For all directions in parameter space
                randPar[j] += EPSILON;                              // Slightly shift the j-th parameter
                for (dim_t k = 0; k < dim; ++k) {                   // Copy current test vector to tmp
                    tmp[k] = vecs[i][k];
                }

                for (uint8_t k = 0; k < 5; ++k) {                   // Calculate the evolution matrices at the shifted
                    cplx_t* blockMat = qbMat[k](qubits);            // parameters
                    pqbMat[k] = zexpm(blockMat, randPar[k], dim);
                    free(blockMat);
                }

                for (int8_t k = 4; k >= 0; --k) {                   // Evolve the temporary reference vector by matrix
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
                printf("Qubits: %d\n", qubits);
                printf("test = ");
                vectorPrint(test, 5);
                printf("ref = ");
                vectorPrint(ref, 5);
                TEST_ASSERT_TRUE(rvectorAlmostEqual(ref, test, 5, APPROXPRECISION));
            }
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
void testMMseq(void) {
    state_t testState;
    for (qubit_t qubits = 2; qubits < MAXQUBITS; ++qubits) {
        const dim_t dim = POW2(qubits, dim_t);
        cplx_t** vecs = generateTestVectors(qubits);

        stateInitEmpty(&testState, qubits);
        const applyQB* blocks = qb + 5 * (qubits - 2);
        const applyPQB* circuit = pqc + 5 * (qubits - 2);

        applyQB u[30];                                              // Concatenation of channels' search unitaries
        for (uint8_t i = 0; i < 5; ++i) {
            for (uint8_t j = 0; j < 5; ++j) {                       // Alternating channels of the defined quantum
                u[i * 10 + j] = blocks[j];                          // blocks
            }
        }
    }
}
/*
 * =====================================================================================================================
 *                                                      main
 * =====================================================================================================================
 */
int main(void) {
    UNITY_BEGIN();
    RUN_TEST(testMeanObs);
    RUN_TEST(testGradPQC1);
    return UNITY_END();
}