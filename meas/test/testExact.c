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
void testGradPQC(void) {
    state_t testState;
    for (qubit_t qubits = 2; qubits < MAXQUBITS; qubits++) {
        const dim_t dim = POW2(qubits, dim_t);
        cplx_t** vecs = generateTestVectors(qubits);

        stateInitEmpty(&testState, qubits);
        const applyQB* blocks = qbs + 4 * (qubits - 2);
        const applyPQB* circuit = pqc + 4 * (qubits - 2);

        cplx_t* pqbMat[4];                                          // Matrix representations of parametrized blocks
        const uint8_t swapc = qubits / 2;                           // Number of adjacent swaps

        /* TESTING */
        for (dim_t i = 0; i < dim + 1; ++i) {
            stateInitVector(&testState, vecs[i]);
            double* test = malloc(4 * sizeof (double));
            gradPQC(&testState, 4, circuit, randPar, blocks, diagObs, test);

            double ref[4];
            cplx_t* tmp = malloc(dim * sizeof (cplx_t));
            if (!tmp) {
                fprintf(stderr, "testGradPQC: tmp allocation failed\n");
                return;
            }
            for (dim_t j = 0; j < dim; ++j) {                       // Copy current test vector to tmp
                tmp[j] = vecs[i][j];
            }

            for (uint8_t j = 0; j < 3; ++j) {                       // Calculate the evolution matrices at the current
                cplx_t* prod = identityMat(qubits);                 // parameters
                for (uint8_t k = 0; k < qubits; k++) {
                    cplx_t* gateMat = singleQubitGateMat(pauliMat[j], qubits, k);
                    cmatMulInPlace(prod, gateMat, dim);
                    free(gateMat);
                }
                pqbMat[j] = zexpm(prod, randPar[j], dim);
                free(prod);
            }
            cplx_t* prod = identityMat(qubits);
            for (uint8_t j = 0; j < swapc; ++j) {
                cplx_t* gateMat = swapGateMat(qubits, 2 * j, 2 * j + 1);
                cmatMulInPlace(prod, gateMat, dim);
                free(gateMat);
            }
            pqbMat[3] = zexpm(prod, randPar[3], dim);
            free(prod);

            for (depth_t j = 0; j < 4; ++j) {                       // Evolve the temporary reference vector by matrix
                cmatVecMulInPlace(pqbMat[j], tmp, dim);             // multiplication
            }
            double mean = 0;                                        // Mean value of the observable at the current
            for (dim_t j = 0; j < dim; ++j) {                       // parameters
                mean += pow(cabs(tmp[j]), 2) * diagObs[j];
            }

            for (uint8_t j = 0; j < 4; ++j) {
                free(pqbMat[j]);
            }

            for (depth_t j = 0; j < 4; ++j) {                       // For all directions in parameter space
                randPar[j] += EPSILON;                              // Slightly shift the j-th parameter
                for (dim_t k = 0; k < dim; ++k) {                   // Copy current test vector to tmp
                    tmp[k] = vecs[i][k];
                }

                for (uint8_t k = 0; k < 3; ++k) {                   // Calculate the evolution matrices at the shifted
                    prod = identityMat(qubits);                     // parameters
                    for (uint8_t l = 0; l < qubits; l++) {
                        cplx_t* gateMat = singleQubitGateMat(pauliMat[k], qubits, l);
                        cmatMulInPlace(prod, gateMat, dim);
                        free(gateMat);
                    }
                    pqbMat[k] = zexpm(prod, randPar[k], dim);
                    free(prod);
                }
                prod = identityMat(qubits);
                for (uint8_t k = 0; k < swapc; ++k) {
                    cplx_t* gateMat = swapGateMat(qubits, 2 * k, 2 * k + 1);
                    cmatMulInPlace(prod, gateMat, dim);
                    free(gateMat);
                }
                pqbMat[3] = zexpm(prod, randPar[3], dim);
                free(prod);

                for (depth_t k = 0; k < 4; ++k) {                   // Evolve the temporary reference vector by matrix
                    cmatVecMulInPlace(pqbMat[k], tmp, dim);         // multiplication
                }
                double val = 0;                                     // Mean value of the observable at the shifted
                for (dim_t k = 0; k < dim; ++k) {                   // parameters
                    val += pow(cabs(tmp[k]), 2) * diagObs[k];
                }

                ref[j] = 1/EPSILON * (val - mean);
                randPar[j] -= EPSILON;
                for (uint8_t k = 0; k < 4; ++k) {
                    free(pqbMat[k]);
                }
            }

            TEST_ASSERT_TRUE(rvectorAlmostEqual(ref, test, 4, APPROXPRECISION));
            free(tmp);
            free(test);
        }
        stateFreeVector(&testState);
        freeTestVectors(vecs, qubits);
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
    RUN_TEST(testGradPQC);
    return UNITY_END();
}