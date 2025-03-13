//
// Created by Timo Ziegler on 13.03.25.
//
/*
 * =====================================================================================================================
 *                                                      includes
 * =====================================================================================================================
 */
#ifndef EXACT_H
#include <exact.h>
#endif

#ifndef QHIPSTER_H
#include "qhipster.h"
#endif

#ifndef GATEMAT_H
#include "gateMat.h"
#endif

#ifndef UNITY_FRAMEWORK_H
#include "unity.h"
#endif

#ifndef UTILS_H
#include "utils.h"
#endif

/*
 * =====================================================================================================================
 *                                                      macros
 * =====================================================================================================================
 */
#define APPROXPRECISION         1e-4
#define EPSILON                 1e-9
#define MAXQUBITS               6
#define PRECISION               1e-8

/*
 * =====================================================================================================================
 *                                                  setUp and tearDown
 * =====================================================================================================================
 */
void setUp(void) {}
void tearDown(void) {}

/*
 * =====================================================================================================================
 *                                                  testMeanObs
 * =====================================================================================================================
 */
double diagObs[64] = {-0.672013803, -1.080744850, 1.989524770, 1.988930799, 1.731153443, 0.219447958, -1.496797918,
    0.153447271, 1.394199856, 1.434125830, -0.525672657, 0.082399361, -0.256029741, 0.729175258, 1.631006491,
    2.747422307, -1.333258664, 0.147877182, 0.678307326, -2.064295279, -2.955271079, 0.834284985, -0.036103279,
    -0.174833204, 2.422058764, -2.213238216, 1.295379069, -2.132723608, 2.890693833, -2.400779788, -2.791631103,
    -0.748354991, -1.296158621, 1.286226791, -1.301539827, 3.019915018, 1.863008808, 0.373589359, -1.390747670,
    -2.016940260, 0.848049481, 1.515307618, -2.638919196, 1.680740517, -2.302014922, 0.864815536, -2.727445042,
    1.404416495, 1.087979959, -0.563509955, -2.881082581, 2.242557345, 1.721171693, 2.040283920, 1.542535211,
    0.520890342, 0.069996072, -3.053401020, 0.427764464, 1.344208689, -2.416954506, -2.499271157, 2.576468900,
    0.289965927
    };

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
double randPar[4] = {0.342377381, 1.416765874, -1.458084103, 1.767723341};

void rx2(state_t* state, const double par) {
    applyRX(state, 0, par);
    applyRX(state, 1, par);
}
void ry2(state_t* state, const double par) {
    applyRY(state, 0, par);
    applyRY(state, 1, par);
}
void rz2(state_t* state, const double par) {
    applyRZ(state, 0, par);
    applyRZ(state, 1, par);
}
void rswap2(state_t* state, const double par) {
    applyRswap(state, 0, 1, par);
}
void rx3(state_t* state, const double par) {
    applyRX(state, 0, par);
    applyRX(state, 1, par);
    applyRX(state, 2, par);
}
void ry3(state_t* state, const double par) {
    applyRY(state, 0, par);
    applyRY(state, 1, par);
    applyRY(state, 2, par);
}
void rz3(state_t* state, const double par) {
    applyRZ(state, 0, par);
    applyRZ(state, 1, par);
    applyRZ(state, 2, par);
}
void rswap3(state_t* state, const double par) {
    applyRswap(state, 0, 1, par);
}
void rx4(state_t* state, const double par) {
    applyRX(state, 0, par);
    applyRX(state, 1, par);
    applyRX(state, 2, par);
    applyRX(state, 3, par);
}
void ry4(state_t* state, const double par) {
    applyRY(state, 0, par);
    applyRY(state, 1, par);
    applyRY(state, 2, par);
    applyRY(state, 3, par);
}
void rz4(state_t* state, const double par) {
    applyRZ(state, 0, par);
    applyRZ(state, 1, par);
    applyRZ(state, 2, par);
    applyRZ(state, 3, par);
}
void rswap4(state_t* state, const double par) {
    applyRswap(state, 0, 1, par);
    applyRswap(state, 2, 3, par);
}
void rx5(state_t* state, const double par) {
    applyRX(state, 0, par);
    applyRX(state, 1, par);
    applyRX(state, 2, par);
    applyRX(state, 3, par);
    applyRX(state, 4, par);
}
void ry5(state_t* state, const double par) {
    applyRY(state, 0, par);
    applyRY(state, 1, par);
    applyRY(state, 2, par);
    applyRY(state, 3, par);
    applyRY(state, 4, par);
}
void rz5(state_t* state, const double par) {
    applyRZ(state, 0, par);
    applyRZ(state, 1, par);
    applyRZ(state, 2, par);
    applyRZ(state, 3, par);
    applyRZ(state, 4, par);
}
void rswap5(state_t* state, const double par) {
    applyRswap(state, 0, 1, par);
    applyRswap(state, 2, 3, par);
}
void rx6(state_t* state, const double par) {
    applyRX(state, 0, par);
    applyRX(state, 1, par);
    applyRX(state, 2, par);
    applyRX(state, 3, par);
    applyRX(state, 4, par);
    applyRX(state, 5, par);
}
void ry6(state_t* state, const double par) {
    applyRY(state, 0, par);
    applyRY(state, 1, par);
    applyRY(state, 2, par);
    applyRY(state, 3, par);
    applyRY(state, 4, par);
    applyRY(state, 5, par);
}
void rz6(state_t* state, const double par) {
    applyRZ(state, 0, par);
    applyRZ(state, 1, par);
    applyRZ(state, 2, par);
    applyRZ(state, 3, par);
    applyRZ(state, 4, par);
    applyRZ(state, 5, par);
}
void rswap6(state_t* state, const double par) {
    applyRswap(state, 0, 1, par);
    applyRswap(state, 2, 3, par);
    applyRswap(state, 4, 5, par);
}

applyPQB pqbs[20] = {rx2, ry2, rz2, rswap2, rx3, ry3, rz3, rswap3, rx4, ry4, rz4, rswap4, rx5, ry5, rz5, rswap5, \
    rx6, ry6, rz6, rswap6};

void x2(state_t* state) {
    applyX(state, 0);
    applyX(state, 1);
}
void y2(state_t* state) {
    applyY(state, 0);
    applyY(state, 1);
}
void z2(state_t* state) {
    applyZ(state, 0);
    applyZ(state, 1);
}
void swap2(state_t* state) {
    applySWAP(state, 0, 1);
}
void x3(state_t* state) {
    applyX(state, 0);
    applyX(state, 1);
    applyX(state, 2);
}
void y3(state_t* state) {
    applyY(state, 0);
    applyY(state, 1);
    applyY(state, 2);
}
void z3(state_t* state) {
    applyZ(state, 0);
    applyZ(state, 1);
    applyZ(state, 2);
}
void swap3(state_t* state) {
    applySWAP(state, 0, 1);
}
void x4(state_t* state) {
    applyX(state, 0);
    applyX(state, 1);
    applyX(state, 2);
    applyX(state, 3);
}
void y4(state_t* state) {
    applyY(state, 0);
    applyY(state, 1);
    applyY(state, 2);
    applyY(state, 3);
}
void z4(state_t* state) {
    applyZ(state, 0);
    applyZ(state, 1);
    applyZ(state, 2);
    applyZ(state, 3);
}
void swap4(state_t* state) {
    applySWAP(state, 0, 1);
    applySWAP(state, 2, 3);
}
void x5(state_t* state) {
    applyX(state, 0);
    applyX(state, 1);
    applyX(state, 2);
    applyX(state, 3);
    applyX(state, 4);
}
void y5(state_t* state) {
    applyY(state, 0);
    applyY(state, 1);
    applyY(state, 2);
    applyY(state, 3);
    applyY(state, 4);
}
void z5(state_t* state) {
    applyZ(state, 0);
    applyZ(state, 1);
    applyZ(state, 2);
    applyZ(state, 3);
    applyZ(state, 4);
}
void swap5(state_t* state) {
    applySWAP(state, 0, 1);
    applySWAP(state, 2, 3);
}
void x6(state_t* state) {
    applyX(state, 0);
    applyX(state, 1);
    applyX(state, 2);
    applyX(state, 3);
    applyX(state, 4);
    applyX(state, 5);
}
void y6(state_t* state) {
    applyY(state, 0);
    applyY(state, 1);
    applyY(state, 2);
    applyY(state, 3);
    applyY(state, 4);
    applyY(state, 5);
}
void z6(state_t* state) {
    applyZ(state, 0);
    applyZ(state, 1);
    applyZ(state, 2);
    applyZ(state, 3);
    applyZ(state, 4);
    applyZ(state, 5);
}
void swap6(state_t* state) {
    applySWAP(state, 0, 1);
    applySWAP(state, 2, 3);
    applySWAP(state, 4, 5);
}

applyQB qbs[20] = {x2, y2, z2, swap2, x3, y3, z3, swap3, x4, y4, z4, swap4, x5, y5, z5, swap5, x6, y6, z6, swap6};

cplx_t* (*rMat[])(qubit_t, qubit_t, double) = {RXGateMat, RYGateMat, RZGateMat};

void testGradPQC(void) {
    state_t testState;
    for (qubit_t qubits = 2; qubits < MAXQUBITS; qubits++) {
        const dim_t dim = POW2(qubits, dim_t);
        stateInitEmpty(&testState, qubits);
        cplx_t** vecs = generateTestVectors(qubits);

        /*
         * Define reference matrices in an array holding the addresses of the PQBs holding the addresses of matrices
         */
        cplx_t** pqbMat[4];
        for (uint8_t i = 0; i < 3; i++) {
            if ((pqbMat[i] = malloc(qubits * sizeof (cplx_t*))) == NULL) {
                fprintf(stderr, "testGradPQC: pqbMat[%d] allocation failed\n", i);
                for (uint8_t j = 0; j < i; j++) {
                    for (uint8_t k = 0; k < qubits; k++) {          // Free all matrix memory allocated up to current i
                        free(pqbMat[j][k]);
                    }
                    free(pqbMat[j]);                                // Free all matrix memory address space allocated
                }                                                   // up to current i
                freeTestVectors(vecs, qubits);
                stateFreeVector(&testState);
                exit(EXIT_FAILURE);
            }
        }
        const uint8_t swapc = qubits / 2;                           // Number of adjacent swaps
        if ((pqbMat[3] = malloc(swapc * sizeof (cplx_t*))) == NULL) {
            fprintf(stderr, "testGradPQC: pqbMat[3] allocation failed\n");
            for (uint8_t j = 0; j < 3; j++) {
                for (uint8_t k = 0; k < qubits; k++) {              // Free all matrix memory allocated for Pauli
                    free(pqbMat[j][k]);                             // rotations
                }
                free(pqbMat[j]);                                    // Free all matrix memory address space allocated
            }                                                       // for Pauli rotations
        }

        /* TESTING */
        for (dim_t i = 0; i < dim + 1; i++) {
            stateInitVector(&testState, vecs[i]);
            double* test = NULL;
            gradPQC(&testState, 4, pqbs + 4 * (qubits - 2), randPar, qbs + 4 * (qubits - 2), diagObs, &test);

            double ref[4];
            cplx_t* tmp = malloc(dim * sizeof (cplx_t));
            if (!tmp) {
                fprintf(stderr, "testGradPQC: tmp allocation failed\n");
                return;
            }
            for (dim_t j = 0; j < dim; j++) {                       // Copy current test vector to tmp
                tmp[j] = vecs[i][j];
            }

            for (uint8_t j = 0; j < 3; j++) {
                for (uint8_t k = 0; k < qubits; k++) {              // Initialize Pauli rotations with current parameter
                    pqbMat[j][k] = rMat[j](qubits, k, randPar[j]);
                }
            }
            for (uint8_t j = 0; j < swapc; j++) {                   // Initialize SWAP rotation with current parameters
                pqbMat[3][j] = RswapGateMat(qubits, 2 * j, 2 * j + 1, randPar[3]);
            }

            for (uint8_t j = 0; j < 3; j++) {
                for (uint8_t k = 0; k < qubits; k++) {              // Multiply tmp by Pauli rotations
                    cmatVecMulInPlace(pqbMat[j][k], tmp, dim);
                }
            }
            for (uint8_t j = 0; j < swapc; j++) {                   // Multiply tmp by SWAP rotations
                cmatVecMulInPlace(pqbMat[3][j], tmp, dim);
            }

            double mean = 0;                                        // Mean value of tmp with the current parameters
            for (dim_t j = 0; j < dim; j++) {
                mean += pow(cabs(tmp[j]), 2) * diagObs[j];
            }

            double ref2[4];
            stateInitVector(&testState, vecs[i]);
            applyPQC(&testState, 4, pqbs, randPar);
            double meanTest = meanObs(&testState, diagObs);

            for (depth_t j = 0; j < 4; ++j) {
                randPar[j] += EPSILON;                              // Slightly shift the j-th parameter value

                stateInitVector(&testState, vecs[i]);
                applyPQC(&testState, 4, pqbs, randPar);
                ref2[j] = 1. / EPSILON * (meanObs(&testState, diagObs) - meanTest);

                for (dim_t k= 0; k < dim; k++) {                    // Copy current test vector to tmp
                    tmp[k] = vecs[i][k];
                }

                for (uint8_t k = 0; k < 3; k++) {
                    for (uint8_t l = 0; l < qubits; l++) {          // Pauli rotations with shifted parameter
                        pqbMat[k][l] = rMat[k](qubits, l, randPar[k]);
                    }
                }
                for (uint8_t k = 0; k < swapc; k++) {               // SWAP rotation with shifted parameters
                    pqbMat[3][k] = RswapGateMat(qubits, 2 * k, 2 * k + 1, randPar[3]);
                }

                for (uint8_t k = 0; k < 3; k++) {
                    for (uint8_t l = 0; l < qubits; l++) {          // Multiply tmp by Pauli rotations
                        cmatVecMulInPlace(pqbMat[k][l], tmp, dim);
                    }
                }
                for (uint8_t k = 0; k < swapc; k++) {               // Multiply tmp by SWAP rotations
                    cmatVecMulInPlace(pqbMat[3][k], tmp, dim);
                }

                double val = 0;                                     // Mean value of tmp with the shifted parameters
                for (dim_t k = 0; k < dim; k++) {
                    val += pow(cabs(tmp[k]), 2) * diagObs[k];
                }
                ref[j] = 1. / EPSILON * (val - mean);
                randPar[j] -= EPSILON;
            }
            free(tmp);

            printf("test = ");
            vectorPrint(test, 4);
            printf("ref = ");
            vectorPrint(ref, 4);
            printf("ref = ");
            vectorPrint(ref2, 4);
            TEST_ASSERT_TRUE(rvectorAlmostEqual(ref, test, 4, APPROXPRECISION));
        }
        for (uint8_t j = 0; j < 3; j++) {
            for (uint8_t k = 0; k < qubits; k++) {              // Free all allocated matrix memory
                free(pqbMat[j][k]);
            }
            free(pqbMat[j]);                                    // Free all allocated matrix memory address space
        }
        for (uint8_t j = 0; j < swapc; j++) {
            free(pqbMat[3][j]);
        }
        free(pqbMat[3]);
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