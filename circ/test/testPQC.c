//
// Created by Timo Ziegler on 13.03.25.
//
/*
 * =====================================================================================================================
 *                                                      includes
 * =====================================================================================================================
 */
#ifndef PQC_H
#include "pqc.h"
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
 *                                                 testApplyPQC
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

cplx_t* (*rMat[])(qubit_t, qubit_t, double) = {RXGateMat, RYGateMat, RZGateMat};

void testApplyPQC(void) {
    state_t testState;
    for (qubit_t qubits = 2; qubits < MAXQUBITS; qubits++) {
        const dim_t dim = POW2(qubits, dim_t);
        stateInitEmpty(&testState, qubits);
        cplx_t** vecs = generateTestVectors(qubits);

        /*
         * Set up the reference matrices in an array holding the addresses of the PQBs holding the addresses of matrices
         */
        cplx_t** pqbMat[4];
        for (uint8_t i = 0; i < 3; i++) {
            if ((pqbMat[i] = malloc(qubits * sizeof (cplx_t*))) == NULL) {
                fprintf(stderr, "testApplyPQC: pqbMat[%d] allocation failed\n", i);
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
            for (uint8_t j = 0; j < qubits; j++) {
                pqbMat[i][j] = rMat[i](qubits, j, randPar[i]);      // i determines whether X,Y, or Z and j the qubit
            }
        }
        const uint8_t swapc = qubits / 2;                           // Number of adjacent swaps
        if ((pqbMat[3] = malloc(swapc * sizeof (cplx_t*))) == NULL) {
            fprintf(stderr, "testApplyPQC: pqbMat[3] allocation failed\n");
            for (uint8_t j = 0; j < 3; j++) {
                for (uint8_t k = 0; k < qubits; k++) {              // Free all matrix memory allocated for Pauli
                    free(pqbMat[j][k]);                             // rotations
                }
                free(pqbMat[j]);                                    // Free all matrix memory address space allocated
            }                                                       // for Pauli rotations
        }
        for (uint8_t j = 0; j < swapc; j++) {
            pqbMat[3][j] = RswapGateMat(qubits, 2 * j, 2 * j + 1, randPar[3]);
        }

        /* TESTING */
        for (uint8_t i = 0; i < dim + 1; ++i) {
            stateInitVector(&testState, vecs[i]);
            applyPQC(&testState, 4, pqbs + 4 * (qubits - 2), randPar);

            for (uint8_t j = 0; j < 3; j++) {
                for (uint8_t k = 0; k < qubits; k++) {
                    cmatVecMulInPlace(pqbMat[j][k], vecs[i], dim);
                }
            }
            for (uint8_t j = 0; j < swapc; j++) {
                cmatVecMulInPlace(pqbMat[3][j], vecs[i], dim);
            }

            TEST_ASSERT_TRUE(cvectorAlmostEqual(vecs[i], testState.vec, dim, PRECISION));
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
        freeTestVectors(vecs, qubits);
        stateFreeVector(&testState);
    }
}

int main(void) {
    UNITY_BEGIN();
    RUN_TEST(testApplyPQC);
    return UNITY_END();
}