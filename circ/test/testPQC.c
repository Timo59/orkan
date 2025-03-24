//
// Created by Timo Ziegler on 13.03.25.
//
/*
 * =====================================================================================================================
 *                                                      includes
 * =====================================================================================================================
 */
#ifndef TESTPQC_H
#include "testPQC.h"
#endif

/*
 * =====================================================================================================================
 *                                                      testEvoQB
 * =====================================================================================================================
 */
void testEvoQB(void) {
    state_t testState;
    for (qubit_t qubits = 2; qubits < MAXQUBITS; ++qubits) {
        const dim_t dim = POW2(qubits, dim_t);
        stateInitEmpty(&testState, qubits);
        cplx_t** vecs = generateTestVectors(qubits);

        cplx_t* pqbMat[4];                                          // Imaginary time evolution in matrix representation
        for (uint8_t i = 0; i < 3; i++) {
            cplx_t* prod = identityMat(qubits);                     // Product of respective Pauli gate on all qubits
            for (uint8_t j = 0; j < qubits; j++) {
                cplx_t* gateMat = singleQubitGateMat(pauliMat[i], qubits, j);
                cmatMulInPlace(prod, gateMat, dim);
                free(gateMat);
            }
            pqbMat[i] = zexpm(prod, randPar[i], dim);               // Imaginary time evolution matrix
            free(prod);
        }
        const uint8_t swapc = qubits / 2;                           // Number of adjacent swaps
        cplx_t* prod = identityMat(qubits);
        for (uint8_t j = 0; j < swapc; j++) {
            cplx_t* gateMat = swapGateMat(qubits, 2 * j, 2 * j + 1);
            cmatMulInPlace(prod, gateMat, dim);
            free(gateMat);
        }
        pqbMat[3] = zexpm(prod, randPar[3], dim);
        free(prod);

        /* TESTING */
        for (uint8_t i = 0; i < dim + 1; ++i) {
            stateInitVector(&testState, vecs[i]);                   // Initialize testState with the i-th test vector
            for (uint8_t j = 0; j < 4; j++) {                       // Evolve testState with all four quantum blocks
                evoQB(&testState, qbs[j + 4 * (qubits - 2)], randPar[j]);
            }

            for (uint8_t j = 0; j < 4; j++) {                       // Evolve the test vector itself by matrix multi-
                cmatVecMulInPlace(pqbMat[j], vecs[i], dim);         // plication
            }
            TEST_ASSERT_TRUE(cvectorAlmostEqual(vecs[i], testState.vec, dim, PRECISION));
        }
        for (uint8_t i = 0; i < 4; ++i) {
            free(pqbMat[i]);
        }
        stateFreeVector(&testState);
        freeTestVectors(vecs, qubits);
    }
}

/*
 * =====================================================================================================================
 *                                                      testLCQB
 * =====================================================================================================================
 */
void testLCQB(void) {
    state_t testState;
    for (qubit_t qubits = 2; qubits < MAXQUBITS; ++qubits) {
        const dim_t dim = POW2(qubits, dim_t);
        cplx_t** vecs = generateTestVectors(qubits);

        stateInitEmpty(&testState, qubits);
        applyQB* blocks = qbs + 4 * (qubits - 2);

        cplx_t* lcqbMat = calloc(dim * dim, sizeof (cplx_t));       // Matrix representation of the linear combination
        if (!lcqbMat) {                                             // of quantum blocks
            fprintf(stderr, "testLCQB: lcqbMat allocation failed\n");
            return;
        }
        for (uint8_t i = 0; i < 3; i++) {
            cplx_t* prod = identityMat(qubits);                     // Product of respective Pauli gates on all qubits
            for (uint8_t j = 0; j < qubits; j++) {
                cplx_t* gateMat = singleQubitGateMat(pauliMat[i], qubits, j);
                cmatMulInPlace(prod, gateMat, dim);
                free(gateMat);
            }
            cscalarMatMulInPlace(coeff[i], prod, dim);
            cmatAddInPlace(lcqbMat, prod, dim);
            free(prod);
        }
        const uint8_t swapc = qubits / 2;                           // Number of adjacent swaps
        cplx_t* prod = identityMat(qubits);
        for (uint8_t j = 0; j < swapc; j++) {
            cplx_t* gateMat = swapGateMat(qubits, 2 * j, 2 * j + 1);
            cmatMulInPlace(prod, gateMat, dim);
            free(gateMat);
        }
        cscalarMatMulInPlace(coeff[3], prod, dim);
        cmatAddInPlace(lcqbMat, prod, dim);
        free(prod);

        for (uint8_t i = 0; i < dim + 1; ++i) {
            stateInitVector(&testState, vecs[i]);
            lcQB(&testState, 4, blocks, coeff);

            cmatVecMulInPlace(lcqbMat, vecs[i], dim);

            TEST_ASSERT_TRUE(cvectorAlmostEqual(vecs[i], testState.vec, dim, PRECISION));
        }
        free(lcqbMat);
        stateFreeVector(&testState);
        freeTestVectors(vecs, qubits);
    }
}

/*
 * =====================================================================================================================
 *                                                 testApplyPQC
 * =====================================================================================================================
 */

// void testApplyPQC(void) {
//     state_t testState;
//     for (qubit_t qubits = 2; qubits < MAXQUBITS; qubits++) {
//         const dim_t dim = POW2(qubits, dim_t);
//         stateInitEmpty(&testState, qubits);
//         cplx_t** vecs = generateTestVectors(qubits);
//
//         /*
//          * Set up the reference matrices in an array holding the addresses of the PQBs holding the addresses of matrices
//          */
//         cplx_t** pqbMat[4];
//         for (uint8_t i = 0; i < 3; i++) {
//             if ((pqbMat[i] = malloc(qubits * sizeof (cplx_t*))) == NULL) {
//                 fprintf(stderr, "testApplyPQC: pqbMat[%d] allocation failed\n", i);
//                 for (uint8_t j = 0; j < i; j++) {
//                     for (uint8_t k = 0; k < qubits; k++) {          // Free all matrix memory allocated up to current i
//                         free(pqbMat[j][k]);
//                     }
//                     free(pqbMat[j]);                                // Free all matrix memory address space allocated
//                 }                                                   // up to current i
//                 freeTestVectors(vecs, qubits);
//                 stateFreeVector(&testState);
//                 exit(EXIT_FAILURE);
//             }
//             for (uint8_t j = 0; j < qubits; j++) {
//                 pqbMat[i][j] = rMat[i](qubits, j, randPar[i]);      // i determines whether X,Y, or Z and j the qubit
//             }
//         }
//         const uint8_t swapc = qubits / 2;                           // Number of adjacent swaps
//         if ((pqbMat[3] = malloc(swapc * sizeof (cplx_t*))) == NULL) {
//             fprintf(stderr, "testApplyPQC: pqbMat[3] allocation failed\n");
//             for (uint8_t j = 0; j < 3; j++) {
//                 for (uint8_t k = 0; k < qubits; k++) {              // Free all matrix memory allocated for Pauli
//                     free(pqbMat[j][k]);                             // rotations
//                 }
//                 free(pqbMat[j]);                                    // Free all matrix memory address space allocated
//             }                                                       // for Pauli rotations
//         }
//         for (uint8_t j = 0; j < swapc; j++) {
//             pqbMat[3][j] = RswapGateMat(qubits, 2 * j, 2 * j + 1, randPar[3]);
//         }
//
//         /* TESTING */
//         for (uint8_t i = 0; i < dim + 1; ++i) {
//             stateInitVector(&testState, vecs[i]);
//             applyPQC(&testState, 4, pqbs + 4 * (qubits - 2), randPar);
//
//             for (uint8_t j = 0; j < 3; j++) {
//                 for (uint8_t k = 0; k < qubits; k++) {
//                     cmatVecMulInPlace(pqbMat[j][k], vecs[i], dim);
//                 }
//             }
//             for (uint8_t j = 0; j < swapc; j++) {
//                 cmatVecMulInPlace(pqbMat[3][j], vecs[i], dim);
//             }
//
//             TEST_ASSERT_TRUE(cvectorAlmostEqual(vecs[i], testState.vec, dim, PRECISION));
//         }
//
//         for (uint8_t j = 0; j < 3; j++) {
//             for (uint8_t k = 0; k < qubits; k++) {              // Free all allocated matrix memory
//                 free(pqbMat[j][k]);
//             }
//             free(pqbMat[j]);                                    // Free all allocated matrix memory address space
//         }
//         for (uint8_t j = 0; j < swapc; j++) {
//             free(pqbMat[3][j]);
//         }
//         free(pqbMat[3]);
//         freeTestVectors(vecs, qubits);
//         stateFreeVector(&testState);
//     }
// }

/*
 * =====================================================================================================================
 *                                                      main
 * =====================================================================================================================
 */
int main(void) {
    UNITY_BEGIN();
    RUN_TEST(testEvoQB);
    RUN_TEST(testLCQB);
    return UNITY_END();
}