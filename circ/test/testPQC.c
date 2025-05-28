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
 *                                                      applyDiag
 * =====================================================================================================================
 */
void testApplyDiag(void) {
    state_t testState;
    for (qubit_t qubits = 2; qubits < MAXQUBITS; ++qubits) {
        const dim_t dim = POW2(qubits, dim_t);
        stateInitEmpty(&testState, qubits);
        cplx_t** vecs = generateTestVectors(qubits);

        // Define the diagonal hermitian operator's matrix representation
        cplx_t* mat = diagMat(qubits);

        /* TESTING */
        for (uint8_t i = 0; i < dim + 1; ++i) {
            stateInitVector(&testState, vecs[i]);
            applyDiag(&testState, diagObs);

            cmatVecMulInPlace(mat, vecs[i], dim);

            TEST_ASSERT_TRUE(cvectorAlmostEqual(vecs[i], testState.vec, dim, PRECISION));
        }

        free(mat);
        freeTestVectors(vecs, qubits);
        stateFreeVector(&testState);
    }
}

/*
 * =====================================================================================================================
 *                                                      applyHerm
 * =====================================================================================================================
 */
void testApplyHerm(void) {
    state_t testState;
    for (qubit_t qubits = 2; qubits < MAXQUBITS; ++qubits) {
        const dim_t dim = POW2(qubits, dim_t);
        stateInitEmpty(&testState, qubits);
        cplx_t** vecs = generateTestVectors(qubits);

        // Define the hermitian operator struct
        herm_t testHerm;
        testHerm.len = 5;                       // xi, yi, zi, swapi, diag for i = qubits defined in testPQC.h
        testHerm.comp = qb + testHerm.len * (qubits - 2);
        testHerm.weight = dcoeff;

        // Define the hermitian operator's matrix representation
        cplx_t* testHermMat = calloc(dim * dim, sizeof(cplx_t));
        if (testHermMat == NULL) {
            fprintf(stderr, "testApplyHerm(): testHermMat allocation failed\n");
            exit(EXIT_FAILURE);
        }
        for (dim_t i = 0; i < 5; ++i) {
            cplx_t* tmp = qbMat[i](qubits);
            const cplx_t alpha = dcoeff[i];
            cscalarMatMulInPlace(alpha, tmp, dim);
            cmatAddInPlace(testHermMat, tmp, dim);
            free(tmp);
        }

        /* TESTING */
        for (uint8_t i = 0; i < dim + 1; ++i) {
            stateInitVector(&testState, vecs[i]);
            applyHerm(&testState, &testHerm);

            cmatVecMulInPlace(testHermMat, vecs[i], dim);

            TEST_ASSERT_TRUE(cvectorAlmostEqual(vecs[i], testState.vec, dim, PRECISION));
        }

        free(testHermMat);
        freeTestVectors(vecs, qubits);
        stateFreeVector(&testState);
    }
}

/*
 * =====================================================================================================================
 *                                                      applyLCQB
 * =====================================================================================================================
 */
void testApplyLCQB(void) {
    state_t testState;
    for (qubit_t qubits = 2; qubits < MAXQUBITS; ++qubits) {
        const dim_t dim = POW2(qubits, dim_t);
        stateInitEmpty(&testState, qubits);
        cplx_t** vecs = generateTestVectors(qubits);

        // Define the hermitian operator struct
        lcqb_t testLCQB;
        testLCQB.len = 5;                       // xi, yi, zi, swapi, diag for i = qubits defined in testPQC.h
        testLCQB.comp = qb + testLCQB.len * (qubits - 2);
        testLCQB.weight = zcoeff;

        // Define the LCQB's matrix representation
        cplx_t* testLCQBMat = calloc(dim * dim, sizeof (cplx_t));
        if (testLCQBMat == NULL) {
            fprintf(stderr, "testApplyLCQB(): testLCQBMat allocation failed\n");
            exit(EXIT_FAILURE);
        }
        for (uint8_t i = 0; i < 5; ++i) {
            cplx_t* tmp = qbMat[i](qubits);
            cscalarMatMulInPlace(zcoeff[i], tmp, dim);
            cmatAddInPlace(testLCQBMat, tmp, dim);
            free(tmp);
        }

        /* TESTING */
        for (dim_t i = 0; i < dim + 1; ++i) {
            stateInitVector(&testState, vecs[i]);
            applyLCQB(&testState, &testLCQB);

            cmatVecMulInPlace(testLCQBMat, vecs[i], dim);

            TEST_ASSERT_TRUE(cvectorAlmostEqual(vecs[i], testState.vec, dim, PRECISION));
        }
        free(testLCQBMat);
        stateFreeVector(&testState);
        freeTestVectors(vecs, qubits);
    }
}

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
        const applyQB* blocks = qb + 4 * (qubits - 2);

        cplx_t* pqbMat[4];                                          // Matrix representations of parametrized blocks
        for (uint8_t i = 0; i < 4; ++i) {
            cplx_t* gateMat = qbMat[i](qubits);
            pqbMat[i] = zexpm(gateMat, randPar[i], dim);
            free(gateMat);
        }

        /* TESTING */
        for (dim_t i = 0; i < dim + 1; ++i) {
            stateInitVector(&testState, vecs[i]);                   // Initialize testState with the i-th test vector
            for (uint8_t j = 0; j < 4; j++) {                       // Evolve testState with all four quantum blocks
                evoQB(&testState, blocks[j], randPar[j]);
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
 *                                                      testEvoDiag
 * =====================================================================================================================
 */

void testEvoDiag(void) {
    state_t testState;
    for (qubit_t qubits = 2; qubits < MAXQUBITS; ++qubits) {
        const dim_t dim = POW2(qubits, dim_t);
        cplx_t** vecs = generateTestVectors(qubits);
        stateInitEmpty(&testState, qubits);

        cplx_t* mat = diagMat(qubits);
        cplx_t* evoMat = zexpm(mat, randPar[0], dim);
        free(mat);

        /* TESTING */
        for (dim_t i = 0; i < dim + 1; ++i) {
            stateInitVector(&testState, vecs[i]);
            evoDiag(&testState, diagObs, randPar[0]);

            cmatVecMulInPlace(evoMat, vecs[i], dim);

            TEST_ASSERT_TRUE(cvectorAlmostEqual(vecs[i], testState.vec, dim, PRECISION));
        }

        free(evoMat);
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
    RUN_TEST(testApplyDiag);
    RUN_TEST(testApplyHerm);
    RUN_TEST(testApplyLCQB);
    // RUN_TEST(testEvoQB);
    // RUN_TEST(testEvoDiag);
    return UNITY_END();
}