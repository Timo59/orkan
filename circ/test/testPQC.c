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
            apply(&testState, diagObs);

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
        testHerm.comp = cg + testHerm.len * (qubits - 2);
        testHerm.weight = dcoeff;

        // Define the hermitian operator's matrix representation
        cplx_t* testHermMat = calloc(dim * dim, sizeof(cplx_t));
        if (testHermMat == NULL) {
            fprintf(stderr, "testApplyHerm(): testHermMat allocation failed\n");
            exit(EXIT_FAILURE);
        }
        for (dim_t i = 0; i < 5; ++i) {
            cplx_t* tmp = cgMat[i](qubits);
            const cplx_t alpha = dcoeff[i];
            cscalarMatMulInPlace(alpha, tmp, dim);
            cmatAddInPlace(testHermMat, tmp, dim);
            free(tmp);
        }

        /* TESTING */
        for (uint8_t i = 0; i < dim + 1; ++i) {
            stateInitVector(&testState, vecs[i]);
            apply(&testState, &testHerm);

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
void testApplyLCCG(void) {
    state_t testState;
    for (qubit_t qubits = 2; qubits < MAXQUBITS; ++qubits) {
        const dim_t dim = POW2(qubits, dim_t);
        stateInitEmpty(&testState, qubits);
        cplx_t** vecs = generateTestVectors(qubits);

        // Define the hermitian operator struct
        lccg_t testLCCG;
        testLCCG.len = 5;                       // xi, yi, zi, swapi, diag for i = qubits defined in testPQC.h
        testLCCG.comp = cg + testLCCG.len * (qubits - 2);
        testLCCG.weight = zcoeff;

        // Define the LCCG's matrix representation
        cplx_t* testLCCGMat = calloc(dim * dim, sizeof (cplx_t));
        if (testLCCGMat == NULL) {
            fprintf(stderr, "testApplyLCQB(): testLCQBMat allocation failed\n");
            exit(EXIT_FAILURE);
        }
        for (uint8_t i = 0; i < 5; ++i) {
            cplx_t* tmp = cgMat[i](qubits);
            cscalarMatMulInPlace(zcoeff[i], tmp, dim);
            cmatAddInPlace(testLCCGMat, tmp, dim);
            free(tmp);
        }

        /* TESTING */
        for (dim_t i = 0; i < dim + 1; ++i) {
            stateInitVector(&testState, vecs[i]);
            apply(&testState, &testLCCG);

            cmatVecMulInPlace(testLCCGMat, vecs[i], dim);

            TEST_ASSERT_TRUE(cvectorAlmostEqual(vecs[i], testState.vec, dim, PRECISION));
        }
        free(testLCCGMat);
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
        cplx_t* evoMat = zexpm(mat, randPar[0] / 2., dim);
        free(mat);

        /* TESTING */
        for (dim_t i = 0; i < dim + 1; ++i) {
            stateInitVector(&testState, vecs[i]);
            evolve(&testState, diagObs, randPar[0]);

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
 *                                                      testEvoQB
 * =====================================================================================================================
 */
void testEvoCG(void) {
    state_t testState;
    for (qubit_t qubits = 2; qubits < MAXQUBITS; ++qubits) {
        const dim_t dim = POW2(qubits, dim_t);
        stateInitEmpty(&testState, qubits);
        cplx_t** vecs = generateTestVectors(qubits);
        const applyCG* testBlock = cg + 5 * (qubits - 2);

        cplx_t* pcgMat[4];                                          // Matrix representations of parametrized composite
        for (uint8_t i = 0; i < 4; ++i) {                           // gates
            cplx_t* gateMat = cgMat[i](qubits);
            pcgMat[i] = zexpm(gateMat, dcoeff[i] / 2., dim);
            free(gateMat);
        }

        /* TESTING */
        for (dim_t i = 0; i < dim + 1; ++i) {
            stateInitVector(&testState, vecs[i]);                   // Initialize testState with the i-th test vector
            for (uint8_t j = 0; j < 4; j++) {                       // Evolve testState with all four composite quantum
                evolve(&testState, testBlock[j], dcoeff[j]);         // gates
            }

            for (uint8_t j = 0; j < 4; j++) {                       // Evolve the test vector itself by matrix multi-
                cmatVecMulInPlace(pcgMat[j], vecs[i], dim);         // plication
            }

            TEST_ASSERT_TRUE(cvectorAlmostEqual(vecs[i], testState.vec, dim, PRECISION));
        }

        for (uint8_t i = 0; i < 4; ++i) {
            free(pcgMat[i]);
        }
        stateFreeVector(&testState);
        freeTestVectors(vecs, qubits);
    }
}

/*
 * =====================================================================================================================
 *                                                      testEvoHerm
 * =====================================================================================================================
 */
void testEvoHerm(void) {
    state_t testState;
    for (qubit_t qubits = 2; qubits < MAXQUBITS; ++qubits) {
        const dim_t dim = POW2(qubits, dim_t);
        stateInitEmpty(&testState, qubits);
        cplx_t** vecs = generateTestVectors(qubits);

        // Define the hermitian operator struct
        herm_t testHerm;
        testHerm.len = 4;                       // xi, yi, zi, swapi, diag for i = qubits defined in testPQC.h
        testHerm.comp = cg + (testHerm.len + 1) * (qubits - 2);
        testHerm.weight = dcoeff;

        // Define the matrix representations of the parametrized composite gates (Since the defined composite quantum
        // gates do not commute pairwise, we have to take the product of the matrix exponentials rather than the
        // exponential of the sum)
        cplx_t* testHermEvoMat = identityMat(qubits);
        for (uint8_t i = 1; i < 5; ++i) {
            cplx_t* gateMat = cgMat[4 - i](qubits);
            cplx_t* pqbMat = zexpm(gateMat, randPar[0] * dcoeff[4 - i] / 2., dim);
            cmatMulInPlace(testHermEvoMat, pqbMat, dim);
            free(gateMat);
            free(pqbMat);
        }

        /* TESTING */
        for (dim_t i = 0; i < dim + 1; ++i) {
            stateInitVector(&testState, vecs[i]);                   // Initialize testState with the i-th test vector
            evolve(&testState, &testHerm, randPar[0]);

            cmatVecMulInPlace(testHermEvoMat, vecs[i], dim);        // Evolve the test vector by matrix multiplication
        }

        free(testHermEvoMat);
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
    printf("sizeof(dim_t) = %zu\n", sizeof(dim_t));
    UNITY_BEGIN();
    RUN_TEST(testApplyDiag);
    RUN_TEST(testApplyHerm);
    RUN_TEST(testApplyLCCG);
    RUN_TEST(testEvoDiag);
    RUN_TEST(testEvoCG);
    RUN_TEST(testEvoHerm);
    return UNITY_END();
}