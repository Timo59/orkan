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
    for (qubit_t qubits = 2; qubits < MAXQUBITS; ++qubits) {
        const dim_t dim = POW2(qubits, dim_t);
        stateInitEmpty(&testState, qubits);
        vecs = generateTestVectors(qubits);

        // Define the diagonal hermitian operator's matrix representation
        testMat = diagMat(qubits);

        /* TESTING */
        for (uint8_t i = 0; i < dim + 1; ++i) {
            stateInitVector(&testState, vecs[i]);
            apply(&testState, diagObs);

            cmatVecMulInPlace(testMat, vecs[i], dim);

            TEST_ASSERT_TRUE(cvectorAlmostEqual(vecs[i], testState.vec, dim, PRECISION));
        }

        free(testMat);
        testMat = NULL;
        freeTestVectors(vecs, qubits);
        vecs = NULL;
        stateFreeVector(&testState);
        testState.vec = NULL;
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

        // Initialize the hermitian operator struct
        testHerm.len = 5;                       // xi, yi, zi, swapi, diag for i = qubits defined in testPQC.h
        testHerm.comp = cg + testHerm.len * (qubits - 2);
        testHerm.weight = dcoeff;

        // Define the hermitian operator's matrix representation
        if ((testMat = calloc(dim * dim, sizeof(cplx_t))) == NULL) {
            fprintf(stderr, "testApplyHerm(): testHermMat allocation failed\n");
            exit(EXIT_FAILURE);
        }
        for (dim_t i = 0; i < 5; ++i) {
            cplx_t* tmp = cgMat[i](qubits);
            const cplx_t alpha = dcoeff[i];
            cscalarMatMulInPlace(alpha, tmp, dim);
            cmatAddInPlace(testMat, tmp, dim);
            free(tmp);
            tmp = NULL;
        }

        /* TESTING */
        for (uint8_t i = 0; i < dim + 1; ++i) {
            stateInitVector(&testState, vecs[i]);
            apply(&testState, &testHerm);

            cmatVecMulInPlace(testMat, vecs[i], dim);

            TEST_ASSERT_TRUE(cvectorAlmostEqual(vecs[i], testState.vec, dim, PRECISION));
        }

        free(testMat);
        testMat = NULL;
        freeTestVectors(vecs, qubits);
        vecs = NULL;
        stateFreeVector(&testState);
        testState.vec = NULL;
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
        if ((testMat = calloc(dim * dim, sizeof (cplx_t))) == NULL) {
            fprintf(stderr, "testApplyLCQB(): testLCQBMat allocation failed\n");
            exit(EXIT_FAILURE);
        }
        for (uint8_t i = 0; i < 5; ++i) {
            cplx_t* tmp = cgMat[i](qubits);
            cscalarMatMulInPlace(zcoeff[i], tmp, dim);
            cmatAddInPlace(testMat, tmp, dim);
            free(tmp);
            tmp = NULL;
        }

        /* TESTING */
        for (dim_t i = 0; i < dim + 1; ++i) {
            stateInitVector(&testState, vecs[i]);
            apply(&testState, &testLCCG);

            cmatVecMulInPlace(testMat, vecs[i], dim);

            TEST_ASSERT_TRUE(cvectorAlmostEqual(vecs[i], testState.vec, dim, PRECISION));
        }
        free(testMat);
        testMat = NULL;
        stateFreeVector(&testState);
        testState.vec = NULL;
        freeTestVectors(vecs, qubits);
        vecs = NULL;
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

        testMat = diagMat(qubits);
        evoMat = zexpm(testMat, randPar[0] / 2., dim);
        free(testMat);
        testMat = NULL;

        /* TESTING */
        for (dim_t i = 0; i < dim + 1; ++i) {
            stateInitVector(&testState, vecs[i]);
            evolve(&testState, diagObs, randPar[0]);

            cmatVecMulInPlace(evoMat, vecs[i], dim);

            TEST_ASSERT_TRUE(cvectorAlmostEqual(vecs[i], testState.vec, dim, PRECISION));
        }

        free(evoMat);
        evoMat = NULL;
        stateFreeVector(&testState);
        testState.vec = NULL;
        freeTestVectors(vecs, qubits);
        vecs = NULL;
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

        evoMat = identityMat(qubits);                               // Matrix representations of parametrized composite
        for (uint8_t i = 0; i < 4; ++i) {                           // gates
            testMat = cgMat[3 - i](qubits);
            cplx_t* pcgMat = zexpm(testMat, dcoeff[3 - i] / 2., dim);
            cmatMulInPlace(evoMat, pcgMat, dim);
            free(testMat);
            testMat = NULL;
            free(pcgMat);
            pcgMat = NULL;
        }

        /* TESTING */
        for (dim_t i = 0; i < dim + 1; ++i) {
            stateInitVector(&testState, vecs[i]);                   // Initialize testState with the i-th test vector
            for (uint8_t j = 0; j < 4; j++) {                       // Evolve testState with all four composite quantum
                evolve(&testState, testBlock[j], dcoeff[j]);         // gates
            }

            cmatVecMulInPlace(evoMat, vecs[i], dim);                // Evolve the test vector itself by matrix multi-
                                                                    // plication
            TEST_ASSERT_TRUE(cvectorAlmostEqual(vecs[i], testState.vec, dim, PRECISION));
        }

        free(evoMat);
        evoMat = NULL;
        stateFreeVector(&testState);
        testState.vec = NULL;
        freeTestVectors(vecs, qubits);
        vecs = NULL;
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
        testHerm.len = 4;                       // xi, yi, zi, swapi, diag for i = qubits defined in testPQC.h
        testHerm.comp = cg + (testHerm.len + 1) * (qubits - 2);
        testHerm.weight = dcoeff;

        // Define the matrix representations of the parametrized composite gates (Since the defined composite quantum
        // gates do not commute pairwise, we have to take the product of the matrix exponentials rather than the
        // exponential of the sum)
        cplx_t* evoMat = identityMat(qubits);
        for (uint8_t i = 1; i < 5; ++i) {
            testMat = cgMat[4 - i](qubits);
            cplx_t* pqbMat = zexpm(testMat, randPar[0] * dcoeff[4 - i] / 2., dim);
            cmatMulInPlace(evoMat, pqbMat, dim);
            free(testMat);
            testMat = NULL;
            free(pqbMat);
            pqbMat = NULL;
        }

        /* TESTING */
        for (dim_t i = 0; i < dim + 1; ++i) {
            stateInitVector(&testState, vecs[i]);                   // Initialize testState with the i-th test vector
            evolve(&testState, &testHerm, randPar[0]);

            cmatVecMulInPlace(evoMat, vecs[i], dim);                // Evolve the test vector by matrix multiplication
        }

        free(evoMat);
        evoMat = NULL;
        stateFreeVector(&testState);
        testState.vec = NULL;
        freeTestVectors(vecs, qubits);
        vecs = NULL;
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
#ifdef OPENBLAS_USE64BITINT
    printf("Using 64-bit OpenBLAS\n");
#else
    printf("Using 32-bit OpenBLAS\n");
#endif
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