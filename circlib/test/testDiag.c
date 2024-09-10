//
// Created by Timo Ziegler on 09.09.24.
//

/*
 * =====================================================================================================================
 *                                                      includes
 * =====================================================================================================================
 */

#include "helpfunc.h"
#include "diag.h"
#include "unity.h"

/*
 * =====================================================================================================================
 *                                                      macros
 * =====================================================================================================================
 */

#define MAXQUBITS               4
#define PRECISION               1e-8
#define STRINGC                 (1 << MAXQUBITS)

cplx_t coeffsI[STRINGC] = {-2.84182 + 0.77871 * I, 2.87720 - 0.06262 * I, -2.34638 + 1.18731 * I, 1.42133 - 2.30075 * I,
                           3.03348 - 0.44808 * I, -0.38861 - 1.34968 * I, 2.66357 - 2.13062 * I, -2.74569 * 1.33655 * I,
                           0.45088 + 0.68697 * I, 2.26591 - 0.37656 * I, -1.44839 - 0.19035 * I, 2.77288 - 2.64995 * I,
                           0.49731 + 1.91852 * I, 0.76366 + 2.66516 * I, -0.39762 - 1.45944 * I, 0.56492 - 1.37514 * I};

double coeffs[STRINGC] = {-1.24296, 0.74044, 2.26267, -1.27728, 1.33188, 3.13913, -0.71095, -2.95816, -0.45472, 1.58387,
                          -2.15871, 1.40738, -0.16989, -2.09528, -2.69691, 3.09517};

double angle[1] = {2.85648};

/*
 * =====================================================================================================================
 *                                                  setUp and tearDown
 * =====================================================================================================================
 */

void setUp(void) {}

void tearDown(void) {}

/*
 * =====================================================================================================================
 *                                              testApplyOpDiag
 * =====================================================================================================================
 */

/*
 * This function checks whether using applyOpDiag coincides with applying the Kronecker product of Pauli matrices on
 * instances ranging from 1 to 4 qubits.
 *
 * Input:
 *      -
 * Output:
 *      There is no return value. Only checks whether the tested state vector and the reference state
 *      vector are equal to some precision after an evolution and when reversing the evolution.
 */

void testApplyOpDiag(void) {
    state_t testState;

    for (qubit_t qubits = 1; qubits <= MAXQUBITS; ++qubits) {
        dim_t dim = POW2(qubits, dim_t);
        cplx_t** testVectors = generateTestVectors(qubits);
        cplx_t** refVectors = generateTestVectors(qubits);

        pauli_t* comps = allPauliStringsDiag(qubits);
        cplx_t* op = pauliObservableDiag(comps, )
        for (dim_t i = 0; i < dim + 1; ++i) {
            stateInitVector(&testState, testVectors[i], qubits);
            applyOpDiag(testState, op);
        }
        freeTestVectors(testVectors, qubits);
        freeTestVectors(refVectors, qubits);
    }
}