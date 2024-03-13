/* 
 * =====================================================================================================================
 *                                                      includes
 * =====================================================================================================================
 */

#include "helpfunc.h"
#include "exact.h"
#include "cplxutil.h"
#include "linalg.h"
#include "unity.h"
//#include <cblas.h>

/*
 * =====================================================================================================================
 *                                                      macros
 * =====================================================================================================================
 */

#define PRECISION               1e-8
#define APPROXPRECISION         1e-5

/*
 * =====================================================================================================================
 *                                                  setUp and tearDown
 * =====================================================================================================================
 */

void setUp(void) {}

void tearDown(void) {}

/* 
 * =====================================================================================================================
 *                                                      test expValObs
 * =====================================================================================================================
 */

void testExpValObs(void) {
    state_t testState;
    cplx_t **testVectors;
    cplx_t **refBra;
    cplx_t **refKet;
    double result;
    double resultDiag;
    double reference;

    /*
                                                        1 qubit
    --------------------------------------------------------------------------------------------------------------------
    */

    qubit_t qubits = 1;
    dim_t dim = POW2(qubits, dim_t);

    /*
     * Initialize the state vectors
     */
    testVectors = generateTestVectors(qubits);
    refBra = generateTestVectors(qubits);
    refKet = generateTestVectors(qubits);

    /*
     * Define the diagonal observable
     */
    complength_t length = dim;
    pauli_t* components = (pauli_t *) malloc(qubits * length * sizeof(pauli_t));
    components[0] = ID;
    components[1] = Z;
    double *coefficients = (double *) malloc(length * sizeof(double));
    coefficients[0] = 1.13126;
    coefficients[1] = 0.24848;

    /*
     * Initialize the diagonal observable
     */
    pauliObs_t pauliObservable;                         // Initialize the Pauli representation of the diagonal
    pauliObservable.components = components;            // hamiltonian
    pauliObservable.coefficients = coefficients;        //
    pauliObservable.length = length;                    //
    pauliObservable.qubits = qubits;                    //
                                                        //
    obs_t observablePauli;                              //
    observablePauli.type = PAULI;                       //
    observablePauli.pauliObs = &pauliObservable;        //

    double* diagObservable = pauliObservableDiag(components, coefficients, length, qubits);
    obs_t observableDiag;                                                       // Initialize the hamiltonian's diagonal
    observableDiag.type = DIAG;                                                 // representation
    observableDiag.diagObs = diagObservable;                                    //

    cplx_t* observableMat = pauliObservableMat(components, coefficients, length, qubits);    // Initialize the reference
                                                                                            // matrix according to the
                                                                                            //  Pauli string

    /*
     * Calculate the expectation values
     */

    for (dim_t i = 0; i < dim + 1; ++i) {
        stateInitVector(&testState, testVectors[i], qubits);
        result = expValObs(&testState, &observablePauli);       // Calculate the expectation value of
                                                                                // the test state using the Pauli
                                                                                // representation

        resultDiag = expValObs(&testState, &observableDiag);    // Calculate the expectation value of
                                                                                // the test state using the diagonal
                                                                                // representation

        cmatVecMulInPlace(observableMat, refKet[i], dim);                 // Calculate the reference expectation
        reference = creal(cinnerProduct(refBra[i], refKet[i], dim));      // value from matrix-vector product

        TEST_ASSERT_TRUE(fabs(result - reference) < PRECISION);
        TEST_ASSERT_TRUE(fabs(resultDiag - reference) < PRECISION);
    }

    /*
     * Free the observables
     */
    free(components);
    free(coefficients);
    free(diagObservable);
    free(observableMat);

    /*
     * Initialize the reference vectors (test vectors are not touched by expVal)
     */
    refKet = generateTestVectors(qubits);

    /*
     * Define the Pauli observable
     */
    length = dim * dim;
    components = (pauli_t *) malloc(qubits * length * sizeof(pauli_t));
    components[0] = ID;
    components[1] = X;
    components[2] = Y;
    components[3] = Z;
    coefficients = (double *) malloc(length * sizeof(double));
    coefficients[0] = -0.16123;
    coefficients[1] = 0.74276;
    coefficients[2] = -0.61566;
    coefficients[3] = 0.96877;

    /*
     * Initialize the Pauli observable
     */
    pauliObservable.components = components;
    pauliObservable.coefficients = coefficients;
    pauliObservable.length = length;

    observableMat = pauliObservableMat(components, coefficients, length, qubits);

    /*
     * Calculate the expectation values
     */

    for (dim_t i = 0; i < dim + 1; ++i) {
        stateInitVector(&testState, testVectors[i], qubits);
        result = expValObs(&testState, &observablePauli);     // Calculate the expectation value of
                                                                              // the test state

        cmatVecMulInPlace(observableMat, refKet[i], dim);               // Calculate the reference expectation
        reference = creal(cinnerProduct(refBra[i], refKet[i], dim));    // value from matrix-vector product

        TEST_ASSERT_TRUE(fabs(result - reference) < PRECISION);
    }

    /*
     * Free test vectors
     */
    freeTestVectors(testVectors, qubits);
    freeTestVectors(refBra, qubits);
    freeTestVectors(refKet, qubits);

    /*
     * Free observables
     */
    free(components);
    free(coefficients);
    free(observableMat);

    /*
                                                        2 qubits
    --------------------------------------------------------------------------------------------------------------------
    */

    qubits = 2;
    dim = POW2(qubits, dim_t);

    /*
     * Initialize the state vectors
     */
    testVectors = generateTestVectors(qubits);
    refBra = generateTestVectors(qubits);
    refKet = generateTestVectors(qubits);

    /*
     * Define the diagonal observable
     */
    length = dim;
    components = (pauli_t *) malloc(qubits * length * sizeof(pauli_t));
    components[0] = ID;
    components[1] = ID;
    components[2] = ID;
    components[3] = Z;
    components[4] = Z;
    components[5] = ID;
    components[6] = Z;
    components[7] = Z;
    coefficients = (double *) malloc(length * sizeof(double));
    coefficients[0] = -2.72482;
    coefficients[1] = -2.81765;
    coefficients[0] = -0.09561;
    coefficients[1] = -2.97345;

    /*
     * Initialize the diagonal observable
     */
    pauliObservable.components = components;            // Initialize the Pauli representation of the diagonal
    pauliObservable.coefficients = coefficients;        // hamiltonian
    pauliObservable.length = length;                    //
    pauliObservable.qubits = qubits;                    //
                                                        //
    observablePauli.type = PAULI;                       //
    observablePauli.pauliObs = &pauliObservable;        //

    diagObservable = pauliObservableDiag(components, coefficients, length, qubits);
    observableDiag.type = DIAG;                                                 // Initialize the hamiltonian's diagonal
    observableDiag.diagObs = diagObservable;                                    // representation

    observableMat = pauliObservableMat(components, coefficients, length, qubits);    // Initialize the reference
                                                                                    // matrix according to the
                                                                                    // Pauli string

    /*
     * Calculate the expectation values
     */

    for (dim_t i = 0; i < dim + 1; ++i) {
        stateInitVector(&testState, testVectors[i], qubits);
        result = expValObs(&testState, &observablePauli);       // Calculate the expectation value of
        // the test state using the Pauli
        // representation

        resultDiag = expValObs(&testState, &observableDiag);    // Calculate the expectation value of
        // the test state using the diagonal
        // representation

        cmatVecMulInPlace(observableMat, refKet[i], dim);                 // Calculate the reference expectation
        reference = creal(cinnerProduct(refBra[i], refKet[i], dim));      // value from matrix-vector product

        TEST_ASSERT_TRUE(fabs(result - reference) < PRECISION);
        TEST_ASSERT_TRUE(fabs(resultDiag - reference) < PRECISION);
    }

    /*
     * Free observables
     */
    free(components);
    free(coefficients);
    free(diagObservable);
    free(observableMat);

    /*
     * Initialize the reference vectors (test vectors are not touched by expVal)
     */
    refKet = generateTestVectors(qubits);

    /*
     * Define the Pauli observable
     */
    length = dim * dim;
    components = (pauli_t *) malloc(qubits * length * sizeof(pauli_t));
    components[0] = ID;
    components[1] = ID;
    components[2] = ID;
    components[3] = X;
    components[4] = ID;
    components[5] = Y;
    components[6] = ID;
    components[7] = Z;
    components[8] = X;
    components[9] = ID;
    components[10] = X;
    components[11] = X;
    components[12] = X;
    components[13] = Y;
    components[14] = X;
    components[15] = Z;
    components[16] = Y;
    components[17] = ID;
    components[18] = Y;
    components[19] = X;
    components[20] = Y;
    components[21] = Y;
    components[22] = Y;
    components[23] = Z;
    components[24] = Z;
    components[25] = ID;
    components[26] = Z;
    components[27] = X;
    components[28] = Z;
    components[29] = Y;
    components[30] = Z;
    components[31] = Z;
    coefficients = (double *) malloc(length * sizeof(double));
    coefficients[0] = -0.67065;
    coefficients[1] = 2.96533;
    coefficients[2] = -1.46769;
    coefficients[3] = -2.03648;
    coefficients[4] = -0.14236;
    coefficients[5] = 2.52121;
    coefficients[6] = -1.19394;
    coefficients[7] = 0.63006;
    coefficients[8] = 2.34530;
    coefficients[9] = -2.38711;
    coefficients[10] = 2.57155;
    coefficients[11] = -1.83499;
    coefficients[12] = 0.20876;
    coefficients[13] = 2.81501;
    coefficients[14] = -0.86782;
    coefficients[15] = -2.98010;

    /*
     * Initialize the Pauli observable
     */
    pauliObservable.components = components;
    pauliObservable.coefficients = coefficients;
    pauliObservable.length = length;

    observableMat = pauliObservableMat(components, coefficients, length, qubits);

    /*
     * Calculate the expectation values
     */

    for (dim_t i = 0; i < dim + 1; ++i) {
        stateInitVector(&testState, testVectors[i], qubits);
        result = expValObs(&testState, &observablePauli);     // Calculate the expectation value of
        // the test state

        cmatVecMulInPlace(observableMat, refKet[i], dim);               // Calculate the reference expectation
        reference = creal(cinnerProduct(refBra[i], refKet[i], dim));    // value from matrix-vector product

        TEST_ASSERT_TRUE(fabs(result - reference) < PRECISION);
    }

    /*
     * Free test vectors
     */
    freeTestVectors(testVectors, qubits);
    freeTestVectors(refBra, qubits);
    freeTestVectors(refKet, qubits);

    /*
     * Free observables
     */
    free(components);
    free(coefficients);
    free(observableMat);

    /*
                                                        3 qubits
    --------------------------------------------------------------------------------------------------------------------
    */

    qubits = 3;
    dim = POW2(qubits, dim_t);

    /*
     * Initialize the state vectors
     */
    testVectors = generateTestVectors(qubits);
    refBra = generateTestVectors(qubits);
    refKet = generateTestVectors(qubits);

    /*
     * Define the diagonal observable
     */
    length = dim;
    components = (pauli_t *) malloc(qubits * length * sizeof(pauli_t));
    components[0] = ID;
    components[1] = ID;
    components[2] = ID;
    components[3] = ID;
    components[4] = ID;
    components[5] = Z;
    components[6] = ID;
    components[7] = Z;
    components[8] = ID;
    components[9] = ID;
    components[10] = Z;
    components[11] = Z;
    components[12] = Z;
    components[13] = ID;
    components[14] = ID;
    components[15] = Z;
    components[16] = ID;
    components[17] = Z;
    components[18] = Z;
    components[19] = Z;
    components[20] = ID;
    components[21] = Z;
    components[22] = Z;
    components[23] = Z;
    coefficients = (double *) malloc(length * sizeof(double));
    coefficients[0] = -0.46007;
    coefficients[1] = -1.00870;
    coefficients[2] = -1.30764;
    coefficients[3] = 1.04407;
    coefficients[4] = 2.48997;
    coefficients[5] = 0.20652;
    coefficients[6] = -0.02780;
    coefficients[7] = 2.37344;

    /*
     * Initialize the diagonal observable
     */
    pauliObservable.components = components;            // Initialize the Pauli representation of the diagonal
    pauliObservable.coefficients = coefficients;        // hamiltonian
    pauliObservable.length = length;                    //
    pauliObservable.qubits = qubits;                    //
                                                        //
    observablePauli.type = PAULI;                       //
    observablePauli.pauliObs = &pauliObservable;        //

    diagObservable = pauliObservableDiag(components, coefficients, length, qubits);
    observableDiag.type = DIAG;                                                 // Initialize the hamiltonian's diagonal
    observableDiag.diagObs = diagObservable;                                    // representation

    observableMat = pauliObservableMat(components, coefficients, length, qubits);    // Initialize the reference
                                                                                    // matrix according to the
                                                                                    // Pauli string

    /*
     * Calculate the expectation values
     */

    for (dim_t i = 0; i < dim + 1; ++i) {
        stateInitVector(&testState, testVectors[i], qubits);
        result = expValObs(&testState, &observablePauli);       // Calculate the expectation value of
                                                                                // the test state using the Pauli
                                                                                // representation

        resultDiag = expValObs(&testState, &observableDiag);    // Calculate the expectation value of
                                                                                // the test state using the diagonal
                                                                                // representation

        cmatVecMulInPlace(observableMat, refKet[i], dim);                 // Calculate the reference expectation
        reference = creal(cinnerProduct(refBra[i], refKet[i], dim));      // value from matrix-vector product

        TEST_ASSERT_TRUE(fabs(result - reference) < PRECISION);
        TEST_ASSERT_TRUE(fabs(resultDiag - reference) < PRECISION);
    }

    /*
     * Free observables
     */
    free(components);
    free(coefficients);
    free(diagObservable);
    free(observableMat);

    /*
     * Initialize the reference vectors (test vectors are not touched by expVal)
     */
    refKet = generateTestVectors(qubits);

    /*
     * Define the Pauli observable
     */
    length = 2 * dim;
    components = (pauli_t *) malloc(qubits * length * sizeof(pauli_t));
    components[0] = ID;
    components[1] = ID;
    components[2] = ID;
    components[3] = X;
    components[4] = ID;
    components[5] = Y;
    components[6] = ID;
    components[7] = Z;
    components[8] = X;
    components[9] = ID;
    components[10] = X;
    components[11] = X;
    components[12] = X;
    components[13] = Y;
    components[14] = X;
    components[15] = Z;
    components[16] = Y;
    components[17] = ID;
    components[18] = Y;
    components[19] = X;
    components[20] = Y;
    components[21] = Y;
    components[22] = Y;
    components[23] = Z;
    components[24] = Z;
    components[25] = ID;
    components[26] = Z;
    components[27] = X;
    components[28] = Z;
    components[29] = Y;
    components[30] = Z;
    components[31] = Z;
    components[32] = ID;
    components[33] = ID;
    components[34] = Y;
    components[35] = X;
    components[36] = ID;
    components[37] = Y;
    components[38] = Y;
    components[39] = ID;
    components[40] = Z;
    components[41] = ID;
    components[42] = Z;
    components[43] = X;
    components[44] = ID;
    components[45] = Y;
    components[46] = Z;
    components[47] = Z;
    coefficients = (double *) malloc(length * sizeof(double));
    coefficients[0] = 2.16540;
    coefficients[1] = 1.69670;
    coefficients[2] = 0.19151;
    coefficients[3] = 1.96392;
    coefficients[4] = -0.10861;
    coefficients[5] = 1.72947;
    coefficients[6] = -0.97590;
    coefficients[7] = 2.04988;
    coefficients[8] = -2.14632;
    coefficients[9] = 1.89464;
    coefficients[10] = 2.16018;
    coefficients[11] = -2.71582;
    coefficients[12] = -0.70224;
    coefficients[13] = 0.54424;
    coefficients[14] = -0.75072;
    coefficients[15] = -0.64382;

    /*
     * Initialize the Pauli observable
     */
    pauliObservable.components = components;
    pauliObservable.coefficients = coefficients;
    pauliObservable.length = length;

    observableMat = pauliObservableMat(components, coefficients, length, qubits);

    /*
     * Calculate the expectation values
     */

    for (dim_t i = 0; i < dim + 1; ++i) {
        stateInitVector(&testState, testVectors[i], qubits);
        result = expValObs(&testState, &observablePauli);     // Calculate the expectation value of
                                                                              // the test state

        cmatVecMulInPlace(observableMat, refKet[i], dim);               // Calculate the reference expectation
        reference = creal(cinnerProduct(refBra[i], refKet[i], dim));    // value from matrix-vector product

        TEST_ASSERT_TRUE(fabs(result - reference) < PRECISION);
    }

    /*
     * Free test vectors
     */
    freeTestVectors(testVectors, qubits);
    freeTestVectors(refBra, qubits);
    freeTestVectors(refKet, qubits);

    /*
     * Free observables
     */
    free(components);
    free(coefficients);
    free(observableMat);

    /*
                                                4 qubit
    ------------------------------------------------------------------------------------------------
    */

    qubits = 4;
    dim = POW2(qubits, dim_t);

    /*
     * Initialize the state vectors
     */
    testVectors = generateTestVectors(qubits);
    refBra = generateTestVectors(qubits);
    refKet = generateTestVectors(qubits);

    /*
     * Define the diagonal observable
     */
    length = dim;
    components = (pauli_t *) malloc(qubits * length * sizeof(pauli_t));
    components[0] = ID;
    components[1] = ID;
    components[2] = ID;
    components[3] = ID;
    components[4] = ID;
    components[5] = ID;
    components[6] = ID;
    components[7] = Z;
    components[8] = ID;
    components[9] = ID;
    components[10] = Z;
    components[11] = ID;
    components[12] = ID;
    components[13] = ID;
    components[14] = Z;
    components[15] = Z;
    components[16] = ID;
    components[17] = Z;
    components[18] = ID;
    components[19] = ID;
    components[20] = ID;
    components[21] = Z;
    components[22] = ID;
    components[23] = Z;
    components[24] = ID;
    components[25] = Z;
    components[26] = Z;
    components[27] = ID;
    components[28] = ID;
    components[29] = Z;
    components[30] = Z;
    components[31] = Z;
    components[32] = Z;
    components[33] = ID;
    components[34] = ID;
    components[35] = ID;
    components[36] = Z;
    components[37] = ID;
    components[38] = ID;
    components[39] = Z;
    components[40] = Z;
    components[41] = ID;
    components[42] = Z;
    components[43] = ID;
    components[44] = Z;
    components[45] = ID;
    components[46] = Z;
    components[47] = Z;
    components[48] = Z;
    components[49] = Z;
    components[50] = ID;
    components[51] = ID;
    components[52] = Z;
    components[53] = Z;
    components[54] = ID;
    components[55] = Z;
    components[56] = Z;
    components[57] = Z;
    components[58] = Z;
    components[59] = ID;
    components[60] = Z;
    components[61] = Z;
    components[62] = Z;
    components[63] = Z;
    coefficients = (double *) malloc(length * sizeof(double));
    coefficients[0] = -1.79196;
    coefficients[1] = 1.46429;
    coefficients[2] = -1.75880;
    coefficients[3] = 0.67855;
    coefficients[4] = -1.82106;
    coefficients[5] = -0.57601;
    coefficients[6] = 2.22127;
    coefficients[7] = 1.42834;
    coefficients[8] = 1.90210;
    coefficients[9] = 1.04138;
    coefficients[10] = 2.25605;
    coefficients[11] = -1.84384;
    coefficients[12] = 1.86982;
    coefficients[13] = -2.34669;
    coefficients[14] = -0.17885;
    coefficients[15] = 2.51645;

    /*
     * Initialize the diagonal observable
     */
    pauliObservable.components = components;            // Initialize the Pauli representation of the diagonal
    pauliObservable.coefficients = coefficients;        // hamiltonian
    pauliObservable.length = length;                    //
    pauliObservable.qubits = qubits;                    //
    //
    observablePauli.type = PAULI;                       //
    observablePauli.pauliObs = &pauliObservable;        //

    diagObservable = pauliObservableDiag(components, coefficients, length, qubits);
    observableDiag.type = DIAG;                                                 // Initialize the hamiltonian's diagonal
    observableDiag.diagObs = diagObservable;                                    // representation

    observableMat = pauliObservableMat(components, coefficients, length, qubits);    // Initialize the reference
    // matrix according to the
    // Pauli string

    /*
     * Calculate the expectation values
     */
    for (dim_t i = 0; i < dim + 1; ++i) {
        stateInitVector(&testState, testVectors[i], qubits);
        result = expValObs(&testState, &observablePauli);       // Calculate the expectation value of
                                                                                // the test state using the Pauli
                                                                                // representation

        resultDiag = expValObs(&testState, &observableDiag);    // Calculate the expectation value of
                                                                                // the test state using the diagonal
                                                                                // representation

        cmatVecMulInPlace(observableMat, refKet[i], dim);                 // Calculate the reference expectation
        reference = creal(cinnerProduct(refBra[i], refKet[i], dim));      // value from matrix-vector product

        TEST_ASSERT_TRUE(fabs(result - reference) < PRECISION);
        TEST_ASSERT_TRUE(fabs(resultDiag - reference) < PRECISION);
    }

    /*
     * Free observables
     */
    free(components);
    free(coefficients);
    free(diagObservable);
    free(observableMat);

    /*
     * Initialize the reference vectors (test vectors are not touched by expVal)
     */
    refKet = generateTestVectors(qubits);

    /*
     * Define the Pauli observable
     */
    length = 12;
    components = (pauli_t *) malloc(qubits * length * sizeof(pauli_t));
    components[0] = ID;
    components[1] = ID;
    components[2] = ID;
    components[3] = X;
    components[4] = ID;
    components[5] = Y;
    components[6] = ID;
    components[7] = Z;
    components[8] = X;
    components[9] = ID;
    components[10] = X;
    components[11] = X;
    components[12] = X;
    components[13] = Y;
    components[14] = X;
    components[15] = Z;
    components[16] = Y;
    components[17] = ID;
    components[18] = Y;
    components[19] = X;
    components[20] = Y;
    components[21] = Y;
    components[22] = Y;
    components[23] = Z;
    components[24] = Z;
    components[25] = ID;
    components[26] = Z;
    components[27] = X;
    components[28] = Z;
    components[29] = Y;
    components[30] = Z;
    components[31] = Z;
    components[32] = ID;
    components[33] = ID;
    components[34] = Y;
    components[35] = X;
    components[36] = ID;
    components[37] = Y;
    components[38] = Y;
    components[39] = ID;
    components[40] = Z;
    components[41] = ID;
    components[42] = Z;
    components[43] = X;
    components[44] = ID;
    components[45] = Y;
    components[46] = Z;
    components[47] = Z;
    coefficients = (double *) malloc(length * sizeof(double));
    coefficients[0] = -1.39395;
    coefficients[1] = 2.08556;
    coefficients[2] = -0.92705;
    coefficients[3] = -1.16536;
    coefficients[4] = -0.89633;
    coefficients[5] = 0.16350;
    coefficients[6] = -0.92537;
    coefficients[7] = -2.49582;
    coefficients[8] = 0.50758;
    coefficients[9] = -1.73637;
    coefficients[10] = -1.27648;
    coefficients[11] = -1.87498;

    /*
     * Initialize the Pauli observable
     */
    pauliObservable.components = components;
    pauliObservable.coefficients = coefficients;
    pauliObservable.length = length;

    observableMat = pauliObservableMat(components, coefficients, length, qubits);

    /*
     * Calculate the expectation values
     */

    for (dim_t i = 0; i < dim + 1; ++i) {
        stateInitVector(&testState, testVectors[i], qubits);
        result = expValObs(&testState, &observablePauli);     // Calculate the expectation value of
                                                                              // the test state

        cmatVecMulInPlace(observableMat, refKet[i], dim);               // Calculate the reference expectation
        reference = creal(cinnerProduct(refBra[i], refKet[i], dim));    // value from matrix-vector product

        TEST_ASSERT_TRUE(fabs(result - reference) < PRECISION);
    }

    /*
     * Free test vectors
     */
    freeTestVectors(testVectors, qubits);
    freeTestVectors(refBra, qubits);
    freeTestVectors(refKet, qubits);

    /*
     * Free observables
     */
    free(components);
    free(coefficients);
    free(observableMat);
}

/* 
 * =====================================================================================================================
 *                                                      test gradientPQC
 * =====================================================================================================================
 */

void testGradientPQC(void) {
    qubit_t qubits;
    dim_t dim;
    depth_t circdepth;
    double* parameters;

    state_t testState;
    cplx_t** testVectors;
    cplx_t** refVectors;

    complength_t lengthObs;
    pauli_t* compObs;
    double* coeffObs;
    pauliObs_t pauliObservable;
    obs_t observable;

    complength_t lengthEvoOps;
    pauli_t* compEvoOps;
    double* coeffEvoOps;
    pauliObs_t* evoOpsPauli;
    obs_t** evoOps;

    /*
                                                3 qubit, circuit depth = 2
    --------------------------------------------------------------------------------------------------------------------
    */
    qubits = 3;
    dim = POW2(qubits, dim_t);

    /*
     * Define the PQC setting
     */
    circdepth = 2;
    parameters = (double*) malloc(circdepth * sizeof(double));
    parameters[0] = -2.62501;
    parameters[1] = 0.81390;

    /*
     * Initialize the state vectors
     */
    testVectors = generateTestVectors(qubits);
    refVectors = generateTestVectors(qubits);

    /*
     * Define the observable
     */
    lengthObs = 2 * dim;
    compObs = (pauli_t *) malloc(qubits * lengthObs * sizeof(pauli_t));
    compObs[0] = ID;
    compObs[1] = ID;
    compObs[2] = ID;
    compObs[3] = X;
    compObs[4] = ID;
    compObs[5] = Y;
    compObs[6] = ID;
    compObs[7] = Z;
    compObs[8] = X;
    compObs[9] = ID;
    compObs[10] = X;
    compObs[11] = X;
    compObs[12] = X;
    compObs[13] = Y;
    compObs[14] = X;
    compObs[15] = Z;
    compObs[16] = Y;
    compObs[17] = ID;
    compObs[18] = Y;
    compObs[19] = X;
    compObs[20] = Y;
    compObs[21] = Y;
    compObs[22] = Y;
    compObs[23] = Z;
    compObs[24] = Z;
    compObs[25] = ID;
    compObs[26] = Z;
    compObs[27] = X;
    compObs[28] = Z;
    compObs[29] = Y;
    compObs[30] = Z;
    compObs[31] = Z;
    compObs[32] = ID;
    compObs[33] = ID;
    compObs[34] = Y;
    compObs[35] = X;
    compObs[36] = ID;
    compObs[37] = Y;
    compObs[38] = Y;
    compObs[39] = ID;
    compObs[40] = Z;
    compObs[41] = ID;
    compObs[42] = Z;
    compObs[43] = X;
    compObs[44] = ID;
    compObs[45] = Y;
    compObs[46] = Z;
    compObs[47] = Z;
    coeffObs = (double *) malloc(lengthObs * sizeof(double));
    coeffObs[0] = 2.16540;
    coeffObs[1] = 1.69670;
    coeffObs[2] = 0.19151;
    coeffObs[3] = 1.96392;
    coeffObs[4] = -0.10861;
    coeffObs[5] = 1.72947;
    coeffObs[6] = -0.97590;
    coeffObs[7] = 2.04988;
    coeffObs[8] = -2.14632;
    coeffObs[9] = 1.89464;
    coeffObs[10] = 2.16018;
    coeffObs[11] = -2.71582;
    coeffObs[12] = -0.70224;
    coeffObs[13] = 0.54424;
    coeffObs[14] = -0.75072;
    coeffObs[15] = -0.64382;

    /*
     * Initialize the Pauli observable
     */
    pauliObservable.components = compObs;
    pauliObservable.coefficients = coeffObs;
    pauliObservable.qubits = qubits;
    pauliObservable.length = lengthObs;

    observable.type = PAULI;
    observable.pauliObs = &pauliObservable;

    /*
     * Define the evolution operators
     */
    lengthEvoOps = dim;
    compEvoOps = (pauli_t*) malloc(qubits * lengthEvoOps * sizeof(pauli_t));
    compEvoOps[0] = ID;
    compEvoOps[1] = ID;
    compEvoOps[2] = ID;
    compEvoOps[3] = ID;
    compEvoOps[4] = ID;
    compEvoOps[5] = Z;
    compEvoOps[6] = ID;
    compEvoOps[7] = Z;
    compEvoOps[8] = ID;
    compEvoOps[9] = ID;
    compEvoOps[10] = Z;
    compEvoOps[11] = Z;
    compEvoOps[12] = Z;
    compEvoOps[13] = ID;
    compEvoOps[14] = ID;
    compEvoOps[15] = Z;
    compEvoOps[16] = ID;
    compEvoOps[17] = Z;
    compEvoOps[18] = Z;
    compEvoOps[19] = Z;
    compEvoOps[20] = ID;
    compEvoOps[21] = Z;
    compEvoOps[22] = Z;
    compEvoOps[23] = Z;
    coeffEvoOps = (double*) malloc(circdepth * lengthEvoOps * sizeof(double));
    coeffEvoOps[0] = -0.39448;                                              // array is a concatenation of coefficients
    coeffEvoOps[1] = 1.70424;                                               // for all evolution operators
    coeffEvoOps[2] = -1.83775;
    coeffEvoOps[3] = 2.39278;
    coeffEvoOps[4] = 1.23840;
    coeffEvoOps[5] = 0.89824;
    coeffEvoOps[6] = -1.86435;
    coeffEvoOps[7] = 2.88865;
    coeffEvoOps[8] = -1.30915;
    coeffEvoOps[9] = 0.15139;
    coeffEvoOps[10] = 0.60963;
    coeffEvoOps[11] = -2.89646;
    coeffEvoOps[12] = -0.18224;
    coeffEvoOps[13] = 0.89809;
    coeffEvoOps[14] = -0.46685;
    coeffEvoOps[15] = 2.87469;

    /*
     * Initialize the evolution operators
     */
    evoOpsPauli = (pauliObs_t*) malloc(circdepth * sizeof(pauliObs_t));
    for (depth_t i = 0; i < circdepth; ++i) {
        evoOpsPauli[i].components = compEvoOps;
        evoOpsPauli[i].coefficients = coeffEvoOps + (i * lengthEvoOps);
        evoOpsPauli[i].qubits = qubits;
        evoOpsPauli[i].length = lengthEvoOps;
    }

    evoOps = (obs_t**) malloc(circdepth * sizeof(obs_t*));
    for (depth_t i = 0; i < circdepth; ++i) {
        evoOps[i] = (obs_t*) malloc(sizeof(obs_t));
        evoOps[i]->type = PAULI;
        evoOps[i]->pauliObs = evoOpsPauli + i;
    }

    /*
     * Calculate the gradients
     */
    for (dim_t i = 0; i < dim + 1; ++i) {
        /*
         * Gradient from gradient PQC
         */
        double* result;
        stateInitVector(&testState, testVectors[i], qubits);
        result = gradientPQC(&testState, parameters, &observable, (const obs_t **) evoOps, circdepth);

        /*
         * Gradient from finite difference method
         */
        double* reference = finiteGradientPQC(refVectors[i], qubits, dim, circdepth, compObs, coeffObs, lengthObs, \
                                              compEvoOps, coeffEvoOps, lengthEvoOps, parameters, 1e-9);

//        double* reference2 = (double*) malloc(circdepth * sizeof(double));
//        for (depth_t j = 0; j < circdepth; ++j) {
//            refEvoOps[j] = expTrotterizedPauliObservableMat(compEvoOps, \
//                                                            coeffEvoOps + (j * lengthEvoOps), \
//                                                            lengthEvoOps, \
//                                                            parameters[j], \
//                                                            qubits);
//        }
//
//        cplx_t** refOps = (cplx_t**) malloc(circdepth * sizeof(cplx_t*));
//        for (depth_t j = 0; j < circdepth; ++j) {
//            refOps[j] = pauliObservableMat(compEvoOps, \
//                                          coeffEvoOps + (j * lengthEvoOps), \
//                                          lengthEvoOps, \
//                                          qubits);
//        }
//
//        refBra = cmatVecMul(refEvoOps[0], refVectors[i], dim);
//        for (depth_t k = 1; k < circdepth; ++k) {
//            cmatVecMulInPlace(refEvoOps[k], refBra, dim);
//        }
//        cmatVecMulInPlace(observableMat, refBra, dim);
//
//        for (depth_t k = 0; k < circdepth; ++k) {
//            refKet = cmatVecMul(refEvoOps[0], refVectors[i], dim);
//            for (depth_t kk = 1; kk <= k; ++kk) {
//                cmatVecMulInPlace(refEvoOps[kk], refKet, dim);
//            }
//            cmatVecMulInPlace(refOps[k], refKet, dim);
//            for (depth_t kk = k + 1; kk < circdepth; ++kk) {
//                cmatVecMulInPlace(refEvoOps[kk], refKet, dim);
//            }
//            reference2[k] = 2 * cimag(cinnerProduct(refBra, refKet, dim));
//        }

        TEST_ASSERT_TRUE(rvectorAlmostEqual(result, reference, circdepth, APPROXPRECISION));

        free(result);
        free(reference);
    }

    /*
    * Free test vectors
    */
    freeTestVectors(testVectors, qubits);
    freeTestVectors(refVectors, qubits);

    /*
     * Free observable and evolution operators
     */
    free(parameters);
    free(compObs);
    free(coeffObs);
    free(compEvoOps);
    free(coeffEvoOps);
    free(evoOpsPauli);
    for (depth_t i = 0; i < circdepth; ++i) {
        free(evoOps[i]);
    }
    free(evoOps);

    /*
                                                3 qubit, circuit depth = 3
    --------------------------------------------------------------------------------------------------------------------
    */
    qubits = 3;
    dim = POW2(qubits, dim_t);

    /*
     * Define the PQC setting
     */
    circdepth = 3;
    parameters = (double*) malloc(circdepth * sizeof(double));
    parameters[0] = 0.33305;
    parameters[1] = 0.24576;
    parameters[2] = 1.27367;

    /*
     * Initialize the state vectors
     */
    testVectors = generateTestVectors(qubits);
    refVectors = generateTestVectors(qubits);

    /*
     * Define the observable
     */
    lengthObs = 2 * dim;
    compObs = (pauli_t *) malloc(qubits * lengthObs * sizeof(pauli_t));
    compObs[0] = ID;
    compObs[1] = ID;
    compObs[2] = ID;
    compObs[3] = X;
    compObs[4] = ID;
    compObs[5] = Y;
    compObs[6] = ID;
    compObs[7] = Z;
    compObs[8] = X;
    compObs[9] = ID;
    compObs[10] = X;
    compObs[11] = X;
    compObs[12] = X;
    compObs[13] = Y;
    compObs[14] = X;
    compObs[15] = Z;
    compObs[16] = Y;
    compObs[17] = ID;
    compObs[18] = Y;
    compObs[19] = X;
    compObs[20] = Y;
    compObs[21] = Y;
    compObs[22] = Y;
    compObs[23] = Z;
    compObs[24] = Z;
    compObs[25] = ID;
    compObs[26] = Z;
    compObs[27] = X;
    compObs[28] = Z;
    compObs[29] = Y;
    compObs[30] = Z;
    compObs[31] = Z;
    compObs[32] = ID;
    compObs[33] = ID;
    compObs[34] = Y;
    compObs[35] = X;
    compObs[36] = ID;
    compObs[37] = Y;
    compObs[38] = Y;
    compObs[39] = ID;
    compObs[40] = Z;
    compObs[41] = ID;
    compObs[42] = Z;
    compObs[43] = X;
    compObs[44] = ID;
    compObs[45] = Y;
    compObs[46] = Z;
    compObs[47] = Z;
    coeffObs = (double *) malloc(lengthObs * sizeof(double));
    coeffObs[0] = -0.10682;
    coeffObs[1] = -0.47524;
    coeffObs[2] = 0.91975;
    coeffObs[3] = 1.58143;
    coeffObs[4] = 0.67805;
    coeffObs[5] = -0.64653;
    coeffObs[6] = -1.81747;
    coeffObs[7] = 0.87833;
    coeffObs[8] = -0.93854;
    coeffObs[9] = 1.09448;
    coeffObs[10] = 1.05518;
    coeffObs[11] = -1.12855;
    coeffObs[12] = -0.57662;
    coeffObs[13] = -1.58328;
    coeffObs[14] = 2.65643;
    coeffObs[15] = -0.03557;

    /*
     * Initialize the Pauli observable
     */
    pauliObservable.components = compObs;
    pauliObservable.coefficients = coeffObs;
    pauliObservable.qubits = qubits;
    pauliObservable.length = lengthObs;

    observable.type = PAULI;
    observable.pauliObs = &pauliObservable;

    /*
     * Define the evolution operators
     */
    lengthEvoOps = dim;
    compEvoOps = (pauli_t*) malloc(qubits * lengthEvoOps * sizeof(pauli_t));
    compEvoOps[0] = ID;
    compEvoOps[1] = ID;
    compEvoOps[2] = ID;
    compEvoOps[3] = ID;
    compEvoOps[4] = ID;
    compEvoOps[5] = Z;
    compEvoOps[6] = ID;
    compEvoOps[7] = Z;
    compEvoOps[8] = ID;
    compEvoOps[9] = ID;
    compEvoOps[10] = Z;
    compEvoOps[11] = Z;
    compEvoOps[12] = Z;
    compEvoOps[13] = ID;
    compEvoOps[14] = ID;
    compEvoOps[15] = Z;
    compEvoOps[16] = ID;
    compEvoOps[17] = Z;
    compEvoOps[18] = Z;
    compEvoOps[19] = Z;
    compEvoOps[20] = ID;
    compEvoOps[21] = Z;
    compEvoOps[22] = Z;
    compEvoOps[23] = Z;
    coeffEvoOps = (double*) malloc(circdepth * lengthEvoOps * sizeof(double));
    coeffEvoOps[0] = -0.75749;                                             // array is a concatenation of coefficients
    coeffEvoOps[1] = 1.14340;                                              // for all evolution operators
    coeffEvoOps[2] = -0.16411;
    coeffEvoOps[3] = -2.59214;
    coeffEvoOps[4] = 1.37578;
    coeffEvoOps[5] = 1.33389;
    coeffEvoOps[6] = 1.87961;
    coeffEvoOps[7] = -0.48781;
    coeffEvoOps[8] = -0.73497;
    coeffEvoOps[9] = 2.10394;
    coeffEvoOps[10] = -1.67602;
    coeffEvoOps[11] = -2.17367;
    coeffEvoOps[12] = 0.48264;
    coeffEvoOps[13] = 2.56591;
    coeffEvoOps[14] = 1.91342;
    coeffEvoOps[15] = -1.19817;
    coeffEvoOps[16] = 2.39568;
    coeffEvoOps[17] = 1.34672;
    coeffEvoOps[18] = 0.28640;
    coeffEvoOps[19] = -0.49067;
    coeffEvoOps[20] = -0.07299;
    coeffEvoOps[21] = 0.34801;
    coeffEvoOps[22] = 2.87215;
    coeffEvoOps[23] = -0.31542;

    /*
     * Initialize the evolution operators
     */
    evoOpsPauli = (pauliObs_t*) malloc(circdepth * sizeof(pauliObs_t));
    for (depth_t i = 0; i < circdepth; ++i) {
        evoOpsPauli[i].components = compEvoOps;
        evoOpsPauli[i].coefficients = coeffEvoOps + (i * lengthEvoOps);
        evoOpsPauli[i].qubits = qubits;
        evoOpsPauli[i].length = lengthEvoOps;
    }

    evoOps = (obs_t**) malloc(circdepth * sizeof(obs_t));
    for (depth_t i = 0; i < circdepth; ++i) {
        evoOps[i] = (obs_t*) malloc(sizeof(obs_t));
        evoOps[i]->type = PAULI;
        evoOps[i]->pauliObs = evoOpsPauli + i;
    }

    /*
     * Calculate the gradients
     */
    for (dim_t i = 0; i < dim + 1; ++i) {
        /*
         * Gradient from gradient PQC
         */
        double* result;
        stateInitVector(&testState, testVectors[i], qubits);
        result = gradientPQC(&testState, parameters, &observable, (const obs_t **) evoOps, circdepth);

        /*
         * Gradient from finite difference method
         */
        double* reference = finiteGradientPQC(refVectors[i], qubits, dim, circdepth, compObs, coeffObs, lengthObs, \
                                              compEvoOps, coeffEvoOps, lengthEvoOps, parameters, 1e-9);

        TEST_ASSERT_TRUE(rvectorAlmostEqual(result, reference, circdepth, APPROXPRECISION));
    }

    /*
    * Free test vectors
    */
    freeTestVectors(testVectors, qubits);
    freeTestVectors(refVectors, qubits);

    /*
     * Free observable and evolution operators
     */
    free(compObs);
    free(coeffObs);
    free(compEvoOps);
    free(coeffEvoOps);
    free(evoOpsPauli);
    for (depth_t i = 0; i < circdepth; ++i) {
        free(evoOps[i]);
    }
    free(evoOps);

    /*
                                                4 qubit, circuit depth = 2
    --------------------------------------------------------------------------------------------------------------------
    */
    qubits = 4;
    dim = POW2(qubits, dim_t);

    /*
     * Define the PQC setting
     */
    circdepth = 2;
    parameters = (double*) malloc(circdepth * sizeof(double));
    parameters[0] = -2.62785;
    parameters[1] = -1.30088;

    /*
     * Initialize the state vectors
     */
    testVectors = generateTestVectors(qubits);
    refVectors = generateTestVectors(qubits);

    /*
     * Define the observable
     */
    lengthObs = 12;
    compObs = (pauli_t *) malloc(qubits * lengthObs * sizeof(pauli_t));
    compObs[0] = ID;
    compObs[1] = ID;
    compObs[2] = ID;
    compObs[3] = X;
    compObs[4] = ID;
    compObs[5] = Y;
    compObs[6] = ID;
    compObs[7] = Z;
    compObs[8] = X;
    compObs[9] = ID;
    compObs[10] = X;
    compObs[11] = X;
    compObs[12] = X;
    compObs[13] = Y;
    compObs[14] = X;
    compObs[15] = Z;
    compObs[16] = Y;
    compObs[17] = ID;
    compObs[18] = Y;
    compObs[19] = X;
    compObs[20] = Y;
    compObs[21] = Y;
    compObs[22] = Y;
    compObs[23] = Z;
    compObs[24] = Z;
    compObs[25] = ID;
    compObs[26] = Z;
    compObs[27] = X;
    compObs[28] = Z;
    compObs[29] = Y;
    compObs[30] = Z;
    compObs[31] = Z;
    compObs[32] = ID;
    compObs[33] = ID;
    compObs[34] = Y;
    compObs[35] = X;
    compObs[36] = ID;
    compObs[37] = Y;
    compObs[38] = Y;
    compObs[39] = ID;
    compObs[40] = Z;
    compObs[41] = ID;
    compObs[42] = Z;
    compObs[43] = X;
    compObs[44] = ID;
    compObs[45] = Y;
    compObs[46] = Z;
    compObs[47] = Z;
    coeffObs = (double *) malloc(lengthObs * sizeof(double));
    coeffObs[0] = 0.83694;
    coeffObs[1] = 2.59855;
    coeffObs[2] = 1.14935;
    coeffObs[3] = 2.34602;
    coeffObs[4] = 2.50542;
    coeffObs[5] = -2.47148;
    coeffObs[6] = 1.12971;
    coeffObs[7] = 0.73640;
    coeffObs[8] = 2.01390;
    coeffObs[9] = 2.44850;
    coeffObs[10] = -1.82246;
    coeffObs[11] = 0.61758;

    /*
     * Initialize the Pauli observable
     */
    pauliObservable.components = compObs;
    pauliObservable.coefficients = coeffObs;
    pauliObservable.qubits = qubits;
    pauliObservable.length = lengthObs;

    observable.type = PAULI;
    observable.pauliObs = &pauliObservable;

    /*
     * Define the evolution operators
     */
    lengthEvoOps = dim;
    compEvoOps = (pauli_t*) malloc(qubits * lengthEvoOps * sizeof(pauli_t));
    compEvoOps[0] = ID;
    compEvoOps[1] = ID;
    compEvoOps[2] = ID;
    compEvoOps[3] = ID;
    compEvoOps[4] = ID;
    compEvoOps[5] = ID;
    compEvoOps[6] = ID;
    compEvoOps[7] = Z;
    compEvoOps[8] = ID;
    compEvoOps[9] = ID;
    compEvoOps[10] = Z;
    compEvoOps[11] = ID;
    compEvoOps[12] = ID;
    compEvoOps[13] = ID;
    compEvoOps[14] = Z;
    compEvoOps[15] = Z;
    compEvoOps[16] = ID;
    compEvoOps[17] = Z;
    compEvoOps[18] = ID;
    compEvoOps[19] = ID;
    compEvoOps[20] = ID;
    compEvoOps[21] = Z;
    compEvoOps[22] = ID;
    compEvoOps[23] = Z;
    compEvoOps[24] = ID;
    compEvoOps[25] = Z;
    compEvoOps[26] = Z;
    compEvoOps[27] = ID;
    compEvoOps[28] = ID;
    compEvoOps[29] = Z;
    compEvoOps[30] = Z;
    compEvoOps[31] = Z;
    compEvoOps[32] = Z;
    compEvoOps[33] = ID;
    compEvoOps[34] = ID;
    compEvoOps[35] = ID;
    compEvoOps[36] = Z;
    compEvoOps[37] = ID;
    compEvoOps[38] = ID;
    compEvoOps[39] = Z;
    compEvoOps[40] = Z;
    compEvoOps[41] = ID;
    compEvoOps[42] = Z;
    compEvoOps[43] = ID;
    compEvoOps[44] = Z;
    compEvoOps[45] = ID;
    compEvoOps[46] = Z;
    compEvoOps[47] = Z;
    compEvoOps[48] = Z;
    compEvoOps[49] = Z;
    compEvoOps[50] = ID;
    compEvoOps[51] = ID;
    compEvoOps[52] = Z;
    compEvoOps[53] = Z;
    compEvoOps[54] = ID;
    compEvoOps[55] = Z;
    compEvoOps[56] = Z;
    compEvoOps[57] = Z;
    compEvoOps[58] = Z;
    compEvoOps[59] = ID;
    compEvoOps[60] = Z;
    compEvoOps[61] = Z;
    compEvoOps[62] = Z;
    compEvoOps[63] = Z;
    coeffEvoOps = (double*) malloc(circdepth * lengthEvoOps * sizeof(double));
    coeffEvoOps[0] = 2.47766;
    coeffEvoOps[1] = -2.76095;
    coeffEvoOps[2] = 0.18111;
    coeffEvoOps[3] = 0.97866;
    coeffEvoOps[4] = -2.76576;
    coeffEvoOps[5] = -2.95371;
    coeffEvoOps[6] = -2.32907;
    coeffEvoOps[7] = -1.24555;
    coeffEvoOps[8] = 2.03087;
    coeffEvoOps[9] = -2.58241;
    coeffEvoOps[10] = 0.97024;
    coeffEvoOps[11] = -1.43111;
    coeffEvoOps[12] = 2.53653;
    coeffEvoOps[13] = 2.22199;
    coeffEvoOps[14] = -1.53345;
    coeffEvoOps[15] = -1.04142;
    coeffEvoOps[16] = -0.40504;
    coeffEvoOps[17] = -0.26551;
    coeffEvoOps[18] = -2.04636;
    coeffEvoOps[19] = 0.2012;
    coeffEvoOps[20] = 2.43244;
    coeffEvoOps[21] = 0.21154;
    coeffEvoOps[22] = 1.68992;
    coeffEvoOps[23] = -0.36902;
    coeffEvoOps[24] = 0.74709;
    coeffEvoOps[25] = -0.56174;
    coeffEvoOps[26] = -2.93947;
    coeffEvoOps[27] = -0.64292;
    coeffEvoOps[28] = -1.3835;
    coeffEvoOps[29] = 2.29929;
    coeffEvoOps[30] = -2.37431;
    coeffEvoOps[31] = 0.7905;

    /*
     * Initialize the evolution operators
     */
    evoOpsPauli = (pauliObs_t*) malloc(circdepth * sizeof(pauliObs_t));
    for (depth_t i = 0; i < circdepth; ++i) {
        evoOpsPauli[i].components = compEvoOps;
        evoOpsPauli[i].coefficients = coeffEvoOps + (i * lengthEvoOps);
        evoOpsPauli[i].qubits = qubits;
        evoOpsPauli[i].length = lengthEvoOps;
    }

    evoOps = (obs_t**) malloc(circdepth * sizeof(obs_t));
    for (depth_t i = 0; i < circdepth; ++i) {
        evoOps[i] = (obs_t*) malloc(sizeof(obs_t));
        evoOps[i]->type = PAULI;
        evoOps[i]->pauliObs = evoOpsPauli + i;
    }

    /*
     * Calculate the gradients
     */
    for (dim_t i = 0; i < dim + 1; ++i) {
        /*
 * Gradient from gradient PQC
 */
        double* result;
        stateInitVector(&testState, testVectors[i], qubits);
        result = gradientPQC(&testState, parameters, &observable, (const obs_t **) evoOps, circdepth);

        /*
         * Gradient from finite difference method
         */
        double* reference = finiteGradientPQC(refVectors[i], qubits, dim, circdepth, compObs, coeffObs, lengthObs, \
                                              compEvoOps, coeffEvoOps, lengthEvoOps, parameters, 1e-9);

        TEST_ASSERT_TRUE(rvectorAlmostEqual(result, reference, circdepth, APPROXPRECISION));
    }

    /*
    * Free test vectors
    */
    freeTestVectors(testVectors, qubits);
    freeTestVectors(refVectors, qubits);

    /*
     * Free observable and evolution operators
     */
    free(compObs);
    free(coeffObs);
    free(compEvoOps);
    free(coeffEvoOps);
    free(evoOpsPauli);
    for (depth_t i = 0; i < circdepth; ++i) {
        free(evoOps[i]);
    }
    free(evoOps);

    /*
                                                4 qubit, circuit depth = 5
    --------------------------------------------------------------------------------------------------------------------
    */
    qubits = 4;
    dim = POW2(qubits, dim_t);

    /*
     * Define the PQC setting
     */
    circdepth = 5;
    parameters = (double*) malloc(circdepth * sizeof(double));
    parameters[0] = -1.31679;
    parameters[1] = 2.91097;
    parameters[2] = -0.06082;
    parameters[3] = 0.13077;
    parameters[4] = -0.24145;

    /*
     * Initialize the state vectors
     */
    testVectors = generateTestVectors(qubits);
    refVectors = generateTestVectors(qubits);

    /*
     * Define the observable
     */
    lengthObs = 12;
    compObs = (pauli_t *) malloc(qubits * lengthObs * sizeof(pauli_t));
    compObs[0] = ID;
    compObs[1] = ID;
    compObs[2] = ID;
    compObs[3] = X;
    compObs[4] = ID;
    compObs[5] = Y;
    compObs[6] = ID;
    compObs[7] = Z;
    compObs[8] = X;
    compObs[9] = ID;
    compObs[10] = X;
    compObs[11] = X;
    compObs[12] = X;
    compObs[13] = Y;
    compObs[14] = X;
    compObs[15] = Z;
    compObs[16] = Y;
    compObs[17] = ID;
    compObs[18] = Y;
    compObs[19] = X;
    compObs[20] = Y;
    compObs[21] = Y;
    compObs[22] = Y;
    compObs[23] = Z;
    compObs[24] = Z;
    compObs[25] = ID;
    compObs[26] = Z;
    compObs[27] = X;
    compObs[28] = Z;
    compObs[29] = Y;
    compObs[30] = Z;
    compObs[31] = Z;
    compObs[32] = ID;
    compObs[33] = ID;
    compObs[34] = Y;
    compObs[35] = X;
    compObs[36] = ID;
    compObs[37] = Y;
    compObs[38] = Y;
    compObs[39] = ID;
    compObs[40] = Z;
    compObs[41] = ID;
    compObs[42] = Z;
    compObs[43] = X;
    compObs[44] = ID;
    compObs[45] = Y;
    compObs[46] = Z;
    compObs[47] = Z;
    coeffObs = (double *) malloc(lengthObs * sizeof(double));
    coeffObs[0] = -2.20587;
    coeffObs[1] = -2.13108;
    coeffObs[2] = -2.55177;
    coeffObs[3] = -0.68342;
    coeffObs[4] = -1.3959;
    coeffObs[5] = 2.4153;
    coeffObs[6] = 1.01669;
    coeffObs[7] = -2.47334;
    coeffObs[8] = 2.16305;
    coeffObs[9] = -0.08817;
    coeffObs[10] = 1.61984;
    coeffObs[11] = -0.65635;

    /*
     * Initialize the Pauli observable
     */
    pauliObservable.components = compObs;
    pauliObservable.coefficients = coeffObs;
    pauliObservable.qubits = qubits;
    pauliObservable.length = lengthObs;

    observable.type = PAULI;
    observable.pauliObs = &pauliObservable;

    /*
     * Define the evolution operators
     */
    lengthEvoOps = dim;
    compEvoOps = (pauli_t*) malloc(qubits * lengthEvoOps * sizeof(pauli_t));
    compEvoOps[0] = ID;
    compEvoOps[1] = ID;
    compEvoOps[2] = ID;
    compEvoOps[3] = ID;
    compEvoOps[4] = ID;
    compEvoOps[5] = ID;
    compEvoOps[6] = ID;
    compEvoOps[7] = Z;
    compEvoOps[8] = ID;
    compEvoOps[9] = ID;
    compEvoOps[10] = Z;
    compEvoOps[11] = ID;
    compEvoOps[12] = ID;
    compEvoOps[13] = ID;
    compEvoOps[14] = Z;
    compEvoOps[15] = Z;
    compEvoOps[16] = ID;
    compEvoOps[17] = Z;
    compEvoOps[18] = ID;
    compEvoOps[19] = ID;
    compEvoOps[20] = ID;
    compEvoOps[21] = Z;
    compEvoOps[22] = ID;
    compEvoOps[23] = Z;
    compEvoOps[24] = ID;
    compEvoOps[25] = Z;
    compEvoOps[26] = Z;
    compEvoOps[27] = ID;
    compEvoOps[28] = ID;
    compEvoOps[29] = Z;
    compEvoOps[30] = Z;
    compEvoOps[31] = Z;
    compEvoOps[32] = Z;
    compEvoOps[33] = ID;
    compEvoOps[34] = ID;
    compEvoOps[35] = ID;
    compEvoOps[36] = Z;
    compEvoOps[37] = ID;
    compEvoOps[38] = ID;
    compEvoOps[39] = Z;
    compEvoOps[40] = Z;
    compEvoOps[41] = ID;
    compEvoOps[42] = Z;
    compEvoOps[43] = ID;
    compEvoOps[44] = Z;
    compEvoOps[45] = ID;
    compEvoOps[46] = Z;
    compEvoOps[47] = Z;
    compEvoOps[48] = Z;
    compEvoOps[49] = Z;
    compEvoOps[50] = ID;
    compEvoOps[51] = ID;
    compEvoOps[52] = Z;
    compEvoOps[53] = Z;
    compEvoOps[54] = ID;
    compEvoOps[55] = Z;
    compEvoOps[56] = Z;
    compEvoOps[57] = Z;
    compEvoOps[58] = Z;
    compEvoOps[59] = ID;
    compEvoOps[60] = Z;
    compEvoOps[61] = Z;
    compEvoOps[62] = Z;
    compEvoOps[63] = Z;
    coeffEvoOps = (double*) malloc(circdepth * lengthEvoOps * sizeof(double));
    coeffEvoOps[0] = 0.22325;
    coeffEvoOps[1] = -1.94488;
    coeffEvoOps[2] = 1.42335;
    coeffEvoOps[3] = -2.07013;
    coeffEvoOps[4] = -1.10257;
    coeffEvoOps[5] = 1.24892;
    coeffEvoOps[6] = 0.85684;
    coeffEvoOps[7] = 1.20935;
    coeffEvoOps[8] = -2.11096;
    coeffEvoOps[9] = -1.0522;
    coeffEvoOps[10] = -0.54858;
    coeffEvoOps[11] = -0.8375;
    coeffEvoOps[12] = -0.70558;
    coeffEvoOps[13] = -1.46657;
    coeffEvoOps[14] = 2.40218;
    coeffEvoOps[15] = 0.88001;
    coeffEvoOps[16] = 0.66932;
    coeffEvoOps[17] = 0.70457;
    coeffEvoOps[18] = 2.84133;
    coeffEvoOps[19] = 2.48495;
    coeffEvoOps[20] = 2.14827;
    coeffEvoOps[21] = 1.32815;
    coeffEvoOps[22] = -1.44737;
    coeffEvoOps[23] = -2.3747;
    coeffEvoOps[24] = -1.47105;
    coeffEvoOps[25] = 1.17588;
    coeffEvoOps[26] = -1.19095;
    coeffEvoOps[27] = -0.93032;
    coeffEvoOps[28] = -0.82618;
    coeffEvoOps[29] = 0.33308;
    coeffEvoOps[30] = -1.29641;
    coeffEvoOps[31] = 1.53211;
    coeffEvoOps[32] = 0.92809;
    coeffEvoOps[33] = 1.50841;
    coeffEvoOps[34] = 1.65796;
    coeffEvoOps[35] = 0.69129;
    coeffEvoOps[36] = -2.65986;
    coeffEvoOps[37] = 0.66674;
    coeffEvoOps[38] = 1.33686;
    coeffEvoOps[39] = 2.25722;
    coeffEvoOps[40] = -0.8822;
    coeffEvoOps[41] = 0.24263;
    coeffEvoOps[42] = 0.48472;
    coeffEvoOps[43] = -1.98169;
    coeffEvoOps[44] = -2.56158;
    coeffEvoOps[45] = -2.0674;
    coeffEvoOps[46] = 1.76359;
    coeffEvoOps[47] = -2.89002;
    coeffEvoOps[48] = 2.09041;
    coeffEvoOps[49] = 1.42077;
    coeffEvoOps[50] = 2.12384;
    coeffEvoOps[51] = -0.42847;
    coeffEvoOps[52] = 2.85773;
    coeffEvoOps[53] = 1.68218;
    coeffEvoOps[54] = 2.54823;
    coeffEvoOps[55] = 2.2562;
    coeffEvoOps[56] = 1.03429;
    coeffEvoOps[57] = -1.8238;
    coeffEvoOps[58] = 1.21566;
    coeffEvoOps[59] = -2.86379;
    coeffEvoOps[60] = -2.20972;
    coeffEvoOps[61] = 2.27518;
    coeffEvoOps[62] = 2.49163;
    coeffEvoOps[63] = 0.32578;
    coeffEvoOps[64] = 1.67526;
    coeffEvoOps[65] = -1.18713;
    coeffEvoOps[66] = -0.12346;
    coeffEvoOps[67] = -0.71438;
    coeffEvoOps[68] = -1.5839;
    coeffEvoOps[69] = -1.10577;
    coeffEvoOps[70] = 2.56117;
    coeffEvoOps[71] = -2.56436;
    coeffEvoOps[72] = 0.42191;
    coeffEvoOps[73] = -2.45698;
    coeffEvoOps[74] = 1.33464;
    coeffEvoOps[75] = -0.67253;
    coeffEvoOps[76] = -2.66351;
    coeffEvoOps[77] = -0.86787;
    coeffEvoOps[78] = 0.807;
    coeffEvoOps[79] = 0.61927;

    /*
     * Initialize the evolution operators
     */
    evoOpsPauli = (pauliObs_t*) malloc(circdepth * sizeof(pauliObs_t));
    for (depth_t i = 0; i < circdepth; ++i) {
        evoOpsPauli[i].components = compEvoOps;
        evoOpsPauli[i].coefficients = coeffEvoOps + (i * lengthEvoOps);
        evoOpsPauli[i].qubits = qubits;
        evoOpsPauli[i].length = lengthEvoOps;
    }

    // evoOps = (obs_t**) malloc(circdepth * sizeof(obs_t));
    for (depth_t i = 0; i < circdepth; ++i) {
        evoOps[i] = (obs_t*) malloc(sizeof(obs_t));
        evoOps[i]->type = PAULI;
        evoOps[i]->pauliObs = evoOpsPauli + i;
    }

    /*
     * Calculate the gradients
     */
    for (dim_t i = 0; i < dim + 1; ++i) {
        /*
 * Gradient from gradient PQC
 */
        double* result;
        stateInitVector(&testState, testVectors[i], qubits);
        result = gradientPQC(&testState, parameters, &observable, (const obs_t **) evoOps, circdepth);

        /*
         * Gradient from finite difference method
         */
        double* reference = finiteGradientPQC(refVectors[i], qubits, dim, circdepth, compObs, coeffObs, lengthObs, \
                                              compEvoOps, coeffEvoOps, lengthEvoOps, parameters, 1e-9);

        TEST_ASSERT_TRUE(rvectorAlmostEqual(result, reference, circdepth, APPROXPRECISION));
    }

    /*
    * Free test vectors
    */
    freeTestVectors(testVectors, qubits);
    freeTestVectors(refVectors, qubits);

    /*
     * Free observable and evolution operators
     */
    free(compObs);
    free(coeffObs);
    free(compEvoOps);
    free(coeffEvoOps);
    free(evoOpsPauli);
    for (depth_t i = 0; i < circdepth; ++i) {
        free(evoOps[i]);
    }
    free(evoOps);
}

/*
 * =====================================================================================================================
 *                                                      test hessianPQC
 * =====================================================================================================================
 */

void testHessianPQC(void) {
    qubit_t qubits;
    dim_t dim;
    depth_t circdepth;
    double *parameters;

    state_t testState;
    cplx_t **testVectors;
    cplx_t **refVectors;

    complength_t lengthObs;
    pauli_t *compObs;
    double *coeffObs;
    pauliObs_t pauliObservable;
    obs_t observable;

    complength_t lengthEvoOps;
    pauli_t *compEvoOps;
    double *coeffEvoOps;
    pauliObs_t *evoOpsPauli;
    obs_t **evoOps;

    /*
                                                3 qubit, circuit depth = 2
    --------------------------------------------------------------------------------------------------------------------
    */
    qubits = 3;
    dim = POW2(qubits, dim_t);

    /*
     * Define the PQC setting
     */
    circdepth = 2;
    parameters = (double *) malloc(circdepth * sizeof(double));
    parameters[0] = 1.39628;
    parameters[1] = -2.83693;

    /*
     * Initialize the state vectors
     */
    testVectors = generateTestVectors(qubits);
    refVectors = generateTestVectors(qubits);

    /*
     * Define the observable
     */
    lengthObs = 2 * dim;
    compObs = (pauli_t *) malloc(qubits * lengthObs * sizeof(pauli_t));
    compObs[0] = ID;
    compObs[1] = ID;
    compObs[2] = ID;
    compObs[3] = X;
    compObs[4] = ID;
    compObs[5] = Y;
    compObs[6] = ID;
    compObs[7] = Z;
    compObs[8] = X;
    compObs[9] = ID;
    compObs[10] = X;
    compObs[11] = X;
    compObs[12] = X;
    compObs[13] = Y;
    compObs[14] = X;
    compObs[15] = Z;
    compObs[16] = Y;
    compObs[17] = ID;
    compObs[18] = Y;
    compObs[19] = X;
    compObs[20] = Y;
    compObs[21] = Y;
    compObs[22] = Y;
    compObs[23] = Z;
    compObs[24] = Z;
    compObs[25] = ID;
    compObs[26] = Z;
    compObs[27] = X;
    compObs[28] = Z;
    compObs[29] = Y;
    compObs[30] = Z;
    compObs[31] = Z;
    compObs[32] = ID;
    compObs[33] = ID;
    compObs[34] = Y;
    compObs[35] = X;
    compObs[36] = ID;
    compObs[37] = Y;
    compObs[38] = Y;
    compObs[39] = ID;
    compObs[40] = Z;
    compObs[41] = ID;
    compObs[42] = Z;
    compObs[43] = X;
    compObs[44] = ID;
    compObs[45] = Y;
    compObs[46] = Z;
    compObs[47] = Z;
    coeffObs = (double *) malloc(lengthObs * sizeof(double));
    coeffObs[0] = 2.42261;
    coeffObs[1] = -0.17700;
    coeffObs[2] = -2.07679;
    coeffObs[3] = -0.45238;
    coeffObs[4] = -1.76552;
    coeffObs[5] = -0.61069;
    coeffObs[6] = -2.43846;
    coeffObs[7] = 2.24760;
    coeffObs[8] = -2.97751;
    coeffObs[9] = -1.10633;
    coeffObs[10] = -0.82939;
    coeffObs[11] = -1.56654;
    coeffObs[12] = -2.66825;
    coeffObs[13] = -0.64442;
    coeffObs[14] = 1.05226;
    coeffObs[15] = 1.38177;

    /*
     * Initialize the Pauli observable
     */
    pauliObservable.components = compObs;
    pauliObservable.coefficients = coeffObs;
    pauliObservable.qubits = qubits;
    pauliObservable.length = lengthObs;

    observable.type = PAULI;
    observable.pauliObs = &pauliObservable;

    /*
     * Define the evolution operators
     */
    lengthEvoOps = dim;
    compEvoOps = (pauli_t *) malloc(qubits * lengthEvoOps * sizeof(pauli_t));
    compEvoOps[0] = ID;
    compEvoOps[1] = ID;
    compEvoOps[2] = ID;
    compEvoOps[3] = ID;
    compEvoOps[4] = ID;
    compEvoOps[5] = Z;
    compEvoOps[6] = ID;
    compEvoOps[7] = Z;
    compEvoOps[8] = ID;
    compEvoOps[9] = ID;
    compEvoOps[10] = Z;
    compEvoOps[11] = Z;
    compEvoOps[12] = Z;
    compEvoOps[13] = ID;
    compEvoOps[14] = ID;
    compEvoOps[15] = Z;
    compEvoOps[16] = ID;
    compEvoOps[17] = Z;
    compEvoOps[18] = Z;
    compEvoOps[19] = Z;
    compEvoOps[20] = ID;
    compEvoOps[21] = Z;
    compEvoOps[22] = Z;
    compEvoOps[23] = Z;
    coeffEvoOps = (double *) malloc(circdepth * lengthEvoOps * sizeof(double));
    coeffEvoOps[0] = 2.50220;                                                // array is a concatenation of coefficients
    coeffEvoOps[1] = -0.81084;                                               // for all evolution operators
    coeffEvoOps[2] = -2.14617;
    coeffEvoOps[3] = 1.33180;
    coeffEvoOps[4] = 1.76367;
    coeffEvoOps[5] = -0.59290;
    coeffEvoOps[6] = -2.32258;
    coeffEvoOps[7] = 1.38857;
    coeffEvoOps[8] = 2.52907;
    coeffEvoOps[9] = -2.59997;
    coeffEvoOps[10] = -1.71381;
    coeffEvoOps[11] = 0.93723;
    coeffEvoOps[12] = -0.58177;
    coeffEvoOps[13] = -1.83702;
    coeffEvoOps[14] = -1.20720;
    coeffEvoOps[15] = -0.39430;

    /*
     * Initialize the evolution operators
     */
    evoOpsPauli = (pauliObs_t *) malloc(circdepth * sizeof(pauliObs_t));
    for (depth_t i = 0; i < circdepth; ++i) {
        evoOpsPauli[i].components = compEvoOps;
        evoOpsPauli[i].coefficients = coeffEvoOps + (i * lengthEvoOps);
        evoOpsPauli[i].qubits = qubits;
        evoOpsPauli[i].length = lengthEvoOps;
    }

    evoOps = (obs_t **) malloc(circdepth * sizeof(obs_t *));
    for (depth_t i = 0; i < circdepth; ++i) {
        evoOps[i] = (obs_t *) malloc(sizeof(obs_t));
        evoOps[i]->type = PAULI;
        evoOps[i]->pauliObs = evoOpsPauli + i;
    }


    /*
     * Test objective:
     * Compare the hessian obtained from hessianPQC with the approximation
     *      H(x).(epsilon * vec{e_{i}}) = (grad f)(x + epsilon * vec{e_{i}}) - (grad f)(x)
     */
    for (dim_t i = 0; i < dim + 1; ++i) {
        printf("Vector: ");
        cvectorPrint(testVectors[i], dim);

        /*
         * 1. Calculate the hessian via hessianPQC
         */
        stateInitVector(&testState, testVectors[i], qubits);
        double* result = hessianPQC(&testState, parameters, &observable, (const obs_t **) evoOps, circdepth);

        /*
         * Calculate the hessian via finite difference method
         */
        double* reference = finiteHessianPQC(refVectors[i],qubits, dim, circdepth, compObs, coeffObs, lengthObs, \
                                             compEvoOps, coeffEvoOps, lengthEvoOps, parameters, 1e-9);

        printf("Result:\n");
        matrixPrint(result, circdepth);
        printf("Reference:\n");
        matrixPrint(reference, circdepth);

        //TEST_ASSERT_TRUE(rvectorAlmostEqual(result, reference, circdepth, APPROXPRECISION));

        free(result);
        free(reference);
    }
    /*
    * Free test vectors
    */
    freeTestVectors(testVectors, qubits);
    freeTestVectors(refVectors, qubits);

    /*
     * Free observable and evolution operators
     */
    free(compObs);
    free(coeffObs);
    free(compEvoOps);
    free(coeffEvoOps);
    free(evoOpsPauli);
    for (depth_t i = 0; i < circdepth; ++i) {
        free(evoOps[i]);
    }
    free(evoOps);
}



/*
 * =================================================================================================
 *                                    test momentMatQAOA
 * =================================================================================================
 */

/*
void testMomentMatQAOA(void) {
    */
/*
     *                                        2 qubits
     *-----------------------------------------------------------------------------------------------------------------
     *//*

    qubit_t qubits = 2;

    */
/*
     * Initialize the state
     *//*

    state_t testState;
    cplx_t** testVectors = generateTestVectors(qubits);
    cplx_t** refVectors = generateTestVectors(qubits);

    */
/*
     * Initialize the Hamiltonian
    *//*

    pauli_t compHamiltonian[8] = {ID, ID, Z, ID, ID, Z, Z, Z};
    double coeffHamiltonian[4] = {1.46803, 2.52930, -1.42054, -2.18498};

    pauliObs_t hamiltonianPauli;
    hamiltonianPauli.components = compHamiltonian;
    hamiltonianPauli.coefficients = coeffHamiltonian;
    hamiltonianPauli.length = 4;
    hamiltonianPauli.qubits = qubits;

    obs_t hamiltonian;
    hamiltonian.type = PAULI;
    hamiltonian.pauliObs = &hamiltonianPauli;
}
*/

/*
 * =================================================================================================
 *                                              main
 * =================================================================================================
 */

int main(void) {
    UNITY_BEGIN();
    RUN_TEST(testExpValObs);
    RUN_TEST(testGradientPQC);
    RUN_TEST(testHessianPQC);
//    RUN_TEST(testMomentMatQAOA);
    return UNITY_END();
}
