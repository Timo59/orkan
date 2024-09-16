/*
 * =====================================================================================================================
 *                                                      includes
 * =====================================================================================================================
 */

#include "helpfunc.h"
#include "pauli.h"
#include "unity.h"


/*
 * =====================================================================================================================
 *                                                      macros
 * =====================================================================================================================
 */

#define MAXQUBITS               4
#define PRECISION               1e-8
#define PSTRINGSC               (1 << (2 * MAXQUBITS))

double coeffs[PSTRINGSC] = {-0.77440, 2.43107, 2.00394, -0.55128, 0.55065, -0.80536, 2.23597, -3.06746, -1.54823,
                            -1.60315, -0.37278, -1.51451, -2.09702, -0.93134, 1.14874, -2.09257, -2.11681, 1.56819,
                            -1.65715, 0.67240, -1.39775, -1.21508, 1.90700, -2.31899, -2.92445, -2.54491, 1.49788,
                            0.30807, 0.54191, -2.17137, -1.97164, -2.83741, 1.60705, 2.11855, 0.02400, -0.16501,
                            -2.57686, -3.09622, -3.14106, 0.78665, -2.04357, -2.54521, 0.91212, -0.51588, 1.41691,
                            -2.04886, 2.29273, 1.65618, -1.85889, 0.64541, 0.52785, -2.45732, 0.50906, 0.41829,
                            -1.42353, -3.07453, 0.17388, -1.35264, 2.67657, 0.61807, -0.81420, 2.30384, -2.55706,
                            0.79168, 0.22730, 0.60121, -0.98412, 0.28509, -1.89957, 1.56349, -0.73497, 2.86948, 0.91245,
                            -2.51270, 1.37985, -0.28207, 1.50461, 0.39119, -1.54438, -1.95145, 0.75922, 2.49509,
                            0.76597, 2.01656, 1.69181, -0.82800, 2.83133, 0.15439, 0.81941, 2.17834, 0.93565,
                            -1.25391, -2.89448, 1.75100, 2.63705, -2.96011, -1.69040, -0.33305, 1.70643, -1.84817,
                            -2.33319, -0.97798, -2.72351, 1.57497, -2.73733, 1.17078, -1.08713, 2.24957, -1.89495,
                            -2.75137, -0.24599, -1.70457, 2.61730, -2.92135, -1.69270, -1.79968, -2.58849, 1.33295,
                            2.94636, 0.68571, 0.89565, -1.01885, -2.42970, 2.88031, 0.87694, 0.99229, 1.92511, 1.33809,
                            -1.67410, -0.18586, 0.41559, -3.03914, -0.58092, 0.29331, -0.13534, -2.23067, -1.32516,
                            -1.92004, -1.49208, -2.76493, -1.84452, 0.47071, 2.72963, -1.11332, -1.30182, 0.91594,
                            3.00805, -0.34837, -0.78294, -1.84118, -2.99507, -3.13418, -2.60612, -2.51759, 2.91219,
                            2.13354, 1.81222, -0.58226, -1.05772, 1.25702, -1.16129, -2.41894, 2.65087, -1.15217,
                            0.24417, 2.59828, 0.23416, 1.61763, -0.15563, -2.80235, 2.68974, 2.50860, 1.37030, 1.44378,
                            1.60758, 0.80128, -0.31447, -1.63624, 0.74842, 1.50803, -1.89638, -0.52538, -1.87867,
                            1.60308, -2.66022, 0.80951, -1.35723, -1.54496, 2.54774, 0.08865, -1.96006, -0.40487,
                            1.36232, -1.85312, 0.32163, 3.10785, -0.56181, 3.12310, -1.28192, 1.88849, -2.01835,
                            -0.09475, 2.08552, 1.05407, -0.67002, 1.74697, 2.56223, 2.69804, 3.04177, 2.42893, -3.02526,
                            0.68351, 0.16499, -2.85656, -1.82865, 1.86240, 1.08765, 2.74518, 0.40135, 2.67996, -0.20404,
                            -2.96269, 0.20928, 3.00434, 0.32514, -1.22278, -2.41959, 2.08024, -0.68699, -1.42931,
                            -1.35932, -0.15992, -1.39280, -1.75033, -2.55720, 1.24729, -1.65536, -1.20058, 1.45968,
                            -2.00608, -1.61730, 1.87437, -0.45750, -2.00305, -0.82357, -2.19217, 1.85407, 0.46075,
                            1.85401, 1.15177, -2.24874, 2.16107, 2.81549, 1.55742, 3.00655, -1.08365};

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
 *                                              test ApplyObsevablePauli
 * =====================================================================================================================
 */

/*
 * This function checks whether using applyObservablePauli coincides with applying the Kronecker product of Pauli matrices
 * on instances ranging from 1 to 4 qubits.
 *
 * Input:
 *      -
 * Output:
 *      There is no return value. Only checks whether the tested state vector and the reference state
 *      vector are equal to some precision after an evolution and when reversing the evolution.
 */

void testApplyObsPauli(void) {
    state_t testState;          // State equipped with the test statevectors
    pauliObs_t obs;

    for (qubit_t qubits = 1; qubits <= MAXQUBITS; ++qubits) {

        dim_t dim = POW2(qubits, dim_t);                        // Hilbert space dimension
        cplx_t** testVectors = generateTestVectors(qubits);     // Statevectors altered by the test function
        cplx_t** refVectors = generateTestVectors(qubits);      // Statevectors altered by matrix multiplication
                                                                // (reference)

        pauli_t* comps = allPauliStrings(qubits);               // Concatenation of all Pauli strings
        obs.comps = comps;
        obs.coeffs = coeffs;
        complength_t length = (1 << (2 * qubits));
        obs.length = length;
        obs.qubits = qubits;

        /* Compute the observable's matrix representation */
        cplx_t* obsMatrix = pauliObsMat(comps, coeffs, length, qubits);

        for (dim_t i = 0; i < dim + 1; ++i) {
            stateInitVector(&testState, testVectors[i], qubits);    // Initialize test state
            applyObsPauliBlas_omp_new(&testState, &obs);                 // Apply test function

            cmatVecMulInPlace(obsMatrix, refVectors[i],dim);        // Multiply reference vector with matrix in
                                                                        // place
            TEST_ASSERT_TRUE(cvectorAlmostEqual(refVectors[i], testState.vec, dim, PRECISION));

            stateFreeVector(&testState);
        }

        freeTestVectors(refVectors, qubits);
        free(testVectors);
        free(obsMatrix);
        free(comps);
    }
}

/* 
 * =====================================================================================================================
 *                                              test evolveWithPauliString
 * =====================================================================================================================
 */

/*
 * This function checks whether applying evolveWithPauliString coincides with applying the Kronecker product of Pauli
 * matrices on instances ranging from 1 to 4 qubits.
 *
 * Input:
 *      -
 * Output:
 *      There is no return value. Only checks whether the tested state vector and the reference state vector are equal
 *      to some precision after an evolution and when reversing the evolution.
 */

void testEvolvePauliStr(void) {
    state_t testState;          // State equipped with the test statevectors

    for (qubit_t qubits = 1; qubits <= MAXQUBITS; ++qubits) {

        dim_t dim = POW2(qubits, dim_t);                        // Hilbert space dimension
        dim_t stringc = (1 << (2 * qubits));                    // Number of Pauli strings for qubits
        pauli_t* strings = allPauliStrings(qubits);             // Concatenation of all Pauli strings

        for (dim_t i = 0; i < stringc; ++i) {
            cplx_t** testVectors = generateTestVectors(qubits);     // Statevectors altered by the test function
            cplx_t** refVectors = generateTestVectors(qubits);      // Statevectors altered by matrix multiplication
                                                                    // (reference)
            /* Compute the observable's matrix representation */
            cplx_t* evoMatrix = expPauliStrMat(strings + i, coeffs[i], qubits);

            for (dim_t j = 0; j < dim + 1; ++j) {
                stateInitVector(&testState, testVectors[j], qubits);        // Initialize test state
                evolvePauliStrOmp(&testState, strings + i, coeffs[i]);  // Apply test function

                /* Multiply evoMatrix and reference statevector */
                cplx_t* reference = cmatVecMul(evoMatrix, refVectors[j], dim);

                TEST_ASSERT_TRUE(cvectorAlmostEqual(reference, testState.vec, dim, PRECISION));

                evolvePauliStrOmp(&testState, strings + i, -coeffs[i]); // Reverse test function

                TEST_ASSERT_TRUE(cvectorAlmostEqual(refVectors[j], testState.vec, dim, PRECISION));

                stateFreeVector(&testState);
            }
            freeTestVectors(refVectors, qubits);
            free(testVectors);

            free(evoMatrix);
            evoMatrix = NULL;
        }

        free(strings);
        strings = NULL;
    }
}

/*
 * =====================================================================================================================
 *                                 test evolveWithTrotterizedObservablePauli
 * =====================================================================================================================
 */

/*
 * This function checks whether applying evolveWithTrotterizedObservablePauli coincides with applying the Kronecker
 * product of Pauli matrices on instances ranging from 1 to MAXQUBITS qubits.
 *
 * Input:
 *      -
 * Output:
 *      There is no return value. Only checks whether the tested state vector and the reference state vector are equal
 *      to some precision after an evolution and when reversing the evolution.
 */

void testEvolveObsPauliTrotter(void) {
    state_t testState;          // State equipped with the test statevectors

    for (qubit_t qubits = 1; qubits <= MAXQUBITS; ++qubits) {

        dim_t dim = POW2(qubits, dim_t);                        // Hilbert space dimension
        cplx_t** testVectors = generateTestVectors(qubits);     // Statevectors altered by the test function
        cplx_t** refVectors = generateTestVectors(qubits);      // Statevectors altered by matrix multiplication
                                                                // (reference)
        pauliObs_t obs;                                         // Declare the observable
        pauli_t* comps = allPauliStrings(qubits);               // Concatenation of all Pauli strings
        obs.comps = comps;
        obs.coeffs = coeffs;
        complength_t length = (1 << (2 * qubits));
        obs.length = length;
        obs.qubits = qubits;

        /* Compute the observable's matrix representation */
        cplx_t* evoMatrix = expPauliObsMatTrotter(comps, coeffs, length, angle[0], qubits);

        for (dim_t i = 0; i < dim + 1; ++i) {
            stateInitVector(&testState, testVectors[i], qubits);    // Initialize test state
            /* Apply test function */
            evolveObsPauliTrotter(&testState, &obs, angle[0]);   // Apply test function

            /* Multiply evoMatrix and reference statevector */
            cplx_t* reference = cmatVecMul(evoMatrix, refVectors[i], dim);

            TEST_ASSERT_TRUE(cvectorAlmostEqual(reference, testState.vec, dim, PRECISION));

            /* Reverse test function */
            evolveObsPauliTrotterDagger(&testState, &obs, angle[0]);

            TEST_ASSERT_TRUE(cvectorAlmostEqual(refVectors[i], testState.vec, dim, PRECISION));

            stateFreeVector(&testState);
        }

        freeTestVectors(refVectors, qubits);
        refVectors = NULL;
        free(testVectors);

        free(comps);
        comps = NULL;
        free(evoMatrix);
    }
}


/*
 * =====================================================================================================================
 *                                                      main
 * =====================================================================================================================
 */

//int main(void) {
//
//    UNITY_BEGIN();
//    RUN_TEST(testApplyObservablePauli);
//    RUN_TEST(testEvolveWithPauliString);
//    RUN_TEST(testEvolveWithTrotterizedObservablePauli);
//    return UNITY_END();
//}
