/*
 * =====================================================================================================================
 *                                                      includes
 * =====================================================================================================================
 */

#include "helpfunc.h"
#include "pauli.h"
#include "cplxutil.h"
#include "linalg.h"
#include "unity.h"


/*
 * =====================================================================================================================
 *                                                      macros
 * =====================================================================================================================
 */

#define PRECISION               1e-8
#define MAXQUBITS               4
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

/*
 * =====================================================================================================================
 *                                                  setUp and tearDown
 * =====================================================================================================================
 */

void setUp(void) {}

void tearDown(void) {}

/*
 * =====================================================================================================================
 *                                              test applyOperatorPauli
 * =====================================================================================================================
 */

void testApplyObservablePauli(void) {
    state_t testState;          // State equipped with the test statevectors

    for (qubit_t qubits = 1; qubits <= MAXQUBITS; ++qubits) {

        dim_t dim = POW2(qubits, dim_t);                        // Hilbert space dimension
        cplx_t** testVectors = generateTestVectors(qubits);     // Statevectors altered by the test function
        cplx_t** refVectors = generateTestVectors(qubits);      // Statevectors altered by matrix multiplication
                                                                // (reference)
        pauliObs_t obs;                                     // Initialize the observable
        pauli_t* comps = allPauliStrings(qubits);           // Concatenation of all Pauli strings
        obs.components = comps;
        obs.coefficients = coeffs;
        complength_t length = (1 << (2 * qubits));
        obs.length = length;
        obs.qubits = qubits;

        for (dim_t i = 0; i < dim + 1; ++i) {
            stateInitVector(&testState, testVectors[i], qubits);    // Initialize test state
            applyObservablePauliOmp(&testState, &obs);                 // Apply test function

            /* Compute the observable's matrix representation */
            cplx_t* obsMatrix = pauliObservableMat(comps, coeffs, length, qubits);
            cmatVecMulInPlace(obsMatrix, refVectors[i],dim);        // Multiply reference vector with matrix in
                                                                    // place

            TEST_ASSERT_TRUE(cvectorAlmostEqual(refVectors[i], testState.vector, dim, PRECISION));
        }
    }
}

/* 
 * =====================================================================================================================
 *                                              test evolveWithPauliString
 * =====================================================================================================================
 */

/*
This function checks whether the evolution induced by evolveWithComponent coincides with the one
induced by the Kronecker product of Pauli matrices.

Input:
    -

Output:
    There is no return value. Only checks whether the tested state vector and the reference state
    vector are equal to some precision after an evolution and when reversing the evolution.
*/
void testEvolveWithPauliString(void) {
    dim_t dim;                  // unsigend integer holding the dimension

    cplx_t** testVectors;       // double complex pointer to an array of pointers defining the set
                                // of reference vectors

    cplx_t** stateVectors;      // double complex pointer to an array of pointers defining the set
                                // of test vectors

    state_t testState;          // custom struct holding the double complex test state vector, the
                                // number of qubits and the Hilbert's space dimension

    cplx_t* evoMatrix;          // double complex pointer to the reference matrix emerging from 
                                // Kronecker products within the Pauli string

    cplx_t* result;             // double complex pointer holding the result of the reference vector
                                // from matrix vector multiplication

    double angle;               // floating variable with double precision holding the angle of the
                                // evolution

    qubit_t index;              // signed integer holding 

    /*
    Starting with one, instances of Pauli strings are created up to the number MAXQUBITS defined in
    the macros.
    */
    for (qubit_t qubits = 1; qubits <= MAXQUBITS; ++qubits) {
        dim = POW2(qubits, dim_t);                                          // dimension is two to
                                                                            // the power of qubits

        testVectors = generateTestVectors(qubits);      // the reference vectors are initialised as
                                                        // all computational basis vectors plus the
                                                        //extra vector

        stateVectors = generateTestVectors(qubits);     // the test vectors are initialised as all
                                                        // computational basis vectors plus the
                                                        // extra vector

        pauli_t* paulistr;                              // define a Pauli string
        paulistr = (pauli_t*) calloc(qubits, sizeof(pauli_t));

        angle = qubits * 0.25;                          // fix the angle to a quarter of the number
                                                        // of qubits
        /*
        Starting at the all identity Pauli string gradually change each position through all Paulis
        until all possible configurations are considered.
        */                                                               
        do {
            evoMatrix = expPauliStringMat(paulistr, angle, qubits); // compute the evolution matrix
                                                                    // of the current Pauli string
                                                                    // and store it in evoMatrix

            /*
            Starting at the first test and reference vector, respectively, apply the multiplication
            with evoMatrix as well as the evolution method from pauli.c and compare the results.
            Then reverse the evolution and compare the resulting test vector with the initial 
            reference vector.
            */
            for (dim_t i = 0; i < dim + 1; ++i) {
                stateInitVector(&testState, stateVectors[i], qubits);       // initialise a pure
                                                                            // state with the
                                                                            // current test vector 
                                                                            // as its state vector

                result = cmatVecMul(evoMatrix, testVectors[i], dim);        // multiply the current
                                                                            // reference vector with
                                                                            // evoMatrix ad store the
                                                                            // result in result

                evolveWithPauliString(&testState, paulistr, angle);         // evolve the state with
                                                                            // the Pauli string by
                                                                            // the specified angle

                TEST_ASSERT_TRUE(cvectorAlmostEqual(result, testState.vector, dim, PRECISION));
                                                                            // check whether the
                                                                            // result and the state
                                                                            // vector coincide up to
                                                                            // a certain threshold
                                                                            // If not return a
                                                                            // warning

                free(result);                                               // free the memory 
                                                                            // allocated to result

                evolveWithPauliString(&testState, paulistr, -angle);         // reverse the evolution
                                                                            // of the state with the
                                                                            // Pauli string

                TEST_ASSERT_TRUE(cvectorAlmostEqual(testVectors[i], testState.vector, dim, \
                                                    PRECISION));            // check whether the
                                                                            // initial reference
                                                                            // vector and the state
                                                                            // coincide up to a
                                                                            // certain threshold
                                                                            // If not return a
                                                                            // warning   
            }

            /*
            Starting at the first position increase the integer assigned to pauli at that position.
            As soon as it reaches four that position is reset to zero and the next index is altered.
            */
            for (index = 0; index < qubits; ++index) {
                ++paulistr[index];              // increase the Pauli at the position index by one
                                                // which coincides with going through all Paulis

                /*
                If the Pauli at position index is smaller than 4, i.e., corresponds to a Pauli,
                break the loop and check the evolution of this Pauli string.
                */
                if (paulistr[index] < 4) {
                    break;
                }
                /*
                As soon as the Pauli at position index reaches 4, reset it to zero and increase the
                index by one. The next index is also increased by one then but one breaks the loop
                when it is still not four. In the next round, the first position is again altered
                until it reaches four and so on.
                */
                paulistr[index] = 0;
            }

            /*
            When all possible combinations are computed, so the index reaches the number of qubits,
            the do while loop is left.
            */
            if (index == qubits) {
                break;
            }

        } while (1);
        freeTestVectors(testVectors, qubits);       // the memory allocated to the reference vectors
                                                    // is freed

        freeTestVectors(stateVectors, qubits);      // the memory allocated to the test vectors is
                                                    // freed
    }
}

/*
 * =================================================================================================
 *                                 test evolveWithTrotterizedObservablePauli
 * =================================================================================================
 */

/*
This function checks whether the evolution of a state by a sum of Pauli strings induced by the 
evolveWithTrotterizedObservablePauli coincides with the evolution from matrix multiplication of
Kronecker products of these Pauli strings.

Input:
    -

Output:
    There is no return value. Only checks whether the tested state vector and the reference state
    vector are equal to some precision after an evolution and when reversing the evolution.
*/
void testEvolveWithTrotterizedObservablePauli(void) {
    qubit_t qubits;  

    dim_t dim;

    cplx_t** testVectors;       // double complex pointer to an array of pointers defining the set
                                // of reference vectors

    cplx_t** stateVectors;      // double complex pointer to an array of pointers defining the set
                                // of test vectors

    state_t testState;          // custom struct holding the double complex test state vector, the
                                // number of qubits and the Hilbert's space dimension

    pauliObs_t observable;      // custom struct holding the components, the coefficients and their
                                // number as well as the number of qubits

    cplx_t* evoMatrix;          // double complex pointer to the reference matrix emerging from 
                                // Kronecker products within the Pauli string

    cplx_t* result;             // double complex pointer holding the result of the reference vector
                                // from matrix vector multiplication

    double angle;               // floating variable with double precision holding the angle of the
                                // evolution

    /* 
                                        1 Qubit
    ------------------------------------------------------------------------------------------------
    */  

    qubits = 1;
    dim = POW2(qubits, dim_t);

    testVectors = generateTestVectors(qubits);  // the reference vectors are initialised as
                                                // all computational basis vectors plus the extra
                                                // vector

    stateVectors = generateTestVectors(qubits); // the test vectors are initialised as all
                                                // computational basis vectors plus the extra
                                                // vector

    pauli_t components1[2] = {ID, Z};           // declare and initialise an array of Pauli strings

    double coefficients1[2] = {101./53., 17./14.};
                                                // declare and initialise the coefficients in an
                                                // array
    observable.components = components1;
    observable.coefficients = coefficients1;
    observable.length = 2;
    observable.qubits = qubits;

    angle = qubits * 0.25;                      // initialise the angle of the evolution

    evoMatrix = expTrotterizedPauliObservableMat(components1, coefficients1, 2, angle, qubits);

    for (dim_t i = 0; i < dim + 1; ++i) {
        stateInitVector(&testState, stateVectors[i], qubits);       // initialise a pure
                                                                    // state with the
                                                                    // current test vector 
                                                                    // as its state vector

        result = cmatVecMul(evoMatrix, testVectors[i], dim);        // multiply the current
                                                                    // reference vector with
                                                                    // evoMatrix ad store the
                                                                    // result in result

        evolveWithTrotterizedObservablePauli(&testState, &observable, angle);
                                                                    // evolve the state with
                                                                    // the Pauli string by
                                                                    // the specified angle

        TEST_ASSERT_TRUE(cvectorAlmostEqual(result, testState.vector, dim, PRECISION));
                                                                    // check whether the
                                                                    // result and the state
                                                                    // vector coincide up to
                                                                    // a certain threshold
                                                                    // If not return a
                                                                    // warning

        free(result);                                               // free the memory 
                                                                    // allocated to result

        evolveWithTrotterizedObservablePauli(&testState, &observable, -angle);        
                                                                    // reverse the evolution
                                                                    // of the state with the
                                                                    // Pauli string

        TEST_ASSERT_TRUE(cvectorAlmostEqual(testVectors[i], testState.vector, dim, \
                                            PRECISION));            // check whether the
                                                                    // initial reference
                                                                    // vector and the state
                                                                    // coincide up to a
                                                                    // certain threshold
                                                                    // If not return a
                                                                    // warning 
    }

    /* 
                                        2 Qubit
    ------------------------------------------------------------------------------------------------
    */  

    qubits = 2;
    dim = POW2(qubits, dim_t);

    testVectors = generateTestVectors(qubits);  // the reference vectors are initialised as
                                                // all computational basis vectors plus the extra
                                                // vector

    stateVectors = generateTestVectors(qubits); // the test vectors are initialised as all
                                                // computational basis vectors plus the extra
                                                // vector

    pauli_t components2[8] = {ID, ID, ID, Z, Z, ID, Z, Z}; 
                                                // declare and initialise an array of Pauli strings

    double coefficients2[4] = {.6532, -1.8123, 2.1386, -3.8998};
                                                // declare and initialise the coefficients in an
                                                // array
    observable.components = components2;
    observable.coefficients = coefficients2;
    observable.length = 4;
    observable.qubits = qubits;

    angle = qubits * 0.25;                      // initialise the angle of the evolution

    evoMatrix = expTrotterizedPauliObservableMat(components2, coefficients2, 4, angle, qubits);

    for (dim_t i = 0; i < dim + 1; ++i) {
        stateInitVector(&testState, stateVectors[i], qubits);       // initialise a pure
                                                                    // state with the
                                                                    // current test vector 
                                                                    // as its state vector
        // printf("Initial vector:\n");
        // cvectorPrint(testState.vector, dim);
        // printf("\n");

        result = cmatVecMul(evoMatrix, testVectors[i], dim);        // multiply the current
                                                                    // reference vector with
                                                                    // evoMatrix ad store the
                                                                    // result in result

        evolveWithTrotterizedObservablePauli(&testState, &observable, angle);
                                                                    // evolve the state with
                                                                    // the Pauli string by
                                                                    // the specified angle


        // printf("State vector:\n");
        // cvectorPrint(testState.vector, dim);
        // printf("\n");
        // printf("Reference vector:\n");
        // cvectorPrint(result, dim);

        TEST_ASSERT_TRUE(cvectorAlmostEqual(result, testState.vector, dim, PRECISION));
                                                                    // check whether the
                                                                    // result and the state
                                                                    // vector coincide up to
                                                                    // a certain threshold
                                                                    // If not return a
                                                                    // warning

        free(result);                                               // free the memory 
                                                                    // allocated to result

        evolveWithTrotterizedObservablePauli(&testState, &observable, -angle);        
                                                                    // reverse the evolution
                                                                    // of the state with the
                                                                    // Pauli string

        TEST_ASSERT_TRUE(cvectorAlmostEqual(testVectors[i], testState.vector, dim, \
                                            PRECISION));            // check whether the
                                                                    // initial reference
                                                                    // vector and the state
                                                                    // coincide up to a
                                                                    // certain threshold
                                                                    // If not return a
                                                                    // warning   
    }

    /* 
                                        3 Qubit
    ------------------------------------------------------------------------------------------------
    */  

    qubits = 3;
    dim = POW2(qubits, dim_t);

    testVectors = generateTestVectors(qubits);  // the reference vectors are initialised as
                                                // all computational basis vectors plus the extra
                                                // vector

    stateVectors = generateTestVectors(qubits); // the test vectors are initialised as all
                                                // computational basis vectors plus the extra
                                                // vector

    pauli_t components3[24] = {ID, ID, ID, X, ID, ID, ID, X, ID, ID, ID, X, X, X, ID, X, ID, X, \
                              ID, X, X, X, X, X}; 
                                                // declare and initialise an array of Pauli strings

    double coefficients3[8] = {1.4593, -0.8641, -0.6432, -3.0139, -2.1424, 1.2383, 0.1050, 2.3081};
                                                // declare and initialise the coefficients in an
                                                // array
    observable.components = components3;
    observable.coefficients = coefficients3;
    observable.length = 8;
    observable.qubits = qubits;

    angle = qubits * 0.25;                      // initialise the angle of the evolution

    evoMatrix = expTrotterizedPauliObservableMat(components3, coefficients3, 8, angle, qubits);

    for (dim_t i = 0; i < dim + 1; ++i) {
        stateInitVector(&testState, stateVectors[i], qubits);       // initialise a pure
                                                                    // state with the
                                                                    // current test vector 
                                                                    // as its state vector
        // printf("Initial vector:\n");
        // cvectorPrint(testState.vector, dim);
        // printf("\n");

        result = cmatVecMul(evoMatrix, testVectors[i], dim);        // multiply the current
                                                                    // reference vector with
                                                                    // evoMatrix ad store the
                                                                    // result in result

        evolveWithTrotterizedObservablePauli(&testState, &observable, angle);
                                                                    // evolve the state with
                                                                    // the Pauli string by
                                                                    // the specified angle


        // printf("State vector:\n");
        // cvectorPrint(testState.vector, dim);
        // printf("\n");
        // printf("Reference vector:\n");
        // cvectorPrint(result, dim);

        TEST_ASSERT_TRUE(cvectorAlmostEqual(result, testState.vector, dim, PRECISION));
                                                                    // check whether the
                                                                    // result and the state
                                                                    // vector coincide up to
                                                                    // a certain threshold
                                                                    // If not return a
                                                                    // warning

        free(result);                                               // free the memory 
                                                                    // allocated to result

        evolveWithTrotterizedObservablePauli(&testState, &observable, -angle);        
                                                                    // reverse the evolution
                                                                    // of the state with the
                                                                    // Pauli string

        TEST_ASSERT_TRUE(cvectorAlmostEqual(testVectors[i], testState.vector, dim, \
                                            PRECISION));            // check whether the
                                                                    // initial reference
                                                                    // vector and the state
                                                                    // coincide up to a
                                                                    // certain threshold
                                                                    // If not return a
                                                                    // warning   
    }

        /* 
                                        4 Qubit
    ------------------------------------------------------------------------------------------------
    */  

    qubits = 4;
    dim = POW2(qubits, dim_t);

    testVectors = generateTestVectors(qubits);  // the reference vectors are initialised as
                                                // all computational basis vectors plus the extra
                                                // vector

    stateVectors = generateTestVectors(qubits); // the test vectors are initialised as all
                                                // computational basis vectors plus the extra
                                                // vector

    pauli_t components4[32] = {ID, ID, ID, ID, ID, ID, ID, Y, ID, Y, ID, ID, Y, Y, ID, ID, \
                               ID, Y, ID, Y, Y, ID, Y, Y, Y, Y, Y, ID, Y, Y, Y, Y}; 
                                                // declare and initialise an array of Pauli strings

    double coefficients4[8] = {-1.4593, 0.8641, 0.6432, 3.0139, 2.1424, -1.2383, -0.1050, -2.3081};
                                                // declare and initialise the coefficients in an
                                                // array
    observable.components = components4;
    observable.coefficients = coefficients4;
    observable.length = 8;
    observable.qubits = qubits;

    angle = qubits * 0.25;                      // initialise the angle of the evolution

    evoMatrix = expTrotterizedPauliObservableMat(components4, coefficients4, 8, angle, qubits);

    for (dim_t i = 0; i < dim + 1; ++i) {
        stateInitVector(&testState, stateVectors[i], qubits);       // initialise a pure
                                                                    // state with the
                                                                    // current test vector 
                                                                    // as its state vector
        // printf("Initial vector:\n");
        // cvectorPrint(testState.vector, dim);
        // printf("\n");

        result = cmatVecMul(evoMatrix, testVectors[i], dim);        // multiply the current
                                                                    // reference vector with
                                                                    // evoMatrix ad store the
                                                                    // result in result

        evolveWithTrotterizedObservablePauli(&testState, &observable, angle);
                                                                    // evolve the state with
                                                                    // the Pauli string by
                                                                    // the specified angle


        // printf("State vector:\n");
        // cvectorPrint(testState.vector, dim);
        // printf("\n");
        // printf("Reference vector:\n");
        // cvectorPrint(result, dim);

        TEST_ASSERT_TRUE(cvectorAlmostEqual(result, testState.vector, dim, PRECISION));
                                                                    // check whether the
                                                                    // result and the state
                                                                    // vector coincide up to
                                                                    // a certain threshold
                                                                    // If not return a
                                                                    // warning

        free(result);                                               // free the memory 
                                                                    // allocated to result

        evolveWithTrotterizedObservablePauli(&testState, &observable, -angle);        
                                                                    // reverse the evolution
                                                                    // of the state with the
                                                                    // Pauli string

        TEST_ASSERT_TRUE(cvectorAlmostEqual(testVectors[i], testState.vector, dim, \
                                            PRECISION));            // check whether the
                                                                    // initial reference
                                                                    // vector and the state
                                                                    // coincide up to a
                                                                    // certain threshold
                                                                    // If not return a
                                                                    // warning   
    }
}

/*
 * =================================================================================================
 *                                              main
 * =================================================================================================
 */

int main(void) {

    UNITY_BEGIN();
    RUN_TEST(testApplyObservablePauli);
    RUN_TEST(testEvolveWithPauliString);
    RUN_TEST(testEvolveWithTrotterizedObservablePauli);
    return UNITY_END();
}
