/*
 * =================================================================================================
 *                                              includes
 * =================================================================================================
 */

#include "helpfunc.h"
#include "pauli.h"
#include "cplxutil.h"
#include "linalg.h"
#include "unity.h"


/*
 * =================================================================================================
 *                                              macros
 * =================================================================================================
 */

#define PRECISION               1e-8
#define MAXQUBITS               4

/*
 * =================================================================================================
 *                                              setUp and tearDown
 * =================================================================================================
 */

void setUp(void) {}

void tearDown(void) {}

/* 
 * =================================================================================================
 *                                     test evolveWithPauliString
 * =================================================================================================
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
    RUN_TEST(testEvolveWithPauliString);
    RUN_TEST(testEvolveWithTrotterizedObservablePauli);
    return UNITY_END();
}
