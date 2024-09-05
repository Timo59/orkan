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

/*
 * =====================================================================================================================
 *                                                      macros
 * =====================================================================================================================
 */

#define APPROXPRECISION         1e-5
#define MAXQUBITS               4
#define PRECISION               1e-8
#define PSTRINGSC               (1 << (2 * MAXQUBITS))

double coeffs[PSTRINGSC] = {-0.91579, -0.25652, -1.51508, 1.05323, -1.23536, 1.49867, -1.53926, -0.49826, 1.22999,
                            -1.92813, 2.85517, -2.47542, -0.38280, -2.09975, 1.00740, 2.35575, -2.45489, -0.30425,
                            1.06238, -2.77212, 2.41397, -1.34010, -0.22039, -2.24888, 2.71889, -1.92035, -1.50087,
                            2.09585, -2.08985, 0.68230, 2.44196, -2.82084, 2.48521, -1.22644, -1.00846, 1.69127,
                            0.06418, -2.35933, -0.19995, 0.31603, -1.73091, 1.14717, 2.33569, 0.57928, -0.25737,
                            -0.45009, 2.78098, 2.95678, 1.91381, 2.54352, -2.15822, 1.72263, 2.05691, 2.39287, 1.58683,
                            -2.90241, 2.59916, 2.73879, 2.95639, 0.56235, 1.44027, -0.81242, 0.38418, 1.18629, -1.17247,
                            -0.12436, 1.62449, -2.43752, 1.66196, -0.47336, 1.84992, 1.95881, -1.00553, 2.13265,
                            0.56086, -1.56981, 1.04719, 2.60808, -1.84578, 1.20584, -2.35521, -1.32020, 2.01014,
                            1.26375, 1.27048, 2.59090, 1.23691, -1.80077, 1.93343, -2.91482, -0.53376, 3.03015, 1.48057,
                            2.76056, -0.43952, 0.21942, 2.04924, -1.52948, 0.75794, 1.35086, -0.38534, 2.40243,
                            -2.48500, -2.20105, -2.82368, -2.47896, -1.40897, -2.55432, 0.71394, 2.17555, 1.98635,
                            0.40811, 0.67655, -0.49530, -2.67357, -1.41395, -1.15755, -0.06063, -2.72028, -1.61941,
                            -2.16315, -1.75902, 0.20458, 0.80771, 0.57488, 1.72225, -1.41010, -0.81181, -2.49148,
                            1.94778, 2.34203, -1.32821, 0.96816, -1.26649, -2.30224, -2.87518, 2.88520, 2.96937,
                            -2.76896, 2.62508, -2.13933, -2.31723, -1.37332, -0.77016, -2.62147, -2.50779, -0.60819,
                            1.91038, 1.60209, 1.67136, 1.46729, -1.07858, -2.11997, -0.13966, 2.05480, 0.97619,
                            -0.11150, 0.30797, 0.83008, -0.52598, -2.93730, 0.21551, -2.75580, 0.49142, 1.47556,
                            2.71840, -2.22514, 1.00911, -2.89986, 1.31044, -0.51764, -2.79606, 1.47126, 1.95153,
                            -0.63209, -1.85194, 2.84156, 2.90295, -1.84033, -0.04310, 1.75775, 0.95017, -2.10064,
                            -2.76770, -2.32057, -2.28410, -2.83722, -2.45986, 1.78606, 0.27618, -0.54108, 2.49172,
                            -2.89120, -0.52169, -2.68745, 0.03343, -0.11790, 3.09736, -3.02081, -1.37294, -1.99044,
                            2.31419, -1.32603, 1.60068, -0.04666, 1.81738, -0.28430, -2.46847, -1.26959, 2.27198,
                            -0.83803, -0.98007, 2.58755, -2.09000, 0.98269, -2.82227, 1.92823, 0.86559, 2.27735,
                            1.34150, -1.16303, -0.63376, -0.07276, 0.39234, 1.33977, 1.88932, -0.83778, 1.88427,
                            0.37126, 1.13776, 2.62312, 0.62600, 1.28199, -0.55797, 2.12048, 0.85703, 1.50861, 1.45908,
                            2.99503, 0.04893, 1.38351, -0.47853, -1.81710, -2.53381, -0.17439, 0.45945, -0.70423,
                            0.10357, 1.14029, 0.34762, 2.70903, 2.34953, 1.29763, 2.76502, -2.88920, -1.20820};


/*
 * =====================================================================================================================
 *                                                  setUp and tearDown
 * =====================================================================================================================
 */

void setUp(void) {}

void tearDown(void) {}

/*
 * =====================================================================================================================
 *                                                  test expValObs
 * =====================================================================================================================
 */

/*
 * This function checks whether calculating the mean using expValObs coincides with calculating it using Kronecker
 * matrix-vector products on 1 to MAXQUBITS qubits.
 * Input:
 *      -
 * Output:
 *      There is no return value. Only checks whether the test and the reference results are equal to some precision for
 *      an observable which is a linear combination of all possible Pauli strings.
*/
void testExpValObs(void) {
    state_t testState;          // State equipped with the test statevectors
    obs_t obsDiag;              // Declare the diagonal observable
    obs_t obsPauli;             // Declare the Pauli observable
    pauliObs_t pauliObs;        // Pauli observable struct

    for (qubit_t qubits = 1; qubits <= MAXQUBITS; ++qubits) {

        dim_t dim = POW2(qubits, dim_t);                        // Hilbert space dimension
        cplx_t** testVectors = generateTestVectors(qubits);     // Statevectors altered by the test function
        cplx_t** refVectors = generateTestVectors(qubits);      // Statevectors altered by matrix multiplication
                                                                // (reference)
        /* Diagonal observable */
        pauli_t* comps = allPauliStringsDiag(qubits);           // Concatenation of all diagonal Pauli strings
        complength_t length = (1 << qubits);                    // Number of components, i.e., number of Pauli strings

        double* diagObs = pauliObservableDiag(comps, coeffs, length, qubits);
        obsDiag.type = DIAG;
        obsDiag.diagObs = diagObs;

        pauliObs.components = comps;
        pauliObs.coefficients = coeffs;
        pauliObs.length = length;
        pauliObs.qubits = qubits;
        obsPauli.type = PAULI;
        obsPauli.pauliObs = &pauliObs;

        /* Compute the observable's matrix representation */
        cplx_t* obsMatrix = pauliObservableMat(comps, coeffs, length, qubits);

        for (dim_t i = 0; i < dim + 1; ++i) {
            stateInitVector(&testState, testVectors[i], qubits);    // Initialize test state
            double result_diag = expValObs(&testState, &obsDiag);   // Apply test functions
            double result_pauli = expValObs(&testState, &obsPauli);

            cplx_t* refKet = cmatVecMul(obsMatrix, refVectors[i], dim);
            double reference = creal(cinnerProduct(refVectors[i], refKet, dim));

            TEST_ASSERT_TRUE(ralmostEqual(result_diag, reference, PRECISION));
            TEST_ASSERT_TRUE(ralmostEqual(result_pauli, reference, PRECISION));

            free(refKet);
        }

        /* Pauli observable */
        free(comps);
        comps = allPauliStrings(qubits);
        length = 1 << (2 * qubits);

        free(obsDiag.diagObs);

        pauliObs.components = comps;
        pauliObs.length = length;

        free(obsMatrix);
        obsMatrix = pauliObservableMat(comps, coeffs, length, qubits);

        for (dim_t i = 0; i < dim + 1; ++i) {
            stateInitVector(&testState, testVectors[i], qubits);    // Initialize test state
            double result = expValObs(&testState, &obsPauli);   // Apply test function

            cplx_t* refKet = cmatVecMul(obsMatrix, refVectors[i], dim);
            double reference = creal(cinnerProduct(refVectors[i], refKet, dim));

            TEST_ASSERT_TRUE(ralmostEqual(result, reference, PRECISION));

            free(refKet);
            stateFreeVector(&testState);
        }
        free(testVectors);
        freeTestVectors(refVectors, qubits);
        free(comps);
        free(obsMatrix);
    }
}


/* 
 * =====================================================================================================================
 *                                                  test expValObsPQC
 * =====================================================================================================================
 */

void testExpValObsPQC(void) {
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

    complength_t lengthEvoOps;
    pauli_t *compEvoOps;
    double *coeffEvoOps;
    pauliObs_t *evoOpsPauli;
    obs_t **evoOps;

    /*
                                                        1 qubit, circdepth = 2
    --------------------------------------------------------------------------------------------------------------------
    */

    qubits = 1;
    dim = POW2(qubits, dim_t);

    circdepth = 2;
    parameters = calloc(circdepth, sizeof(double));
    parameters[0] = -2.82317;
    parameters[1] = 0.24576;

    /*
     * Initialize the state vectors
     */
    testVectors = generateTestVectors(qubits);
    refVectors = generateTestVectors(qubits);

    /*
     * Define the diagonal observable
     */
    lengthObs = dim;
    compObs = calloc(qubits * lengthObs, sizeof(pauli_t));
    compObs[0] = ID;
    compObs[1] = Z;
    coeffObs = calloc(lengthObs, sizeof(double));
    coeffObs[0] = 1.13126;
    coeffObs[1] = 0.24848;

    /*
     * Initialize the diagonal observable
     */
    pauliObservable.components = compObs;               // Initialize the Pauli representation of the diagonal
    pauliObservable.coefficients = coeffObs;            // hamiltonian
    pauliObservable.length = lengthObs;                 //
    pauliObservable.qubits = qubits;                    //
    //
    obs_t observablePauli;                              //
    observablePauli.type = PAULI;                       //
    observablePauli.pauliObs = &pauliObservable;        //

    double *diagObservable = pauliObservableDiag(compObs, coeffObs, lengthObs, qubits);
    obs_t observableDiag;                                                       // Initialize the hamiltonian's diagonal
    observableDiag.type = DIAG;                                                 // representation
    observableDiag.diagObs = diagObservable;                                    //

    /*
     * Define the evolution operators
     */
    lengthEvoOps = dim;
    compEvoOps = calloc(dim * qubits, sizeof(pauli_t));
    compEvoOps[0] = ID;
    compEvoOps[1] = Z;
    coeffEvoOps = calloc(lengthEvoOps, sizeof(double));
    coeffEvoOps[0] = 0.53858;
    coeffEvoOps[1] = -1.09597;

    /*
     * Initialize the evolution operators
     */

    evoOpsPauli = calloc(circdepth, sizeof(pauliObs_t));
    for (depth_t i = 0; i < circdepth; ++i) {
        evoOpsPauli[i].components = compEvoOps;
        evoOpsPauli[i].coefficients = coeffEvoOps + (i * lengthEvoOps);
        evoOpsPauli[i].qubits = qubits;
        evoOpsPauli[i].length = lengthEvoOps;
    }

    evoOps = calloc(circdepth, sizeof(obs_t*));
    for (depth_t i = 0; i < circdepth; ++i) {
        evoOps[i] = calloc(1, sizeof(obs_t));
        evoOps[i]->type = PAULI;
        evoOps[i]->pauliObs = evoOpsPauli + i;
    }

    /*
     * Calculate the expectation values
     */

    for (dim_t i = 0; i < dim + 1; ++i) {
        stateInitVector(&testState, testVectors[i], qubits);
        double result = expValObsPQC(&testState, \
                                     parameters, \
                                     &observablePauli, \
                                     (const obs_t **) evoOps, \
                                     circdepth);

        double resultDiag = expValObsPQC(&testState, \
                                         parameters, \
                                         &observablePauli, \
                                         (const obs_t **) evoOps, \
                                         circdepth);

        double reference = finiteExpValPQC(refVectors[i], qubits, dim, circdepth, compObs, coeffObs, lengthObs, \
                                           compEvoOps, coeffEvoOps, lengthEvoOps, parameters);

        TEST_ASSERT_TRUE(fabs(result - reference) < PRECISION);
        TEST_ASSERT_TRUE(fabs(resultDiag - reference) < PRECISION);
    }

    /*
     * Free observable and evolution operators
     */
    free(compObs);
    free(coeffObs);

    /*
     * Define the Pauli observable
     */
    lengthObs = dim * dim;
    compObs = calloc(qubits * lengthObs, sizeof(pauli_t));
    compObs[0] = ID;
    compObs[1] = X;
    compObs[2] = Y;
    compObs[3] = Z;
    coeffObs = calloc(lengthObs, sizeof(double));
    coeffObs[0] = -0.16123;
    coeffObs[1] = 0.74276;
    coeffObs[2] = -0.61566;
    coeffObs[3] = 0.96877;

    /*
     * Initialize the Pauli observable
     */
    pauliObservable.components = compObs;
    pauliObservable.coefficients = coeffObs;
    pauliObservable.length = lengthObs;

    /*
     * Calculate the expectation values
     */

    for (dim_t i = 0; i < dim + 1; ++i) {
        stateInitVector(&testState, testVectors[i], qubits);
        double result = expValObsPQC(&testState, \
                                     parameters, \
                                     &observablePauli, \
                                     (const obs_t **) evoOps, \
                                     circdepth);

        double reference = finiteExpValPQC(refVectors[i], qubits, dim, circdepth, compObs, coeffObs, lengthObs, \
                                           compEvoOps, coeffEvoOps, lengthEvoOps, parameters);

        TEST_ASSERT_TRUE(fabs(result - reference) < PRECISION);
    }

    /*
     * Free test vectors
     */
    freeTestVectors(testVectors, qubits);
    freeTestVectors(refVectors, qubits);

    /*
     * Free observables and evolution operators
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
    parameters = calloc(circdepth, sizeof(double));
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
    compObs = calloc(qubits * lengthObs, sizeof(pauli_t));
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
    coeffObs = calloc(lengthObs, sizeof(double));
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
    compEvoOps = calloc(qubits * lengthEvoOps, sizeof(pauli_t));
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
    coeffEvoOps = calloc(circdepth * lengthEvoOps, sizeof(double));
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
    evoOpsPauli = calloc(circdepth, sizeof(pauliObs_t));
    for (depth_t i = 0; i < circdepth; ++i) {
        evoOpsPauli[i].components = compEvoOps;
        evoOpsPauli[i].coefficients = coeffEvoOps + (i * lengthEvoOps);
        evoOpsPauli[i].qubits = qubits;
        evoOpsPauli[i].length = lengthEvoOps;
    }

    evoOps = calloc(circdepth, sizeof(obs_t*));
    for (depth_t i = 0; i < circdepth; ++i) {
        evoOps[i] = calloc(1, sizeof(obs_t));
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
        stateInitVector(&testState, testVectors[i], qubits);
        double* result = gradientPQC(&testState, parameters, &observable, (const obs_t **) evoOps, circdepth);

        /*
         * Gradient from approxGradientPQC
         */
        double* result2 = approxGradientPQC(&testState, \
                                            parameters, \
                                            &observable, \
                                            (const obs_t**) evoOps, \
                                            circdepth, \
                                            1e-9);

        /*
         * Gradient from finite difference method
         */
        double* reference = finiteGradientPQC(refVectors[i], qubits, dim, circdepth, compObs, coeffObs, lengthObs, \
                                              compEvoOps, coeffEvoOps, lengthEvoOps, parameters, 1e-9);

        TEST_ASSERT_TRUE(rvectorAlmostEqual(result, reference, circdepth, APPROXPRECISION));
        TEST_ASSERT_TRUE(rvectorAlmostEqual(result2, reference, circdepth, APPROXPRECISION));

        free(result);
        free(result2);
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
    parameters = calloc(circdepth, sizeof(double));
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
    compObs = calloc(qubits * lengthObs, sizeof(pauli_t));
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
    coeffObs = calloc(lengthObs, sizeof(double));
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
    compEvoOps = calloc(qubits * lengthEvoOps, sizeof(pauli_t));
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
    coeffEvoOps = calloc(circdepth * lengthEvoOps, sizeof(double));
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
    evoOpsPauli = calloc(circdepth, sizeof(pauliObs_t));
    for (depth_t i = 0; i < circdepth; ++i) {
        evoOpsPauli[i].components = compEvoOps;
        evoOpsPauli[i].coefficients = coeffEvoOps + (i * lengthEvoOps);
        evoOpsPauli[i].qubits = qubits;
        evoOpsPauli[i].length = lengthEvoOps;
    }

    evoOps = calloc(circdepth, sizeof(obs_t*));
    for (depth_t i = 0; i < circdepth; ++i) {
        evoOps[i] = calloc(1, sizeof(obs_t));
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
        stateInitVector(&testState, testVectors[i], qubits);
        double* result = gradientPQC(&testState, parameters, &observable, (const obs_t **) evoOps, circdepth);

        /*
         * Gradient from approxGradientPQC
         */
        double* result2 = approxGradientPQC(&testState, \
                                            parameters, \
                                            &observable, \
                                            (const obs_t**) evoOps, \
                                            circdepth, \
                                            1e-9);

        /*
         * Gradient from finite difference method
         */
        double* reference = finiteGradientPQC(refVectors[i], qubits, dim, circdepth, compObs, coeffObs, lengthObs, \
                                              compEvoOps, coeffEvoOps, lengthEvoOps, parameters, 1e-9);

        TEST_ASSERT_TRUE(rvectorAlmostEqual(result, reference, circdepth, APPROXPRECISION));
        TEST_ASSERT_TRUE(rvectorAlmostEqual(result2, reference, circdepth, APPROXPRECISION));

        free(result);
        free(result2);
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
    parameters = calloc(circdepth, sizeof(double));
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
    compObs = calloc(qubits * lengthObs, sizeof(pauli_t));
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
    coeffObs = calloc(lengthObs, sizeof(double));
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
    compEvoOps = calloc(qubits * lengthEvoOps, sizeof(pauli_t));
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
    coeffEvoOps = calloc(circdepth * lengthEvoOps, sizeof(double));
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
    evoOpsPauli = calloc(circdepth, sizeof(pauliObs_t));
    for (depth_t i = 0; i < circdepth; ++i) {
        evoOpsPauli[i].components = compEvoOps;
        evoOpsPauli[i].coefficients = coeffEvoOps + (i * lengthEvoOps);
        evoOpsPauli[i].qubits = qubits;
        evoOpsPauli[i].length = lengthEvoOps;
    }

    evoOps = calloc(circdepth, sizeof(obs_t*));
    for (depth_t i = 0; i < circdepth; ++i) {
        evoOps[i] = calloc(1, sizeof(obs_t));
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
        stateInitVector(&testState, testVectors[i], qubits);
        double* result = gradientPQC(&testState, parameters, &observable, (const obs_t **) evoOps, circdepth);

        /*
         * Gradient from approxGradientPQC
         */
        double* result2 = approxGradientPQC(&testState, \
                                            parameters, \
                                            &observable, \
                                            (const obs_t**) evoOps, \
                                            circdepth, \
                                            1e-9);

        /*
         * Gradient from finite difference method
         */
        double* reference = finiteGradientPQC(refVectors[i], qubits, dim, circdepth, compObs, coeffObs, lengthObs, \
                                              compEvoOps, coeffEvoOps, lengthEvoOps, parameters, 1e-9);

        TEST_ASSERT_TRUE(rvectorAlmostEqual(result, reference, circdepth, APPROXPRECISION));
        TEST_ASSERT_TRUE(rvectorAlmostEqual(result2, reference, circdepth, APPROXPRECISION));

        free(result);
        free(result2);
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
    compObs = calloc(qubits * lengthObs, sizeof(pauli_t));
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
    coeffObs = calloc(lengthObs, sizeof(double));
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
    compEvoOps = calloc(qubits * lengthEvoOps, sizeof(pauli_t));
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
    coeffEvoOps = calloc(circdepth * lengthEvoOps, sizeof(double));
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
    evoOpsPauli = calloc(circdepth, sizeof(pauliObs_t));
    for (depth_t i = 0; i < circdepth; ++i) {
        evoOpsPauli[i].components = compEvoOps;
        evoOpsPauli[i].coefficients = coeffEvoOps + (i * lengthEvoOps);
        evoOpsPauli[i].qubits = qubits;
        evoOpsPauli[i].length = lengthEvoOps;
    }

    for (depth_t i = 0; i < circdepth; ++i) {
        evoOps[i] = calloc(1, sizeof(obs_t));
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
        stateInitVector(&testState, testVectors[i], qubits);
        double* result = gradientPQC(&testState, parameters, &observable, (const obs_t **) evoOps, circdepth);

        /*
         * Gradient from approxGradientPQC
         */
        double* result2 = approxGradientPQC(&testState, \
                                            parameters, \
                                            &observable, \
                                            (const obs_t**) evoOps, \
                                            circdepth, \
                                            1e-9);

        /*
         * Gradient from finite difference method
         */
        double* reference = finiteGradientPQC(refVectors[i], qubits, dim, circdepth, compObs, coeffObs, lengthObs, \
                                              compEvoOps, coeffEvoOps, lengthEvoOps, parameters, 1e-9);

        TEST_ASSERT_TRUE(rvectorAlmostEqual(result, reference, circdepth, APPROXPRECISION));
        TEST_ASSERT_TRUE(rvectorAlmostEqual(result2, reference, circdepth, APPROXPRECISION));

        free(result);
        free(result2);
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
 * =====================================================================================================================
 *                                                      test hessianPQC
 * =====================================================================================================================
 */

double* grad(const state_t* state, const double params[], const obs_t* observable, \
                    const obs_t* evoOps[], depth_t circdepth) {
    return gradientPQC(state, params, observable, evoOps, circdepth);
}

void testHessianPQC(void) {
    qubit_t qubits;
    dim_t dim;
    depth_t circdepth;
    double* parameters;
    double epsilon = 1e-6;

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
    parameters = calloc(circdepth, sizeof(double));
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
    lengthObs = 16;
    compObs = calloc(qubits * lengthObs, sizeof(pauli_t));
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
    coeffObs = calloc(lengthObs, sizeof(double));
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
    lengthEvoOps = 8;
    compEvoOps = calloc(qubits * lengthEvoOps, sizeof(pauli_t));
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
    coeffEvoOps = calloc(circdepth * lengthEvoOps, sizeof(double));
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
    evoOpsPauli = calloc(circdepth, sizeof(pauliObs_t));
    for (depth_t i = 0; i < circdepth; ++i) {
        evoOpsPauli[i].components = compEvoOps;
        evoOpsPauli[i].coefficients = coeffEvoOps + (i * lengthEvoOps);
        evoOpsPauli[i].qubits = qubits;
        evoOpsPauli[i].length = lengthEvoOps;
    }

    evoOps = calloc(circdepth, sizeof(obs_t*));
    for (depth_t i = 0; i < circdepth; ++i) {
        evoOps[i] = calloc(1, sizeof(obs_t));
        evoOps[i]->type = PAULI;
        evoOps[i]->pauliObs = evoOpsPauli + i;
    }

    /*
     * Test objective:
     * Compare the hessian obtained from hessianPQC with the approximation
     *      H(x).(epsilon * vec{e_{i}}) = (grad f)(x + epsilon * vec{e_{i}}) - (grad f)(x)
     */
    for (dim_t i = 0; i < dim + 1; ++i) {
        stateInitVector(&testState, testVectors[i], qubits);

        /*
         * 1. Calculate the hessian via hessianPQC
         */
        double* result = hessianPQC(&testState, \
                                    parameters, \
                                    &observable, \
                                    (const obs_t **) evoOps, \
                                    circdepth);

        /*
         * 2. Calculate the hessian via approximateHessianPQC
         */
        double* result2 = approxHessianPQC(&testState, \
                                           parameters, \
                                           &observable, \
                                           (const obs_t**) evoOps, \
                                           circdepth, \
                                           grad, \
                                           epsilon);

        /*
         * 3. Calculate the hessian via finite difference method
         */
        double* reference = finiteHessianPQC(refVectors[i],qubits, dim, circdepth, compObs, coeffObs, lengthObs, \
                                             compEvoOps, coeffEvoOps, lengthEvoOps, parameters, epsilon);

        TEST_ASSERT_TRUE(rvectorAlmostEqual(result, reference, circdepth, 1e-1));
        TEST_ASSERT_TRUE(rvectorAlmostEqual(result2, reference, circdepth, 1e-1));

        free(result);
        free(result2);
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
 * =====================================================================================================================
 *                                                  test momentMatQAOA
 * =====================================================================================================================
 */

void applyU(state_t* state, const obs_t* srchOp) {
    return (evolveWithTrotterizedObservable(state, srchOp, 0.05));
}

void testMomMat(void) {
    qubit_t qubits;
    dim_t dim;
    state_t testState;
    cplx_t** testVectors;

    complength_t lengthObs;
    depth_t obsc;
    pauli_t* compObs;
    double* coeffObs;
    pauliObs_t* obsPauli;
    obs_t** observables;

    complength_t lengthSrchObs;
    depth_t dimMat;
    pauli_t* srchComps;
    double* srchCoeffs;
    pauliObs_t* srchObsPauli;
    obs_t** srchObs;

    /*
     *                                                  4 qubits
     *------------------------------------------------------------------------------------------------------------------
     */
    qubits = 4;
    dim = POW2(4, dim_t);
    testVectors = generateTestVectors(qubits);

    /* Initialize two observables to be diagonal in the computational basis */

    lengthObs = dim;
    obsc = 2;
    compObs = (pauli_t*) calloc(qubits * lengthObs, sizeof(pauli_t));
    compObs[0] = ID;
    compObs[1] = ID;
    compObs[2] = ID;
    compObs[3] = ID;
    compObs[4] = ID;
    compObs[5] = ID;
    compObs[6] = ID;
    compObs[7] = Z;
    compObs[8] = ID;
    compObs[9] = ID;
    compObs[10] = Z;
    compObs[11] = ID;
    compObs[12] = ID;
    compObs[13] = ID;
    compObs[14] = Z;
    compObs[15] = Z;
    compObs[16] = ID;
    compObs[17] = Z;
    compObs[18] = ID;
    compObs[19] = ID;
    compObs[20] = ID;
    compObs[21] = Z;
    compObs[22] = ID;
    compObs[23] = Z;
    compObs[24] = ID;
    compObs[25] = Z;
    compObs[26] = Z;
    compObs[27] = ID;
    compObs[28] = ID;
    compObs[29] = Z;
    compObs[30] = Z;
    compObs[31] = Z;
    compObs[32] = Z;
    compObs[33] = ID;
    compObs[34] = ID;
    compObs[35] = ID;
    compObs[36] = Z;
    compObs[37] = ID;
    compObs[38] = ID;
    compObs[39] = Z;
    compObs[40] = Z;
    compObs[41] = ID;
    compObs[42] = Z;
    compObs[43] = ID;
    compObs[44] = Z;
    compObs[45] = ID;
    compObs[46] = Z;
    compObs[47] = Z;
    compObs[48] = Z;
    compObs[49] = Z;
    compObs[50] = ID;
    compObs[51] = ID;
    compObs[52] = Z;
    compObs[53] = Z;
    compObs[54] = ID;
    compObs[55] = Z;
    compObs[56] = Z;
    compObs[57] = Z;
    compObs[58] = Z;
    compObs[59] = ID;
    compObs[60] = Z;
    compObs[61] = Z;
    compObs[62] = Z;
    compObs[63] = Z;

    coeffObs = calloc(obsc * lengthObs, sizeof(double));
    coeffObs[0] = -3.02593;
    coeffObs[1] = -0.16359;
    coeffObs[2] = -1.10004;
    coeffObs[3] = 1.56744;
    coeffObs[4] = -0.60557;
    coeffObs[5] = -1.43826;
    coeffObs[6] = 0.57069;
    coeffObs[7] = -0.30491;
    coeffObs[8] = -1.90849;
    coeffObs[9] = -0.26485;
    coeffObs[10] = -1.23567;
    coeffObs[11] = -0.84204;
    coeffObs[12] = -0.43410;
    coeffObs[13] = -0.27031;
    coeffObs[14] = 3.07625;
    coeffObs[15] = -1.91561;
    coeffObs[16] = -0.20932;
    coeffObs[17] = 2.15950;
    coeffObs[18] = -1.72770;
    coeffObs[19] = 1.57711;
    coeffObs[20] = 0.07107;
    coeffObs[21] = -1.22595;
    coeffObs[22] = 1.26687;
    coeffObs[23] = 2.54369;
    coeffObs[24] = -2.93828;
    coeffObs[25] = 0.08806;
    coeffObs[26] = -2.63439;
    coeffObs[27] = -0.63076;
    coeffObs[28] = -0.81555;
    coeffObs[29] = 0.80927;
    coeffObs[30] = -0.20716;
    coeffObs[31] = 3.04246;

    obsPauli = calloc(obsc, sizeof(pauliObs_t));
    for (depth_t i = 0; i < obsc; ++i) {
        obsPauli[i].components = compObs;
        obsPauli[i].coefficients = coeffObs + (i * lengthObs);
        obsPauli[i].length = lengthObs;
        obsPauli[i].qubits = qubits;
    }

    observables = calloc(obsc, sizeof(obs_t*));
    for (depth_t i = 0; i < obsc; ++i) {
        observables[i] = calloc(1, sizeof(obs_t));
        observables[i]->type = PAULI;
        observables[i]->pauliObs = obsPauli + i;
    }

    /* Initialize the rotations around the three axes and the identity as search unitaries */
    lengthSrchObs = qubits;
    dimMat = 4;
    srchComps = calloc(dimMat * lengthSrchObs * qubits, sizeof(pauli_t));
    srchComps[0] = ID;
    srchComps[1] = ID;
    srchComps[2] = ID;
    srchComps[3] = ID;
    srchComps[4] = ID;
    srchComps[5] = ID;
    srchComps[6] = ID;
    srchComps[7] = ID;
    srchComps[8] = ID;
    srchComps[9] = ID;
    srchComps[10] = ID;
    srchComps[11] = ID;
    srchComps[12] = ID;
    srchComps[13] = ID;
    srchComps[14] = ID;
    srchComps[15] = ID;
    srchComps[16] = X;
    srchComps[17] = ID;
    srchComps[18] = ID;
    srchComps[19] = ID;
    srchComps[20] = ID;
    srchComps[21] = X;
    srchComps[22] = ID;
    srchComps[23] = ID;
    srchComps[24] = ID;
    srchComps[25] = ID;
    srchComps[26] = X;
    srchComps[27] = ID;
    srchComps[28] = ID;
    srchComps[29] = ID;
    srchComps[30] = ID;
    srchComps[31] = X;
    srchComps[32] = Y;
    srchComps[33] = ID;
    srchComps[34] = ID;
    srchComps[35] = ID;
    srchComps[36] = ID;
    srchComps[37] = Y;
    srchComps[38] = ID;
    srchComps[39] = ID;
    srchComps[40] = ID;
    srchComps[41] = ID;
    srchComps[42] = Y;
    srchComps[43] = ID;
    srchComps[44] = ID;
    srchComps[45] = ID;
    srchComps[46] = ID;
    srchComps[47] = Y;
    srchComps[48] = Z;
    srchComps[49] = ID;
    srchComps[50] = ID;
    srchComps[51] = ID;
    srchComps[52] = ID;
    srchComps[53] = Z;
    srchComps[54] = ID;
    srchComps[55] = ID;
    srchComps[56] = ID;
    srchComps[57] = ID;
    srchComps[58] = Z;
    srchComps[59] = ID;
    srchComps[60] = ID;
    srchComps[61] = ID;
    srchComps[62] = ID;
    srchComps[63] = Z;

    srchCoeffs = calloc(lengthObs, sizeof(double));
    srchCoeffs[0] = 1.;
    srchCoeffs[1] = 1.;
    srchCoeffs[2] = 1.;
    srchCoeffs[3] = 1.;

    srchObsPauli = calloc(dimMat, sizeof(pauliObs_t));
    for (depth_t i = 0; i < dimMat; ++i) {
        srchObsPauli[i].components = srchComps + (i * lengthSrchObs * qubits);
        srchObsPauli[i].coefficients = srchCoeffs;
        srchObsPauli[i].length = lengthSrchObs;
        srchObsPauli[i].qubits = qubits;
    }

    srchObs = calloc(dimMat, sizeof(obs_t*));
    for (depth_t i = 0; i < dimMat; ++i) {
        srchObs[i] = calloc(1, sizeof(obs_t));
        srchObs[i]->type = PAULI;
        srchObs[i]->pauliObs = srchObsPauli + i;
    }

    /* Define and allocate the reference matrix */
    cplx_t** reference = calloc(obsc, sizeof(cplx_t*));
    for (depth_t i = 0; i < obsc; ++i) {
        reference[i] = calloc(dimMat * dimMat, sizeof(cplx_t));
    }

    /* Define and allocate the observables' matrices */
    cplx_t** obsMat = calloc(obsc, sizeof(cplx_t*));

    /* Compute the observables' matrices */
    for (depth_t i = 0; i < obsc; ++i) {
        obsMat[i] = pauliObservableMat(compObs, \
                                       coeffObs + (i * lengthObs), \
                                       lengthObs, \
                                       qubits);
    }

    /* Define and allocate the search unitaries' matrices */
    cplx_t** srchObsMat = calloc(dimMat, sizeof(cplx_t*));

    /* Compute the search unitary matrices from the search operators */
    for (depth_t i = 0; i < dimMat; ++i) {
        srchObsMat[i] = expTrotterizedPauliObservableMat(srchComps + (i * lengthSrchObs * qubits), \
                                                         srchCoeffs, \
                                                         lengthSrchObs , \
                                                         0.05, qubits);
    }

    for (dim_t i = 0; i < dim + 1; ++i) {

        /* Compute the moment matrices using the exact function */
        stateInitVector(&testState, testVectors[i], qubits);
        cplx_t** result = momMat(&testState, \
                                 (const obs_t **) observables, \
                                 obsc, \
                                 (const obs_t **) srchObs, \
                                 dimMat, \
                                 applyU);

        /* Compute the reference moment matrices in column major form using matrix multiplication */
        for (depth_t j = 0; j < dimMat; ++j) {
            cplx_t* refKet = cmatVecMul(srchObsMat[j], testState.vector, dim);

            for (depth_t k = 0; k < obsc; ++k) {
                cplx_t* refBra = cmatVecMul(obsMat[k], refKet, dim);
                reference[k][j * dimMat + j] = creal(cinnerProduct(refBra, refKet, dim));
                free(refBra);
            }

            for (depth_t k = 0; k < dimMat; ++k) {
                cplx_t* tmpBra = cmatVecMul(srchObsMat[k], testState.vector, dim);
                for (depth_t l = 0; l < obsc; ++l) {
                    cplx_t* refBra = cmatVecMul(obsMat[l], tmpBra, dim);
                    reference[l][j * dimMat + k] = cinnerProduct(refBra, refKet, dim);
                    reference[l][k * dimMat + j] = conj(reference[l][j * dimMat + k]);
                    free(refBra);
                }
                free(tmpBra);
            }
            free(refKet);
        }

        for (depth_t j = 0; j < obsc; ++j) {
            TEST_ASSERT_TRUE(cvectorAlmostEqual(result[j], reference[j], dimMat * dimMat, 1e-6));
        }

        free(result);
        stateFreeVector(&testState);
    }

    /* Free all allocated memory */
    for (depth_t i = 0; i < obsc; ++i) {
        free(reference[i]);
    }
    free(reference);
    free(testVectors);
    free(compObs);
    free(coeffObs);
    free(obsPauli);
    for (depth_t i = 0; i < obsc; ++i) {
        free(observables[i]);
    }
    free(observables);
    for (depth_t i = 0; i < obsc; ++i) {
        free(obsMat[i]);
    }
    free(obsMat);

    free(srchComps);
    free(srchCoeffs);
    free(srchObsPauli);
    for (depth_t i = 0; i < dimMat; ++i) {
        free(srchObs[i]);
    }
    free(srchObs);
    for (depth_t i = 0; i < dimMat; ++i) {
        free(srchObsMat[i]);
    }
    free(srchObsMat);
}

void testMomMatPQC(void) {
    qubit_t qubits;
    dim_t dim;
    state_t testState;
    cplx_t** testVectors;

    complength_t lengthObs;
    depth_t obsc;
    pauli_t* compObs;
    double* coeffObs;
    pauliObs_t* obsPauli;
    obs_t** observables;

    complength_t lengthSrchObs;
    depth_t dimMat;
    pauli_t* srchComps;
    double* srchCoeffs;
    pauliObs_t* srchObsPauli;
    obs_t** srchObs;

    /*
     *                                                  4 qubits
     *------------------------------------------------------------------------------------------------------------------
     */
    qubits = 4;
    dim = POW2(4, dim_t);
    testVectors = generateTestVectors(qubits);


    /* Initialize two observables to be diagonal in the computational basis */

    lengthObs = dim;
    obsc = 2;
    compObs = calloc(qubits * lengthObs, sizeof(pauli_t));
    compObs[0] = ID;
    compObs[1] = ID;
    compObs[2] = ID;
    compObs[3] = ID;
    compObs[4] = ID;
    compObs[5] = ID;
    compObs[6] = ID;
    compObs[7] = Z;
    compObs[8] = ID;
    compObs[9] = ID;
    compObs[10] = Z;
    compObs[11] = ID;
    compObs[12] = ID;
    compObs[13] = ID;
    compObs[14] = Z;
    compObs[15] = Z;
    compObs[16] = ID;
    compObs[17] = Z;
    compObs[18] = ID;
    compObs[19] = ID;
    compObs[20] = ID;
    compObs[21] = Z;
    compObs[22] = ID;
    compObs[23] = Z;
    compObs[24] = ID;
    compObs[25] = Z;
    compObs[26] = Z;
    compObs[27] = ID;
    compObs[28] = ID;
    compObs[29] = Z;
    compObs[30] = Z;
    compObs[31] = Z;
    compObs[32] = Z;
    compObs[33] = ID;
    compObs[34] = ID;
    compObs[35] = ID;
    compObs[36] = Z;
    compObs[37] = ID;
    compObs[38] = ID;
    compObs[39] = Z;
    compObs[40] = Z;
    compObs[41] = ID;
    compObs[42] = Z;
    compObs[43] = ID;
    compObs[44] = Z;
    compObs[45] = ID;
    compObs[46] = Z;
    compObs[47] = Z;
    compObs[48] = Z;
    compObs[49] = Z;
    compObs[50] = ID;
    compObs[51] = ID;
    compObs[52] = Z;
    compObs[53] = Z;
    compObs[54] = ID;
    compObs[55] = Z;
    compObs[56] = Z;
    compObs[57] = Z;
    compObs[58] = Z;
    compObs[59] = ID;
    compObs[60] = Z;
    compObs[61] = Z;
    compObs[62] = Z;
    compObs[63] = Z;

    coeffObs = calloc(obsc * lengthObs, sizeof(double));
    coeffObs[0] = -3.02593;
    coeffObs[1] = -0.16359;
    coeffObs[2] = -1.10004;
    coeffObs[3] = 1.56744;
    coeffObs[4] = -0.60557;
    coeffObs[5] = -1.43826;
    coeffObs[6] = 0.57069;
    coeffObs[7] = -0.30491;
    coeffObs[8] = -1.90849;
    coeffObs[9] = -0.26485;
    coeffObs[10] = -1.23567;
    coeffObs[11] = -0.84204;
    coeffObs[12] = -0.43410;
    coeffObs[13] = -0.27031;
    coeffObs[14] = 3.07625;
    coeffObs[15] = -1.91561;
    coeffObs[16] = -0.20932;
    coeffObs[17] = 2.15950;
    coeffObs[18] = -1.72770;
    coeffObs[19] = 1.57711;
    coeffObs[20] = 0.07107;
    coeffObs[21] = -1.22595;
    coeffObs[22] = 1.26687;
    coeffObs[23] = 2.54369;
    coeffObs[24] = -2.93828;
    coeffObs[25] = 0.08806;
    coeffObs[26] = -2.63439;
    coeffObs[27] = -0.63076;
    coeffObs[28] = -0.81555;
    coeffObs[29] = 0.80927;
    coeffObs[30] = -0.20716;
    coeffObs[31] = 3.04246;

    obsPauli = (pauliObs_t*) malloc(obsc * sizeof(pauliObs_t));
    for (depth_t i = 0; i < obsc; ++i) {
        obsPauli[i].components = compObs;
        obsPauli[i].coefficients = coeffObs + (i * lengthObs);
        obsPauli[i].length = lengthObs;
        obsPauli[i].qubits = qubits;
    }

    observables = calloc(obsc, sizeof(obs_t*));
    for (depth_t i = 0; i < obsc; ++i) {
        observables[i] = calloc(1, sizeof(obs_t));
        observables[i]->type = PAULI;
        observables[i]->pauliObs = obsPauli + i;
    }

    /* Initialize the rotations around the three axes and the identity as search unitaries */
    lengthSrchObs = qubits;
    dimMat = 4;
    srchComps = calloc(dimMat * lengthSrchObs * qubits, sizeof(pauli_t));
    srchComps[0] = ID;
    srchComps[1] = ID;
    srchComps[2] = ID;
    srchComps[3] = ID;
    srchComps[4] = ID;
    srchComps[5] = ID;
    srchComps[6] = ID;
    srchComps[7] = ID;
    srchComps[8] = ID;
    srchComps[9] = ID;
    srchComps[10] = ID;
    srchComps[11] = ID;
    srchComps[12] = ID;
    srchComps[13] = ID;
    srchComps[14] = ID;
    srchComps[15] = ID;
    srchComps[16] = X;
    srchComps[17] = ID;
    srchComps[18] = ID;
    srchComps[19] = ID;
    srchComps[20] = ID;
    srchComps[21] = X;
    srchComps[22] = ID;
    srchComps[23] = ID;
    srchComps[24] = ID;
    srchComps[25] = ID;
    srchComps[26] = X;
    srchComps[27] = ID;
    srchComps[28] = ID;
    srchComps[29] = ID;
    srchComps[30] = ID;
    srchComps[31] = X;
    srchComps[32] = Y;
    srchComps[33] = ID;
    srchComps[34] = ID;
    srchComps[35] = ID;
    srchComps[36] = ID;
    srchComps[37] = Y;
    srchComps[38] = ID;
    srchComps[39] = ID;
    srchComps[40] = ID;
    srchComps[41] = ID;
    srchComps[42] = Y;
    srchComps[43] = ID;
    srchComps[44] = ID;
    srchComps[45] = ID;
    srchComps[46] = ID;
    srchComps[47] = Y;
    srchComps[48] = Z;
    srchComps[49] = ID;
    srchComps[50] = ID;
    srchComps[51] = ID;
    srchComps[52] = ID;
    srchComps[53] = Z;
    srchComps[54] = ID;
    srchComps[55] = ID;
    srchComps[56] = ID;
    srchComps[57] = ID;
    srchComps[58] = Z;
    srchComps[59] = ID;
    srchComps[60] = ID;
    srchComps[61] = ID;
    srchComps[62] = ID;
    srchComps[63] = Z;

    srchCoeffs = calloc(lengthSrchObs, sizeof(double));
    srchCoeffs[0] = 1.;
    srchCoeffs[1] = 1.;
    srchCoeffs[2] = 1.;
    srchCoeffs[3] = 1.;

    double* angles = calloc(dimMat, sizeof(double));
    angles[0] = 0.05;
    angles[1] = 0.05;
    angles[2] = 0.05;
    angles[3] = 0.05;

    srchObsPauli = calloc(dimMat, sizeof(pauliObs_t));
    for (depth_t i = 0; i < dimMat; ++i) {
        srchObsPauli[i].components = srchComps + (i * lengthSrchObs * qubits);
        srchObsPauli[i].coefficients = srchCoeffs;
        srchObsPauli[i].length = lengthSrchObs;
        srchObsPauli[i].qubits = qubits;
    }

    srchObs = calloc(dimMat, sizeof(obs_t*));
    for (depth_t i = 0; i < dimMat; ++i) {
        srchObs[i] = calloc(1, sizeof(obs_t));
        srchObs[i]->type = PAULI;
        srchObs[i]->pauliObs = srchObsPauli + i;
    }

    /* Define and allocate the reference matrix */
    cplx_t** reference = calloc(obsc, sizeof(cplx_t*));
    for (depth_t i = 0; i < obsc; ++i) {
        reference[i] = (cplx_t*) malloc((dimMat + 1) * (dimMat + 1) * sizeof(cplx_t));
    }

    /* Define and allocate the observables' matrices */
    cplx_t** obsMat = calloc(obsc, sizeof(cplx_t*));

    /* Compute the observables' matrices */
    for (depth_t i = 0; i < obsc; ++i) {
        obsMat[i] = pauliObservableMat(compObs, \
                                       coeffObs + (i * lengthObs), \
                                       lengthObs, \
                                       qubits);
    }

    /* Define and allocate the search unitaries' matrices */
    cplx_t** srchObsMat = (cplx_t**) malloc((dimMat + 1) * sizeof(cplx_t*));

    srchObsMat[0] = singleQubitGateMat(IDMAT, qubits, 0);


    /* Compute the search unitary matrices from the search operators */
    for (depth_t i = 0; i < dimMat; ++i) {
        srchObsMat[i + 1] = expTrotterizedPauliObservableMat(srchComps + (i * lengthSrchObs * qubits), \
                                                         srchCoeffs, \
                                                         lengthSrchObs , \
                                                         0.05, qubits);
    }

    for (dim_t i = 0; i < dim + 1; ++i) {

        /* Compute the moment matrices using the exact function */
        stateInitVector(&testState, testVectors[i], qubits);
        cplx_t** result = momMatPQC(&testState, \
                                    (const obs_t **) observables, \
                                    obsc, \
                                    (const obs_t **) srchObs, \
                                    dimMat, \
                                    angles);
        printf("Result: \n");
        for (depth_t j = 0; j < obsc; ++j) {
            matrixPrint(result[j], dimMat + 1);
        }

        /* Compute the reference moment matrices in column major form using matrix multiplication */
        for (depth_t j = 0; j < (dimMat + 1); ++j) {
            cplx_t* refKet = cmatVecMul(srchObsMat[j], testState.vector, dim);

            for (depth_t k = 0; k < obsc; ++k) {
                cplx_t* refBra = cmatVecMul(obsMat[k], refKet, dim);
                reference[k][j * (dimMat + 1) + j] = creal(cinnerProduct(refBra, refKet, dim));
                free(refBra);
            }

            for (depth_t k = j + 1; k < (dimMat + 1); ++k) {
                cplx_t* tmpBra = cmatVecMul(srchObsMat[k], testState.vector, dim);
                for (depth_t l = 0; l < obsc; ++l) {
                    cplx_t* refBra = cmatVecMul(obsMat[l], tmpBra, dim);
                    reference[l][j * (dimMat + 1) + k] = cinnerProduct(refBra, refKet, dim);
                    reference[l][k * (dimMat + 1) + j] = conj(reference[l][j * (dimMat + 1) + k]);
                    free(refBra);
                }
                free(tmpBra);
            }
            free(refKet);
        }

        printf("Reference: \n");
        for (depth_t j = 0; j < obsc; ++j) {
            matrixPrint(reference[j], dimMat + 1);
        }

        for (depth_t j = 0; j < obsc; ++j) {
            TEST_ASSERT_TRUE(cvectorAlmostEqual(result[j], reference[j], (dimMat + 1) * (dimMat + 1), 1e-6));
        }

        free(result);
        stateFreeVector(&testState);
    }

    /* Free all allocated memory */
    for (depth_t i = 0; i < obsc; ++i) {
        free(reference[i]);
    }
    free(reference);
    free(testVectors);
    free(compObs);
    free(coeffObs);
    free(obsPauli);
    for (depth_t i = 0; i < obsc; ++i) {
        free(observables[i]);
    }
    free(observables);
    for (depth_t i = 0; i < obsc; ++i) {
        free(obsMat[i]);
    }
    free(obsMat);

    free(srchComps);
    free(srchCoeffs);
    free(srchObsPauli);
    for (depth_t i = 0; i < dimMat; ++i) {
        free(srchObs[i]);
    }
    free(srchObs);
    for (depth_t i = 0; i < dimMat; ++i) {
        free(srchObsMat[i]);
    }
    free(srchObsMat);
}

/*
 * =====================================================================================================================
 *                                                      main
 * =====================================================================================================================
 */

int main(void) {
    UNITY_BEGIN();
    RUN_TEST(testExpValObs);
    RUN_TEST(testExpValObsPQC);
    RUN_TEST(testGradientPQC);
//    RUN_TEST(testHessianPQC);
    RUN_TEST(testMomMat);
    RUN_TEST(testMomMatPQC);
    return UNITY_END();
}
