/* 
 * =====================================================================================================================
 *                                                      includes
 * =====================================================================================================================
 */

#ifndef EXACT_H
#include "exact.h"
#endif

#ifndef FINDIFF_H
#include "finDiff.h"
#endif

#ifndef HELPFUNC_H
#include "utils.h"
#endif

#ifndef UNITY_FRAMEWORK_H
#include "unity.h"
#endif


/*
 * =====================================================================================================================
 *                                                      macros
 * =====================================================================================================================
 */

#define APPROXPRECISION         1e-4
#define DELTA                   1e-6
#define EPSILON                 1e-9
#define MAXCIRCDEPTH            4
#define MAXQUBITS               4
#define PRECISION               1e-8
#define PSTRINGSC               ((1 << (2 * MAXQUBITS)) * MAXCIRCDEPTH)

double coeffs[PSTRINGSC] = {0.51139, 1.95053, 0.66727, -2.11905, 2.26233, 1.35227, -1.73031, -1.15815, -1.91041,
                            1.87147, 3.09033, 1.02791, -0.99938, -1.88385, -1.39707, 1.65929, -3.09955, -2.17565,
                            2.48871, 0.51081, -2.58661, 1.49061, 1.88940, -1.34715, 3.03955, 1.23211, 2.80910, -0.60036,
                            -2.64360, 1.40809, -1.96627, -2.87907, -2.88355, 0.63188, -0.61143, -2.24444, -0.62184,
                            -2.90745, 2.90393, -1.50515, -0.04786, 2.02768, -2.43060, 1.83989, 0.95332, -1.52473,
                            2.97657, 2.86420, 2.49007, 2.12637, -2.58924, 0.25564, 0.34811, -1.73592, 1.56375, 2.82349,
                            -1.82513, 0.76837, -0.18511, -0.16624, -2.45121, -1.75851, -3.10641, 2.35299, 1.43762,
                            2.76983, -1.96671, -0.73667, 1.69139, -1.07372, -2.92576, 0.02534, -0.71348, 0.07071,
                            -3.02836, -0.51941, 1.61372, -2.78892, -3.06266, -0.37482, 3.12958, 1.41753, -0.43756,
                            3.10822, -0.40606, -0.29795, 0.96884, -2.28602, 1.27814, -2.80878, 1.25143, 0.38066,
                            -1.68643, 0.36258, -0.68315, 0.12629, 2.49414, 2.26605, -0.85737, -1.42797, -2.06949,
                            -2.08796, 2.20261, -1.43226, -1.69188, -0.90615, -2.51329, -0.69795, -2.28935, -0.04953,
                            0.59151, -1.83018, -1.44849, 1.23337, 2.59623, 1.12234, 1.39195, -0.51915, 2.38677, 0.76239,
                            1.59011, 2.62509, 1.60008, 2.66190, 1.60987, 0.95967, 1.16755, 0.04462, -1.88234, 1.64222,
                            2.82077, 0.33606, 0.94486, -0.32558, -1.90228, 0.62116, 3.00294, 3.06838, 1.21646, 0.22309,
                            1.15924, -1.70447, -1.86835, -0.13829, -1.78180, -2.50856, -3.12674, -2.99372, 0.39346,
                            -0.39182, -0.16177, 0.28365, -2.20857, 0.06332, -2.38298, -1.91661, -1.54321, 1.54211,
                            1.76292, 2.46719, 1.07745, -1.85114, 2.76741, -2.39311, 2.22944, -2.14911, 2.11646,
                            -0.37847, -1.10131, 2.40082, 0.84214, 0.96206, 1.01660, -0.78703, -0.11848, 1.36139,
                            -0.71425, -0.93528, -1.76758, 1.21787, 0.22640, -0.74087, -2.81168, -0.51407, -2.84095,
                            2.13100, 1.01672, 0.89057, 0.72767, 2.13528, 1.81048, 0.09652, -1.69665, 0.65844, 2.41809,
                            -1.36448, -1.25854, -2.11551, -1.28757, -0.19859, -0.74596, 1.46630, 0.33067, 0.05223,
                            2.26652, -2.88458, -0.69098, 2.84389, -2.86971, -2.53665, -0.02970, -0.18473, -1.77616,
                            -2.61001, 0.83981, -1.92443, 2.53637, 2.20789, 1.56268, 1.66118, -0.37489, 2.34612,
                            -1.70441, -1.41267, -2.70352, 3.09487, 2.36502, -1.29958, -2.39139, 2.75134, 3.08195,
                            -1.86908, -3.07748, -1.66036, -0.41642, -1.87230, 0.17552, -1.56418, -2.35762, -0.43139,
                            2.73425, -3.04271, 1.29312, -2.20116, -0.94523, -0.54829, -2.07539, -2.35920, 2.19664,
                            -2.53442, -0.93800, 1.44565, -2.56374, 1.98230, -2.86949, -2.49495, 0.60522, -2.28480,
                            -1.79339, 3.03862, -0.88420, 1.74822, -0.17979, -1.35466, 0.02520, -2.67470, -2.07880,
                            2.44546, -0.09106, 1.59381, -0.96150, 0.35742, 2.88014, -1.27507, -0.47398, -1.22673,
                            -1.68050, 0.28014, -0.06485, 0.12599, -2.91194, 2.59218, 2.54377, -1.74164, 1.50924,
                            -2.43269, -0.71455, 0.87536, 0.36279, 1.66694, 0.40870, -1.37135, 1.77564, 1.19161,
                            -1.25992, -0.38115, 2.76202, -0.64145, -1.07195, -0.24300, 1.24704, -0.52678, -1.53777,
                            -1.16209, -1.85578, -3.09099, 2.11460, -2.78222, -3.08001, 1.00542, 0.48727, -2.53355,
                            -2.48320, -0.99121, 0.44693, 2.86058, 1.24224, 2.64067, -2.42872, -0.09577, -2.29889,
                            -2.03724, -2.17104, 0.47927, 3.03569, -1.72780, -0.07751, 1.55025, -0.16935, -0.48073,
                            -0.83875, 2.87060, -0.95392, -0.78647, 3.12557, -2.23086, 1.89060, -0.47157, -2.97527,
                            -2.43154, 0.73102, 0.59786, 2.51815, 2.03577, -2.50672, -1.48970, 2.37404, 0.10188,
                            -2.23695, -3.05172, 1.71418, 1.84550, -0.47959, 2.84604, -1.03134, 0.83133, 1.51637,
                            0.72312, -1.67281, -2.23513, -1.12363, -2.60774, -2.75475, 0.44491, -0.44406, 1.81590,
                            1.25662, -1.44510, 1.40271, -2.82312, -2.95329, -2.08122, 2.66993, 0.68513, -0.83546,
                            -2.98846, -2.51786, 1.83460, 1.03644, 1.31368, 3.07523, -1.23877, -0.25951, -2.75627,
                            2.30077, 2.84525, 2.81171, -1.85460, 0.68757, -0.45970, -1.62156, -3.08966, 1.68482,
                            -2.47231, -2.08778, -1.17875, 0.58635, -3.05420, 2.53887, 2.19525, 0.38194, 2.05804,
                            -1.31547, 2.76878, 2.30501, 1.27691, -1.48389, 3.05643, -0.00621, 1.01591, 2.19113, 0.82058,
                            0.67704, 2.53850, 0.13009, 0.04800, 0.73761, 2.34690, 0.01942, 1.87070, 0.51564, 2.91166,
                            -0.93513, 0.54447, -1.31699, -1.04136, 2.48883, 3.10615, -3.13205, 1.83197, -0.63320,
                            1.76861, -0.27807, 0.55417, 2.13834, -0.29319, 2.81429, 2.42045, -0.15453, -0.63348,
                            0.69694, -1.64308, 1.68140, 1.17120, -1.99091, 1.10963, 2.20759, -2.85056, -2.47543,
                            3.07784, -2.19641, -0.13718, -2.32498, -1.89471, 3.05367, -2.12697, -2.61817, 2.99187,
                            -0.93662, 0.35261, 2.20407, 1.36422, 0.05334, -1.15646, 2.61469, -0.02390, -1.74802,
                            2.09471, 2.48034, -2.67864, -2.89619, -2.74894, 0.97649, 0.99258, 1.21938, 2.96121,
                            -1.28247, -0.35908, 0.15602, -0.76313, -2.69766, -1.19959, -0.00175, 2.78713, 2.08171,
                            -0.62763, 2.04446, -0.49511, -1.17784, -2.15424, 1.80462, 0.45639, 0.31280, -3.11294,
                            -1.35746, 1.15297, 2.37797, -0.43380, -2.11951, -2.62371, 2.28956, -1.34802, 1.32173,
                            0.85083, 1.85287, 0.37455, -2.06590, -1.99114, -2.90804, 2.57062, 0.06647, -1.59878,
                            -2.69219, 1.03864, 0.19101, -0.01343, -1.55901, -1.22445, 0.06973, -1.27328, -0.84423,
                            -0.28907, -2.28840, -0.42904, -1.17445, 0.55826, -1.86550, -0.81152, 0.94297, -2.43928,
                            2.73043, 0.43485, 1.73703, 1.02383, 2.00566, -1.93909, 1.31345, 0.21445, 1.82756, -0.92249,
                            1.27349, 1.78288, -0.69996, 1.48529, 0.34694, -2.96057, 0.53774, 2.47748, 0.30927, -2.49610,
                            1.06931, 0.69663, -2.77373, 2.55121, -2.64046, -1.68346, 2.24356, 0.62434, 0.50341, 0.77080,
                            2.79216, 2.56958, -1.84312, -1.93601, 2.07832, 1.70025, -0.30175, 1.91275, -1.47344,
                            1.46447, 0.27533, 2.14766, 1.47790, 2.62504, -2.87978, -0.69158, 0.19543, 1.97901,
                            -0.49932, -2.16049, -2.98013, -1.24106, -1.07173, -2.63159, -0.44106, 0.21118, -0.33968,
                            -2.63518, 2.08051, 0.17271, 0.89791, 1.44897, 1.74611, -0.80467, -1.04297, -1.75427,
                            0.19448, -2.86978, -0.37598, 1.29803, -1.03437, -1.17823, 0.64045, 1.12886, 2.33493,
                            -1.92424, -0.64312, -2.56065, -1.35198, 1.18736, 0.10395, 2.98741, -2.04651, -1.05485,
                            2.48144, 0.41553, 2.75063, -2.41112, -0.19789, -2.09808, 3.12485, 0.34573, -2.22634,
                            -2.81426, -2.37454, 3.12539, -2.53899, -1.19324, 2.92840, 0.47964, -2.53061, 0.31097,
                            0.83954, -0.80138, -1.03792, 0.51992, -2.55064, 1.55128, -1.57775, -1.19823, 0.15297,
                            -2.38643, 0.47734, -2.24217, -2.81552, -1.58253, 0.18684, 2.97040, 2.78973, -0.56667,
                            -0.50373, -2.76300, -0.85822, -0.19070, -0.11773, 2.65887, 0.29674, -1.20039, -2.18116,
                            -2.58898, -2.46817, 2.11521, -2.22228, -1.12308, 1.52800, -2.70596, 1.71603, 2.60767,
                            1.12966, -0.57960, 1.56008, -1.91953, -1.12685, -2.91488, -0.82677, 2.29515, -2.67772,
                            -1.23257, 0.32259, 2.39709, -0.37969, -2.38664, -0.11685, -1.63726, 2.61409, -2.73499,
                            0.39288, 0.95628, -1.41029, 2.05420, -2.44205, 1.46135, -2.95291, 0.43444, 1.24553,
                            -1.87879, -2.53255, 2.80932, 2.28554, 0.12711, 2.59821, -1.82030, -1.12240, -0.39970,
                            1.32370, 2.49092, 2.82877, -0.24857, 0.14047, 1.52383, 2.47708, 0.08759, -1.87276, -1.44702,
                            -0.91764, -2.88249, 0.26762, -2.53352, 1.57858, -1.13262, -2.25379, -0.71158, -2.30312,
                            0.70888, 0.90930, 2.31687, -0.42580, 2.34600, 1.82194, -2.40389, 2.30429, -2.00021, 0.67602,
                            -2.98567, -2.24551, -0.88780, 1.76980, 2.37586, -2.66286, 0.36882, 1.16104, 0.37849,
                            -1.58317, -2.33493, -1.31233, -2.96242, -1.13997, -1.84332, -1.26202, 2.57051, 2.44813,
                            -1.90624, 2.21008, -1.92623, -1.19501, -2.17406, 1.39127, 2.96543, 2.54791, 2.37220,
                            2.04891, -2.25405, -3.10350, -2.50543, -0.29790, -2.57823, 2.69895, -0.52789, -0.40107,
                            -1.73099, 2.99384, -2.69949, -0.14827, 2.26248, -1.27668, -0.73729, -2.47530, -2.26229,
                            -0.82134, 3.03771, 0.21818, 0.73202, -2.97780, 2.75373, -1.74176, -1.91304, -1.05452,
                            0.36436, 0.24049, 2.02641, 1.18546, 1.40324, 2.68570, 1.65134, 0.30198, -2.67718, 1.36649,
                            0.58377, -2.43663, -3.09451, 2.28995, -2.44309, 2.10250, 1.63199, 1.71574, 2.01835,
                            -2.43199, 1.45209, 0.39788, 0.30147, 1.75722, -0.13659, -0.10783, 0.06229, -2.44538,
                            1.44817, 1.68584, -2.56047, 3.04902, 1.97895, 2.96370, 2.76029, -2.37310, -2.49095,
                            -0.06960, -0.85387, 1.16268, 1.53600, -0.12371, 2.88537, -1.86084, 2.03242, -1.79394,
                            3.11615, -2.27119, 1.56889, -1.03037, -2.54125, -1.37225, -2.26581, -1.72696, -1.21156,
                            1.52030, 1.43661, -2.84004, -3.03293, 0.88135, -1.81340, -0.86574, 2.39843, 0.08795,
                            -0.55270, 0.42344, 0.02875, 2.48646, -3.10335, 0.90102, 1.18978, -2.88776, 2.14599, 0.76452,
                            -1.46900, -1.25653, -1.50358, -0.70193, 1.64639, -1.13185, -2.56013, 2.67554, 1.10358,
                            2.75207, 2.12486, -2.13852, -0.16209, 0.31660, -0.99960, 0.29775, 1.02614, 3.06464,
                            -0.11420, 0.82198, 2.95167, 1.76196, 1.89234, -2.70912, -2.73842, -0.73353, 2.39541,
                            -1.13518, 2.55999, -2.81618, -1.10073, 0.66192, 3.01326, -0.01354, 2.65666, -2.40225,
                            2.48566, 2.79502, 0.77968, -1.36742, -2.30674, 1.76693, -2.41793, -0.37606, -2.60649,
                            1.76364, 0.62243, -2.61839, 0.88229, -0.79452, 1.12988, 1.48016, 1.39267, -0.42901,
                            -2.54922, -2.84138, 3.09998, 1.86872, -0.43229, -0.45152, -2.27467, 0.76977, -3.12753,
                            -1.01664, -1.84180, 2.26142, 1.76876, 1.03932, 1.40369, 0.19098, -2.17759, 1.59644,
                            -2.04547, 1.04845, -0.97393, -2.84350, -1.11257, 0.87445, 0.63051, 0.25802, -2.07518,
                            -0.49396, -1.40497, 1.36673, -1.98342, -2.00055, -0.34673, 2.74126, 1.58725, 2.20844,
                            -1.36463, -1.97730, 0.86022, 0.31392, -2.76651, -2.38860, -2.26204, 2.82309, -0.04134,
                            3.02479, -3.05624, -1.81008, 1.92298, -0.47135, -2.56661, 2.27788, -0.74818, -1.63017,
                            -1.64943, -2.25413, 1.60952, -2.46003, 3.12651, -1.87918, 2.01741, -2.37515, 0.30874,
                            0.53719, -2.08100, -2.64575, 0.29875, -2.58131, -2.36022, 0.98075, 0.94413, 0.54718,
                            0.42442, 0.21579, 2.82013, 2.15529, -1.61807, -0.62961, -3.02648, 0.75909, -1.28681,
                            -1.62035, -2.76664, -1.32404, -1.61353, 1.58994, 1.46776, -0.74882, 1.42237, -1.55085,
                            1.83755, 1.97871, -1.51677, -0.17833, -1.52515, -2.17454, -1.33803, -0.19802, 1.52326,
                            0.64931, 0.41097, 3.02335, -0.27681, 2.91757, 3.09209, -1.41274, -1.41490, -1.25075,
                            -1.37539};

double par[MAXCIRCDEPTH] = {1.46487, 0.67381, 2.79234, 2.91183};



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

        double* diagObs = pauliObsMatDiag(comps, coeffs, length, qubits);
        obsDiag.type = DIAG;
        obsDiag.diagObs = diagObs;

        pauliObs.comps = comps;
        pauliObs.coeffs = coeffs;
        pauliObs.length = length;
        pauliObs.qubits = qubits;
        obsPauli.type = PAULI;
        obsPauli.pauliObs = &pauliObs;

        /* Compute the observable's matrix representation */
        cplx_t* obsMatrix = pauliObsMat(comps, coeffs, length, qubits);

        for (dim_t i = 0; i < dim + 1; ++i) {
            stateInitVector(&testState, testVectors[i], qubits);    // Initialize test state
            double result_diag = expValObs(&testState, &obsDiag);   // Apply test functions
            double result_pauli = expValObs(&testState, &obsPauli);

            cplx_t* refKet = cmatVecMul(obsMatrix, refVectors[i], dim);
            double reference = creal(cInner(refVectors[i], refKet, dim));

            TEST_ASSERT_TRUE(ralmostEqual(result_diag, reference, PRECISION));
            TEST_ASSERT_TRUE(ralmostEqual(result_pauli, reference, PRECISION));

            free(refKet);
        }
        /* Pauli observable */
        free(comps);
        comps = allPauliStrings(qubits);
        length = 1 << (2 * qubits);

        free(diagObs);

        pauliObs.comps = comps;
        pauliObs.length = length;

        free(obsMatrix);
        obsMatrix = pauliObsMat(comps, coeffs, length, qubits);

        for (dim_t i = 0; i < dim + 1; ++i) {
            stateInitVector(&testState, testVectors[i], qubits);    // Initialize test state
            double result = expValObs(&testState, &obsPauli);   // Apply test function

            cplx_t* refKet = cmatVecMul(obsMatrix, refVectors[i], dim);
            double reference = creal(cInner(refVectors[i], refKet, dim));

            TEST_ASSERT_TRUE(ralmostEqual(result, reference, PRECISION));

            free(refKet);
        }
        freeTestVectors(testVectors, qubits);
        freeTestVectors(refVectors, qubits);
        free(comps);
        free(obsMatrix);
    }
}


/* 
 * =====================================================================================================================
 *                                              test expValObsPQC
 * =====================================================================================================================
 */

/*
 * This function checks whether calculating the mean value of an observable with respect to a PQC using expValObsPQC
 * coincides with calculating it using Kronecker matrix-vector products on 1 to MAXQUBITS qubits.
 * It checks all possible combinations of types for the observables and for the evolution operators
 * Input:
 *      -
 * Output:
 *      There is no return value. Only checks whether the test and the reference results are equal to some precision for
 *      an observable which is a linear combination of all possible Pauli strings.
*/
void testExpValObsPQC(void) {
    state_t testState;          // State equipped with the test statevectors
    obs_t obsDiag;              // Declare the diagonal observable
    obs_t obsPauli;             // Declare the Pauli observable
    pauliObs_t pauliObs;        // Pauli observable struct

    for (qubit_t qubits = 1; qubits <= MAXQUBITS; ++qubits) {

        dim_t dim = POW2(qubits, dim_t);                        // Hilbert space dimension
        cplx_t** testVectors = generateTestVectors(qubits);     // Statevectors altered by the test function
        cplx_t** refVectors = generateTestVectors(qubits);      // Statevectors altered by matrix multiplication
                                                                // (reference)

        pauli_t* comps_diag = allPauliStringsDiag(qubits);      // Concatenation of all diagonal Pauli strings
        complength_t length_diag = (1 << qubits);               // Number of components, i.e., number of Pauli strings
        pauli_t* comps_pauli = allPauliStrings(qubits);         // Concatenation of all diagonal Pauli strings
        complength_t length_pauli = (1 << (2 * qubits));        // Number of Pauli components

        double* diagObs = pauliObsMatDiag(comps_diag, coeffs, length_diag, qubits);
        obsDiag.type = DIAG;                                    // Initialize the diagonal observable with all possible
        obsDiag.diagObs = diagObs;                              // diagonal Pauli strings and random coefficients

        pauliObs.comps = comps_pauli;                      // Initialize the Pauli observable with all possible
        pauliObs.coeffs = coeffs;                         // Pauli strings and random coefficients
        pauliObs.length = length_pauli;
        pauliObs.qubits = qubits;
        obsPauli.type = PAULI;
        obsPauli.pauliObs = &pauliObs;

        for (depth_t circdepth = 1; circdepth <= MAXCIRCDEPTH; ++circdepth) {
            obs_t* evoOps[circdepth];                           // Declare the pointer to the evolution operators

            /* Initialize the evolution operators to diagonal operators */
            for (depth_t i = 0; i < circdepth; ++i) {
                evoOps[i] = malloc(sizeof(*(evoOps[i])));
                double* evoOpsDiag = pauliObsMatDiag(comps_diag, coeffs + (i * length_diag), length_diag, qubits);
                evoOps[i]->type = DIAG;
                evoOps[i]->diagObs = evoOpsDiag;
            }

            for (dim_t i = 0; i < dim + 1; ++i) {
                stateInitVector(&testState, testVectors[i], qubits);    // Initialize test state

                /* Test diagonal observable and diagonal evolution operators */
                double result = expValObsPQC(&testState, par, &obsDiag, (const obs_t**) evoOps, circdepth);
                double reference = finiteExpValPQC(refVectors[i],
                                                   qubits,
                                                   dim,
                                                   circdepth,
                                                   comps_diag,
                                                   coeffs,
                                                   length_diag,
                                                   comps_diag,
                                                   coeffs,
                                                   length_diag,
                                                   par);
                TEST_ASSERT_TRUE(ralmostEqual(result, reference, PRECISION));

                /* Test Pauli observable and diagonal evolution operators */
                result = expValObsPQC(&testState, par, &obsPauli, (const obs_t**) evoOps, circdepth);
                reference = finiteExpValPQC(refVectors[i],
                                                   qubits,
                                                   dim,
                                                   circdepth,
                                                   comps_pauli,
                                                   coeffs,
                                                   length_pauli,
                                                   comps_diag,
                                                   coeffs,
                                                   length_diag,
                                                   par);
                TEST_ASSERT_TRUE(ralmostEqual(result, reference, PRECISION));
            }

            /* Initialize the evolution operators to Pauli operators */
            pauliObs_t* evoOpsPauli[circdepth];
            for (depth_t i = 0; i < circdepth; ++i) {
                free(evoOps[i]->diagObs);                                  // free diagonal observable array
                evoOpsPauli[i] = malloc(sizeof(*(evoOpsPauli[i])));

                evoOpsPauli[i]->comps = comps_pauli;
                evoOpsPauli[i]->coeffs = coeffs + (i * length_pauli);
                evoOpsPauli[i]->length = length_pauli;
                evoOpsPauli[i]->qubits = qubits;

                evoOps[i]->type = PAULI;
                evoOps[i]->pauliObs = evoOpsPauli[i];
            }

            for (dim_t i = 0; i < dim + 1; ++i) {
                stateInitVector(&testState, testVectors[i], qubits);    // Initialize test state
                /* Test diagonal observable and Pauli evolution operators */
                double result = expValObsPQC(&testState, par, &obsDiag, (const obs_t**) evoOps, circdepth);
                double reference = finiteExpValPQC(refVectors[i],
                                                   qubits,
                                                   dim,
                                                   circdepth,
                                                   comps_diag,
                                                   coeffs,
                                                   length_diag,
                                                   comps_pauli,
                                                   coeffs,
                                                   length_pauli,
                                                   par);
                TEST_ASSERT_TRUE(ralmostEqual(result, reference, PRECISION));

                /* Test Pauli observable and diagonal Pauli operators */
                result = expValObsPQC(&testState, par, &obsPauli, (const obs_t**) evoOps, circdepth);
                reference = finiteExpValPQC(refVectors[i],
                                            qubits,
                                            dim,
                                            circdepth,
                                            comps_pauli,
                                            coeffs,
                                            length_pauli,
                                            comps_pauli,
                                            coeffs,
                                            length_pauli,
                                            par);
                TEST_ASSERT_TRUE(ralmostEqual(result, reference, PRECISION));
            }

            for (depth_t i = 0; i < circdepth; ++i) {
                free(evoOpsPauli[i]);
                free(evoOps[i]);
            }
        }
        freeTestVectors(testVectors, qubits);
        freeTestVectors(refVectors, qubits);
        free(comps_diag);
        free(comps_pauli);
        free(diagObs);
    }
}


/* 
 * =====================================================================================================================
 *                                                      test gradientPQC
 * =====================================================================================================================
 */

/*
 * This function checks whether calculating the gradient of an observable's mean value with respect to a PQC using
 * gradientPQC and approximating it using approxGradientPQC coincide with calculating it using Kronecker matrix-vector
 * products on 1 to MAXQUBITS qubits.
 * It checks all possible combinations of types for the observables and for the evolution operators
 * Input:
 *      -
 * Output:
 *      There is no return value. Only checks whether the test and the reference results are equal to some precision for
 *      an observable which is a linear combination of all possible Pauli strings.
*/
void testGradientPQC(void) {
    state_t testState;          // State equipped with the test statevectors
    obs_t obsDiag;              // Declare the diagonal observable
    obs_t obsPauli;             // Declare the Pauli observable
    pauliObs_t pauliObs;        // Pauli observable struct

    for (qubit_t qubits = 1; qubits <= MAXQUBITS; ++qubits) {

        dim_t dim = POW2(qubits, dim_t);                        // Hilbert space dimension
        cplx_t** testVectors = generateTestVectors(qubits);     // Statevectors altered by the test function
        cplx_t** refVectors = generateTestVectors(qubits);      // Statevectors altered by matrix multiplication
                                                                // (reference)

        pauli_t* comps_diag = allPauliStringsDiag(qubits);      // Concatenation of all diagonal Pauli strings
        complength_t length_diag = (1 << qubits);               // Number of components, i.e., number of Pauli strings
        pauli_t* comps_pauli = allPauliStrings(qubits);         // Concatenation of all diagonal Pauli strings
        complength_t length_pauli = (1 << (2 * qubits));        // Number of Pauli components
        pauli_t* compsY = allPauliStringsY(qubits);

        double* diagObs = pauliObsMatDiag(comps_diag, coeffs, length_diag, qubits);
        obsDiag.type = DIAG;                                    // Initialize the diagonal observable with all possible
        obsDiag.diagObs = diagObs;                              // diagonal Pauli strings and random coefficients

        pauliObs.comps = comps_pauli;                      // Initialize the Pauli observable with all possible
        pauliObs.coeffs = coeffs;                         // Pauli strings and random coefficients
        pauliObs.length = length_pauli;
        pauliObs.qubits = qubits;
        obsPauli.type = PAULI;
        obsPauli.pauliObs = &pauliObs;

        for (depth_t circdepth = 1; circdepth <= MAXCIRCDEPTH; ++circdepth) {
            obs_t* evoOps[circdepth];                           // Declare the pointer to the evolution operators
            double* evoOpsDiag[circdepth];

            /* Initialize the evolution operators to diagonal operators */
            for (depth_t i = 0; i < circdepth; ++i) {
                evoOps[i] = malloc(sizeof(*(evoOps[i])));
                evoOpsDiag[i] = pauliObsMatDiag(comps_diag, coeffs + (i * length_diag), length_diag, qubits);
                evoOps[i]->type = DIAG;
                evoOps[i]->diagObs = evoOpsDiag[i];
            }

            for (dim_t i = 0; i < dim + 1; ++i) {
                stateInitVector(&testState, testVectors[i], qubits);    // Initialize test state

                /* Test diagonal observable and diagonal evolution operators */
                double* result = gradPQC(&testState, par, &obsDiag, (const obs_t **) evoOps, circdepth);
                double* result2 = approxGradientPQC(&testState,     // Test gradient using finite difference
                                                    par,
                                                    &obsDiag,
                                                    (const obs_t**) evoOps,
                                                    circdepth,
                                                    EPSILON);
                double* reference = finiteGradientPQC(refVectors[i],
                                                      qubits,
                                                      dim,
                                                      circdepth,
                                                      comps_diag,
                                                      coeffs,
                                                      length_diag,
                                                      comps_diag,
                                                      coeffs,
                                                      length_diag,
                                                      par,
                                                      EPSILON);
//                printf("Diag-diag\n");
//                printf("Result: ");
//                vectorPrint(result, circdepth);
//                printf("Reference: ");
//                vectorPrint(reference, circdepth);
                TEST_ASSERT_TRUE(rvectorAlmostEqual(result, reference, circdepth, APPROXPRECISION));
                TEST_ASSERT_TRUE(rvectorAlmostEqual(result2, reference, circdepth, APPROXPRECISION));
                free(result);
                free(result2);
                free(reference);

                /* Test Pauli observable and diagonal evolution operators */
                result = gradPQC(&testState, par, &obsPauli, (const obs_t **) evoOps, circdepth);
                result2 = approxGradientPQC(&testState,     // Test gradient using finite difference
                                                    par,
                                                    &obsPauli,
                                                    (const obs_t**) evoOps,
                                                    circdepth,
                                                    EPSILON);
                reference = finiteGradientPQC(refVectors[i],
                                                      qubits,
                                                      dim,
                                                      circdepth,
                                                      comps_pauli,
                                                      coeffs,
                                                      length_pauli,
                                                      comps_diag,
                                                      coeffs,
                                                      length_diag,
                                                      par,
                                                      EPSILON);
//                printf("Diag-pauli\n");
//                printf("Result: ");
//                vectorPrint(result, circdepth);
//                printf("Reference: ");
//                vectorPrint(reference, circdepth);
                TEST_ASSERT_TRUE(rvectorAlmostEqual(result, reference, circdepth, APPROXPRECISION));
                TEST_ASSERT_TRUE(rvectorAlmostEqual(result2, reference, circdepth, APPROXPRECISION));
                free(result);
                free(result2);
                free(reference);
            }

            /* Initialize the evolution operators to Pauli operators */
            pauliObs_t* evoOpsPauli[circdepth];
            for (depth_t i = 0; i < circdepth; ++i) {
                free(evoOps[i]->diagObs);                                   // free diagonal observable array
                evoOpsPauli[i] = malloc(sizeof(*(evoOpsPauli[i])));

                evoOpsPauli[i]->comps = compsY;
                evoOpsPauli[i]->coeffs = coeffs + (i * length_diag);
                evoOpsPauli[i]->length = length_diag;
                evoOpsPauli[i]->qubits = qubits;

                evoOps[i]->type = PAULI;
                evoOps[i]->pauliObs = evoOpsPauli[i];
            }

            for (dim_t i = 0; i < dim + 1; ++i) {
                stateInitVector(&testState, testVectors[i], qubits);    // Initialize test state
                /* Test diagonal observable and Pauli evolution operators */
                double* result = gradPQC(&testState, par, &obsDiag, (const obs_t **) evoOps, circdepth);
                double* result2 = approxGradientPQC(&testState,     // Test gradient using finite difference
                                            par,
                                            &obsDiag,
                                            (const obs_t**) evoOps,
                                            circdepth,
                                            EPSILON);
                double* reference = finiteGradientPQC(refVectors[i],
                                                      qubits,
                                                      dim,
                                                      circdepth,
                                                      comps_diag,
                                                      coeffs,
                                                      length_diag,
                                                      compsY,
                                                      coeffs,
                                                      length_diag,
                                                      par,
                                                      EPSILON);
//                printf("Pauli-diag\n");
//                printf("Result: ");
//                vectorPrint(result, circdepth);
//                printf("Reference: ");
//                vectorPrint(reference, circdepth);
                TEST_ASSERT_TRUE(rvectorAlmostEqual(result, reference, circdepth, APPROXPRECISION));
                TEST_ASSERT_TRUE(rvectorAlmostEqual(result2, reference, circdepth, APPROXPRECISION));
                free(result);
                free(result2);
                free(reference);

                /* Test Pauli observable and diagonal Pauli operators */
                result = gradPQC(&testState, par, &obsPauli, (const obs_t **) evoOps, circdepth);
                result2 = approxGradientPQC(&testState,     // Test gradient using finite difference
                                            par,
                                            &obsPauli,
                                            (const obs_t**) evoOps,
                                            circdepth,
                                            EPSILON);
                reference = finiteGradientPQC(refVectors[i],
                                              qubits,
                                              dim,
                                              circdepth,
                                              comps_pauli,
                                              coeffs,
                                              length_pauli,
                                              compsY,
                                              coeffs,
                                              length_diag,
                                              par,
                                              EPSILON);
//                printf("Pauli-pauli\n");
//                printf("Result: ");
//                vectorPrint(result, circdepth);
//                printf("Reference: ");
//                vectorPrint(reference, circdepth);
                TEST_ASSERT_TRUE(rvectorAlmostEqual(result, reference, circdepth, APPROXPRECISION));
                TEST_ASSERT_TRUE(rvectorAlmostEqual(result2, reference, circdepth, APPROXPRECISION));
                free(result);
                free(result2);
                free(reference);
            }

            for (depth_t i = 0; i < circdepth; ++i) {
                free(evoOpsPauli[i]);
                free(evoOps[i]);
            }
        }
        freeTestVectors(testVectors, qubits);
        freeTestVectors(refVectors, qubits);
        free(comps_diag);
        free(comps_pauli);
        free(compsY);
        free(diagObs);
    }
}

/*
 * =====================================================================================================================
 *                                                  test hessianPQC
 * =====================================================================================================================
 */

double* grad(const state_t* state,
             const double params[],
             const obs_t* observable,
             const obs_t* evoOps[],
             depth_t circdepth)
{
    return gradPQC(state, params, observable, evoOps, circdepth);
}

/*
 * This function checks whether calculating the hessian of an observable's mean value with respect to a PQC using
 * hessianPQC and approximating it using approxHessianPQC coincide with calculating it using Kronecker matrix-vector
 * products on 1 to MAXQUBITS qubits.
 * It checks all possible combinations of types for the observables and for the evolution operators
 * Input:
 *      -
 * Output:
 *      There is no return value. Only checks whether the test and the reference results are equal to some precision for
 *      an observable which is a linear combination of all possible Pauli strings.
*/
void testHessianPQC(void) {
    state_t testState;          // State equipped with the test statevectors
    obs_t obsDiag;              // Declare the diagonal observable
    obs_t obsPauli;             // Declare the Pauli observable
    pauliObs_t pauliObs;        // Pauli observable struct

    for (qubit_t qubits = 1; qubits <= MAXQUBITS; ++qubits) {

        dim_t dim = POW2(qubits, dim_t);                        // Hilbert space dimension
        cplx_t** testVectors = generateTestVectors(qubits);     // Statevectors altered by the test function
        cplx_t** refVectors = generateTestVectors(qubits);      // Statevectors altered by matrix multiplication
                                                                // (reference)
        double precision = pow(10, qubits - 3);                 // Precision for comparing hessians

        pauli_t* comps_diag = allPauliStringsDiag(qubits);      // Concatenation of all diagonal Pauli strings
        complength_t length_diag = (1 << qubits);               // Number of components, i.e., number of Pauli strings
        pauli_t* comps_pauli = allPauliStrings(qubits);         // Concatenation of all diagonal Pauli strings
        complength_t length_pauli = (1 << (2 * qubits));        // Number of Pauli components
        pauli_t* compsY = allPauliStringsY(qubits);

        double* diagObs = pauliObsMatDiag(comps_diag, coeffs, length_diag, qubits);
        obsDiag.type = DIAG;                                    // Initialize the diagonal observable with all possible
        obsDiag.diagObs = diagObs;                              // diagonal Pauli strings and random coefficients

        pauliObs.comps = comps_pauli;                      // Initialize the Pauli observable with all possible
        pauliObs.coeffs = coeffs;                         // Pauli strings and random coefficients
        pauliObs.length = length_pauli;
        pauliObs.qubits = qubits;
        obsPauli.type = PAULI;
        obsPauli.pauliObs = &pauliObs;

        for (depth_t circdepth = 1; circdepth <= MAXCIRCDEPTH; ++circdepth) {
            obs_t* evoOps[circdepth];                           // Declare the pointer to the evolution operators
            double* evoOpsDiag[circdepth];

            /* Initialize the evolution operators to diagonal operators */
            for (depth_t i = 0; i < circdepth; ++i) {
                evoOps[i] = malloc(sizeof(*(evoOps[i])));
                evoOpsDiag[i] = pauliObsMatDiag(comps_diag, coeffs + (i * length_diag), length_diag, qubits);
                evoOps[i]->type = DIAG;
                evoOps[i]->diagObs = evoOpsDiag[i];
            }

            for (dim_t i = 0; i < dim + 1; ++i) {
                stateInitVector(&testState, testVectors[i], qubits);    // Initialize test state

                /* Test diagonal observable and diagonal evolution operators */
                double* result = hessianPQC(&testState, par, &obsDiag, (const obs_t**) evoOps, circdepth);
                double* result2 = approxHessianPQC(&testState,     // Test gradient using finite difference
                                                    par,
                                                    &obsDiag,
                                                    (const obs_t**) evoOps,
                                                    circdepth,
                                                    grad,
                                                    DELTA);
                double* reference = finiteHessianPQC(refVectors[i],
                                                      qubits,
                                                      dim,
                                                      circdepth,
                                                      comps_diag,
                                                      coeffs,
                                                      length_diag,
                                                      comps_diag,
                                                      coeffs,
                                                      length_diag,
                                                      par,
                                                      DELTA);
                TEST_ASSERT_TRUE(rvectorAlmostEqual(result, reference, circdepth, precision));
                TEST_ASSERT_TRUE(rvectorAlmostEqual(result2, reference, circdepth, precision));
                free(result);
                free(result2);
                free(reference);

                /* Test Pauli observable and diagonal evolution operators */
                result = hessianPQC(&testState, par, &obsPauli, (const obs_t**) evoOps, circdepth);
                result2 = approxHessianPQC(&testState,     // Test gradient using finite difference
                                            par,
                                            &obsPauli,
                                            (const obs_t**) evoOps,
                                            circdepth,
                                            grad,
                                            DELTA);
                reference = finiteHessianPQC(refVectors[i],
                                              qubits,
                                              dim,
                                              circdepth,
                                              comps_pauli,
                                              coeffs,
                                              length_pauli,
                                              comps_diag,
                                              coeffs,
                                              length_diag,
                                              par,
                                              DELTA);
                TEST_ASSERT_TRUE(rvectorAlmostEqual(result, reference, circdepth, precision));
                TEST_ASSERT_TRUE(rvectorAlmostEqual(result2, reference, circdepth, precision));
                free(result);
                free(result2);
                free(reference);
            }

            /* Initialize the evolution operators to Pauli operators */
            pauliObs_t* evoOpsPauli[circdepth];
            for (depth_t i = 0; i < circdepth; ++i) {
                free(evoOps[i]->diagObs);                                   // free diagonal observable array
                evoOpsPauli[i] = malloc(sizeof(*(evoOpsPauli[i])));

                evoOpsPauli[i]->comps = compsY;
                evoOpsPauli[i]->coeffs = coeffs + (i * length_diag);
                evoOpsPauli[i]->length = length_diag;
                evoOpsPauli[i]->qubits = qubits;

                evoOps[i]->type = PAULI;
                evoOps[i]->pauliObs = evoOpsPauli[i];
            }

            for (dim_t i = 0; i < dim + 1; ++i) {
                stateInitVector(&testState, testVectors[i], qubits);    // Initialize test state
                /* Test diagonal observable and Pauli evolution operators */
                double* result = hessianPQC(&testState, par, &obsDiag, (const obs_t**) evoOps, circdepth);
                double* result2 = approxHessianPQC(&testState,     // Test gradient using finite difference
                                                    par,
                                                    &obsDiag,
                                                    (const obs_t**) evoOps,
                                                    circdepth,
                                                    grad,
                                                    DELTA);
                double* reference = finiteHessianPQC(refVectors[i],
                                                      qubits,
                                                      dim,
                                                      circdepth,
                                                      comps_diag,
                                                      coeffs,
                                                      length_diag,
                                                      compsY,
                                                      coeffs,
                                                      length_diag,
                                                      par,
                                                      DELTA);
                TEST_ASSERT_TRUE(rvectorAlmostEqual(result, reference, circdepth, precision));
                TEST_ASSERT_TRUE(rvectorAlmostEqual(result2, reference, circdepth, precision));
                free(result);
                free(result2);
                free(reference);

                /* Test Pauli observable and diagonal Pauli operators */
                result = hessianPQC(&testState, par, &obsPauli, (const obs_t**) evoOps, circdepth);
                result2 = approxHessianPQC(&testState,     // Test gradient using finite difference
                                            par,
                                            &obsPauli,
                                            (const obs_t**) evoOps,
                                            circdepth,
                                            grad,
                                            DELTA);
                reference = finiteHessianPQC(refVectors[i],
                                              qubits,
                                              dim,
                                              circdepth,
                                              comps_pauli,
                                              coeffs,
                                              length_pauli,
                                              compsY,
                                              coeffs,
                                              length_diag,
                                              par,
                                              DELTA);
                TEST_ASSERT_TRUE(rvectorAlmostEqual(result, reference, circdepth, precision));
                TEST_ASSERT_TRUE(rvectorAlmostEqual(result2, reference, circdepth, precision));
                free(result);
                free(result2);
                free(reference);
            }

            for (depth_t i = 0; i < circdepth; ++i) {
                free(evoOpsPauli[i]);
                free(evoOps[i]);
            }
        }
        freeTestVectors(testVectors, qubits);
        freeTestVectors(refVectors, qubits);
        free(comps_diag);
        free(comps_pauli);
        free(compsY);
        free(diagObs);
    }
}

/*
 * =====================================================================================================================
 *                                                  test momentMatQAOA
 * =====================================================================================================================
 */

void applyU(state_t* state, const obs_t* srchOp) {
    return (evolveObsTrotter(state, srchOp, 0.05));
}

void testMomMat(void) {
    state_t testState;          // State equipped with the test statevectors
    depth_t obsc = 5;
    obs_t* obs[obsc];           // Array of observables with one diagonal followed by 4 Pauli observables

    for (qubit_t qubits = 1; qubits <= MAXQUBITS; ++qubits) {

        dim_t dim = POW2(qubits, dim_t);                        // Hilbert space dimension
        cplx_t** testVectors = generateTestVectors(qubits);     // Statevectors altered by the test function

        pauli_t* comps_diag = allPauliStringsDiag(qubits);      // Concatenation of all diagonal Pauli strings
        complength_t length_diag = (1 << qubits);               // Number of components, i.e., number of Pauli strings
        pauli_t* comps_pauli = allPauliStrings(qubits);         // Concatenation of all diagonal Pauli strings
        complength_t length_pauli = (1 << (2 * qubits));        // Number of Pauli components

        double* diagObs = pauliObsMatDiag(comps_diag, coeffs, length_diag, qubits);
        obs[0] = malloc(sizeof(*(obs[0])));
        obs[0]->type = DIAG;                                    // Initialize the diagonal observable with all possible
        obs[0]->diagObs = diagObs;                              // diagonal Pauli strings and random coefficients

        pauliObs_t* obsPauli[obsc];
        for (int i = 0; i < obsc - 1; ++i) {
            obsPauli[i] = malloc(sizeof(*(obsPauli[i])));
            obsPauli[i]->comps = comps_pauli;                      // Initialize 4 Pauli observables with all
            obsPauli[i]->coeffs = coeffs + (i * length_pauli);    // possible Pauli strings and the 4 sets of
            obsPauli[i]->length = length_pauli;                         // coefficients
            obsPauli[i]->qubits = qubits;

            obs[i + 1] = malloc(sizeof(*(obs[i])));
            obs[i + 1]->type = PAULI;
            obs[i + 1]->pauliObs = obsPauli[i];
        }

        cplx_t* obsMat[obsc];                                                   // Array of observable's matrix
                                                                                // representation
        obsMat[0] = pauliObsMat(comps_diag, coeffs, length_diag, qubits);
        for (depth_t i = 0; i < obsc - 1; ++i) {
            obsMat[i + 1] = pauliObsMat(comps_pauli,
                                        coeffs + (i * length_pauli),
                                        length_pauli,
                                        qubits);
        }

        cplx_t* srchOpsMat[obsc];                                               // Array of search operator's matrix
                                                                                // representation
        srchOpsMat[0] = expPauliObsMatTrotter(comps_diag, coeffs, length_diag, 0.05, qubits);
        for (depth_t i = 0; i < obsc - 1; ++i) {
            srchOpsMat[i + 1] = expPauliObsMatTrotter(comps_pauli,
                                                      coeffs + (i * length_pauli),
                                                      length_pauli,
                                                      0.05,
                                                      qubits);
        }

        for (dim_t i = 0; i < dim + 1; ++i) {

            /* Compute the moment matrices using the test function */
            stateInitVector(&testState, testVectors[i], qubits);
            cplx_t** result = momMat(&testState, (const obs_t**) obs, obsc, (const obs_t**) obs, obsc, applyU);

            cplx_t** reference = calloc(obsc, sizeof(cplx_t*));         // Declare and allocate the reference
            for (depth_t j = 0; j < obsc; ++j) {                                    // matrix
                reference[j] = calloc(obsc * obsc, sizeof(cplx_t));
            }

            /* Compute the reference moment matrices in column major form using matrix multiplication */
            for (depth_t j = 0; j < obsc; ++j) {                                            // Iterate columns; apply
                cplx_t* refKet = cmatVecMul(srchOpsMat[j], testState.vec, dim);  // search unitary

                for (depth_t k = 0; k < obsc; ++k) {                                        // Iterate moment matrices;
                    cplx_t* refBra = cmatVecMul(obsMat[k], refKet, dim);            // compute diagonals
                    reference[k][j * obsc + j] = creal(cInner(refBra, refKet, dim));
                    free(refBra);
                }

                for (depth_t k = 0; k < obsc; ++k) {                                            // Iterate rows; apply
                    cplx_t* tmpBra = cmatVecMul(srchOpsMat[k], testState.vec, dim);  // search unitaries
                    for (depth_t l = 0; l < obsc; ++l) {                                    // Iterate moment matrices;
                        cplx_t* refBra = cmatVecMul(obsMat[l], tmpBra, dim);        // compute entries
                        reference[l][j * obsc + k] = cInner(refBra, refKet, dim);
                        reference[l][k * obsc + j] = conj(reference[l][j * obsc + k]);
                        free(refBra);
                    }
                    free(tmpBra);
                }
                free(refKet);
            }

            for (depth_t j = 0; j < obsc; ++j) {
                TEST_ASSERT_TRUE(cvectorAlmostEqual(result[j], reference[j], obsc * obsc, 1e-6));
                free(result[j]);
                free(reference[j]);
            }
            stateFreeVector(&testState);
            free(result);
            free(reference);
        }
        free(obs[0]->diagObs);
        for (depth_t i = 0; i < obsc - 1; ++i) {
            free(obsPauli[i]);
        }
        for (depth_t i = 0; i < obsc; ++i) {
            free(obsMat[i]);
            free(srchOpsMat[i]);
            free(obs[i]);
        }
        free(comps_diag);
        free(comps_pauli);
        free(testVectors);
    }
}

void testMomMatPQC(void) {
    state_t testState;          // State equipped with the test statevectors
    depth_t obsc = 5;
    obs_t* obs[obsc];           // Array of observables with one diagonal followed by 4 Pauli observables

    for (qubit_t qubits = 1; qubits <= MAXQUBITS; ++qubits) {

        dim_t dim = POW2(qubits, dim_t);                        // Hilbert space dimension
        cplx_t** testVectors = generateTestVectors(qubits);     // Statevectors altered by the test function

        pauli_t* comps_diag = allPauliStringsDiag(qubits);      // Concatenation of all diagonal Pauli strings
        complength_t length_diag = (1 << qubits);               // Number of components, i.e., number of Pauli strings
        pauli_t* comps_pauli = allPauliStrings(qubits);         // Concatenation of all Pauli strings
        complength_t length_pauli = (1 << (2 * qubits));        // Number of Pauli components

        double* diagObs = pauliObsMatDiag(comps_diag, coeffs, length_diag, qubits);
        obs[0] = malloc(sizeof(*(obs[0])));
        obs[0]->type = DIAG;                                    // Initialize the diagonal observable with all possible
        obs[0]->diagObs = diagObs;                              // diagonal Pauli strings and random coefficients

        pauliObs_t* obsPauli[obsc];
        for (int i = 0; i < obsc - 1; ++i) {
            obsPauli[i] = malloc(sizeof(*(obsPauli[i])));
            obsPauli[i]->comps = comps_pauli;                      // Initialize 4 Pauli observables with all
            obsPauli[i]->coeffs = coeffs + (i * length_pauli);    // possible Pauli strings and the 4 sets of
            obsPauli[i]->length = length_pauli;                         // coefficients
            obsPauli[i]->qubits = qubits;

            obs[i + 1] = malloc(sizeof(*(obs[i])));
            obs[i + 1]->type = PAULI;
            obs[i + 1]->pauliObs = obsPauli[i];
        }

        cplx_t* obsMat[obsc];                                                       // Array of observable's matrix
                                                                                    // representation
        obsMat[0] = pauliObsMat(comps_diag, coeffs, length_diag, qubits);
        for (depth_t i = 0; i < obsc - 1; ++i) {
            obsMat[i + 1] = pauliObsMat(comps_pauli,
                                        coeffs + (i * length_pauli),
                                        length_pauli,
                                        qubits);
        }

        cplx_t* srchOpsMat[obsc + 1];                                           // Array of search operator's matrix
                                                                                // representation preceded by id
        srchOpsMat[0] = identityMat(qubits);
        srchOpsMat[1] = expPauliObsMatTrotter(comps_diag, coeffs, length_diag, coeffs[0], qubits);
        for (depth_t i = 0; i < obsc - 1; ++i) {
            srchOpsMat[i + 2] = expPauliObsMatTrotter(comps_pauli,
                                                      coeffs + (i * length_pauli),
                                                      length_pauli,
                                                      coeffs[i + 1],
                                                      qubits);
        }

        for (dim_t i = 0; i < dim + 1; ++i) {

            /* Compute the moment matrices using the test function */
            stateInitVector(&testState, testVectors[i], qubits);
            cplx_t** result = momMatPQC(&testState, (const obs_t**) obs, obsc, (const obs_t**) obs, obsc, coeffs);

            cplx_t** reference = calloc(obsc, sizeof(cplx_t*));         // Declare and allocate the reference
            for (depth_t j = 0; j < obsc; ++j) {                                    // matrix
                reference[j] = calloc((obsc + 1) * (obsc + 1), sizeof(cplx_t));
            }

            /* Compute the reference moment matrices in column major form using matrix multiplication */
            for (depth_t j = 0; j < obsc + 1; ++j) {                                       // Iterate columns; apply
                cplx_t* refKet = cmatVecMul(srchOpsMat[j], testState.vec, dim);  // search unitary

                for (depth_t k = 0; k < obsc; ++k) {                                        // Iterate moment matrices;
                    cplx_t* refBra = cmatVecMul(obsMat[k], refKet, dim);            // compute diagonals
                    reference[k][j * (obsc + 1) + j] = creal(cInner(refBra, refKet, dim));
                    free(refBra);
                }

                for (depth_t k = 0; k < obsc + 1; ++k) {                                        // Iterate rows; apply
                    cplx_t* tmpBra = cmatVecMul(srchOpsMat[k], testState.vec, dim);  // search unitaries
                    for (depth_t l = 0; l < obsc; ++l) {                                    // Iterate moment matrices;
                        cplx_t* refBra = cmatVecMul(obsMat[l], tmpBra, dim);        // compute entries
                        reference[l][j * (obsc + 1) + k] = cInner(refBra, refKet, dim);
                        reference[l][k * (obsc + 1) + j] = conj(reference[l][j * (obsc + 1) + k]);
                        free(refBra);
                    }
                    free(tmpBra);
                }
                free(refKet);
            }

            for (depth_t j = 0; j < obsc; ++j) {
                TEST_ASSERT_TRUE(cvectorAlmostEqual(result[j], reference[j], obsc * obsc, 1e-6));
                free(result[j]);
                free(reference[j]);
            }
            stateFreeVector(&testState);
            free(result);
            free(reference);
        }
        free(obs[0]->diagObs);
        for (depth_t i = 0; i < obsc - 1; ++i) {
            free(obsPauli[i]);
        }
        for (depth_t i = 0; i < obsc; ++i) {
            free(obsMat[i]);
            free(srchOpsMat[i]);
            free(obs[i]);
        }
        free(srchOpsMat[obsc]);
        free(comps_diag);
        free(comps_pauli);
        free(testVectors);
    }
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
    RUN_TEST(testHessianPQC);
    RUN_TEST(testMomMat);
    RUN_TEST(testMomMatPQC);
    return UNITY_END();
}
