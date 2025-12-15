// test.h - Functions and objects to check library
#ifndef Q_TEST_H
#define Q_TEST_H

/*
 * =====================================================================================================================
 * includes
 * =====================================================================================================================
 */

#ifndef _COMPLEX_H_
#include <complex.h>
#endif

#ifndef QTYPES_H
#include  "q_types.h"
#endif

#ifndef _STDIO_H_
#include <stdio.h>
#endif

#ifndef _STDLIB_H_
#include <stdlib.h>
#endif

#ifndef UNITY_FRAMEWORK_H
#include "unity.h"
#endif

/*
 * =====================================================================================================================
 * C++ check
 * =====================================================================================================================
 */

#ifdef __cplusplus
extern "C" {
#endif

/*
 * =====================================================================================================================
 * Macros
 * =====================================================================================================================
 */

#define MAXQUBITS   4
#define SQRT2       1.4142135623730951
#define INVSQRT2    0.7071067811865475
#define INVSQRT4    0.5
#define INVSQRT8    0.3535533905932738
#define INVSQRT16   0.25


/*
 * =====================================================================================================================
 * Check functions (Unity extensions)
 * =====================================================================================================================
 */

/*
 * @brief   Check whether all elements in two DOUBLE PRECISION arrays coincide up to tolerance
 */
#define TEST_ASSERT_EQUAL_DOUBLE_ARRAY_TOL(expected, actual, num_elements, tol) \
    do { \
        for (index_t _i = 0; _i < (num_elements); ++_i) { \
            int _fail = fabs((expected)[_i] - (actual)[_i]) > (tol); \
            if (_fail) { \
                char _msg[256]; \
                snprintf(_msg, sizeof(_msg), \
                    "Element %u Expected %.7g, Was %.7g", _i, (expected)[_i], (actual)[_i]); \
                    UNITY_TEST_FAIL(__LINE__, _msg); \
            } \
        } \
    } while (0)


/*
 * @brief   Check whether all elements' absolute values in two DOUBLE PRECISION arrays coincide up to tolerance
 */
#define TEST_ASSERT_EQUAL_DOUBLE_ARRAY_ABS(expected, actual, num_elements) \
    do { \
        for (index_t _i = 0; _i < (num_elements); ++_i) { \
            double _expected = fabs((expected)[_i]); \
            double _actual = fabs((actual)[_i]); \
            int _fail = fabs(_expected - _actual) > (1e-7); \
            if (_fail) { \
                char _msg[256]; \
                snprintf(_msg, sizeof(_msg), \
                    "Element %zu Expected %.7g, Was %.7g"; _i, _expected, _actual); \
                    UNITY_TEST_FAIL(__LINE__, _msg); \
            } \
        } \
    } while (0)


/*
 * @brief   Check whether two complex numbers coincide up to the tolerance
 */
#define TEST_ASSERT_COMPLEX_WITHIN(expected, actual, tol) \
    do { \
        double _re_expected = creal((expected)); \
        double _im_expected = cimag((expected)); \
        double _re_actual = creal((actual)); \
        double _im_actual = cimag((actual)); \
        \
        int _re_fail = fabs(_re_expected - _re_actual) > (tol); \
        int _im_fail = fabs(_im_expected - _im_actual) > (tol); \
        if (_re_fail || _im_fail) { \
            char _msg[256]; \
            snprintf(_msg, sizeof(_msg), \
                " Element mismatch: Expected %.7g + %.7g*I, Was %.7g + %.7g*I", \
                _re_expected, _im_expected, _re_actual, _im_actual); \
                UNITY_TEST_FAIL(__LINE__, _msg); \
        } \
    } while(0)


/*
 * @brief   Check whether all elements in two DOUBLE COMPLEX arrays coincide up to tolerance
 */
#define TEST_ASSERT_EQUAL_COMPLEX_ARRAY_TOL(expected, actual, num_elements, tol) \
    do { \
        for (size_t _i = 0; _i < (num_elements); ++_i) { \
            double _re_expected = creal((expected)[_i]); \
            double _im_expected = cimag((expected)[_i]); \
            double _re_actual   = creal((actual)[_i]); \
            double _im_actual   = cimag((actual)[_i]); \
            \
            int _re_fail = fabs(_re_expected - _re_actual) > (tol); \
            int _im_fail = fabs(_im_expected - _im_actual) > (tol); \
            if (_re_fail || _im_fail) { \
                char _msg[256]; \
                snprintf(_msg, sizeof(_msg), \
                    "Element %zu mismatch: Expected %.7g + %.7g*I, Was %.7g + %.7g*I", \
                    _i, _re_expected, _im_expected, _re_actual, _im_actual); \
                    UNITY_TEST_FAIL(__LINE__, _msg); \
            } \
        } \
    } while (0)

/*
 * =====================================================================================================================
 * Pure quantum test states
 * =====================================================================================================================
 */

// 1-qubit pure states (size = 2 amplitudes)
const static cplx_t state_pure_one_qubit_zero[2] = {1.0 + I*0.0, 0.0 + I*0.0};
const static cplx_t state_pure_one_qubit_one[2] = {0.0 + I*0.0, 0.0 + I*1.0};
const static cplx_t state_pure_one_qubit_plus[2] = {INVSQRT2 + I*0.0, INVSQRT2 + I*0.0};
const static cplx_t *states_pure_one_qubit[3] = {
    state_pure_one_qubit_zero, state_pure_one_qubit_one, state_pure_one_qubit_zero
};


// 2-qubit pure states (size = 4 amplitudes)
const static cplx_t state_pure_two_qubits_zero[4] = {1.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0};
const static cplx_t state_pure_two_qubits_one[4] = {0.0 + I*0.0, 1.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0};
const static cplx_t state_pure_two_qubits_two[4] = {0.0 + I*0.0, 0.0 + I*0.0, 1.0 + I*0.0, 0.0 + I*0.0};
const static cplx_t state_pure_two_qubits_three[4] = {0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 1.0 + I*0.0};
const static cplx_t state_pure_two_qubits_plus[4] = {
    INVSQRT4 + I*0.0, INVSQRT4 + I*0.0, INVSQRT4 + I*0.0, INVSQRT4 + I*0.0
};
const static cplx_t state_pure_two_qubits_ghz[4] = {INVSQRT2 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, INVSQRT2 + I*0.0};

const static cplx_t *states_pure_two_qubits[6] = {
    state_pure_two_qubits_zero,
    state_pure_two_qubits_one,
    state_pure_two_qubits_two,
    state_pure_two_qubits_three,
    state_pure_two_qubits_plus,
    state_pure_two_qubits_ghz,
};


// 3-qubit pure states (size = 8 amplitudes)
const static cplx_t state_pure_three_qubits_zero[8] = {
    1.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0
};
const static cplx_t state_pure_three_qubits_one[8] = {
    0.0 + I*0.0, 1.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0
};
const static cplx_t state_pure_three_qubits_two[8] = {
    0.0 + I*0.0, 0.0 + I*0.0, 1.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0
};
const static cplx_t state_pure_three_qubits_three[8] = {
    0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 1.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0
};
const static cplx_t state_pure_three_qubits_four[8] = {
    0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 1.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0
};
const static cplx_t state_pure_three_qubits_five[8] = {
    0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 1.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0
};
const static cplx_t state_pure_three_qubits_six[8] = {
    0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 1.0 + I*0.0, 0.0 + I*0.0
};
const static cplx_t state_pure_three_qubits_seven[8] = {
    0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 1.0 + I*0.0
};
const static cplx_t state_pure_three_qubits_plus[8] = {
    INVSQRT8 + I*0.0, INVSQRT8 + I*0.0, INVSQRT8 + I*0.0, INVSQRT8 + I*0.0,
    INVSQRT8 + I*0.0, INVSQRT8 + I*0.0, INVSQRT8 + I*0.0, INVSQRT8 + I*0.0
};
const static cplx_t state_pure_three_qubits_ghz[8] = {
    INVSQRT2 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, INVSQRT2 + I*0.0
};
const static cplx_t *states_pure_three_qubits[10] = {
    state_pure_three_qubits_zero,
    state_pure_three_qubits_one,
    state_pure_three_qubits_two,
    state_pure_three_qubits_three,
    state_pure_three_qubits_four,
    state_pure_three_qubits_five,
    state_pure_three_qubits_six,
    state_pure_three_qubits_seven,
    state_pure_three_qubits_plus,
    state_pure_three_qubits_ghz,
};


// 4-qubit pure states (size = 16 amplitudes)
// 4-qubit pure states (size = 16 amplitudes)
const static cplx_t state_pure_four_qubits_zero[16] = {
    1.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0,
    0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0
};

const static cplx_t state_pure_four_qubits_one[16] = {
    0.0 + I*0.0, 1.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0,
    0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0
};

const static cplx_t state_pure_four_qubits_two[16] = {
    0.0 + I*0.0, 0.0 + I*0.0, 1.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0,
    0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0
};

const static cplx_t state_pure_four_qubits_three[16] = {
    0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 1.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0,
    0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0
};

const static cplx_t state_pure_four_qubits_four[16] = {
    0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 1.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0,
    0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0
};

const static cplx_t state_pure_four_qubits_five[16] = {
    0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 1.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0,
    0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0
};

const static cplx_t state_pure_four_qubits_six[16] = {
    0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 1.0 + I*0.0, 0.0 + I*0.0,
    0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0
};

const static cplx_t state_pure_four_qubits_seven[16] = {
    0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 1.0 + I*0.0,
    0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0
};

const static cplx_t state_pure_four_qubits_eight[16] = {
    0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0,
    1.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0
};

const static cplx_t state_pure_four_qubits_nine[16] = {
    0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0,
    0.0 + I*0.0, 1.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0
};

const static cplx_t state_pure_four_qubits_ten[16] = {
    0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0,
    0.0 + I*0.0, 0.0 + I*0.0, 1.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0
};

const static cplx_t state_pure_four_qubits_eleven[16] = {
    0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0,
    0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 1.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0
};

const static cplx_t state_pure_four_qubits_twelve[16] = {
    0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0,
    0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 1.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0
};

const static cplx_t state_pure_four_qubits_thirteen[16] = {
    0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0,
    0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 1.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0
};

const static cplx_t state_pure_four_qubits_fourteen[16] = {
    0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0,
    0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 1.0 + I*0.0, 0.0 + I*0.0
};

const static cplx_t state_pure_four_qubits_fifteen[16] = {
    0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0,
    0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 1.0 + I*0.0
};

const static cplx_t state_pure_four_qubits_plus[16] = {
    INVSQRT16 + I*0.0, INVSQRT16 + I*0.0, INVSQRT16 + I*0.0, INVSQRT16 + I*0.0,
    INVSQRT16 + I*0.0, INVSQRT16 + I*0.0, INVSQRT16 + I*0.0, INVSQRT16 + I*0.0,
    INVSQRT16 + I*0.0, INVSQRT16 + I*0.0, INVSQRT16 + I*0.0, INVSQRT16 + I*0.0,
    INVSQRT16 + I*0.0, INVSQRT16 + I*0.0, INVSQRT16 + I*0.0, INVSQRT16 + I*0.0
};

const static cplx_t state_pure_four_qubits_ghz[16] = {
    INVSQRT2 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0,
    0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, INVSQRT2 + I*0.0
};

const static cplx_t* states_pure_four_qubits[18] = {
    state_pure_four_qubits_zero,
    state_pure_four_qubits_one,
    state_pure_four_qubits_two,
    state_pure_four_qubits_three,
    state_pure_four_qubits_four,
    state_pure_four_qubits_five,
    state_pure_four_qubits_six,
    state_pure_four_qubits_seven,
    state_pure_four_qubits_eight,
    state_pure_four_qubits_nine,
    state_pure_four_qubits_ten,
    state_pure_four_qubits_eleven,
    state_pure_four_qubits_twelve,
    state_pure_four_qubits_thirteen,
    state_pure_four_qubits_fourteen,
    state_pure_four_qubits_fifteen,
    state_pure_four_qubits_plus,
    state_pure_four_qubits_ghz,
};

// Array of the pure state instances per qubit count and their respective number
const cplx_t * const * states_pure[4] = {
    states_pure_one_qubit, states_pure_two_qubits, states_pure_three_qubits, states_pure_four_qubits
};

unsigned nstates_pure[4] = {3, 6, 10, 18};

#ifdef __cplusplus
}
#endif

#endif //Q_TEST_H