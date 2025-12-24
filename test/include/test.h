// test.h - Functions and objects to check library
#ifndef Q_TEST_H
#define Q_TEST_H

/*
 * =====================================================================================================================
 * includes
 * =====================================================================================================================
 */

#include <complex.h>

#ifndef QTYPES_H
#include  "q_types.h"
#endif

#include <stdio.h>

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
#define PRECISION   1e-8
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
                    "Element %zu Expected %.7g, Was %.7g", _i, _expected, _actual); \
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
 * @brief   Check whether two pointers are not equal
 */
#define TEST_ASSERT_NOT_EQUAL_PTR(expected, actual) \
    do { \
        if ((expected) == (actual)) { \
            char _msg[256]; \
            snprintf(_msg, sizeof(_msg), \
                "Pointers should not be equal. Both are %p", (void*)(actual)); \
            UNITY_TEST_FAIL(__LINE__, _msg); \
        } \
    } while (0)

/*
 * =====================================================================================================================
 * Pure quantum test states
 * =====================================================================================================================
 */

/*
 * @brief   Returns all computational basis states (eigenstates of Pauli-Z) of a system with <nqubits> qubits
 *
 * @param[in]   nqubits     Number of qubits
 *
 * @returns Array of size 2**<nqubits> of COMPLEX DOUBLE arrays; the state vectors of all computational basis states.
 */
cplx_t** test_cb_pure(unsigned nqubits);


/*
 * @brief   Returns all Hadamard basis states (eigenstates of Pauli-X) of a system with <nqubits> qubits in the
 *          computational basis
 *
 * @param[in]   nqubits     Number of qubits
 *
 * @returns Array of size 2**<nqubits> of COMPLEX DOUBLE arrays; the state vectors of all Hadamard basis states.
 */
cplx_t** test_xb_pure(unsigned nqubits);


/*
 * @brief   Returns all circular basis states (eigenstates of Pauli-Y) of a system with <nqubits> qubits in the
 *          computational basis
 *
 * @param[in]   nqubits     Number of qubits
 *
 * @returns Array of size 2**<nqubits> of COMPLEX DOUBLE arrays; the state vectors of all circular basis states.
 */
cplx_t** test_yb_pure(unsigned nqubits);


/*
 * @brief   Returns all 4 Bell states for 2 qubits
 *
 * @returns Array of size 4 of COMPLEX DOUBLE arrays; the state vectors of all Bell states (|Φ⁺⟩, |Φ⁻⟩, |Ψ⁺⟩, |Ψ⁻⟩).
 */
cplx_t** test_bell_states(void);


/*
 * @brief   Returns the GHZ state for n qubits: |GHZ⟩ = (|0...0⟩ + |1...1⟩)/√2
 *
 * @param[in]   nqubits     Number of qubits (must be >= 2)
 *
 * @returns Array of size 1 containing the GHZ state vector, or NULL if nqubits < 2.
 */
cplx_t** test_ghz_state(unsigned nqubits);


/*
 * @brief   Returns the W state for n qubits: |W⟩ = (|100...0⟩ + |010...0⟩ + ... + |00...01⟩)/√n
 *
 * @param[in]   nqubits     Number of qubits (must be >= 3)
 *
 * @returns Array of size 1 containing the W state vector, or NULL if nqubits < 3.
 */
cplx_t** test_w_state(unsigned nqubits);


/*
 * @brief   Returns a variety of state vectors with <nqubits> qubits
 *
 * @param[in]       nqubits     Number of qubits
 * @param[in,out]   nvecs       Unsigned integer passed by reference. On exit, the number of state vectors
 *
 * @returns Array of size <...> of COMPLEX DOUBLE arrays; the state vectors of quantum test states.
 */
cplx_t** test_gen_states_pure(unsigned nqubits, unsigned *nvecs);


/*
 * @brief   Frees all memory allocated for quantum test states with <nqubits> qubits
 *
 * @param[in]   nqubits     Number of qubits
 */
void test_rm_states_pure(unsigned nqubits, cplx_t **states);


#ifdef __cplusplus
}
#endif

#endif //Q_TEST_H