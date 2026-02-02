// test_gates.h  - Unit tests for quantum gate operations on pure and mixed states

#ifndef TEST_GATE_H
#define TEST_GATE_H

/*
 * =====================================================================================================================
 * includes
 * =====================================================================================================================
 */

#ifndef GATE_H
#include "gate.h"
#endif

#ifndef Q_TEST_H
#include "test.h"
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
 * Type definitions
 * =====================================================================================================================
 */

// Function representing a qubit gate that acts locally on one qubit
typedef qs_error_t (*single_qubit_gate)(state_t *state, qubit_t target);

// Function representing a parameterized single-qubit rotation gate
typedef qs_error_t (*rotation_gate)(state_t *state, qubit_t target, double theta);

/*
 * =====================================================================================================================
 * Single-qubit gate matrices
 * =====================================================================================================================
 */

static const cplx_t IDMAT[4] = {1.0 + 0.0 * I, \
                                0.0 + 0.0 * I, \
                                0.0 + 0.0 * I, \
                                1.0 + 0.0 * I};

static const cplx_t XMAT[4] = {0.0 + 0.0 * I, \
                               1.0 + 0.0 * I, \
                               1.0 + 0.0 * I, \
                               0.0 + 0.0 * I};

static const cplx_t YMAT[4] = {0.0 + 0.0 * I, \
                               0.0 + 1.0 * I, \
                               0.0 - 1.0 * I, \
                               0.0 + 0.0 * I};

static const cplx_t ZMAT[4] = {1.0 + 0.0 * I, \
                               0.0 + 0.0 * I, \
                               0.0 + 0.0 * I, \
                               -1.0 + 0.0 * I};

static const cplx_t SMAT[4] = {1.0 + 0.0 * I, \
                               0.0 + 0.0 * I, \
                               0.0 + 0.0 * I, \
                               0.0 + 1.0 * I};

static const cplx_t HMAT[4] = {INVSQRT2 + 0.0 * I, \
                               INVSQRT2 + 0.0 * I, \
                               INVSQRT2 + 0.0 * I, \
                               -INVSQRT2 + 0.0 * I};

static const cplx_t TMAT[4] = {1.0 + 0.0 * I, \
                               0.0 + 0.0 * I, \
                               0.0 + 0.0 * I, \
                               INVSQRT2 + INVSQRT2 * I};

static const cplx_t SDGMAT[4] = {1.0 + 0.0 * I, \
                                 0.0 + 0.0 * I, \
                                 0.0 + 0.0 * I, \
                                 0.0 - 1.0 * I};

static const cplx_t TDGMAT[4] = {1.0 + 0.0 * I, \
                                 0.0 + 0.0 * I, \
                                 0.0 + 0.0 * I, \
                                 INVSQRT2 - INVSQRT2 * I};

static const cplx_t P0MAT[4] = {1.0 + 0.0 * I, \
                                0.0 + 0.0 * I, \
                                0.0 + 0.0 * I, \
                                0.0 + 0.0 * I};

static const cplx_t P1MAT[4] = {0.0 + 0.0 * I, \
                                0.0 + 0.0 * I, \
                                0.0 + 0.0 * I, \
                                1.0 + 0.0 * I};

/*
 * =====================================================================================================================
 * Two-qubit gate matrices
 * =====================================================================================================================
 */

static const cplx_t SWAPMAT[16] = {1.0 + 0.0 * I, 0.0 + 0.0 * I, 0.0 + 0.0 * I, 0.0 + 0.0 * I, \
                                   0.0 + 0.0 * I, 0.0 + 0.0 * I, 1.0 + 0.0 * I, 0.0 + 0.0 * I, \
                                   0.0 + 0.0 * I, 1.0 + 0.0 * I, 0.0 + 0.0 * I, 0.0 + 0.0 * I, \
                                   0.0 + 0.0 * I, 0.0 + 0.0 * I, 0.0 + 0.0 * I, 1.0 + 0.0 * I};

/*
 * =====================================================================================================================
 * Test harness function declarations
 * =====================================================================================================================
 */

/*
 * @brief   Unit test of the single qubit gate on pure states for each qubit of a multi-qubit state
 *
 * @param[in]   gate    Function representing the single qubit gate
 * @param[in]   mat     Matrix representation of the single qubit operation (2×2)
 */
void testSingleQubitGate(single_qubit_gate gate, const cplx_t *mat);

/*
 * @brief   Unit test of the single qubit gate on mixed states for each qubit of a multi-qubit state
 *
 * @param[in]   gate    Function representing the single qubit gate
 * @param[in]   mat     Matrix representation of the single qubit operation (2×2)
 */
void testSingleQubitGateMixed(single_qubit_gate gate, const cplx_t *mat);

/*
 * @brief   Unit test of a parameterized rotation gate on pure states
 *
 * @param[in]   gate    Function representing the rotation gate
 * @param[in]   mat_fn  Function that builds the 2×2 matrix for a given theta
 */
void testRotationGate(rotation_gate gate, void (*mat_fn)(double theta, cplx_t *mat));

/*
 * @brief   Unit test of a parameterized rotation gate on mixed states
 *
 * @param[in]   gate    Function representing the rotation gate
 * @param[in]   mat_fn  Function that builds the 2×2 matrix for a given theta
 */
void testRotationGateMixed(rotation_gate gate, void (*mat_fn)(double theta, cplx_t *mat));

/*
 * @brief   Build Rx(θ) matrix: [[cos(θ/2), -i·sin(θ/2)], [-i·sin(θ/2), cos(θ/2)]]
 */
void mat_rx(double theta, cplx_t *mat);

/*
 * @brief   Build Ry(θ) matrix: [[cos(θ/2), -sin(θ/2)], [sin(θ/2), cos(θ/2)]]
 */
void mat_ry(double theta, cplx_t *mat);

/*
 * @brief   Build Rz(θ) matrix: [[e^(-iθ/2), 0], [0, e^(iθ/2)]]
 */
void mat_rz(double theta, cplx_t *mat);

/*
 * =====================================================================================================================
 * Multiqubit gates
 * =====================================================================================================================
 */




#ifdef __cplusplus
}
#endif

#endif // TEST_GATE_H
