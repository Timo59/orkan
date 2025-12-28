// test_qhipster.h  - Unit tests for representations of quantum gates on pure states. The reference is given by matrix
//                    vector multiplication of the state vector with the matrix representation of the quantum gate.

#ifndef TEST_QHIPSTER_H
#define TEST_QHIPSTER_H

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
typedef void (*single_qubit_gate)(state_t *state, qubit_t target);

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
 * Single qubit gates
 * =====================================================================================================================
 */

/*
 * @brief   Unit test of the single qubit gate on each qubit of a multi-qubit state of up to MAXQUBITS qubits
 *
 * @param[in]   gate    Function representing the single qubit gate
 * @param[in]   mat     Matrix representation of the single qubit operation
 */
void testSingleQubitGate(single_qubit_gate gate, const cplx_t *mat);

/*
 * =====================================================================================================================
 * Multiqubit gates
 * =====================================================================================================================
 */




#ifdef __cplusplus
}
#endif

#endif // TEST_QHIPSTER_H