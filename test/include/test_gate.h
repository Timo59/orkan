// test_gates.h  - Unit tests for quantum gate operations on pure, packed, and tiled states

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
typedef void (*single_qubit_gate)(state_t *state, qubit_t target);

// Function representing a parameterized single-qubit rotation gate
typedef void (*rotation_gate)(state_t *state, qubit_t target, double theta);

// Function representing a two-qubit gate (e.g., CNOT)
typedef void (*two_qubit_gate)(state_t *state, qubit_t q1, qubit_t q2);

// Function representing a parameterized two-qubit gate (e.g., CP)
typedef void (*two_qubit_rotation_gate)(state_t *state, qubit_t q1, qubit_t q2, double theta);

// Function representing a three-qubit gate (e.g., Toffoli)
typedef void (*three_qubit_gate)(state_t *state, qubit_t q1, qubit_t q2, qubit_t q3);

/*
 * =====================================================================================================================
 * Single-qubit gate matrices (column-major: [col0_row0, col0_row1, col1_row0, col1_row1])
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

// Hy (Hadamard-Y) matrix: [[1,-1],[1,1]]/sqrt(2)
// Column-major: col0=[1,1]/sqrt(2), col1=[-1,1]/sqrt(2)
static const cplx_t HYMAT[4] = {INVSQRT2 + 0.0 * I, \
                                INVSQRT2 + 0.0 * I, \
                                -INVSQRT2 + 0.0 * I, \
                                INVSQRT2 + 0.0 * I};

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

// CNOT (Controlled-X) matrix in column-major format
// Basis order: |00>, |01>, |10>, |11> (target is rightmost/LSB)
// |00> -> |00>, |01> -> |01>, |10> -> |11>, |11> -> |10>
static const cplx_t CXMAT[16] = {1.0 + 0.0 * I, 0.0 + 0.0 * I, 0.0 + 0.0 * I, 0.0 + 0.0 * I, \
                                 0.0 + 0.0 * I, 1.0 + 0.0 * I, 0.0 + 0.0 * I, 0.0 + 0.0 * I, \
                                 0.0 + 0.0 * I, 0.0 + 0.0 * I, 0.0 + 0.0 * I, 1.0 + 0.0 * I, \
                                 0.0 + 0.0 * I, 0.0 + 0.0 * I, 1.0 + 0.0 * I, 0.0 + 0.0 * I};

// CY (Controlled-Y) matrix: diag(I, Y) in control-target basis
// |00> -> |00>, |01> -> |01>, |10> -> i|11>, |11> -> -i|10>
static const cplx_t CYMAT[16] = {1.0 + 0.0 * I, 0.0 + 0.0 * I, 0.0 + 0.0 * I, 0.0 + 0.0 * I, \
                                 0.0 + 0.0 * I, 1.0 + 0.0 * I, 0.0 + 0.0 * I, 0.0 + 0.0 * I, \
                                 0.0 + 0.0 * I, 0.0 + 0.0 * I, 0.0 + 0.0 * I, 0.0 + 1.0 * I, \
                                 0.0 + 0.0 * I, 0.0 + 0.0 * I, 0.0 - 1.0 * I, 0.0 + 0.0 * I};

// CZ (Controlled-Z) matrix: diag(1, 1, 1, -1)
static const cplx_t CZMAT[16] = {1.0 + 0.0 * I, 0.0 + 0.0 * I, 0.0 + 0.0 * I, 0.0 + 0.0 * I, \
                                 0.0 + 0.0 * I, 1.0 + 0.0 * I, 0.0 + 0.0 * I, 0.0 + 0.0 * I, \
                                 0.0 + 0.0 * I, 0.0 + 0.0 * I, 1.0 + 0.0 * I, 0.0 + 0.0 * I, \
                                 0.0 + 0.0 * I, 0.0 + 0.0 * I, 0.0 + 0.0 * I, -1.0 + 0.0 * I};

// CS (Controlled-S) matrix: diag(1, 1, 1, i)
static const cplx_t CSMAT[16] = {1.0 + 0.0 * I, 0.0 + 0.0 * I, 0.0 + 0.0 * I, 0.0 + 0.0 * I, \
                                 0.0 + 0.0 * I, 1.0 + 0.0 * I, 0.0 + 0.0 * I, 0.0 + 0.0 * I, \
                                 0.0 + 0.0 * I, 0.0 + 0.0 * I, 1.0 + 0.0 * I, 0.0 + 0.0 * I, \
                                 0.0 + 0.0 * I, 0.0 + 0.0 * I, 0.0 + 0.0 * I, 0.0 + 1.0 * I};

// CSdg (Controlled-S-dagger) matrix: diag(1, 1, 1, -i)
static const cplx_t CSDGMAT[16] = {1.0 + 0.0 * I, 0.0 + 0.0 * I, 0.0 + 0.0 * I, 0.0 + 0.0 * I, \
                                   0.0 + 0.0 * I, 1.0 + 0.0 * I, 0.0 + 0.0 * I, 0.0 + 0.0 * I, \
                                   0.0 + 0.0 * I, 0.0 + 0.0 * I, 1.0 + 0.0 * I, 0.0 + 0.0 * I, \
                                   0.0 + 0.0 * I, 0.0 + 0.0 * I, 0.0 + 0.0 * I, 0.0 - 1.0 * I};

// CT (Controlled-T) matrix: diag(1, 1, 1, e^(iπ/4))
static const cplx_t CTMAT[16] = {1.0 + 0.0 * I, 0.0 + 0.0 * I, 0.0 + 0.0 * I, 0.0 + 0.0 * I, \
                                 0.0 + 0.0 * I, 1.0 + 0.0 * I, 0.0 + 0.0 * I, 0.0 + 0.0 * I, \
                                 0.0 + 0.0 * I, 0.0 + 0.0 * I, 1.0 + 0.0 * I, 0.0 + 0.0 * I, \
                                 0.0 + 0.0 * I, 0.0 + 0.0 * I, 0.0 + 0.0 * I, INVSQRT2 + INVSQRT2 * I};

// CTdg (Controlled-T-dagger) matrix: diag(1, 1, 1, e^(-iπ/4))
static const cplx_t CTDGMAT[16] = {1.0 + 0.0 * I, 0.0 + 0.0 * I, 0.0 + 0.0 * I, 0.0 + 0.0 * I, \
                                   0.0 + 0.0 * I, 1.0 + 0.0 * I, 0.0 + 0.0 * I, 0.0 + 0.0 * I, \
                                   0.0 + 0.0 * I, 0.0 + 0.0 * I, 1.0 + 0.0 * I, 0.0 + 0.0 * I, \
                                   0.0 + 0.0 * I, 0.0 + 0.0 * I, 0.0 + 0.0 * I, INVSQRT2 - INVSQRT2 * I};

/*
 * =====================================================================================================================
 * Three-qubit gate matrices
 * =====================================================================================================================
 */

// CCX (Toffoli/CCNOT) matrix in column-major format
// Basis order: |q1,q2,q3> = |000>, |001>, |010>, |011>, |100>, |101>, |110>, |111>
// Only flips q3 (target) when both q1 and q2 (controls) are 1: |110> <-> |111>
static const cplx_t CCXMAT[64] = {
    1,0,0,0,0,0,0,0,  /* col 0: |000> -> |000> */
    0,1,0,0,0,0,0,0,  /* col 1: |001> -> |001> */
    0,0,1,0,0,0,0,0,  /* col 2: |010> -> |010> */
    0,0,0,1,0,0,0,0,  /* col 3: |011> -> |011> */
    0,0,0,0,1,0,0,0,  /* col 4: |100> -> |100> */
    0,0,0,0,0,1,0,0,  /* col 5: |101> -> |101> */
    0,0,0,0,0,0,0,1,  /* col 6: |110> -> |111> */
    0,0,0,0,0,0,1,0   /* col 7: |111> -> |110> */
};

/*
 * =====================================================================================================================
 * Test harness function declarations -- pure states
 * =====================================================================================================================
 */

void testSingleQubitGate(single_qubit_gate gate, const cplx_t *mat);
void testRotationGate(rotation_gate gate, void (*mat_fn)(double theta, cplx_t *mat));
void testTwoQubitGate(two_qubit_gate gate, const cplx_t *mat);
void testTwoQubitRotationGate(two_qubit_rotation_gate gate, void (*mat_fn)(double theta, cplx_t *mat));
void testThreeQubitGate(three_qubit_gate gate, const cplx_t *mat);

/*
 * =====================================================================================================================
 * Test harness function declarations -- packed mixed states
 * =====================================================================================================================
 */

void testSingleQubitGateMixed(single_qubit_gate gate, const cplx_t *mat);
void testRotationGateMixed(rotation_gate gate, void (*mat_fn)(double theta, cplx_t *mat));
void testTwoQubitGateMixed(two_qubit_gate gate, const cplx_t *mat);
void testTwoQubitRotationGateMixed(two_qubit_rotation_gate gate, void (*mat_fn)(double theta, cplx_t *mat));
void testThreeQubitGateMixed(three_qubit_gate gate, const cplx_t *mat);

/*
 * =====================================================================================================================
 * Test harness function declarations -- tiled mixed states
 * =====================================================================================================================
 */

void testSingleQubitGateTiled(single_qubit_gate gate, const cplx_t *mat);
void testRotationGateTiled(rotation_gate gate, void (*mat_fn)(double theta, cplx_t *mat));
void testTwoQubitGateTiled(two_qubit_gate gate, const cplx_t *mat);
void testTwoQubitRotationGateTiled(two_qubit_rotation_gate gate, void (*mat_fn)(double theta, cplx_t *mat));
void testThreeQubitGateTiled(three_qubit_gate gate, const cplx_t *mat);

/*
 * =====================================================================================================================
 * Rotation gate matrix builders
 * =====================================================================================================================
 */

void mat_rx(double theta, cplx_t *mat);
void mat_ry(double theta, cplx_t *mat);
void mat_rz(double theta, cplx_t *mat);
void mat_p(double theta, cplx_t *mat);
void mat_cp(double theta, cplx_t *mat);
void mat_cpdg(double theta, cplx_t *mat);

/*
 * =====================================================================================================================
 * CH matrix builder (4x4, used for CH controlled gate test)
 * =====================================================================================================================
 */

// CH (Controlled-Hadamard) -- built dynamically since it contains irrational entries
void mat_ch_build(cplx_t *mat);

// CHy (Controlled-Hadamard-Y) -- built dynamically
void mat_chy_build(cplx_t *mat);


#ifdef __cplusplus
}
#endif

#endif // TEST_GATE_H
