#ifndef Q_TEST_QHIPSTER_H
#define Q_TEST_QHIPSTER_H

/*
 * =====================================================================================================================
 * includes
 * =====================================================================================================================
 */

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
                               0.0 - 1.0 * I, \
                               0.0 + 1.0 * I, \
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
 * Linear algebra
 * =====================================================================================================================
 */

/*
 * @brief   Returns the Kronecker product of the k-by-l matrix A and the m-by-n matrix B
 *
 * @param[in]   k   Unsigned integer; the number of rows in A
 * @param[in]   l   Unsigned integer; the number of columns in A
 * @param[in]   A   Double complex array of size k*l; represents the matrix A in column major format
 * @param[in]   m   Unsigned integer; the number of rows in B
 * @param[in]   n   Unsigned integer; the number of columns in B in column major format
 *
 * @returns Double complex array of size (k*m)*(l*n); The (k*m)-by-(l*n) matrix which is the Kronecker product of the
 *          matrices A and B in column major format
 */
cplx_t* kron(unsigned k, unsigned l, const cplx_t A[], unsigned m, unsigned n, const cplx_t B[]);

/*
 * =====================================================================================================================
 * Multiqubit gates
 * =====================================================================================================================
 */

/*
 * @brief   Returns the matrix representing the application of the identity gate to a multiqubit state
 *
 * @param[in]   nqubits Number of qubits in the multiqubit state
 */
cplx_t* mat_id(unsigned nqubits);

/*
 * @brief   Returns the matrix representing the application of a single-qubit gate to a multiqubit statevector
 *
 * @param[in]   nqubits Number of qubits in the multiqubit state
 * @param[in]   gate    Double complex array of size 2-by-2; matrix representation of the single qubit gate
 * @param[in]   target  Targeted qubit; counted from the RIGHT
 *
 * @returns Double complex array of size 4**(nqubits); the 2**(nqubits)-by-2**(nqubits) matrix which is the Kronecker
 *          product of identities for all inactive qubits with the single-qubit gate matrix on the target.
 */
cplx_t* mat_single_qubit_gate(unsigned nqubits, const cplx_t gate[], unsigned target);

/*
 * =====================================================================================================================
 * Post-X pure quantum states
 * =====================================================================================================================
 */



#ifdef __cplusplus
}
#endif

#endif // Q_TEST_QHIPSTER_H