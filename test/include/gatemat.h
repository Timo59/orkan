// gatemat.h    - Matrix representation of quantum gates

#ifndef GATEMAT_H
#define GATEMAT_H

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
 * Function definitions
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
cplx_t* mat_single_qubit_gate(unsigned nqubits, const cplx_t *gate, unsigned target);

#ifdef __cplusplus
}
#endif

#endif // GATEMAT_H