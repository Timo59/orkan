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

/*
 * @brief   Returns the matrix representing the application of a two-qubit gate to a multiqubit statevector
 *
 * @param[in]   nqubits Number of qubits in the multiqubit state (must be >= 2)
 * @param[in]   gate    Double complex array of size 4-by-4; matrix representation of the two-qubit gate
 *                      in column-major format with basis order |q1=0,q2=0>, |q1=0,q2=1>, |q1=1,q2=0>, |q1=1,q2=1>
 * @param[in]   q1      First qubit index (e.g., control for CNOT); counted from the RIGHT
 * @param[in]   q2      Second qubit index (e.g., target for CNOT); counted from the RIGHT
 *
 * @returns Double complex array of size 4**(nqubits); the 2**(nqubits)-by-2**(nqubits) matrix representing
 *          the gate action. Returns NULL on error.
 */
cplx_t* mat_two_qubit_gate(unsigned nqubits, const cplx_t *gate, unsigned q1, unsigned q2);

/*
 * @brief   Returns the matrix representing the application of a three-qubit gate to a multiqubit statevector
 *
 * @param[in]   nqubits Number of qubits in the multiqubit state (must be >= 3)
 * @param[in]   gate    Double complex array of size 8-by-8; matrix representation of the three-qubit gate
 *                      in column-major format with basis order |q1,q2,q3>:
 *                      |000>, |001>, |010>, |011>, |100>, |101>, |110>, |111>
 * @param[in]   q1      First qubit index; counted from the RIGHT
 * @param[in]   q2      Second qubit index; counted from the RIGHT
 * @param[in]   q3      Third qubit index; counted from the RIGHT
 *
 * @returns Double complex array of size 4**(nqubits); the 2**(nqubits)-by-2**(nqubits) matrix representing
 *          the gate action. Returns NULL on error.
 */
cplx_t* mat_three_qubit_gate(unsigned nqubits, const cplx_t *gate, unsigned q1, unsigned q2, unsigned q3);

#ifdef __cplusplus
}
#endif

#endif // GATEMAT_H