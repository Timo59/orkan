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

/**
 * @brief Build the full 2^n × 2^n matrix for a Pauli string via Kronecker products
 *
 * Constructs M = M(s[0]) ⊗ M(s[1]) ⊗ ... ⊗ M(s[n-1]) where s[0] is qubit n-1
 * (big-endian, matching pauli_from_str convention) and M('I'), M('X'), M('Y'), M('Z')
 * are the standard 2×2 Pauli matrices (column-major). Completely independent of the
 * bitmask representation in pauli_t.
 *
 * @param[in] s         Big-endian string of 'I','X','Y','Z' characters (length n_qubits)
 * @param[in] n_qubits  Number of qubits; 0 returns a 1×1 matrix with value 1
 * @return Heap-allocated column-major matrix of size (2^n_qubits)^2, or NULL on OOM
 */
cplx_t *mat_pauli_str(const char *s, unsigned n_qubits);

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

/*
 * Rotation gate matrix builders — fill a 4-element column-major cplx_t[4] with the
 * 2×2 matrix for the given gate/angle.
 */
void mat_rx(double theta, cplx_t *mat);
void mat_ry(double theta, cplx_t *mat);
void mat_rz(double theta, cplx_t *mat);
void mat_p(double theta, cplx_t *mat);
void mat_rx_pi3(cplx_t *mat);

#ifdef __cplusplus
}
#endif

#endif // GATEMAT_H