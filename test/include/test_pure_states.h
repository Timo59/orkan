// test_pure.c  - Create and destroy state vectors for testing quantum circuits

#ifndef TEST_PURE_STATES_H
#define TEST_PURE_STATES_H

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
cplx_t** test_bell_pure(unsigned nqubits);


/*
 * @brief   Returns the GHZ state for n qubits: |GHZ⟩ = (|0...0⟩ + |1...1⟩)/√2
 *
 * @param[in]   nqubits     Number of qubits (must be >= 2)
 *
 * @returns Array of size 1 containing the GHZ state vector, or NULL if nqubits < 2.
 */
cplx_t** test_ghz_pure(unsigned nqubits);


/*
 * @brief   Returns the W state for n qubits: |W⟩ = (|100...0⟩ + |010...0⟩ + ... + |00...01⟩)/√n
 *
 * @param[in]   nqubits     Number of qubits (must be >= 3)
 *
 * @returns Array of size 1 containing the W state vector, or NULL if nqubits < 3.
 */
cplx_t** test_w_pure(unsigned nqubits);


/*
 * @brief   Returns a variety of state vectors with <nqubits> qubits
 *
 * @param[in]       nqubits     Number of qubits
 * @param[in,out]   nvecs       Unsigned integer passed by reference. On exit, the number of state vectors
 *
 * @returns Array of size <...> of COMPLEX DOUBLE arrays; the state vectors of quantum test states.
 */
cplx_t** test_mk_states_pure(unsigned nqubits, unsigned *nvecs);


/*
 * @brief   Frees all memory allocated for quantum test states with <nqubits> qubits
 *
 * @param[in]   nqubits     Number of qubits
 */
void test_rm_states_pure(unsigned nqubits, cplx_t **states);

#ifdef __cplusplus
}
#endif

#endif // TEST_PURE_STATES_H