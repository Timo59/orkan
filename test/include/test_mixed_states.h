// test_mixed_states.h - Create and destroy density matrices for testing quantum circuits

#ifndef TEST_MIXED_STATES_H
#define TEST_MIXED_STATES_H

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
 * @brief   Returns all computational basis states as density matrices (ρ = |i⟩⟨i|)
 *
 * @param[in]   nqubits     Number of qubits
 *
 * @returns Array of size 2^nqubits of density matrices in packed lower-triangular storage
 */
cplx_t** test_cb_mixed(unsigned nqubits);


/*
 * @brief   Returns all Hadamard basis states as density matrices (ρ = |ψ⟩⟨ψ| where |ψ⟩ = H|i⟩)
 *
 * @param[in]   nqubits     Number of qubits
 *
 * @returns Array of size 2^nqubits of density matrices in packed lower-triangular storage
 */
cplx_t** test_xb_mixed(unsigned nqubits);


/*
 * @brief   Returns all circular basis states as density matrices (ρ = |ψ⟩⟨ψ| where |ψ⟩ = Hy|i⟩)
 *
 * @param[in]   nqubits     Number of qubits
 *
 * @returns Array of size 2^nqubits of density matrices in packed lower-triangular storage
 */
cplx_t** test_yb_mixed(unsigned nqubits);


/*
 * @brief   Returns all 4 Bell states as density matrices
 *
 * @param[in]   nqubits     Number of qubits (must be 2)
 *
 * @returns Array of size 4 of density matrices in packed lower-triangular storage
 */
cplx_t** test_bell_mixed(unsigned nqubits);


/*
 * @brief   Returns the GHZ state as a density matrix: ρ = |GHZ⟩⟨GHZ| where |GHZ⟩ = (|0...0⟩ + |1...1⟩)/√2
 *
 * @param[in]   nqubits     Number of qubits (must be >= 2)
 *
 * @returns Array of size 1 containing the GHZ density matrix, or NULL if nqubits < 2
 */
cplx_t** test_ghz_mixed(unsigned nqubits);


/*
 * @brief   Returns the W state as a density matrix: ρ = |W⟩⟨W| where |W⟩ = (|100...0⟩ + ... + |00...01⟩)/√n
 *
 * @param[in]   nqubits     Number of qubits (must be >= 3)
 *
 * @returns Array of size 1 containing the W density matrix, or NULL if nqubits < 3
 */
cplx_t** test_w_mixed(unsigned nqubits);


/*
 * @brief   Generates the maximally mixed state I/2^n for n qubits
 *
 * @param[in]   nqubits     Number of qubits
 *
 * @returns Array of size 1 containing the maximally mixed state in packed format, or NULL on failure
 */
cplx_t** test_maximally_mixed(unsigned nqubits);


/*
 * @brief   Generates random statistical mixture ρ = Σᵢ pᵢ|ψᵢ⟩⟨ψᵢ| with random probabilities
 *
 * @param[in]   nqubits     Number of qubits
 *
 * @returns Array of size 1 containing the random mixture density matrix, or NULL on failure
 *
 * @note    Number of components and random seed are hardcoded internally for reproducibility
 */
cplx_t** test_random_mixture(unsigned nqubits);


/*
 * @brief   Generates a comprehensive set of mixed state test cases
 *
 * @param[in]       nqubits     Number of qubits
 * @param[in,out]   nmixed      Pointer to store the number of generated mixed states
 *
 * @returns Array of mixed states (density matrices in packed format) or NULL on failure
 *
 * @note    Generated states include:
 *          - Computational basis states as density matrices
 *          - Hadamard basis states as density matrices
 *          - Circular basis states as density matrices
 *          - Maximally mixed state
 *          - Random statistical mixtures
 *          - Entangled states: 4 Bell states (n=2) or GHZ + W states (n>2)
 */
cplx_t** test_mk_states_mixed(unsigned nqubits, unsigned* nmixed);


/*
 * @brief   Frees all memory allocated for mixed state test cases
 *
 * @param[in]   nqubits     Number of qubits
 * @param[in]   states      Array of mixed states to free
 */
void test_rm_states_mixed(unsigned nqubits, cplx_t** states);

#ifdef __cplusplus
}
#endif

#endif // TEST_MIXED_STATES_H
