#ifndef CHANNEL_H
#define CHANNEL_H

/*
 * =====================================================================================================================
 * Includes
 * =====================================================================================================================
 */
#include "state.h"

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
 * Kraus representation type
 * =====================================================================================================================
 *
 * kraus_t stores the Kraus representation of a quantum channel acting on n_qubits qubits
 * as a concatenation of the n_terms Kraus matrices in row-major format.
 */
typedef struct {
  qubit_t n_qubits;
  uint64_t n_terms;
  cplx_t* data;
} kraus_t;

/*
 * =====================================================================================================================
 * Superoperator matrix representation type
 * =====================================================================================================================
 *
 * superop_t stores the superoperator matrix of a quantum channel acting on n_qubits qubits
 * in row-major format.
 */
typedef struct {
  qubit_t n_qubits;
  cplx_t* data;
} superop_t;

/*
 * =====================================================================================================================
 * Compute superoperator from Kraus representation
 * =====================================================================================================================
 */

/*
 * @brief Returns the superoperator from a Kraus representation
 */
superop_t kraus_to_superop(const kraus_t *kraus);

/*
 * =====================================================================================================================
 * Quantum channels
 * =====================================================================================================================
 */

/*
 * @brief Apply the 1-local quantum channel given in terms of its superoperator to the quantum state
 */
void channel_1q(state_t *state, const superop_t *sop, qubit_t target);

#ifdef __cplusplus
}
#endif

#endif // CHANNEL_H
