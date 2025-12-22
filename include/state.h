#ifndef STATE_H
#define STATE_H

/*
 * =====================================================================================================================
 * includes
 * =====================================================================================================================
 */

#ifndef QTYPES_H
#include "q_types.h"
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

/*
 * @brief   Constant indicating whether the quantum state is a pure or mixed state
 */
typedef enum state_type {
    PURE,   // Pure quantum state, i.e., represented by a statevector
    MIXED   // Mixed quantum state, i.e., represented by the lower triangle of a density matrix
}state_type_t;

/*
 * @brief   Structure for multi-qubit quantum state
 */
typedef struct state {
    state_type_t type;  // Constant indicating whether the quantum state is a pure or a mixed state
    cplx_t *data;   // Representation; Statevector for pure, lower triangle of density matrix for mixed state
    qubit_t qubits; // Number of qubits of the quantum system
} state_t;

/*
 * =====================================================================================================================
 * State manipulation
 * =====================================================================================================================
 */

/*
 * @brief   Reset the quantum state and free allocated memory for the statevector
 *
 * @param[in,out]   state   Quantum state
 */
void state_free(state_t *state);

/*
 * @brief   Returns the size of the array representing the quantum state
 *
 * @param[in]       state   Quantum state (ATTENTION: Number of qubits must be initialized; otherwise returns zero)
 */
dim_t state_len(const state_t *state);

/*
 * @brief   Initialize a quantum state with the specified number of qubits and the passed representation
 *
 * @param[in,out]	state	Quantum state
 * @param[in]		qubits	Number of qubits
 * @param[in]       data    Representation of the initial state (if NULL; state is set to zero)
 */
void state_init(state_t *state, qubit_t qubits, cplx_t **data);

/*
 * @brief   Initialize a quantum state to the uniform superposition of all computational basis states (plus state)
 *
 * @param[in,out]	state	Quantum state
 * @param[in]		qubits	Number of qubits
 */
void state_plus(state_t *state, qubit_t qubits);

/*
 * @brief	Returns a deep copy of the quantum state
 *
 * @param[in,out]	state		Address of the state
 *
 * @returns	Deep copy of the quantum state
 */
state_t state_cp(const state_t* state);


#ifdef __cplusplus
}
#endif

#endif  // STATE_H
