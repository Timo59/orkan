// gate.h   - Functions representing the action of local quantum gates to a general quantum state

/*
 * =====================================================================================================================
 * includes
 * =====================================================================================================================
 */

#ifndef GATE_H
#include "gate.h"
#endif

#include <stdio.h>

/*
 * =====================================================================================================================
 * Pauli gates
 * =====================================================================================================================
 */

/*
 * @brief   Applies Pauli-X gate to pure state (internal implementation)
 * @param   state   Quantum state (must be non-NULL with valid data)
 * @param   target  Target qubit index (must be in range [0, state->qubits))
 * @note    Input validation is performed by dispatcher x(). This function trusts its inputs.
 */
extern void x_pure(state_t *state, qubit_t target);

/*
 * @brief   Applies Pauli-X gate to mixed state in packed storage (internal implementation)
 * @param   state   Quantum state with packed lower-triangular density matrix
 * @param   target  Target qubit index
 * @note    Operates directly on packed format via qHiPSTER stride indexing
 * @note    Input validation is performed by dispatcher x(). This function trusts its inputs.
 */
extern void x_mixed(state_t *state, qubit_t target);

void x(state_t *state, const qubit_t target) {
    // Input validation: state points to valid address and is initialized, target is in range
    if (!state) {
        fprintf(stderr, "x(): state is NULL\n");
        return;
    }
    if (!state->data) {
        fprintf(stderr, "x(): state data is NULL\n");
        return;
    }
    if (target >= state->qubits) {
        fprintf(stderr, "x(): target is out of range; Expected [0,%u], Was %u\n", state->qubits, target);
        return;
    }

    if (state->type == PURE) x_pure(state, target);
    else x_mixed(state, target);
}