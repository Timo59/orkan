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

/*
 * @brief   Applies Pauli-Y gate to pure state (internal implementation)
 */
extern void y_pure(state_t *state, qubit_t target);

/*
 * @brief   Applies Pauli-Y gate to mixed state in packed storage (internal implementation)
 */
extern void y_mixed(state_t *state, qubit_t target);

/*
 * @brief   Applies Pauli-Z gate to pure state (internal implementation)
 */
extern void z_pure(state_t *state, qubit_t target);

/*
 * @brief   Applies Pauli-Z gate to mixed state in packed storage (internal implementation)
 */
extern void z_mixed(state_t *state, qubit_t target);

qs_error_t x(state_t *state, const qubit_t target) {
    // Input validation: state points to valid address and is initialized, target is in range
    if (!state) {
        return QS_ERR_NULL;
    }
    if (!state->data) {
        return QS_ERR_NULL;
    }
    if (target >= state->qubits) {
        return QS_ERR_QUBIT;
    }

    if (state->type == PURE) x_pure(state, target);
    else x_mixed(state, target);

    return QS_OK;
}

qs_error_t y(state_t *state, const qubit_t target) {
    if (!state) {
        return QS_ERR_NULL;
    }
    if (!state->data) {
        return QS_ERR_NULL;
    }
    if (target >= state->qubits) {
        return QS_ERR_QUBIT;
    }

    if (state->type == PURE) y_pure(state, target);
    else y_mixed(state, target);

    return QS_OK;
}

qs_error_t z(state_t *state, const qubit_t target) {
    if (!state) {
        return QS_ERR_NULL;
    }
    if (!state->data) {
        return QS_ERR_NULL;
    }
    if (target >= state->qubits) {
        return QS_ERR_QUBIT;
    }

    if (state->type == PURE) z_pure(state, target);
    else z_mixed(state, target);

    return QS_OK;
}