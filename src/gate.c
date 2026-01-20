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

/*
 * =====================================================================================================================
 * Clifford gates
 * =====================================================================================================================
 */

extern void h_pure(state_t *state, qubit_t target);
extern void h_mixed(state_t *state, qubit_t target);

extern void s_pure(state_t *state, qubit_t target);
extern void s_mixed(state_t *state, qubit_t target);

extern void sdg_pure(state_t *state, qubit_t target);
extern void sdg_mixed(state_t *state, qubit_t target);

extern void t_pure(state_t *state, qubit_t target);
extern void t_mixed(state_t *state, qubit_t target);

extern void tdg_pure(state_t *state, qubit_t target);
extern void tdg_mixed(state_t *state, qubit_t target);

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

qs_error_t h(state_t *state, const qubit_t target) {
    if (!state) {
        return QS_ERR_NULL;
    }
    if (!state->data) {
        return QS_ERR_NULL;
    }
    if (target >= state->qubits) {
        return QS_ERR_QUBIT;
    }

    if (state->type == PURE) h_pure(state, target);
    else h_mixed(state, target);

    return QS_OK;
}

qs_error_t s(state_t *state, const qubit_t target) {
    if (!state) {
        return QS_ERR_NULL;
    }
    if (!state->data) {
        return QS_ERR_NULL;
    }
    if (target >= state->qubits) {
        return QS_ERR_QUBIT;
    }

    if (state->type == PURE) s_pure(state, target);
    else s_mixed(state, target);

    return QS_OK;
}

qs_error_t sdg(state_t *state, const qubit_t target) {
    if (!state) {
        return QS_ERR_NULL;
    }
    if (!state->data) {
        return QS_ERR_NULL;
    }
    if (target >= state->qubits) {
        return QS_ERR_QUBIT;
    }

    if (state->type == PURE) sdg_pure(state, target);
    else sdg_mixed(state, target);

    return QS_OK;
}

qs_error_t t(state_t *state, const qubit_t target) {
    if (!state) {
        return QS_ERR_NULL;
    }
    if (!state->data) {
        return QS_ERR_NULL;
    }
    if (target >= state->qubits) {
        return QS_ERR_QUBIT;
    }

    if (state->type == PURE) t_pure(state, target);
    else t_mixed(state, target);

    return QS_OK;
}

qs_error_t tdg(state_t *state, const qubit_t target) {
    if (!state) {
        return QS_ERR_NULL;
    }
    if (!state->data) {
        return QS_ERR_NULL;
    }
    if (target >= state->qubits) {
        return QS_ERR_QUBIT;
    }

    if (state->type == PURE) tdg_pure(state, target);
    else tdg_mixed(state, target);

    return QS_OK;
}