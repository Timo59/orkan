// gate.h   - Functions representing the action of local quantum gates to a general quantum state

/*
 * =====================================================================================================================
 * includes
 * =====================================================================================================================
 */

#ifndef GATE_H
#include "gate.h"
#endif

/*
 * =====================================================================================================================
 * Pauli gates
 * =====================================================================================================================
 */

extern void x_pure(state_t *state, qubit_t target);
extern void x_mixed(state_t *state, qubit_t target);
void x(state_t *state, const qubit_t target) {
    if (state->type == PURE) x_pure(state, target);
    else x_mixed(state, target);
}