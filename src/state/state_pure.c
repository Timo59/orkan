/**
 * @file state_pure.c
 * @brief Pure quantum state implementation
 *
 * Pure states are represented as statevectors of dimension 2^n where n is the
 * number of qubits. The amplitude for basis state |i> is stored at index i.
 *
 * Storage layout: [a_0, a_1, a_2, ..., a_{2^n-1}]
 * where a_i = <i|psi> is the amplitude for computational basis state |i>
 */

#include "state.h"
#include "utils.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>

/*
 * =====================================================================================================================
 * Pure state implementations
 * =====================================================================================================================
 */

idx_t state_pure_len(qubit_t qubits) {
    return POW2(qubits, idx_t);
}


void state_pure_plus(state_t *state, qubit_t qubits) {
    /* Initialize empty state (allocates zero-initialized storage) */
    state_init(state, qubits, NULL);
    if (!state->data) return;

    /* Compute uniform amplitude: 1/sqrt(2^n) */
    const idx_t len = state_pure_len(qubits);
    const cplx_t amplitude = 1.0 / sqrt((double)len);

    /* Set all amplitudes to the uniform value */
    for (idx_t i = 0; i < len; ++i) {
        state->data[i] = amplitude;
    }
}


cplx_t state_pure_get(const state_t *state, idx_t row) {
    assert(row < POW2(state->qubits, idx_t));
    return state->data[row];
}


void state_pure_set(state_t *state, idx_t row, cplx_t val) {
    assert(row < POW2(state->qubits, idx_t));
    state->data[row] = val;
}


void state_pure_print(const state_t *state) {
    printf("Type: PURE\nQubits: %u\n", state->qubits);
    vprint(state->data, state_pure_len(state->qubits));
}
