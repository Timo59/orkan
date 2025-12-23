#include <math.h>

#ifndef STATE_H
#include "state.h"
#endif

#include <stdio.h>

#include <stdlib.h>

#if defined(__APPLE__)
    #include <vecLib/cblas_new.h>
#elif defined(__linux__)
    #include <cblas.h>
#endif

#ifndef UTILS_H
#include "utils.h"
#endif

/*
 * =====================================================================================================================
 * State manipulation
 * =====================================================================================================================
 */

void state_free(state_t *state) {
    free(state->data);
    state->data = NULL;

    state->qubits = 0;
}


dim_t state_len(const state_t *state) {
    // Hilbert space dimension
    const dim_t dim = POW2(state->qubits, dim_t);

    // Hilbert space dimension is the length of the statevector
    if (state->type == PURE) return dim;

    // Length of lower triangle for a density matrix on the Hilbert space
    return dim * (dim + 1) / 2;
}


void state_init(state_t *state, const qubit_t qubits, cplx_t **data) {
    // If state has been initialized before, reset it
    if (state->data) {
        free(state->data);
        state->data = NULL;
    }

	state->qubits = qubits;

    // Initialize state representation to the passed representation or to all zero if the latter is NULL
    if (data) {
        state->data = *data;
        *data = NULL;
    } else {
        const dim_t len = state_len(state); // Size of the array representing the quantum state
        if(!((state->data = calloc(len, sizeof(*state->data))))) {
            fprintf(stderr, "state_init(): state->data allocation failed\n");
            state->qubits = 0;
        }
    }
}


void state_plus(state_t *state, const qubit_t qubits) {
    // Normalization constant
    cplx_t prefactor;
    if (state->type == PURE) {
        prefactor = 1.0 / sqrt(1 << qubits) + I*0.0;
    } else {
        prefactor = 1.0 / (1 << qubits) + I*0.0;
    }

    // Initialize empty quantum state
    state_init(state, qubits, NULL);
    if (!state->data) {
        fprintf(stderr, "state_plus(): state_init() failed\n");
        return;
    }

    // Fill array with normalization constant
    const dim_t len = state_len(state); // Size of the array representing the quantum state
    for (dim_t i = 0; i < len; ++i) {
        state->data[i] = prefactor;
    }
}


state_t state_cp(const state_t* state) {
    // Shallow copy of the values
	state_t out = *state;

    // Deep copy of the representation
    const dim_t len = state_len(&out);  // Size of the array representing the quantum state
    if (!((out.data = malloc(len * sizeof(*out.data))))) {
        fprintf(stderr, "state_cp(): out.data allocation failed\n");
        state_free(&out);
        return out;
    }
    cblas_zcopy(len, state->data, 1, out.data, 1);

    return out;
}
