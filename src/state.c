/*
 * =====================================================================================================================
 * includes
 * =====================================================================================================================
 */
#ifndef __MATH__
#include <math.h>
#endif

#ifndef STATE_H
#include "state.h"
#endif

#ifndef _STDIO_H_
#include <stdio.h>
#endif

#ifndef _STDLIB_H_
#include <stdlib.h>
#endif

#ifdef MACOS
#include <vecLib/cblas_new.h>
#else
#include <cblas.h>
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


void state_init(state_t *state, const qubit_t qubits, cplx_t *data) {
    // If state has been initialized before, reset it
    if (state->data) {
        free(state->data);
    }

	state->qubits = qubits;

    // Initialize state representation to the passed representation or to all zero if the latter is NULL
    if (data) {
        state->data = data;
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
        prefactor = (cplx_t) sqrt(1 << qubits);
    } else {
        prefactor = (cplx_t) (1 << qubits);
    }

    // Initialize the array representing the quantum state
    const dim_t len = state_len(state); // Size of the array representing the quantum state
    cplx_t *data = malloc(len * sizeof(*state->data));
    if (!data) {
        fprintf(stderr, "state_plus(): data allocation failed\n");
        state_free(state);
        return;
    }
    for (dim_t i = 0; i < len; ++i) {
        data[i] = prefactor;
    }

    // Pass the data array to state_init() to initialize the quantum state
    state_init(state, qubits, data);
}


state_t state_cp(const state_t* state) {
    // Shallow copy of the values
	state_t out = *state;

    // Deep copy of the representation
    const dim_t len = state_len(&out);  // Size of the array representing the quantum state
    if (!((out.data = malloc(len * sizeof(cplx_t))))) {
        fprintf(stderr, "state_cp(): out.data allocation failed\n");
        state_free(&out);
        return out;
    }
    cblas_zcopy(len, state->data, 1, out.data, 1);

    return out;
}
