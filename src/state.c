#include "state.h"
#include "utils.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#if defined(__APPLE__)
    #include <vecLib/cblas_new.h>
#elif defined(__linux__)
    #include <cblas.h>
#endif

/*
 * =====================================================================================================================
 * State manipulation
 * =====================================================================================================================
 */

void state_free(state_t *state) {
    if (!state) return;

    free(state->data);
    state->data = NULL;

    state->qubits = 0;
}


dim_t state_len(const state_t *state) {
    // Guard against uninitialized or invalid state type
    if (state->type != PURE && state->type != MIXED) {
        fprintf(stderr, "state_len(): invalid state type\n");
        return 0;
    }

    // Hilbert space dimension
    const dim_t dim = POW2(state->qubits, dim_t);

    // Hilbert space dimension is the length of the statevector
    if (state->type == PURE) return dim;

    // Length of lower triangle for a density matrix on the Hilbert space
    return dim * (dim + 1) / 2;
}


void state_print(const state_t *state) {
    printf("Type: ");
    if (state->type == PURE) printf("PURE");
    else if (state->type == MIXED) printf("MIXED");
    else {
        printf("INVALID\n");
        return;
    }

    printf("\nQubits: %u\n", state->qubits);

    dim_t dim = POW2(state->qubits, dim_t);

    if (state->type == PURE) vprint(state->data, dim);
    else if (state->type == MIXED) mprint_packed(state->data, dim);
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
        state->data = *data;    // Transfer ownership of state vector
        *data = NULL;           // Nullify source
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
        prefactor = 1.0 / sqrt((double)(1ULL << qubits));
    } else {
        prefactor = 1.0 / (double)(1ULL << qubits);
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
