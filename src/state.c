/**
 * @file state.c
 * @brief Dispatch layer for quantum state operations
 *
 * This file implements the public API for quantum states, dispatching to
 * type-specific implementations based on state->type.
 */

#include "state.h"
#include "utils.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#if defined(__APPLE__)
    #include <vecLib/cblas_new.h>
#elif defined(__linux__)
    #include <cblas.h>
#endif

/*
 * =====================================================================================================================
 * Public API - Dispatch implementations
 * =====================================================================================================================
 */

dim_t state_len(const state_t *state) {
    assert(state != NULL);

    switch (state->type) {
        case PURE:
            return state_pure_len(state->qubits);
        case MIXED_PACKED:
            return state_packed_len(state->qubits);
        case MIXED_TILED:
            return state_tiled_len(state->qubits);
        default:
            fprintf(stderr, "state_len(): invalid state type %d\n", state->type);
            return 0;
    }
}


void state_init(state_t *state, const qubit_t qubits, cplx_t **data) {
    assert(state != NULL);

    /* Free existing data if present */
    if (state->data) {
        free(state->data);
        state->data = NULL;
    }

    state->qubits = qubits;

    if (data && *data) {
        /* Transfer ownership */
        state->data = *data;
        *data = NULL;
    } else {
        /* Allocate zero-initialized storage */
        const dim_t len = state_len(state);
        state->data = calloc(len, sizeof(*state->data));
        if (!state->data) {
            fprintf(stderr, "state_init(): allocation failed for %d elements\n", len);
            state->qubits = 0;
        }
    }
}


void state_free(state_t *state) {
    if (!state) return;

    free(state->data);
    state->data = NULL;
    state->qubits = 0;
}


void state_plus(state_t *state, const qubit_t qubits) {
    assert(state != NULL);

    switch (state->type) {
        case PURE:
            state_pure_plus(state, qubits);
            break;
        case MIXED_PACKED:
            state_packed_plus(state, qubits);
            break;
        case MIXED_TILED:
            state_tiled_plus(state, qubits);
            break;
        default:
            fprintf(stderr, "state_plus(): invalid state type %d\n", state->type);
    }
}


state_t state_cp(const state_t *state) {
    assert(state != NULL && state->data != NULL);

    /* Shallow copy of metadata */
    state_t out = *state;

    /* Allocate new storage for the copy */
    const dim_t len = state_len(state);
    out.data = malloc(len * sizeof(*out.data));
    if (!out.data) {
        fprintf(stderr, "state_cp(): allocation failed\n");
        out.qubits = 0;
        return out;
    }

    /* Deep copy of data using BLAS */
    cblas_zcopy(len, state->data, 1, out.data, 1);

    return out;
}


void state_print(const state_t *state) {
    assert(state != NULL);

    switch (state->type) {
        case PURE:
            state_pure_print(state);
            break;
        case MIXED_PACKED:
            state_packed_print(state);
            break;
        case MIXED_TILED:
            state_tiled_print(state);
            break;
        default:
            printf("Type: INVALID (%d)\n", state->type);
    }
}


cplx_t state_get(const state_t *state, dim_t row, dim_t col) {
    assert(state != NULL && state->data != NULL);

    switch (state->type) {
        case PURE:
            return state_pure_get(state, row);
        case MIXED_PACKED:
            return state_packed_get(state, row, col);
        case MIXED_TILED:
            return state_tiled_get(state, row, col);
        default:
            fprintf(stderr, "state_get(): invalid state type %d\n", state->type);
            return 0.0;
    }
}


void state_set(state_t *state, dim_t row, dim_t col, cplx_t val) {
    assert(state != NULL && state->data != NULL);

    switch (state->type) {
        case PURE:
            state_pure_set(state, row, val);
            break;
        case MIXED_PACKED:
            state_packed_set(state, row, col, val);
            break;
        case MIXED_TILED:
            state_tiled_set(state, row, col, val);
            break;
        default:
            fprintf(stderr, "state_set(): invalid state type %d\n", state->type);
    }
}
