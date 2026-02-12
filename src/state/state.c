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
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* 64-byte alignment for AVX-512 and cache line alignment */
#define STATE_ALIGNMENT 64

/**
 * @brief Allocate zero-initialized, aligned memory for state data
 *
 * @param len Number of complex elements
 * @return Aligned, zero-initialized memory, or NULL on failure
 */
static cplx_t *state_alloc_aligned(dim_t len) {
    size_t size = len * sizeof(cplx_t);
    /* Round up to alignment boundary (required by aligned_alloc) */
    size = (size + STATE_ALIGNMENT - 1) & ~(size_t)(STATE_ALIGNMENT - 1);
    cplx_t *ptr = aligned_alloc(STATE_ALIGNMENT, size);
    if (ptr) {
        memset(ptr, 0, size);
    }
    return ptr;
}

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
            assert(0 && "state_len(): invalid state type");
            return 0;  /* unreachable, silences compiler warning */
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
        /* Transfer ownership: if buffer is already aligned, take it directly;
         * otherwise copy into aligned storage and free the original. */
        if (((uintptr_t)*data & (STATE_ALIGNMENT - 1)) == 0) {
            state->data = *data;
        } else {
            const dim_t len = state_len(state);
            state->data = state_alloc_aligned(len);
            if (!state->data) {
                fprintf(stderr, "state_init(): aligned realloc failed for %zu elements\n", (size_t)len);
                state->qubits = 0;
                free(*data);
                *data = NULL;
                return;
            }
            memcpy(state->data, *data, len * sizeof(cplx_t));
            free(*data);
        }
        *data = NULL;
    } else {
        /* Allocate zero-initialized, aligned storage */
        const dim_t len = state_len(state);
        state->data = state_alloc_aligned(len);
        if (!state->data) {
            fprintf(stderr, "state_init(): allocation failed for %zu elements\n", (size_t)len);
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
            assert(0 && "state_plus(): invalid state type");
    }
}


state_t state_cp(const state_t *state) {
    assert(state != NULL && state->data != NULL);

    /* Shallow copy of metadata */
    state_t out = *state;

    /* Allocate aligned storage for the copy */
    const dim_t len = state_len(state);
    out.data = state_alloc_aligned(len);
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
            assert(0 && "state_print(): invalid state type");
    }
}


cplx_t state_get(const state_t *state, dim_t row, dim_t col) {
    assert(state != NULL && state->data != NULL);

    switch (state->type) {
        case PURE:
            assert(col == 0 && "col must be 0 for PURE states");
            return state_pure_get(state, row);
        case MIXED_PACKED:
            return state_packed_get(state, row, col);
        case MIXED_TILED:
            return state_tiled_get(state, row, col);
        default:
            assert(0 && "state_get(): invalid state type");
            return 0.0;  /* unreachable, silences compiler warning */
    }
}


void state_set(state_t *state, dim_t row, dim_t col, cplx_t val) {
    assert(state != NULL && state->data != NULL);

    switch (state->type) {
        case PURE:
            assert(col == 0 && "col must be 0 for PURE states");
            state_pure_set(state, row, val);
            break;
        case MIXED_PACKED:
            state_packed_set(state, row, col, val);
            break;
        case MIXED_TILED:
            state_tiled_set(state, row, col, val);
            break;
        default:
            assert(0 && "state_set(): invalid state type");
    }
}
