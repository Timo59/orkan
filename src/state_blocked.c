/**
 * @file state_blocked.c
 * @brief Implementation of blocked state management for density matrices
 */

#include "state_blocked.h"
#include <stdlib.h>
#include <stdio.h>
#include <complex.h>

#if defined(__APPLE__)
    #include <vecLib/cblas_new.h>
#elif defined(__linux__)
    #include <cblas.h>
#endif

/*
 * =====================================================================================================================
 * Index helpers
 * =====================================================================================================================
 */

dim_t tile_offset(dim_t tile_row, dim_t tile_col, dim_t n_blocks) {
    /*
     * Column-major lower triangular tile storage:
     * Column 0: tiles (0,0), (1,0), (2,0), ... (n_blocks-1, 0)  -> n_blocks tiles
     * Column 1: tiles (1,1), (2,1), ... (n_blocks-1, 1)        -> n_blocks-1 tiles
     * ...
     *
     * Offset to column tile_col = sum of column heights before it
     *                           = n_blocks + (n_blocks-1) + ... + (n_blocks - tile_col + 1)
     *                           = tile_col * n_blocks - tile_col*(tile_col-1)/2
     *                           = tile_col * (2*n_blocks - tile_col + 1) / 2
     *
     * Offset within column = tile_row - tile_col
     */
    dim_t col_offset = tile_col * (2 * n_blocks - tile_col + 1) / 2;
    return col_offset + (tile_row - tile_col);
}


cplx_t* get_tile(const state_blocked_t *state, dim_t tile_row, dim_t tile_col) {
    dim_t tile_idx = tile_offset(tile_row, tile_col, state->n_blocks);
    return state->data + tile_idx * BLOCK_DIM * BLOCK_DIM;
}


/*
 * =====================================================================================================================
 * State management
 * =====================================================================================================================
 */

dim_t state_blocked_len(qubit_t qubits) {
    /* Hilbert space dimension */
    const dim_t dim = (dim_t)1 << qubits;

    /* Number of tiles per dimension (rounded up) */
    const dim_t n_blocks = (dim + BLOCK_DIM - 1) / BLOCK_DIM;

    /* Total tiles in lower triangle: n_blocks*(n_blocks+1)/2 */
    const dim_t n_tiles = n_blocks * (n_blocks + 1) / 2;

    /* Each tile is BLOCK_DIM x BLOCK_DIM */
    return n_tiles * BLOCK_DIM * BLOCK_DIM;
}


void state_blocked_init(state_blocked_t *state, qubit_t qubits) {
    /* Free existing data if any */
    if (state->data) {
        free(state->data);
        state->data = NULL;
    }

    state->qubits = qubits;

    /* Hilbert space dimension */
    const dim_t dim = (dim_t)1 << qubits;

    /* Number of tiles per dimension (rounded up) */
    state->n_blocks = (dim + BLOCK_DIM - 1) / BLOCK_DIM;

    /* Padded dimension */
    state->padded_dim = state->n_blocks * BLOCK_DIM;

    /* Allocate and zero-initialize */
    const dim_t len = state_blocked_len(qubits);
    state->data = calloc(len, sizeof(*state->data));
    if (!state->data) {
        fprintf(stderr, "state_blocked_init(): allocation failed\n");
        state->qubits = 0;
        state->n_blocks = 0;
        state->padded_dim = 0;
    }
}


void state_blocked_free(state_blocked_t *state) {
    free(state->data);
    state->data = NULL;
    state->qubits = 0;
    state->n_blocks = 0;
    state->padded_dim = 0;
}


void state_blocked_plus(state_blocked_t *state, qubit_t qubits) {
    /* Initialize empty state */
    state_blocked_init(state, qubits);
    if (!state->data) {
        fprintf(stderr, "state_blocked_plus(): state_blocked_init() failed\n");
        return;
    }

    /* Hilbert space dimension */
    const dim_t dim = (dim_t)1 << qubits;

    /* Normalization: 1/dim for density matrix of |+⟩^n */
    const cplx_t prefactor = 1.0 / dim + 0.0 * I;

    /*
     * Fill all valid elements with prefactor.
     * Valid elements are those where both row and col are < dim.
     * In blocked format, we need to handle padding carefully.
     */
    const dim_t n_blocks = state->n_blocks;

    for (dim_t tb_col = 0; tb_col < n_blocks; ++tb_col) {
        for (dim_t tb_row = tb_col; tb_row < n_blocks; ++tb_row) {
            cplx_t *tile = get_tile(state, tb_row, tb_col);

            /* Determine valid range within this tile */
            const dim_t row_start = tb_row * BLOCK_DIM;
            const dim_t col_start = tb_col * BLOCK_DIM;
            const dim_t row_end = MIN(row_start + BLOCK_DIM, dim);
            const dim_t col_end = MIN(col_start + BLOCK_DIM, dim);

            /* Fill valid elements in row-major order within tile */
            for (dim_t local_row = 0; local_row < BLOCK_DIM; ++local_row) {
                const dim_t global_row = row_start + local_row;
                if (global_row >= dim) break;

                for (dim_t local_col = 0; local_col < BLOCK_DIM; ++local_col) {
                    const dim_t global_col = col_start + local_col;
                    if (global_col >= dim) break;

                    /* Only fill lower triangle (global_row >= global_col) */
                    if (global_row >= global_col) {
                        tile[local_row * BLOCK_DIM + local_col] = prefactor;
                    }
                }
            }
        }
    }
}


/*
 * =====================================================================================================================
 * Element access
 * =====================================================================================================================
 */

cplx_t state_blocked_get(const state_blocked_t *state, dim_t row, dim_t col) {
    /* Handle Hermitian symmetry: upper triangle accessed via conjugate */
    int conjugate = 0;
    if (row < col) {
        dim_t tmp = row;
        row = col;
        col = tmp;
        conjugate = 1;
    }

    /* Find tile coordinates */
    const dim_t tb_row = row / BLOCK_DIM;
    const dim_t tb_col = col / BLOCK_DIM;
    const dim_t local_row = row % BLOCK_DIM;
    const dim_t local_col = col % BLOCK_DIM;

    /* Get tile and read element */
    const cplx_t *tile = get_tile(state, tb_row, tb_col);
    cplx_t val = tile[local_row * BLOCK_DIM + local_col];

    return conjugate ? conj(val) : val;
}


void state_blocked_set(state_blocked_t *state, dim_t row, dim_t col, cplx_t val) {
    /* Only set lower triangle (row >= col) */
    if (row < col) {
        fprintf(stderr, "state_blocked_set(): row must be >= col\n");
        return;
    }

    /* Find tile coordinates */
    const dim_t tb_row = row / BLOCK_DIM;
    const dim_t tb_col = col / BLOCK_DIM;
    const dim_t local_row = row % BLOCK_DIM;
    const dim_t local_col = col % BLOCK_DIM;

    /* Get tile and write element */
    cplx_t *tile = get_tile(state, tb_row, tb_col);
    tile[local_row * BLOCK_DIM + local_col] = val;
}


double state_blocked_prob(const state_blocked_t *state, dim_t basis_state) {
    /* Probability is diagonal element ρ[i,i] */
    cplx_t val = state_blocked_get(state, basis_state, basis_state);
    return creal(val);
}


double state_blocked_trace(const state_blocked_t *state) {
    const dim_t dim = (dim_t)1 << state->qubits;
    double trace = 0.0;

    for (dim_t i = 0; i < dim; ++i) {
        trace += state_blocked_prob(state, i);
    }

    return trace;
}
