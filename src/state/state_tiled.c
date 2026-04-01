/**
 * @file state_tiled.c
 * @brief Tiled storage for mixed quantum states (density matrices)
 *
 * This implementation stores density matrices in a tiled (blocked) format
 * optimized for cache locality during quantum gate operations. The key insight
 * is that for within-tile operations (target qubit < LOG_TILE_DIM), no conjugation
 * is needed because all 4 butterfly elements are directly accessible within the tile.
 *
 * Memory layout:
 * - Density matrix partitioned into TILE_DIM x TILE_DIM tiles
 * - Only lower-triangular tiles stored (at tile level)
 * - Tiles stored in row-major order: T(0,0), T(1,0), T(1,1), T(2,0), ...
 * - Within tiles: row-major element order (NOT column-major like LAPACK packed)
 * - Tile size: TILE_DIM x TILE_DIM (default 32x32 = 16KB for complex double)
 *
 * Index formulas:
 * - Tile offset for tile (tr, tc): tr*(tr+1)/2 + tc
 * - Element offset within tile (lr, lc): lr*TILE_DIM + lc
 * - Global index: tile_offset * TILE_SIZE + elem_offset
 *
 * WARNING: Intra-tile layout is row-major. Do NOT pass tile pointers to LAPACK
 * routines (zhpmv, zhpr, etc.) which expect column-major packed storage.
 *
 * Reference: thesis sec4.tex
 */

#include "state.h"
#include "index.h"
#include "utils.h"

#include <assert.h>
#include <complex.h>
#include <stdio.h>

/*
 * =====================================================================================================================
 * Helper functions
 * =====================================================================================================================
 */

/**
 * @brief Compute absolute storage index for a lower-triangle element
 *
 * For element (row, col) where row >= col:
 * 1. Compute tile coordinates: tr = row / TILE_DIM, tc = col / TILE_DIM
 * 2. Compute local coordinates: lr = row % TILE_DIM, lc = col % TILE_DIM
 * 3. index = tile_off(tr, tc) + elem_off(lr, lc)
 *
 * @param row Global row index (row >= col)
 * @param col Global column index
 * @return Absolute storage index
 */
static inline idx_t tiled_index(idx_t row, idx_t col) {
    const idx_t tr = row >> LOG_TILE_DIM;
    const idx_t tc = col >> LOG_TILE_DIM;
    const idx_t lr = row & (TILE_DIM - 1);
    const idx_t lc = col & (TILE_DIM - 1);
    return tile_off(tr, tc) + elem_off(lr, lc);
}

/*
 * =====================================================================================================================
 * Tiled state implementations
 * =====================================================================================================================
 */

idx_t state_tiled_len(qubit_t qubits) {
    /*
     * Storage length for tiled format: n_tiles*(n_tiles+1)/2 * TILE_SIZE
     *
     * With 64-bit idx_t, memory exhaustion occurs long before overflow.
     * The assert below guards against dimension overflow for the matrix itself.
     */
    assert(qubits < sizeof(idx_t) * 8 - 1);  /* Prevent dim overflow */

    const idx_t dim = POW2(qubits, idx_t);
    const idx_t n_tiles = (dim + TILE_DIM - 1) >> LOG_TILE_DIM;  /* ceil(dim / TILE_DIM) */

    return (n_tiles * (n_tiles + 1) >> 1) * TILE_SIZE;
}


void state_tiled_plus(state_t *state, qubit_t qubits) {
    /* Initialize empty state (allocates zero-initialized storage) */
    state_init(state, qubits, NULL);
    if (!state->data) return;

    /*
     * The plus state |+>^n = (1/sqrt(2^n)) * sum_i |i> has density matrix:
     *   rho = |+><+| = (1/2^n) * sum_{i,j} |i><j|
     *
     * All matrix elements are equal to 1/2^n.
     * This is consistent with state_packed_plus.
     */
    const idx_t dim = POW2(state->qubits, idx_t);
    const idx_t len = state_tiled_len(state->qubits);
    const cplx_t prefactor = 1.0 / (double)dim;

    if (state->qubits < LOG_TILE_DIM) {
        /* Single tile with padding: only write the physical dim x dim corner.
         * Padding stays zero from state_init's memset. */
        for (idx_t row = 0; row < dim; ++row) {
            for (idx_t col = 0; col < dim; ++col) {
                state->data[row * TILE_DIM + col] = prefactor;
            }
        }
    } else {
        /* dim is a multiple of TILE_DIM, so no padding exists. Bulk fill. */
        for (idx_t i = 0; i < len; ++i) {
            state->data[i] = prefactor;
        }
    }
}


cplx_t state_tiled_get(const state_t *state, idx_t row, idx_t col) {
    const idx_t dim = POW2(state->qubits, idx_t);
    assert(row < dim && col < dim);

    if (row >= col) {
        /* Lower triangle: direct access */
        return state->data[tiled_index(row, col)];
    } else {
        /* Upper triangle: access via conjugate of (col, row) */
        return conj(state->data[tiled_index(col, row)]);
    }
}


void state_tiled_set(state_t *state, idx_t row, idx_t col, cplx_t val) {
    const idx_t dim = POW2(state->qubits, idx_t);
    assert(row < dim && col < dim);

    if (row >= col) {
        /* Lower triangle: direct storage */
        state->data[tiled_index(row, col)] = val;
        /* Diagonal tiles store full TILE_DIM × TILE_DIM: mirror to upper triangle */
        if ((row >> LOG_TILE_DIM) == (col >> LOG_TILE_DIM) && row != col) {
            state->data[tiled_index(col, row)] = conj(val);
        }
    } else {
        /* Upper triangle: store conjugate at transposed lower-triangle position */
        state->data[tiled_index(col, row)] = conj(val);
        /* Diagonal tiles: also store directly at (row, col) */
        if ((row >> LOG_TILE_DIM) == (col >> LOG_TILE_DIM)) {
            state->data[tiled_index(row, col)] = val;
        }
    }
}


void state_tiled_print(const state_t *state) {
    printf("Type: MIXED_TILED\nQubits: %u\n", state->qubits);
    mprint_tiled(state->data, POW2(state->qubits, idx_t));
}
