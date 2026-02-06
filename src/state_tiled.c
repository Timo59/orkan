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
 * - Within tiles: row-major element order
 * - Tile size: TILE_DIM x TILE_DIM (default 32x32 = 16KB for complex double)
 *
 * Index formulas:
 * - Tile offset for tile (tr, tc): tr*(tr+1)/2 + tc
 * - Element offset within tile (lr, lc): lr*TILE_DIM + lc
 * - Global index: tile_offset * TILE_SIZE + elem_offset
 *
 * Reference: thesis sec4.tex
 */

#include "state.h"
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
 * @brief Compute tile offset in storage (lower-triangular, row-major tile ordering)
 *
 * Tiles in row tr: T(tr,0), T(tr,1), ..., T(tr,tr)
 * Tiles before row tr: 0 + 1 + ... + (tr-1) = tr*(tr+1)/2
 *
 * @param tr Tile row index
 * @param tc Tile column index (tc <= tr)
 * @return Tile offset (in units of tiles)
 */
static inline dim_t tile_offset(dim_t tr, dim_t tc) {
    return tr * (tr + 1) / 2 + tc;
}

/**
 * @brief Compute element offset within a tile (row-major)
 *
 * @param lr Local row within tile (0 to TILE_DIM-1)
 * @param lc Local column within tile (0 to TILE_DIM-1)
 * @return Element offset within tile
 */
static inline dim_t elem_offset(dim_t lr, dim_t lc) {
    return lr * TILE_DIM + lc;
}

/**
 * @brief Compute absolute storage index for a lower-triangle element
 *
 * For element (row, col) where row >= col:
 * 1. Compute tile coordinates: tr = row / TILE_DIM, tc = col / TILE_DIM
 * 2. Compute local coordinates: lr = row % TILE_DIM, lc = col % TILE_DIM
 * 3. index = tile_offset(tr, tc) * TILE_SIZE + elem_offset(lr, lc)
 *
 * @param row Global row index (row >= col)
 * @param col Global column index
 * @return Absolute storage index
 */
static inline dim_t tiled_index(dim_t row, dim_t col) {
    const dim_t tr = row >> LOG_TILE_DIM;  /* row / TILE_DIM */
    const dim_t tc = col >> LOG_TILE_DIM;  /* col / TILE_DIM */
    const dim_t lr = row & (TILE_DIM - 1); /* row % TILE_DIM */
    const dim_t lc = col & (TILE_DIM - 1); /* col % TILE_DIM */
    return tile_offset(tr, tc) * TILE_SIZE + elem_offset(lr, lc);
}

/*
 * =====================================================================================================================
 * Tiled state implementations
 * =====================================================================================================================
 */

dim_t state_tiled_len(qubit_t qubits) {
    /*
     * Storage length for tiled format: n_tiles*(n_tiles+1)/2 * TILE_SIZE
     *
     * With 64-bit dim_t, memory exhaustion occurs long before overflow.
     * The assert below guards against dimension overflow for the matrix itself.
     */
    assert(qubits < sizeof(dim_t) * 8 - 1);  /* Prevent dim overflow */

    const dim_t dim = POW2(qubits, dim_t);
    const dim_t n_tiles = (dim + TILE_DIM - 1) >> LOG_TILE_DIM;  /* ceil(dim / TILE_DIM) */

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
    const dim_t dim = POW2(state->qubits, dim_t);
    const dim_t len = state_tiled_len(state->qubits);
    const cplx_t prefactor = 1.0 / (double)dim;

    /* Write all elements including tile padding. Simpler than selective writes
     * and acceptable overhead for one-time initialization. */
    for (dim_t i = 0; i < len; ++i) {
        state->data[i] = prefactor;
    }
}


cplx_t state_tiled_get(const state_t *state, dim_t row, dim_t col) {
    const dim_t dim = POW2(state->qubits, dim_t);
    assert(row < dim && col < dim);

    if (row >= col) {
        /* Lower triangle: direct access */
        return state->data[tiled_index(row, col)];
    } else {
        /* Upper triangle: access via conjugate of (col, row) */
        return conj(state->data[tiled_index(col, row)]);
    }
}


void state_tiled_set(state_t *state, dim_t row, dim_t col, cplx_t val) {
    const dim_t dim = POW2(state->qubits, dim_t);
    assert(row < dim && col < dim);

    if (row >= col) {
        /* Lower triangle: direct storage */
        state->data[tiled_index(row, col)] = val;
    } else {
        /* Upper triangle: store conjugate at (col, row) */
        state->data[tiled_index(col, row)] = conj(val);
    }
}


void state_tiled_print(const state_t *state) {
    printf("Type: MIXED_TILED\nQubits: %u\n", state->qubits);
    mprint_tiled(state->data, POW2(state->qubits, dim_t));
}
