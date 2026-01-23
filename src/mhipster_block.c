/**
 * @file mhipster_block.c
 * @brief Quantum gates for blocked density matrices
 *
 * Implements single-qubit gates using tile-based processing for cache efficiency.
 * Each gate transforms ρ → UρU† by processing tiles independently.
 */

#include "mhipster_block.h"
#include <complex.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#if defined(__APPLE__)
    #include <vecLib/cblas_new.h>
#elif defined(__linux__)
    #include <cblas.h>
#endif

/* Minimum dimension to enable OpenMP parallelization */
#define OMP_THRESHOLD 64

/*
 * =====================================================================================================================
 * Helper functions
 * =====================================================================================================================
 */

/**
 * @brief Process Z gate within a single tile
 *
 * Z gate negates elements where the target bit differs between row and column.
 *
 * @param tile Pointer to tile data (BLOCK_DIM x BLOCK_DIM in row-major)
 * @param tb_row Tile row index
 * @param tb_col Tile column index
 * @param incr Bit increment for target qubit (2^target)
 * @param dim Hilbert space dimension
 */
static void z_tile(cplx_t *tile, dim_t tb_row, dim_t tb_col, dim_t incr, dim_t dim) {
    const dim_t row_start = tb_row * BLOCK_DIM;
    const dim_t col_start = tb_col * BLOCK_DIM;
    const dim_t row_end = MIN(row_start + BLOCK_DIM, dim);
    const dim_t col_end = MIN(col_start + BLOCK_DIM, dim);

    /*
     * Z gate: negate element if target bit differs between row and col.
     * For density matrix ρ → ZρZ, element ρ[i,j] is multiplied by:
     *   (+1)(+1) = +1  if both i,j have target bit = 0
     *   (-1)(-1) = +1  if both i,j have target bit = 1
     *   (+1)(-1) = -1  if i has bit 0 and j has bit 1
     *   (-1)(+1) = -1  if i has bit 1 and j has bit 0
     */
    for (dim_t local_row = 0; local_row < BLOCK_DIM; ++local_row) {
        const dim_t global_row = row_start + local_row;
        if (global_row >= dim) break;

        const int row_bit = (global_row & incr) ? 1 : 0;

        for (dim_t local_col = 0; local_col < BLOCK_DIM; ++local_col) {
            const dim_t global_col = col_start + local_col;
            if (global_col >= dim) break;

            /* Only process lower triangle */
            if (global_row < global_col) continue;

            const int col_bit = (global_col & incr) ? 1 : 0;

            /* Negate if bits differ */
            if (row_bit != col_bit) {
                const dim_t idx = local_row * BLOCK_DIM + local_col;
                tile[idx] = -tile[idx];
            }
        }
    }
}


/*
 * =====================================================================================================================
 * Pauli-Z gate: ρ → ZρZ
 * =====================================================================================================================
 */

void z_blocked(state_blocked_t *state, qubit_t target) {
    const dim_t dim = (dim_t)1 << state->qubits;
    const dim_t n_blocks = state->n_blocks;
    const dim_t n_tiles = n_blocks * (n_blocks + 1) / 2;
    const dim_t incr = (dim_t)1 << target;

    const int do_parallel = (dim >= OMP_THRESHOLD) && (n_tiles > 1);

    /*
     * Process all tiles in lower triangle.
     * Using linearized loop for OpenMP compatibility (avoids collapse with triangular space).
     */
    #pragma omp parallel for schedule(dynamic) if(do_parallel)
    for (dim_t t = 0; t < n_tiles; ++t) {
        /*
         * Derive tile coordinates from linear index.
         * Column-major lower triangular: tile t is at (tb_row, tb_col) where
         * tb_col is the largest c such that c*(2*n-c+1)/2 <= t
         */
        dim_t tb_col = 0;
        dim_t col_start_idx = 0;
        while (tb_col < n_blocks) {
            dim_t next_col_start = (tb_col + 1) * (2 * n_blocks - tb_col) / 2;
            if (next_col_start > t) break;
            col_start_idx = next_col_start;
            ++tb_col;
        }
        dim_t tb_row = tb_col + (t - col_start_idx);

        cplx_t *tile = get_tile(state, tb_row, tb_col);
        z_tile(tile, tb_row, tb_col, incr, dim);
    }
}


/*
 * =====================================================================================================================
 * Stub implementations for other gates (to be implemented)
 * =====================================================================================================================
 */

void x_blocked(state_blocked_t *state, qubit_t target) {
    (void)state; (void)target;
    /* TODO: Implement X gate for blocked format */
}

void y_blocked(state_blocked_t *state, qubit_t target) {
    (void)state; (void)target;
    /* TODO: Implement Y gate for blocked format */
}

void h_blocked(state_blocked_t *state, qubit_t target) {
    (void)state; (void)target;
    /* TODO: Implement H gate for blocked format */
}

void s_blocked(state_blocked_t *state, qubit_t target) {
    (void)state; (void)target;
    /* TODO: Implement S gate for blocked format */
}

void sdg_blocked(state_blocked_t *state, qubit_t target) {
    (void)state; (void)target;
    /* TODO: Implement Sdg gate for blocked format */
}

void t_blocked(state_blocked_t *state, qubit_t target) {
    (void)state; (void)target;
    /* TODO: Implement T gate for blocked format */
}

void tdg_blocked(state_blocked_t *state, qubit_t target) {
    (void)state; (void)target;
    /* TODO: Implement Tdg gate for blocked format */
}
