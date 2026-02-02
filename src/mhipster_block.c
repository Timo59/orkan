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

/* SIMD support detection */
#if defined(__ARM_NEON) || defined(__ARM_NEON__)
    #include <arm_neon.h>
    #define USE_NEON 1
#elif defined(__AVX2__)
    #include <immintrin.h>
    #define USE_AVX2 1
#endif

/*
 * OpenMP tuning parameters
 */
/* Minimum dimension to enable OpenMP parallelization */
#define OMP_DIM_THRESHOLD 64

/* Minimum number of tiles to justify parallelization overhead */
#define OMP_TILE_THRESHOLD 4

/* Minimum tile rows to enable nested parallelism within tiles */
#define OMP_NESTED_ROW_THRESHOLD 128

/*
 * =====================================================================================================================
 * Helper functions
 * =====================================================================================================================
 */

/**
 * @brief Compute tile coordinates from linear index (O(1) closed-form)
 *
 * For column-major lower triangular tile storage, tile t maps to (row, col) where:
 * - Column c contains tiles c, c+1, ..., n-1 (i.e., n-c tiles)
 * - Tiles 0..n-1 are in col 0, tiles n..2n-2 are in col 1, etc.
 *
 * Given linear index t, solve for col c:
 *   sum_{i=0}^{c-1} (n-i) = cn - c(c-1)/2 <= t
 *   c^2 - (2n+1)c + 2t >= 0
 *   c = floor((2n+1 - sqrt((2n+1)^2 - 8t)) / 2)
 *
 * @param t Linear tile index
 * @param n_blocks Number of tile rows/cols
 * @param[out] tb_row Tile row index
 * @param[out] tb_col Tile column index
 */
static inline void tile_coords_from_index(dim_t t, dim_t n_blocks, dim_t *tb_row, dim_t *tb_col) {
    /* Closed-form solution using quadratic formula */
    const double n = (double)n_blocks;
    const double discriminant = (2.0 * n + 1.0) * (2.0 * n + 1.0) - 8.0 * (double)t;
    const dim_t c = (dim_t)((2.0 * n + 1.0 - sqrt(discriminant)) / 2.0);

    /* Column start index: c * n - c*(c-1)/2 = c * (2n - c + 1) / 2 */
    const dim_t col_start = c * (2 * n_blocks - c + 1) / 2;

    *tb_col = c;
    *tb_row = c + (t - col_start);
}

/**
 * @brief Swap n contiguous complex elements between two arrays using SIMD
 *
 * Swaps a[0:n] with b[0:n]. Uses NEON when available for better throughput.
 */
static inline void swap_run(cplx_t *a, cplx_t *b, dim_t n) {
#if USE_NEON
    double *pa = (double *)a;
    double *pb = (double *)b;
    dim_t total_doubles = n * 2;
    dim_t i = 0;

    /* Process 4 complex numbers (8 doubles) per iteration */
    for (; i + 8 <= total_doubles; i += 8) {
        float64x2_t va0 = vld1q_f64(pa + i);
        float64x2_t va1 = vld1q_f64(pa + i + 2);
        float64x2_t va2 = vld1q_f64(pa + i + 4);
        float64x2_t va3 = vld1q_f64(pa + i + 6);
        float64x2_t vb0 = vld1q_f64(pb + i);
        float64x2_t vb1 = vld1q_f64(pb + i + 2);
        float64x2_t vb2 = vld1q_f64(pb + i + 4);
        float64x2_t vb3 = vld1q_f64(pb + i + 6);
        vst1q_f64(pa + i,     vb0);
        vst1q_f64(pa + i + 2, vb1);
        vst1q_f64(pa + i + 4, vb2);
        vst1q_f64(pa + i + 6, vb3);
        vst1q_f64(pb + i,     va0);
        vst1q_f64(pb + i + 2, va1);
        vst1q_f64(pb + i + 4, va2);
        vst1q_f64(pb + i + 6, va3);
    }

    /* Process 2 complex numbers (4 doubles) */
    for (; i + 4 <= total_doubles; i += 4) {
        float64x2_t va0 = vld1q_f64(pa + i);
        float64x2_t va1 = vld1q_f64(pa + i + 2);
        float64x2_t vb0 = vld1q_f64(pb + i);
        float64x2_t vb1 = vld1q_f64(pb + i + 2);
        vst1q_f64(pa + i,     vb0);
        vst1q_f64(pa + i + 2, vb1);
        vst1q_f64(pb + i,     va0);
        vst1q_f64(pb + i + 2, va1);
    }

    /* Process 1 complex number (2 doubles) */
    for (; i + 2 <= total_doubles; i += 2) {
        float64x2_t va = vld1q_f64(pa + i);
        float64x2_t vb = vld1q_f64(pb + i);
        vst1q_f64(pa + i, vb);
        vst1q_f64(pb + i, va);
    }
#else
    /* Fallback to BLAS */
    cblas_zswap((int)n, a, 1, b, 1);
#endif
}

/**
 * @brief Negate n contiguous complex elements using SIMD when available
 *
 * Each complex number is stored as two doubles (real, imag).
 * NEON processes 2 doubles at a time (1 complex number per vector op).
 * For better throughput, unroll to process 4 complex numbers per iteration.
 */
static inline void negate_run(cplx_t *data, dim_t n) {
#if USE_NEON
    double *ptr = (double *)data;
    dim_t total_doubles = n * 2;
    dim_t i = 0;

    /* Process 4 complex numbers (8 doubles) per iteration for better throughput */
    for (; i + 8 <= total_doubles; i += 8) {
        float64x2_t v0 = vld1q_f64(ptr + i);
        float64x2_t v1 = vld1q_f64(ptr + i + 2);
        float64x2_t v2 = vld1q_f64(ptr + i + 4);
        float64x2_t v3 = vld1q_f64(ptr + i + 6);
        vst1q_f64(ptr + i,     vnegq_f64(v0));
        vst1q_f64(ptr + i + 2, vnegq_f64(v1));
        vst1q_f64(ptr + i + 4, vnegq_f64(v2));
        vst1q_f64(ptr + i + 6, vnegq_f64(v3));
    }

    /* Process 2 complex numbers (4 doubles) */
    for (; i + 4 <= total_doubles; i += 4) {
        float64x2_t v0 = vld1q_f64(ptr + i);
        float64x2_t v1 = vld1q_f64(ptr + i + 2);
        vst1q_f64(ptr + i,     vnegq_f64(v0));
        vst1q_f64(ptr + i + 2, vnegq_f64(v1));
    }

    /* Process 1 complex number (2 doubles) */
    for (; i + 2 <= total_doubles; i += 2) {
        float64x2_t v = vld1q_f64(ptr + i);
        vst1q_f64(ptr + i, vnegq_f64(v));
    }
#else
    /* Fallback to BLAS */
    cblas_zdscal((int)n, -1.0, data, 1);
#endif
}

/**
 * @brief Process a single row of z_tile
 *
 * Extracted for potential nested parallelism in large tiles.
 */
static inline void z_tile_row(cplx_t *tile, dim_t local_row, dim_t row_start, dim_t col_start,
                               dim_t tile_cols, dim_t incr, int is_diagonal_tile) {
    const dim_t global_row = row_start + local_row;
    const int row_bit = (global_row & incr) ? 1 : 0;

    /* Determine column range for this row (lower triangle constraint) */
    dim_t col_limit = tile_cols;
    if (is_diagonal_tile) {
        col_limit = MIN(tile_cols, local_row + 1);
    }

    if (col_limit == 0) return;

    /* Process column runs */
    dim_t local_col = 0;
    while (local_col < col_limit) {
        const dim_t global_col = col_start + local_col;
        const int col_bit = (global_col & incr) ? 1 : 0;

        const dim_t next_bit_boundary = (global_col | (incr - 1)) + 1;
        const dim_t run_end_global = MIN(next_bit_boundary, col_start + col_limit);
        const dim_t run_end_local = run_end_global - col_start;
        const dim_t run_len = run_end_local - local_col;

        if (col_bit != row_bit && run_len > 0) {
            negate_run(&tile[local_row * BLOCK_DIM + local_col], run_len);
        }

        local_col = run_end_local;
    }
}

static void z_tile(cplx_t *tile, dim_t tb_row, dim_t tb_col, dim_t incr, dim_t dim) {
    const dim_t row_start = tb_row * BLOCK_DIM;
    const dim_t col_start = tb_col * BLOCK_DIM;
    const dim_t tile_rows = MIN(BLOCK_DIM, dim - row_start);
    const dim_t tile_cols = MIN(BLOCK_DIM, dim - col_start);
    const int is_diagonal_tile = (tb_row == tb_col);

    /*
     * Z gate: negate element if target bit differs between row and col.
     * For density matrix ρ → ZρZ, element ρ[i,j] is multiplied by:
     *   (+1)(+1) = +1  if both i,j have target bit = 0
     *   (-1)(-1) = +1  if both i,j have target bit = 1
     *   (+1)(-1) = -1  if i has bit 0 and j has bit 1
     *   (-1)(+1) = -1  if i has bit 1 and j has bit 0
     *
     * Optimization: Process columns in runs of size 'incr' where all elements
     * in a run have the same target bit value. Use SIMD for bulk negation.
     *
     * For very large tiles, enable nested parallelism on the row loop.
     */
#ifdef _OPENMP
    const int do_nested = (tile_rows >= OMP_NESTED_ROW_THRESHOLD);
    #pragma omp parallel for if(do_nested)
#endif
    for (dim_t local_row = 0; local_row < tile_rows; ++local_row) {
        z_tile_row(tile, local_row, row_start, col_start, tile_cols, incr, is_diagonal_tile);
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

    /*
     * OpenMP parallelization decision:
     * - Need enough tiles to justify thread overhead
     * - Need enough work per tile (dimension threshold)
     * - Use guided scheduling: starts with large chunks, shrinks for load balancing
     */
    const int do_parallel = (dim >= OMP_DIM_THRESHOLD) && (n_tiles >= OMP_TILE_THRESHOLD);

    /*
     * Process all tiles in lower triangle.
     * Using linearized loop for OpenMP compatibility (avoids collapse with triangular space).
     * Schedule(guided) provides good balance between static efficiency and dynamic load balancing.
     */
    #pragma omp parallel for schedule(guided) if(do_parallel)
    for (dim_t t = 0; t < n_tiles; ++t) {
        /* O(1) tile coordinate calculation using closed-form formula */
        dim_t tb_row, tb_col;
        tile_coords_from_index(t, n_blocks, &tb_row, &tb_col);

        cplx_t *tile = get_tile(state, tb_row, tb_col);
        z_tile(tile, tb_row, tb_col, incr, dim);
    }
}


/*
 * =====================================================================================================================
 * Pauli-X gate: ρ → XρX
 * =====================================================================================================================
 *
 * X = |0⟩⟨1| + |1⟩⟨0| swaps amplitudes where target bit = 0 with target bit = 1.
 * For density matrix: element ρ[i,j] becomes ρ[i^incr, j^incr].
 *
 * Swap pairs:
 *   - (row_bit=0, col_bit=0) ↔ (row_bit=1, col_bit=1): Both in lower triangle
 *   - (row_bit=1, col_bit=0) ↔ (row_bit=0, col_bit=1): May need conjugate access
 */

/**
 * @brief Process X gate within a single tile (for target < log2(BLOCK_DIM))
 *
 * Handles swaps where both elements are within the same tile.
 */
static void x_tile_internal(cplx_t *tile, dim_t tb_row, dim_t tb_col, dim_t incr, dim_t dim) {
    const dim_t row_start = tb_row * BLOCK_DIM;
    const dim_t col_start = tb_col * BLOCK_DIM;
    const dim_t tile_rows = MIN(BLOCK_DIM, dim - row_start);
    const dim_t tile_cols = MIN(BLOCK_DIM, dim - col_start);
    const int is_diagonal_tile = (tb_row == tb_col);

    /*
     * Process (0,0) ↔ (1,1) swaps.
     * For each row with row_bit=0, swap columns with col_bit=0 to their (1,1) partners.
     * SIMD optimization: swap runs of 'incr' elements at a time.
     */
    for (dim_t local_row = 0; local_row < tile_rows; ++local_row) {
        const dim_t global_row = row_start + local_row;
        if (global_row & incr) continue;  /* Only process rows with bit=0 */

        const dim_t partner_local_row = local_row + incr;
        if (partner_local_row >= tile_rows) continue;  /* Partner must be in same tile */

        /* Determine column range */
        dim_t col_limit = tile_cols;
        if (is_diagonal_tile) {
            /* For diagonal tile, only process lower triangle of (0,0) block */
            col_limit = MIN(tile_cols, local_row + 1);
        }

        /* Process column runs where col_bit=0 */
        for (dim_t local_col = 0; local_col < col_limit; local_col += 2 * incr) {
            const dim_t global_col = col_start + local_col;
            if (global_col & incr) continue;

            /* Run length: up to 'incr' elements with col_bit=0 */
            dim_t run_len = MIN(incr, col_limit - local_col);

            /* Swap (local_row, local_col:local_col+run_len) with (partner_local_row, partner_local_col:...) */
            const dim_t partner_local_col = local_col + incr;
            if (partner_local_col >= tile_cols) continue;

            dim_t partner_run_len = MIN(run_len, tile_cols - partner_local_col);
            run_len = MIN(run_len, partner_run_len);

            if (run_len > 0) {
                cplx_t *src = &tile[local_row * BLOCK_DIM + local_col];
                cplx_t *dst = &tile[partner_local_row * BLOCK_DIM + partner_local_col];
                swap_run(src, dst, run_len);
            }
        }
    }

    /*
     * Process (1,0) ↔ (0,1) swaps.
     * Element (row, col) with row_bit=1, col_bit=0 swaps with (row-incr, col+incr).
     * Partner in lower triangle iff (row - incr) >= (col + incr), i.e., row - col >= 2*incr.
     */
    for (dim_t local_row = 0; local_row < tile_rows; ++local_row) {
        const dim_t global_row = row_start + local_row;
        if (!(global_row & incr)) continue;  /* Only process rows with bit=1 */

        const dim_t partner_local_row = local_row - incr;
        /* Partner row must be in same tile (guaranteed since incr < BLOCK_DIM and local_row has bit set) */

        /* Determine column range for lower triangle */
        dim_t col_limit = tile_cols;
        if (is_diagonal_tile) {
            col_limit = MIN(tile_cols, local_row + 1);
        }

        /* Process columns with col_bit=0 */
        for (dim_t local_col = 0; local_col < col_limit; local_col += 2 * incr) {
            const dim_t global_col = col_start + local_col;
            if (global_col & incr) continue;

            /* For each element in this run */
            dim_t run_end = MIN(local_col + incr, col_limit);
            for (dim_t c = local_col; c < run_end; ++c) {
                const dim_t gc = col_start + c;
                const dim_t partner_local_col = c + incr;
                const dim_t partner_global_col = gc + incr;

                if (partner_local_col >= tile_cols) continue;

                /* Check if partner is in lower triangle */
                const dim_t partner_global_row = global_row - incr;
                if (partner_global_row >= partner_global_col) {
                    /* Both in lower triangle - direct swap */
                    cplx_t *src = &tile[local_row * BLOCK_DIM + c];
                    cplx_t *dst = &tile[partner_local_row * BLOCK_DIM + partner_local_col];
                    cplx_t tmp = *src;
                    *src = *dst;
                    *dst = tmp;
                } else {
                    /* Partner in upper triangle - swap with conjugate */
                    /* Upper triangle element (r,c) where r<c is stored as conj at (c,r) */
                    /* But both elements are in the same tile, so we can handle this */
                    /* The (0,1) element at (partner_row, partner_col) with partner_row < partner_col */
                    /* is accessed as conj(tile[partner_col][partner_row]) */
                    /* But partner_col might be outside the diagonal constraint... */

                    /* For diagonal tiles: if partner_global_row < partner_global_col, */
                    /* the (0,1) element is at upper triangle position. */
                    /* We store it at the transposed lower triangle position. */
                    if (is_diagonal_tile) {
                        /* Partner at (partner_local_row, partner_local_col) where row < col */
                        /* Access via (partner_local_col, partner_local_row) with conjugate */
                        cplx_t *src = &tile[local_row * BLOCK_DIM + c];
                        cplx_t *dst = &tile[partner_local_col * BLOCK_DIM + partner_local_row];
                        cplx_t tmp = *src;
                        *src = conj(*dst);
                        *dst = conj(tmp);
                    }
                    /* For non-diagonal tiles, the (0,1) partner would be in a different tile */
                    /* which is handled by the cross-tile case */
                }
            }
        }
    }
}

/**
 * @brief Swap elements between two tiles for X gate (cross-tile case)
 *
 * Swaps elements between tile (tb_row, tb_col) and tile (tb_row ^ tile_incr, tb_col ^ tile_incr).
 */
static void x_swap_tiles(state_blocked_t *state, dim_t tb_row, dim_t tb_col,
                          dim_t partner_tb_row, dim_t partner_tb_col, dim_t dim) {
    const dim_t tile_rows = MIN(BLOCK_DIM, dim - tb_row * BLOCK_DIM);
    const dim_t tile_cols = MIN(BLOCK_DIM, dim - tb_col * BLOCK_DIM);
    const dim_t partner_tile_rows = MIN(BLOCK_DIM, dim - partner_tb_row * BLOCK_DIM);
    const dim_t partner_tile_cols = MIN(BLOCK_DIM, dim - partner_tb_col * BLOCK_DIM);
    const int is_diagonal = (tb_row == tb_col);

    cplx_t *tile = get_tile(state, tb_row, tb_col);
    cplx_t *partner_tile = get_tile(state, partner_tb_row, partner_tb_col);

    /* Swap corresponding elements */
    for (dim_t local_row = 0; local_row < MIN(tile_rows, partner_tile_rows); ++local_row) {
        dim_t col_limit = MIN(tile_cols, partner_tile_cols);
        if (is_diagonal) {
            col_limit = MIN(col_limit, local_row + 1);
        }

        swap_run(&tile[local_row * BLOCK_DIM], &partner_tile[local_row * BLOCK_DIM], col_limit);
    }
}

/**
 * @brief Handle (1,0) ↔ (0,1) cross-tile swaps for X gate
 *
 * For elements where (row_bit=1, col_bit=0) in one tile pair with
 * (row_bit=0, col_bit=1) which may be in a different tile.
 */
static void x_swap_tiles_cross(state_blocked_t *state, dim_t tb_row_10, dim_t tb_col_10,
                                 dim_t tb_row_01, dim_t tb_col_01, dim_t dim, dim_t incr) {
    const dim_t row_start_10 = tb_row_10 * BLOCK_DIM;
    const dim_t col_start_10 = tb_col_10 * BLOCK_DIM;
    const dim_t row_start_01 = tb_row_01 * BLOCK_DIM;
    const dim_t col_start_01 = tb_col_01 * BLOCK_DIM;

    const dim_t tile_rows_10 = MIN(BLOCK_DIM, dim - row_start_10);
    const dim_t tile_cols_10 = MIN(BLOCK_DIM, dim - col_start_10);
    const dim_t tile_rows_01 = MIN(BLOCK_DIM, dim - row_start_01);
    const dim_t tile_cols_01 = MIN(BLOCK_DIM, dim - col_start_01);

    const int is_diag_10 = (tb_row_10 == tb_col_10);
    const int is_diag_01 = (tb_row_01 == tb_col_01);

    cplx_t *tile_10 = get_tile(state, tb_row_10, tb_col_10);

    /* Check if (0,1) tile is in lower triangle of tile grid */
    int tile_01_in_lower = (tb_row_01 >= tb_col_01);

    for (dim_t local_row = 0; local_row < tile_rows_10; ++local_row) {
        const dim_t global_row = row_start_10 + local_row;

        dim_t col_limit = tile_cols_10;
        if (is_diag_10) {
            col_limit = MIN(col_limit, local_row + 1);
        }

        for (dim_t local_col = 0; local_col < col_limit; ++local_col) {
            const dim_t global_col = col_start_10 + local_col;

            /* Partner position */
            const dim_t partner_global_row = global_row ^ incr;
            const dim_t partner_global_col = global_col ^ incr;

            /* Check if partner is in the expected (0,1) tile */
            if (partner_global_row / BLOCK_DIM != tb_row_01) continue;
            if (partner_global_col / BLOCK_DIM != tb_col_01) continue;

            const dim_t partner_local_row = partner_global_row - row_start_01;
            const dim_t partner_local_col = partner_global_col - col_start_01;

            if (partner_local_row >= tile_rows_01 || partner_local_col >= tile_cols_01) continue;

            cplx_t *src = &tile_10[local_row * BLOCK_DIM + local_col];

            if (tile_01_in_lower) {
                /* (0,1) tile is in lower triangle storage */
                if (is_diag_01 && partner_local_row < partner_local_col) {
                    /* Partner element in upper part of diagonal tile - use conjugate */
                    cplx_t *tile_01 = get_tile(state, tb_row_01, tb_col_01);
                    cplx_t *dst = &tile_01[partner_local_col * BLOCK_DIM + partner_local_row];
                    cplx_t tmp = *src;
                    *src = conj(*dst);
                    *dst = conj(tmp);
                } else {
                    cplx_t *tile_01 = get_tile(state, tb_row_01, tb_col_01);
                    cplx_t *dst = &tile_01[partner_local_row * BLOCK_DIM + partner_local_col];
                    cplx_t tmp = *src;
                    *src = *dst;
                    *dst = tmp;
                }
            } else {
                /* (0,1) tile is in upper triangle - access via transposed tile */
                cplx_t *tile_01_trans = get_tile(state, tb_col_01, tb_row_01);
                cplx_t *dst = &tile_01_trans[partner_local_col * BLOCK_DIM + partner_local_row];
                cplx_t tmp = *src;
                *src = conj(*dst);
                *dst = conj(tmp);
            }
        }
    }
}

void x_blocked(state_blocked_t *state, qubit_t target) {
    const dim_t dim = (dim_t)1 << state->qubits;
    const dim_t n_blocks = state->n_blocks;
    const dim_t n_tiles = n_blocks * (n_blocks + 1) / 2;
    const dim_t incr = (dim_t)1 << target;

    const int do_parallel = (dim >= OMP_DIM_THRESHOLD) && (n_tiles >= OMP_TILE_THRESHOLD);

    if (incr < BLOCK_DIM) {
        /*
         * Within-tile swaps: target qubit affects positions within tiles.
         * Each tile can be processed independently.
         */
        #pragma omp parallel for schedule(guided) if(do_parallel)
        for (dim_t t = 0; t < n_tiles; ++t) {
            dim_t tb_row, tb_col;
            tile_coords_from_index(t, n_blocks, &tb_row, &tb_col);

            cplx_t *tile = get_tile(state, tb_row, tb_col);
            x_tile_internal(tile, tb_row, tb_col, incr, dim);
        }
    } else {
        /*
         * Cross-tile swaps: target qubit affects which tile elements are in.
         * tile_incr = incr / BLOCK_DIM = 2^(target - log2(BLOCK_DIM))
         */
        const dim_t tile_incr = incr / BLOCK_DIM;

        /*
         * Process (0,0) ↔ (1,1) tile swaps.
         * Tile (tb_row, tb_col) with both having bit=0 swaps with (tb_row|tile_incr, tb_col|tile_incr).
         */
        #pragma omp parallel for schedule(guided) if(do_parallel)
        for (dim_t tb_row = 0; tb_row < n_blocks; ++tb_row) {
            if (tb_row & tile_incr) continue;  /* Only process tiles with row bit=0 */

            for (dim_t tb_col = 0; tb_col <= tb_row; ++tb_col) {
                if (tb_col & tile_incr) continue;  /* Only process tiles with col bit=0 */

                const dim_t partner_tb_row = tb_row | tile_incr;
                const dim_t partner_tb_col = tb_col | tile_incr;

                if (partner_tb_row >= n_blocks || partner_tb_col >= n_blocks) continue;
                if (partner_tb_row < partner_tb_col) continue;  /* Partner must be in lower triangle */

                x_swap_tiles(state, tb_row, tb_col, partner_tb_row, partner_tb_col, dim);
            }
        }

        /*
         * Process (1,0) ↔ (0,1) tile swaps.
         * Tile (tb_row, tb_col) with row_bit=1, col_bit=0 pairs with tile having row_bit=0, col_bit=1.
         */
        #pragma omp parallel for schedule(guided) if(do_parallel)
        for (dim_t tb_row = 0; tb_row < n_blocks; ++tb_row) {
            if (!(tb_row & tile_incr)) continue;  /* Only process tiles with row bit=1 */

            for (dim_t tb_col = 0; tb_col <= tb_row; ++tb_col) {
                if (tb_col & tile_incr) continue;  /* Only process tiles with col bit=0 */

                const dim_t partner_tb_row = tb_row ^ tile_incr;  /* Clear bit: row - tile_incr */
                const dim_t partner_tb_col = tb_col ^ tile_incr;  /* Set bit: col + tile_incr */

                if (partner_tb_col >= n_blocks) continue;

                x_swap_tiles_cross(state, tb_row, tb_col, partner_tb_row, partner_tb_col, dim, incr);
            }
        }
    }
}

/*
 * =====================================================================================================================
 * Pauli-Y gate: ρ → YρY†
 * =====================================================================================================================
 *
 * Y = [[0, -i], [i, 0]]
 *
 * The transformation for density matrix elements (derived from ρ'[i,j] = Σ Y[i,k]ρ[k,l]Y†[l,j]):
 *   ρ'[i0,j0] = ρ[i1,j1]   (swap with (1,1), NO negate)
 *   ρ'[i0,j1] = -ρ[i1,j0]  (swap with (1,0), WITH negate)
 *   ρ'[i1,j0] = -ρ[i0,j1]  (swap with (0,1), WITH negate)
 *   ρ'[i1,j1] = ρ[i0,j0]   (swap with (0,0), NO negate)
 *
 * So Y = X (swap all pairs) + negate only the (0,1)↔(1,0) off-diagonal blocks.
 */

/**
 * @brief Swap and negate n contiguous complex elements between two arrays
 *
 * After: a[i] = -old_b[i], b[i] = -old_a[i] for i in [0,n)
 */
static inline void swap_negate_run(cplx_t *a, cplx_t *b, dim_t n) {
#if USE_NEON
    double *pa = (double *)a;
    double *pb = (double *)b;
    dim_t total_doubles = n * 2;
    dim_t i = 0;

    /* Process 4 complex numbers (8 doubles) per iteration */
    for (; i + 8 <= total_doubles; i += 8) {
        float64x2_t va0 = vld1q_f64(pa + i);
        float64x2_t va1 = vld1q_f64(pa + i + 2);
        float64x2_t va2 = vld1q_f64(pa + i + 4);
        float64x2_t va3 = vld1q_f64(pa + i + 6);
        float64x2_t vb0 = vld1q_f64(pb + i);
        float64x2_t vb1 = vld1q_f64(pb + i + 2);
        float64x2_t vb2 = vld1q_f64(pb + i + 4);
        float64x2_t vb3 = vld1q_f64(pb + i + 6);
        /* Swap and negate */
        vst1q_f64(pa + i,     vnegq_f64(vb0));
        vst1q_f64(pa + i + 2, vnegq_f64(vb1));
        vst1q_f64(pa + i + 4, vnegq_f64(vb2));
        vst1q_f64(pa + i + 6, vnegq_f64(vb3));
        vst1q_f64(pb + i,     vnegq_f64(va0));
        vst1q_f64(pb + i + 2, vnegq_f64(va1));
        vst1q_f64(pb + i + 4, vnegq_f64(va2));
        vst1q_f64(pb + i + 6, vnegq_f64(va3));
    }

    /* Process 2 complex numbers (4 doubles) */
    for (; i + 4 <= total_doubles; i += 4) {
        float64x2_t va0 = vld1q_f64(pa + i);
        float64x2_t va1 = vld1q_f64(pa + i + 2);
        float64x2_t vb0 = vld1q_f64(pb + i);
        float64x2_t vb1 = vld1q_f64(pb + i + 2);
        vst1q_f64(pa + i,     vnegq_f64(vb0));
        vst1q_f64(pa + i + 2, vnegq_f64(vb1));
        vst1q_f64(pb + i,     vnegq_f64(va0));
        vst1q_f64(pb + i + 2, vnegq_f64(va1));
    }

    /* Process 1 complex number (2 doubles) */
    for (; i + 2 <= total_doubles; i += 2) {
        float64x2_t va = vld1q_f64(pa + i);
        float64x2_t vb = vld1q_f64(pb + i);
        vst1q_f64(pa + i, vnegq_f64(vb));
        vst1q_f64(pb + i, vnegq_f64(va));
    }
#else
    /* Fallback: swap then negate both */
    cblas_zswap((int)n, a, 1, b, 1);
    cblas_zdscal((int)n, -1.0, a, 1);
    cblas_zdscal((int)n, -1.0, b, 1);
#endif
}

/**
 * @brief Process Y gate within a single tile (for target < log2(BLOCK_DIM))
 *
 * Handles swaps where both elements are within the same tile.
 * Similar to X but with negation for (0,1)↔(1,0) swaps (off-diagonal blocks).
 */
static void y_tile_internal(cplx_t *tile, dim_t tb_row, dim_t tb_col, dim_t incr, dim_t dim) {
    const dim_t row_start = tb_row * BLOCK_DIM;
    const dim_t col_start = tb_col * BLOCK_DIM;
    const dim_t tile_rows = MIN(BLOCK_DIM, dim - row_start);
    const dim_t tile_cols = MIN(BLOCK_DIM, dim - col_start);
    const int is_diagonal_tile = (tb_row == tb_col);

    /*
     * Process (0,0) ↔ (1,1) swaps WITHOUT negation.
     * For each row with row_bit=0, swap columns with col_bit=0 to their (1,1) partners.
     */
    for (dim_t local_row = 0; local_row < tile_rows; ++local_row) {
        const dim_t global_row = row_start + local_row;
        if (global_row & incr) continue;  /* Only process rows with bit=0 */

        const dim_t partner_local_row = local_row + incr;
        if (partner_local_row >= tile_rows) continue;

        /* Determine column range */
        dim_t col_limit = tile_cols;
        if (is_diagonal_tile) {
            col_limit = MIN(tile_cols, local_row + 1);
        }

        /* Process column runs where col_bit=0 */
        for (dim_t local_col = 0; local_col < col_limit; local_col += 2 * incr) {
            const dim_t global_col = col_start + local_col;
            if (global_col & incr) continue;

            dim_t run_len = MIN(incr, col_limit - local_col);
            const dim_t partner_local_col = local_col + incr;
            if (partner_local_col >= tile_cols) continue;

            dim_t partner_run_len = MIN(run_len, tile_cols - partner_local_col);
            run_len = MIN(run_len, partner_run_len);

            if (run_len > 0) {
                cplx_t *src = &tile[local_row * BLOCK_DIM + local_col];
                cplx_t *dst = &tile[partner_local_row * BLOCK_DIM + partner_local_col];
                swap_run(src, dst, run_len);  /* Swap WITHOUT negate */
            }
        }
    }

    /*
     * Process (1,0) ↔ (0,1) swaps WITH negation.
     * Similar to X gate but with negation added.
     *
     * Special case: When the partner's upper triangle position maps to the same
     * storage location as the source (due to Hermitian storage), we just apply
     * negate+conjugate in place rather than swapping.
     */
    for (dim_t local_row = 0; local_row < tile_rows; ++local_row) {
        const dim_t global_row = row_start + local_row;
        if (!(global_row & incr)) continue;  /* Only process rows with bit=1 */

        const dim_t partner_local_row = local_row - incr;

        dim_t col_limit = tile_cols;
        if (is_diagonal_tile) {
            col_limit = MIN(tile_cols, local_row + 1);
        }

        /* Process columns with col_bit=0 */
        for (dim_t local_col = 0; local_col < col_limit; local_col += 2 * incr) {
            const dim_t global_col = col_start + local_col;
            if (global_col & incr) continue;

            dim_t run_end = MIN(local_col + incr, col_limit);
            for (dim_t c = local_col; c < run_end; ++c) {
                const dim_t gc = col_start + c;
                const dim_t partner_local_col = c + incr;
                const dim_t partner_global_col = gc + incr;

                if (partner_local_col >= tile_cols) continue;

                const dim_t partner_global_row = global_row - incr;
                if (partner_global_row >= partner_global_col) {
                    /* Both in lower triangle - swap with negate */
                    cplx_t *src = &tile[local_row * BLOCK_DIM + c];
                    cplx_t *dst = &tile[partner_local_row * BLOCK_DIM + partner_local_col];
                    cplx_t tmp = *src;
                    *src = -(*dst);
                    *dst = -tmp;
                } else if (is_diagonal_tile) {
                    /*
                     * Partner in upper triangle of diagonal tile.
                     * Upper (partner_global_row, partner_global_col) is stored as
                     * conj(lower(partner_global_col, partner_global_row)).
                     *
                     * Partner storage position: (partner_local_col, partner_local_row)
                     *
                     * Case 1: Same storage (partner_local_col, partner_local_row) == (local_row, c)
                     *         Just negate+conjugate in place.
                     *
                     * Case 2: Different storage, but partner storage may also be a (1,0) position
                     *         that will try to swap with us. To avoid double-swap, only process
                     *         if our position comes "before" partner storage position.
                     */
                    const dim_t partner_storage_row = partner_local_col;
                    const dim_t partner_storage_col = partner_local_row;

                    if (local_row == partner_storage_row && c == partner_storage_col) {
                        /* Same storage location - apply negate+conjugate in place */
                        cplx_t *elem = &tile[local_row * BLOCK_DIM + c];
                        *elem = -conj(*elem);
                    } else {
                        /*
                         * Check if partner storage (partner_storage_row, partner_storage_col)
                         * is also a (1,0) position that will process us.
                         * A position (r,c) is (1,0) if (r & incr) && !(c & incr).
                         * Only process if we're the "canonical" side of the cycle.
                         */
                        const dim_t psr = row_start + partner_storage_row;
                        const dim_t psc = col_start + partner_storage_col;
                        const int partner_is_10 = (psr & incr) && !(psc & incr);

                        if (partner_is_10) {
                            /* Both form a 2-cycle. Only process if source < partner storage. */
                            const dim_t src_idx = global_row * dim + (col_start + c);
                            const dim_t pst_idx = psr * dim + psc;
                            if (src_idx > pst_idx) {
                                /* Partner will handle this swap */
                                continue;
                            }
                        }

                        /* Different storage locations - swap with conjugate AND negate */
                        cplx_t *src = &tile[local_row * BLOCK_DIM + c];
                        cplx_t *dst = &tile[partner_storage_row * BLOCK_DIM + partner_storage_col];
                        cplx_t tmp = *src;
                        *src = -conj(*dst);
                        *dst = -conj(tmp);
                    }
                }
            }
        }
    }
}

/**
 * @brief Swap elements between two tiles for Y gate (cross-tile case)
 *
 * For (0,0)↔(1,1) swaps that span different tiles - NO negation.
 */
static void y_swap_tiles(state_blocked_t *state, dim_t tb_row, dim_t tb_col,
                          dim_t partner_tb_row, dim_t partner_tb_col, dim_t dim) {
    const dim_t tile_rows = MIN(BLOCK_DIM, dim - tb_row * BLOCK_DIM);
    const dim_t tile_cols = MIN(BLOCK_DIM, dim - tb_col * BLOCK_DIM);
    const dim_t partner_tile_rows = MIN(BLOCK_DIM, dim - partner_tb_row * BLOCK_DIM);
    const dim_t partner_tile_cols = MIN(BLOCK_DIM, dim - partner_tb_col * BLOCK_DIM);
    const int is_diagonal = (tb_row == tb_col);

    cplx_t *tile = get_tile(state, tb_row, tb_col);
    cplx_t *partner_tile = get_tile(state, partner_tb_row, partner_tb_col);

    /* Swap corresponding elements (no negation) */
    for (dim_t local_row = 0; local_row < MIN(tile_rows, partner_tile_rows); ++local_row) {
        dim_t col_limit = MIN(tile_cols, partner_tile_cols);
        if (is_diagonal) {
            col_limit = MIN(col_limit, local_row + 1);
        }

        swap_run(&tile[local_row * BLOCK_DIM], &partner_tile[local_row * BLOCK_DIM], col_limit);
    }
}

/**
 * @brief Handle (1,0) ↔ (0,1) cross-tile swaps for Y gate (WITH negation)
 *
 * Similar to X gate's cross-tile swap but with negation added.
 */
static void y_swap_negate_tiles_cross(state_blocked_t *state, dim_t tb_row_10, dim_t tb_col_10,
                                       dim_t tb_row_01, dim_t tb_col_01, dim_t dim, dim_t incr) {
    const dim_t row_start_10 = tb_row_10 * BLOCK_DIM;
    const dim_t col_start_10 = tb_col_10 * BLOCK_DIM;
    const dim_t row_start_01 = tb_row_01 * BLOCK_DIM;
    const dim_t col_start_01 = tb_col_01 * BLOCK_DIM;

    const dim_t tile_rows_10 = MIN(BLOCK_DIM, dim - row_start_10);
    const dim_t tile_cols_10 = MIN(BLOCK_DIM, dim - col_start_10);
    const dim_t tile_rows_01 = MIN(BLOCK_DIM, dim - row_start_01);
    const dim_t tile_cols_01 = MIN(BLOCK_DIM, dim - col_start_01);

    const int is_diag_10 = (tb_row_10 == tb_col_10);
    const int is_diag_01 = (tb_row_01 == tb_col_01);

    cplx_t *tile_10 = get_tile(state, tb_row_10, tb_col_10);

    int tile_01_in_lower = (tb_row_01 >= tb_col_01);

    for (dim_t local_row = 0; local_row < tile_rows_10; ++local_row) {
        const dim_t global_row = row_start_10 + local_row;

        dim_t col_limit = tile_cols_10;
        if (is_diag_10) {
            col_limit = MIN(col_limit, local_row + 1);
        }

        for (dim_t local_col = 0; local_col < col_limit; ++local_col) {
            const dim_t global_col = col_start_10 + local_col;

            const dim_t partner_global_row = global_row ^ incr;
            const dim_t partner_global_col = global_col ^ incr;

            if (partner_global_row / BLOCK_DIM != tb_row_01) continue;
            if (partner_global_col / BLOCK_DIM != tb_col_01) continue;

            const dim_t partner_local_row = partner_global_row - row_start_01;
            const dim_t partner_local_col = partner_global_col - col_start_01;

            if (partner_local_row >= tile_rows_01 || partner_local_col >= tile_cols_01) continue;

            cplx_t *src = &tile_10[local_row * BLOCK_DIM + local_col];

            if (tile_01_in_lower) {
                if (is_diag_01 && partner_local_row < partner_local_col) {
                    /* Partner in upper triangle - swap with conjugate AND negate */
                    cplx_t *tile_01 = get_tile(state, tb_row_01, tb_col_01);
                    cplx_t *dst = &tile_01[partner_local_col * BLOCK_DIM + partner_local_row];
                    cplx_t tmp = *src;
                    *src = -conj(*dst);
                    *dst = -conj(tmp);
                } else {
                    /* Both in lower triangle - swap with negate */
                    cplx_t *tile_01 = get_tile(state, tb_row_01, tb_col_01);
                    cplx_t *dst = &tile_01[partner_local_row * BLOCK_DIM + partner_local_col];
                    cplx_t tmp = *src;
                    *src = -(*dst);
                    *dst = -tmp;
                }
            } else {
                /* Upper triangle tile - swap with conjugate AND negate */
                cplx_t *tile_01_trans = get_tile(state, tb_col_01, tb_row_01);
                cplx_t *dst = &tile_01_trans[partner_local_col * BLOCK_DIM + partner_local_row];
                cplx_t tmp = *src;
                *src = -conj(*dst);
                *dst = -conj(tmp);
            }
        }
    }
}

void y_blocked(state_blocked_t *state, qubit_t target) {
    const dim_t dim = (dim_t)1 << state->qubits;
    const dim_t n_blocks = state->n_blocks;
    const dim_t n_tiles = n_blocks * (n_blocks + 1) / 2;
    const dim_t incr = (dim_t)1 << target;

    const int do_parallel = (dim >= OMP_DIM_THRESHOLD) && (n_tiles >= OMP_TILE_THRESHOLD);

    if (incr < BLOCK_DIM) {
        /*
         * Within-tile case: target qubit affects positions within tiles.
         * Each tile can be processed independently.
         */
        #pragma omp parallel for schedule(guided) if(do_parallel)
        for (dim_t t = 0; t < n_tiles; ++t) {
            dim_t tb_row, tb_col;
            tile_coords_from_index(t, n_blocks, &tb_row, &tb_col);

            cplx_t *tile = get_tile(state, tb_row, tb_col);
            y_tile_internal(tile, tb_row, tb_col, incr, dim);
        }
    } else {
        /*
         * Cross-tile case: target qubit affects which tile elements are in.
         */
        const dim_t tile_incr = incr / BLOCK_DIM;

        /*
         * Process (0,0) ↔ (1,1) tile swaps WITHOUT negation.
         */
        #pragma omp parallel for schedule(guided) if(do_parallel)
        for (dim_t tb_row = 0; tb_row < n_blocks; ++tb_row) {
            if (tb_row & tile_incr) continue;

            for (dim_t tb_col = 0; tb_col <= tb_row; ++tb_col) {
                if (tb_col & tile_incr) continue;

                const dim_t partner_tb_row = tb_row | tile_incr;
                const dim_t partner_tb_col = tb_col | tile_incr;

                if (partner_tb_row >= n_blocks || partner_tb_col >= n_blocks) continue;
                if (partner_tb_row < partner_tb_col) continue;

                y_swap_tiles(state, tb_row, tb_col, partner_tb_row, partner_tb_col, dim);
            }
        }

        /*
         * Process (1,0) ↔ (0,1) tile swaps WITH negation.
         */
        #pragma omp parallel for schedule(guided) if(do_parallel)
        for (dim_t tb_row = 0; tb_row < n_blocks; ++tb_row) {
            if (!(tb_row & tile_incr)) continue;

            for (dim_t tb_col = 0; tb_col <= tb_row; ++tb_col) {
                if (tb_col & tile_incr) continue;

                const dim_t partner_tb_row = tb_row ^ tile_incr;
                const dim_t partner_tb_col = tb_col ^ tile_incr;

                if (partner_tb_col >= n_blocks) continue;

                y_swap_negate_tiles_cross(state, tb_row, tb_col, partner_tb_row, partner_tb_col, dim, incr);
            }
        }
    }
}

/*
 * =====================================================================================================================
 * Hadamard gate: ρ → HρH†
 * =====================================================================================================================
 *
 * H = 1/√2 * [[1, 1], [1, -1]]  (Hermitian: H† = H)
 *
 * The transformation for each 2x2 block is:
 *   ρ'[i0,j0] = (ρ[i0,j0] + ρ[i0,j1] + ρ[i1,j0] + ρ[i1,j1]) / 2
 *   ρ'[i0,j1] = (ρ[i0,j0] - ρ[i0,j1] + ρ[i1,j0] - ρ[i1,j1]) / 2
 *   ρ'[i1,j0] = (ρ[i0,j0] + ρ[i0,j1] - ρ[i1,j0] - ρ[i1,j1]) / 2
 *   ρ'[i1,j1] = (ρ[i0,j0] - ρ[i0,j1] - ρ[i1,j0] + ρ[i1,j1]) / 2
 *
 * where i0,i1 differ by target bit, and j0,j1 differ by target bit.
 */

/**
 * @brief Apply 2x2 Hadamard butterfly to four elements in place
 *
 * Scalar version for single butterfly.
 */
static inline void h_butterfly_4(cplx_t *v00, cplx_t *v01, cplx_t *v10, cplx_t *v11) {
    const cplx_t a = *v00;
    const cplx_t b = *v01;
    const cplx_t c = *v10;
    const cplx_t d = *v11;

    *v00 = 0.5 * (a + b + c + d);
    *v01 = 0.5 * (a - b + c - d);
    *v10 = 0.5 * (a + b - c - d);
    *v11 = 0.5 * (a - b - c + d);
}

#if USE_NEON
/**
 * @brief Apply H butterfly to n contiguous butterflies using NEON SIMD
 *
 * Processes multiple butterflies where elements are contiguous in memory.
 * Each butterfly processes 4 elements: (row0, col), (row0, col+incr), (row1, col), (row1, col+incr)
 * where row1 = row0 + incr.
 *
 * For incr=1 with contiguous columns, we can vectorize across columns.
 * Memory layout in row-major tile:
 *   row0: [v00_0, v01_0, v00_1, v01_1, v00_2, v01_2, ...]  (v00 and v01 interleaved)
 *   row1: [v10_0, v11_0, v10_1, v11_1, v10_2, v11_2, ...]  (v10 and v11 interleaved)
 *
 * @param row0_ptr Pointer to start of row0 data (n*2 complex values)
 * @param row1_ptr Pointer to start of row1 data (n*2 complex values)
 * @param n Number of butterflies to process
 */
static inline void h_butterfly_run_simd(cplx_t *row0_ptr, cplx_t *row1_ptr, dim_t n) {
    double *r0 = (double *)row0_ptr;
    double *r1 = (double *)row1_ptr;

    /* Each complex is 2 doubles, each butterfly needs 2 complex from each row */
    /* So 4 doubles from row0 (v00, v01) and 4 doubles from row1 (v10, v11) per butterfly */
    const dim_t total_doubles = n * 4;  /* n butterflies * 2 complex * 2 doubles */

    /* NEON processes 2 doubles at a time = 1 complex number */
    /* Process 2 butterflies (8 doubles from each row) per iteration for better throughput */
    dim_t i = 0;

    /* Main loop: 2 butterflies at a time */
    for (; i + 8 <= total_doubles; i += 8) {
        /* Load 4 complex from row0: v00_0, v01_0, v00_1, v01_1 */
        float64x2_t r0_0 = vld1q_f64(r0 + i);      /* v00_0 */
        float64x2_t r0_1 = vld1q_f64(r0 + i + 2);  /* v01_0 */
        float64x2_t r0_2 = vld1q_f64(r0 + i + 4);  /* v00_1 */
        float64x2_t r0_3 = vld1q_f64(r0 + i + 6);  /* v01_1 */

        /* Load 4 complex from row1: v10_0, v11_0, v10_1, v11_1 */
        float64x2_t r1_0 = vld1q_f64(r1 + i);      /* v10_0 */
        float64x2_t r1_1 = vld1q_f64(r1 + i + 2);  /* v11_0 */
        float64x2_t r1_2 = vld1q_f64(r1 + i + 4);  /* v10_1 */
        float64x2_t r1_3 = vld1q_f64(r1 + i + 6);  /* v11_1 */

        /* Compute a+b, a-b, c+d, c-d for each butterfly */
        /* Butterfly 0 */
        float64x2_t ab_sum_0 = vaddq_f64(r0_0, r0_1);   /* a+b */
        float64x2_t ab_diff_0 = vsubq_f64(r0_0, r0_1);  /* a-b */
        float64x2_t cd_sum_0 = vaddq_f64(r1_0, r1_1);   /* c+d */
        float64x2_t cd_diff_0 = vsubq_f64(r1_0, r1_1);  /* c-d */

        /* Butterfly 1 */
        float64x2_t ab_sum_1 = vaddq_f64(r0_2, r0_3);
        float64x2_t ab_diff_1 = vsubq_f64(r0_2, r0_3);
        float64x2_t cd_sum_1 = vaddq_f64(r1_2, r1_3);
        float64x2_t cd_diff_1 = vsubq_f64(r1_2, r1_3);

        /* Scale factor 0.5 */
        float64x2_t half = vdupq_n_f64(0.5);

        /* Compute outputs: v00' = 0.5*(a+b+c+d), v01' = 0.5*(a-b+c-d), etc. */
        /* Butterfly 0 */
        float64x2_t out00_0 = vmulq_f64(vaddq_f64(ab_sum_0, cd_sum_0), half);    /* 0.5*(a+b+c+d) */
        float64x2_t out01_0 = vmulq_f64(vaddq_f64(ab_diff_0, cd_diff_0), half);  /* 0.5*(a-b+c-d) */
        float64x2_t out10_0 = vmulq_f64(vsubq_f64(ab_sum_0, cd_sum_0), half);    /* 0.5*(a+b-c-d) */
        float64x2_t out11_0 = vmulq_f64(vsubq_f64(ab_diff_0, cd_diff_0), half);  /* 0.5*(a-b-c+d) */

        /* Butterfly 1 */
        float64x2_t out00_1 = vmulq_f64(vaddq_f64(ab_sum_1, cd_sum_1), half);
        float64x2_t out01_1 = vmulq_f64(vaddq_f64(ab_diff_1, cd_diff_1), half);
        float64x2_t out10_1 = vmulq_f64(vsubq_f64(ab_sum_1, cd_sum_1), half);
        float64x2_t out11_1 = vmulq_f64(vsubq_f64(ab_diff_1, cd_diff_1), half);

        /* Store results */
        vst1q_f64(r0 + i, out00_0);
        vst1q_f64(r0 + i + 2, out01_0);
        vst1q_f64(r0 + i + 4, out00_1);
        vst1q_f64(r0 + i + 6, out01_1);

        vst1q_f64(r1 + i, out10_0);
        vst1q_f64(r1 + i + 2, out11_0);
        vst1q_f64(r1 + i + 4, out10_1);
        vst1q_f64(r1 + i + 6, out11_1);
    }

    /* Handle remaining butterflies (1 at a time) */
    for (; i < total_doubles; i += 4) {
        float64x2_t v00 = vld1q_f64(r0 + i);
        float64x2_t v01 = vld1q_f64(r0 + i + 2);
        float64x2_t v10 = vld1q_f64(r1 + i);
        float64x2_t v11 = vld1q_f64(r1 + i + 2);

        float64x2_t ab_sum = vaddq_f64(v00, v01);
        float64x2_t ab_diff = vsubq_f64(v00, v01);
        float64x2_t cd_sum = vaddq_f64(v10, v11);
        float64x2_t cd_diff = vsubq_f64(v10, v11);

        float64x2_t half = vdupq_n_f64(0.5);

        vst1q_f64(r0 + i, vmulq_f64(vaddq_f64(ab_sum, cd_sum), half));
        vst1q_f64(r0 + i + 2, vmulq_f64(vaddq_f64(ab_diff, cd_diff), half));
        vst1q_f64(r1 + i, vmulq_f64(vsubq_f64(ab_sum, cd_sum), half));
        vst1q_f64(r1 + i + 2, vmulq_f64(vsubq_f64(ab_diff, cd_diff), half));
    }
}
#endif

/**
 * @brief Process H gate within a single tile (for target < log2(BLOCK_DIM))
 *
 * Handles all 2x2 butterfly transformations where elements are within the same tile.
 * Uses Hermitian symmetry: only processes lower triangle, reads upper via conjugate.
 *
 * Key insight: For each unique 2x2 butterfly block defined by (row0, col0) where both
 * have target bit = 0, we only need to process it once. We choose to process each
 * butterfly when row0 >= col0 (the (0,0) element is in lower triangle).
 */
static void h_tile_internal(cplx_t *tile, dim_t tb_row, dim_t tb_col, dim_t incr, dim_t dim) {
    const dim_t row_start = tb_row * BLOCK_DIM;
    const dim_t col_start = tb_col * BLOCK_DIM;
    const dim_t tile_rows = MIN(BLOCK_DIM, dim - row_start);
    const dim_t tile_cols = MIN(BLOCK_DIM, dim - col_start);
    const int is_diagonal_tile = (tb_row == tb_col);

#if USE_NEON
    /*
     * Fast path for incr=1: Butterflies are contiguous in memory.
     * For non-diagonal tiles, we can process entire rows with SIMD.
     * For diagonal tiles, we need to handle the lower triangle constraint.
     */
    if (incr == 1 && !is_diagonal_tile) {
        /* Non-diagonal tile with incr=1: process rows in pairs */
        for (dim_t local_row = 0; local_row < tile_rows; local_row += 2) {
            if (local_row + 1 >= tile_rows) break;

            /* Number of butterflies = tile_cols / 2 (pairs of columns) */
            dim_t n_butterflies = tile_cols / 2;
            if (n_butterflies > 0) {
                cplx_t *row0_ptr = &tile[local_row * BLOCK_DIM];
                cplx_t *row1_ptr = &tile[(local_row + 1) * BLOCK_DIM];
                h_butterfly_run_simd(row0_ptr, row1_ptr, n_butterflies);
            }
        }
        return;
    }

    if (incr == 1 && is_diagonal_tile) {
        /* Diagonal tile with incr=1: partially SIMD, with special handling for diagonal */
        for (dim_t local_row = 0; local_row < tile_rows; local_row += 2) {
            if (local_row + 1 >= tile_rows) break;

            const dim_t global_row = row_start + local_row;
            const dim_t global_row1 = global_row + 1;

            /* For diagonal tile, columns must satisfy col <= row for lower triangle */
            /* For row pair (local_row, local_row+1), max col in lower triangle is local_row */
            /* But we need col pairs (col, col+1), so max col0 is local_row - 1 (if local_row is odd)
             * or local_row (if local_row is even) */

            /* Actually, for the (row, row+1) pair:
             * - (row, col) is lower triangle iff row >= col
             * - (row, col+1) is lower triangle iff row >= col+1, i.e., col <= row-1
             * - (row+1, col) is lower triangle iff row+1 >= col (always true if col <= row)
             * - (row+1, col+1) is lower triangle iff row+1 >= col+1, i.e., col <= row
             *
             * For SIMD we need all 4 elements in lower triangle, so col <= row-1.
             * After that, we handle col = row separately (upper triangle access needed).
             */

            /* SIMD part: col ranges from 0 to local_row-1 (exclusive), in pairs */
            dim_t n_full_butterflies = local_row / 2;  /* Number of complete butterfly pairs */
            if (n_full_butterflies > 0) {
                cplx_t *row0_ptr = &tile[local_row * BLOCK_DIM];
                cplx_t *row1_ptr = &tile[(local_row + 1) * BLOCK_DIM];
                h_butterfly_run_simd(row0_ptr, row1_ptr, n_full_butterflies);
            }

            /* Handle remaining columns that need Hermitian access */
            dim_t remaining_start = n_full_butterflies * 2;
            for (dim_t local_col = remaining_start; local_col <= local_row; local_col += 2) {
                if (local_col + 1 >= tile_cols) break;

                const dim_t local_col1 = local_col + 1;
                const dim_t global_col = col_start + local_col;
                const dim_t global_col1 = global_col + 1;

                cplx_t v00, v01, v10, v11;

                v00 = tile[local_row * BLOCK_DIM + local_col];

                if (global_row >= global_col1) {
                    v01 = tile[local_row * BLOCK_DIM + local_col1];
                } else {
                    v01 = conj(tile[local_col1 * BLOCK_DIM + local_row]);
                }

                v10 = tile[(local_row + 1) * BLOCK_DIM + local_col];

                if (global_row1 >= global_col1) {
                    v11 = tile[(local_row + 1) * BLOCK_DIM + local_col1];
                } else {
                    v11 = conj(tile[local_col1 * BLOCK_DIM + (local_row + 1)]);
                }

                h_butterfly_4(&v00, &v01, &v10, &v11);

                tile[local_row * BLOCK_DIM + local_col] = v00;

                if (global_row >= global_col1) {
                    tile[local_row * BLOCK_DIM + local_col1] = v01;
                } else {
                    tile[local_col1 * BLOCK_DIM + local_row] = conj(v01);
                }

                tile[(local_row + 1) * BLOCK_DIM + local_col] = v10;

                if (global_row1 >= global_col1) {
                    tile[(local_row + 1) * BLOCK_DIM + local_col1] = v11;
                } else {
                    tile[local_col1 * BLOCK_DIM + (local_row + 1)] = conj(v11);
                }
            }
        }
        return;
    }
#endif

    /*
     * General case: Process each unique 2x2 block individually.
     * Iterate over (row0, col0) where both have bit=0 AND row0 >= col0.
     */
    for (dim_t local_row = 0; local_row < tile_rows; ++local_row) {
        const dim_t global_row = row_start + local_row;

        /* Skip if row has target bit = 1 */
        if (global_row & incr) continue;

        /* Partner row (with target bit = 1) */
        const dim_t local_row1 = local_row + incr;
        if (local_row1 >= tile_rows) continue;
        const dim_t global_row1 = global_row + incr;

        /* Determine column range based on lower triangle */
        dim_t col_end = tile_cols;
        if (is_diagonal_tile) {
            /* For diagonal tiles, limit to lower triangle: col <= row */
            col_end = MIN(tile_cols, local_row + 1);
        }

        for (dim_t local_col = 0; local_col < col_end; ++local_col) {
            const dim_t global_col = col_start + local_col;

            /* Skip if col has target bit = 1 */
            if (global_col & incr) continue;

            /* Partner column (with target bit = 1) */
            const dim_t local_col1 = local_col + incr;
            if (local_col1 >= tile_cols) continue;
            const dim_t global_col1 = global_col + incr;

            /* Now we process the butterfly for (row, col), (row, col1), (row1, col), (row1, col1) */
            /* where row has bit=0, row1 has bit=1, col has bit=0, col1 has bit=1 */

            /* Read all 4 elements, handling Hermitian symmetry */
            cplx_t v00, v01, v10, v11;

            /* (row, col): in lower triangle since row >= col (ensured by loop bounds) */
            v00 = tile[local_row * BLOCK_DIM + local_col];

            /* (row, col1): lower triangle iff row >= col1 */
            if (global_row >= global_col1) {
                v01 = tile[local_row * BLOCK_DIM + local_col1];
            } else if (is_diagonal_tile) {
                /* Upper triangle in same tile: access (col1, row) with conjugate */
                v01 = conj(tile[local_col1 * BLOCK_DIM + local_row]);
            } else {
                /* Upper triangle in different tile - this shouldn't happen for within-tile case */
                v01 = conj(tile[local_col1 * BLOCK_DIM + local_row]);
            }

            /* (row1, col): lower triangle iff row1 >= col, which is true since row1 > row >= col */
            v10 = tile[local_row1 * BLOCK_DIM + local_col];

            /* (row1, col1): lower triangle iff row1 >= col1 */
            if (global_row1 >= global_col1) {
                v11 = tile[local_row1 * BLOCK_DIM + local_col1];
            } else if (is_diagonal_tile) {
                /* Upper triangle: access (col1, row1) with conjugate */
                v11 = conj(tile[local_col1 * BLOCK_DIM + local_row1]);
            } else {
                v11 = conj(tile[local_col1 * BLOCK_DIM + local_row1]);
            }

            /* Apply butterfly transformation */
            h_butterfly_4(&v00, &v01, &v10, &v11);

            /* Write back to lower triangle positions only */

            /* (row, col): always in lower triangle */
            tile[local_row * BLOCK_DIM + local_col] = v00;

            /* (row, col1): write directly if in lower triangle, else write conjugate to (col1, row) */
            if (global_row >= global_col1) {
                tile[local_row * BLOCK_DIM + local_col1] = v01;
            } else if (is_diagonal_tile) {
                tile[local_col1 * BLOCK_DIM + local_row] = conj(v01);
            }

            /* (row1, col): always in lower triangle since row1 > row >= col */
            tile[local_row1 * BLOCK_DIM + local_col] = v10;

            /* (row1, col1): write directly if in lower triangle, else write conjugate to (col1, row1) */
            if (global_row1 >= global_col1) {
                tile[local_row1 * BLOCK_DIM + local_col1] = v11;
            } else if (is_diagonal_tile) {
                tile[local_col1 * BLOCK_DIM + local_row1] = conj(v11);
            }
        }
    }
}

/**
 * @brief Compute group coordinates from linear index (for parallelization)
 *
 * Groups are enumerated in column-major order within the lower triangle.
 * Similar to tile_coords_from_index but for the reduced group grid.
 *
 * @param g Linear group index
 * @param n_base Number of base rows/cols (those with tile_bit = 0)
 * @param[out] br Base row index
 * @param[out] bc Base col index
 */
static inline void group_coords_from_index(dim_t g, dim_t n_base, dim_t *br, dim_t *bc) {
    /* Same formula as tile_coords_from_index but for the group grid */
    const double n = (double)n_base;
    const double discriminant = (2.0 * n + 1.0) * (2.0 * n + 1.0) - 8.0 * (double)g;
    const dim_t c = (dim_t)((2.0 * n + 1.0 - sqrt(discriminant)) / 2.0);
    const dim_t col_start = c * (2 * n_base - c + 1) / 2;

    *bc = c;
    *br = c + (g - col_start);
}

#if USE_NEON
/**
 * @brief Process a row of cross-tile H butterflies using NEON SIMD
 *
 * Processes multiple butterflies where each butterfly spans 4 tiles.
 * Elements at position (lr, lc) in each of the 4 tiles form one butterfly.
 *
 * @param row0_00 Pointer to row lr in tile_00
 * @param row0_01 Pointer to row lr (or col lr if transposed) in tile_01
 * @param row0_10 Pointer to row lr in tile_10
 * @param row0_11 Pointer to row lr in tile_11
 * @param n Number of elements (butterflies) to process
 * @param is_01_transposed Whether tile_01 access is transposed (upper triangle)
 * @param is_11_upper_diag Whether writing to tile_11 needs upper-diagonal handling
 * @param lr Current local row (needed for diagonal constraint)
 */
static void h_cross_tile_row_simd(cplx_t *row_00, cplx_t *tile_01, dim_t tile_01_stride,
                                   cplx_t *row_10, cplx_t *row_11, dim_t n,
                                   int is_01_transposed, int is_11_upper_diag, dim_t lr) {
    double *r00 = (double *)row_00;
    double *r10 = (double *)row_10;
    double *r11 = (double *)row_11;

    float64x2_t half = vdupq_n_f64(0.5);

    dim_t lc = 0;

    /* Process 2 butterflies at a time when possible */
    for (; lc + 2 <= n; lc += 2) {
        /* Check diagonal constraint for tile_11 */
        int skip_11_0 = is_11_upper_diag && (lr < lc);
        int skip_11_1 = is_11_upper_diag && (lr < lc + 1);

        /* Load from tile_00 and tile_10 (contiguous) */
        float64x2_t v00_0 = vld1q_f64(r00 + lc * 2);
        float64x2_t v00_1 = vld1q_f64(r00 + lc * 2 + 2);
        float64x2_t v10_0 = vld1q_f64(r10 + lc * 2);
        float64x2_t v10_1 = vld1q_f64(r10 + lc * 2 + 2);

        /* Load from tile_01 (may be transposed) */
        float64x2_t v01_0, v01_1;
        if (is_01_transposed) {
            /* Access tile_01[lc][lr] and tile_01[lc+1][lr] with conjugate */
            double *t01_0 = (double *)&tile_01[lc * tile_01_stride + lr];
            double *t01_1 = (double *)&tile_01[(lc + 1) * tile_01_stride + lr];
            v01_0 = vld1q_f64(t01_0);
            v01_1 = vld1q_f64(t01_1);
            /* Conjugate: negate imaginary part (index 1) */
            v01_0 = vsetq_lane_f64(-vgetq_lane_f64(v01_0, 1), v01_0, 1);
            v01_1 = vsetq_lane_f64(-vgetq_lane_f64(v01_1, 1), v01_1, 1);
        } else {
            /* Direct access: tile_01[lr][lc] */
            double *t01 = (double *)&tile_01[lr * tile_01_stride];
            v01_0 = vld1q_f64(t01 + lc * 2);
            v01_1 = vld1q_f64(t01 + lc * 2 + 2);
        }

        /* Load from tile_11 */
        float64x2_t v11_0 = vld1q_f64(r11 + lc * 2);
        float64x2_t v11_1 = vld1q_f64(r11 + lc * 2 + 2);

        /* Butterfly 0 */
        float64x2_t ab_sum_0 = vaddq_f64(v00_0, v01_0);
        float64x2_t ab_diff_0 = vsubq_f64(v00_0, v01_0);
        float64x2_t cd_sum_0 = vaddq_f64(v10_0, v11_0);
        float64x2_t cd_diff_0 = vsubq_f64(v10_0, v11_0);

        float64x2_t out00_0 = vmulq_f64(vaddq_f64(ab_sum_0, cd_sum_0), half);
        float64x2_t out01_0 = vmulq_f64(vaddq_f64(ab_diff_0, cd_diff_0), half);
        float64x2_t out10_0 = vmulq_f64(vsubq_f64(ab_sum_0, cd_sum_0), half);
        float64x2_t out11_0 = vmulq_f64(vsubq_f64(ab_diff_0, cd_diff_0), half);

        /* Butterfly 1 */
        float64x2_t ab_sum_1 = vaddq_f64(v00_1, v01_1);
        float64x2_t ab_diff_1 = vsubq_f64(v00_1, v01_1);
        float64x2_t cd_sum_1 = vaddq_f64(v10_1, v11_1);
        float64x2_t cd_diff_1 = vsubq_f64(v10_1, v11_1);

        float64x2_t out00_1 = vmulq_f64(vaddq_f64(ab_sum_1, cd_sum_1), half);
        float64x2_t out01_1 = vmulq_f64(vaddq_f64(ab_diff_1, cd_diff_1), half);
        float64x2_t out10_1 = vmulq_f64(vsubq_f64(ab_sum_1, cd_sum_1), half);
        float64x2_t out11_1 = vmulq_f64(vsubq_f64(ab_diff_1, cd_diff_1), half);

        /* Store to tile_00 and tile_10 */
        vst1q_f64(r00 + lc * 2, out00_0);
        vst1q_f64(r00 + lc * 2 + 2, out00_1);
        vst1q_f64(r10 + lc * 2, out10_0);
        vst1q_f64(r10 + lc * 2 + 2, out10_1);

        /* Store to tile_01 */
        if (is_01_transposed) {
            /* Conjugate before storing */
            out01_0 = vsetq_lane_f64(-vgetq_lane_f64(out01_0, 1), out01_0, 1);
            out01_1 = vsetq_lane_f64(-vgetq_lane_f64(out01_1, 1), out01_1, 1);
            double *t01_0 = (double *)&tile_01[lc * tile_01_stride + lr];
            double *t01_1 = (double *)&tile_01[(lc + 1) * tile_01_stride + lr];
            vst1q_f64(t01_0, out01_0);
            vst1q_f64(t01_1, out01_1);
        } else {
            double *t01 = (double *)&tile_01[lr * tile_01_stride];
            vst1q_f64(t01 + lc * 2, out01_0);
            vst1q_f64(t01 + lc * 2 + 2, out01_1);
        }

        /* Store to tile_11 with diagonal handling */
        if (!skip_11_0) {
            vst1q_f64(r11 + lc * 2, out11_0);
        } else {
            /* Store conjugate at transposed position */
            out11_0 = vsetq_lane_f64(-vgetq_lane_f64(out11_0, 1), out11_0, 1);
            double *t11_trans = (double *)&((cplx_t*)r11)[lc * BLOCK_DIM + lr - lc * BLOCK_DIM];
            /* Recalculate: we need tile_11[lc][lr] */
            cplx_t *tile_11_base = (cplx_t *)(r11 - lr * 2);  /* Back to tile start */
            double *dst = (double *)&tile_11_base[lc * BLOCK_DIM + lr];
            vst1q_f64(dst, out11_0);
        }
        if (!skip_11_1) {
            vst1q_f64(r11 + lc * 2 + 2, out11_1);
        } else {
            out11_1 = vsetq_lane_f64(-vgetq_lane_f64(out11_1, 1), out11_1, 1);
            cplx_t *tile_11_base = (cplx_t *)(r11 - lr * 2);
            double *dst = (double *)&tile_11_base[(lc + 1) * BLOCK_DIM + lr];
            vst1q_f64(dst, out11_1);
        }
    }

    /* Handle remaining element */
    if (lc < n) {
        cplx_t v00 = row_00[lc];
        cplx_t v10 = row_10[lc];
        cplx_t v11 = row_11[lc];

        cplx_t v01;
        if (is_01_transposed) {
            v01 = conj(tile_01[lc * tile_01_stride + lr]);
        } else {
            v01 = tile_01[lr * tile_01_stride + lc];
        }

        h_butterfly_4(&v00, &v01, &v10, &v11);

        row_00[lc] = v00;
        row_10[lc] = v10;

        if (is_01_transposed) {
            tile_01[lc * tile_01_stride + lr] = conj(v01);
        } else {
            tile_01[lr * tile_01_stride + lc] = v01;
        }

        if (!is_11_upper_diag || lr >= lc) {
            row_11[lc] = v11;
        } else {
            cplx_t *tile_11_base = row_11 - lr * BLOCK_DIM;
            tile_11_base[lc * BLOCK_DIM + lr] = conj(v11);
        }
    }
}
#endif /* USE_NEON */

/**
 * @brief Process one tile group for cross-tile H gate
 *
 * A tile group consists of 4 tiles that share butterflies:
 * - (base_row, base_col) - the (0,0) tile
 * - (base_row, base_col | tile_incr) - the (0,1) tile
 * - (base_row | tile_incr, base_col) - the (1,0) tile
 * - (base_row | tile_incr, base_col | tile_incr) - the (1,1) tile
 *
 * Elements at the same (local_row, local_col) across these 4 tiles form one butterfly.
 */
static void h_process_tile_group(state_blocked_t *state,
                                  dim_t base_row, dim_t base_col,
                                  dim_t tile_incr, dim_t dim) {
    const dim_t n_blocks = state->n_blocks;

    /* Compute the 4 tile coordinates */
    const dim_t tb_row_00 = base_row;
    const dim_t tb_col_00 = base_col;
    const dim_t tb_row_01 = base_row;
    const dim_t tb_col_01 = base_col | tile_incr;
    const dim_t tb_row_10 = base_row | tile_incr;
    const dim_t tb_col_10 = base_col;
    const dim_t tb_row_11 = base_row | tile_incr;
    const dim_t tb_col_11 = base_col | tile_incr;

    /* Check bounds */
    if (tb_row_11 >= n_blocks || tb_col_11 >= n_blocks) return;
    if (tb_row_11 < tb_col_11) return;  /* (1,1) tile must be in lower triangle */

    /* Compute tile dimensions */
    const dim_t row_start_00 = tb_row_00 * BLOCK_DIM;
    const dim_t col_start_00 = tb_col_00 * BLOCK_DIM;
    const dim_t tile_rows = MIN(BLOCK_DIM, dim - row_start_00);
    const dim_t tile_cols = MIN(BLOCK_DIM, dim - col_start_00);

    const dim_t row_start_11 = tb_row_11 * BLOCK_DIM;
    const dim_t col_start_11 = tb_col_11 * BLOCK_DIM;
    const dim_t tile_rows_11 = MIN(BLOCK_DIM, dim - row_start_11);
    const dim_t tile_cols_11 = MIN(BLOCK_DIM, dim - col_start_11);

    /* Check if tiles are diagonal (need special handling for Hermitian storage) */
    const int is_diagonal_00 = (tb_row_00 == tb_col_00);
    const int is_diagonal_11 = (tb_row_11 == tb_col_11);

    /* Get tile pointers */
    cplx_t *tile_00 = get_tile(state, tb_row_00, tb_col_00);
    cplx_t *tile_11 = get_tile(state, tb_row_11, tb_col_11);

    /* Handle (0,1) tile - may be in upper triangle of tile grid */
    const int tile_01_in_lower = (tb_row_01 >= tb_col_01);
    cplx_t *tile_01;
    dim_t tile_01_stride;
    int is_01_transposed;

    if (tile_01_in_lower) {
        tile_01 = get_tile(state, tb_row_01, tb_col_01);
        tile_01_stride = BLOCK_DIM;
        is_01_transposed = 0;
    } else {
        /* Access via transposed tile with conjugate */
        tile_01 = get_tile(state, tb_col_01, tb_row_01);
        tile_01_stride = BLOCK_DIM;
        is_01_transposed = 1;
    }

    /* (1,0) tile is always in lower triangle since tb_row_10 > tb_row_00 >= tb_col_00 = tb_col_10 */
    cplx_t *tile_10 = get_tile(state, tb_row_10, tb_col_10);

    /* Process all elements in this tile group */
    const dim_t max_rows = MIN(tile_rows, tile_rows_11);
    const dim_t max_cols = MIN(tile_cols, tile_cols_11);

#if USE_NEON
    /* Use SIMD-optimized row processing */
    for (dim_t lr = 0; lr < max_rows; ++lr) {
        dim_t col_limit = max_cols;
        if (is_diagonal_00) col_limit = MIN(col_limit, lr + 1);

        if (col_limit > 0) {
            h_cross_tile_row_simd(&tile_00[lr * BLOCK_DIM],
                                   tile_01, tile_01_stride,
                                   &tile_10[lr * BLOCK_DIM],
                                   &tile_11[lr * BLOCK_DIM],
                                   col_limit,
                                   is_01_transposed,
                                   is_diagonal_11,
                                   lr);
        }
    }
#else
    /* Scalar fallback */
    for (dim_t lr = 0; lr < max_rows; ++lr) {
        dim_t col_limit = max_cols;
        if (is_diagonal_00) col_limit = MIN(col_limit, lr + 1);

        for (dim_t lc = 0; lc < col_limit; ++lc) {
            cplx_t v00 = tile_00[lr * BLOCK_DIM + lc];
            cplx_t v11 = tile_11[lr * BLOCK_DIM + lc];

            /* Handle (0,1) element */
            cplx_t v01;
            if (is_01_transposed) {
                v01 = conj(tile_01[lc * BLOCK_DIM + lr]);
            } else {
                v01 = tile_01[lr * BLOCK_DIM + lc];
            }

            /* Handle (1,0) element - always in lower triangle */
            cplx_t v10 = tile_10[lr * BLOCK_DIM + lc];

            /* Apply butterfly */
            h_butterfly_4(&v00, &v01, &v10, &v11);

            /* Write back */
            tile_00[lr * BLOCK_DIM + lc] = v00;

            if (is_01_transposed) {
                tile_01[lc * BLOCK_DIM + lr] = conj(v01);
            } else {
                tile_01[lr * BLOCK_DIM + lc] = v01;
            }

            tile_10[lr * BLOCK_DIM + lc] = v10;

            /* Handle diagonal constraint for (1,1) tile */
            if (!is_diagonal_11 || lr >= lc) {
                tile_11[lr * BLOCK_DIM + lc] = v11;
            } else {
                /* (1,1) element is in upper part of diagonal tile */
                tile_11[lc * BLOCK_DIM + lr] = conj(v11);
            }
        }
    }
#endif
}

/**
 * @brief Process H gate for cross-tile case (target >= log2(BLOCK_DIM))
 *
 * Uses tile-group parallelization: tiles partition into independent groups of 4,
 * where each group can be processed in parallel without data races.
 *
 * Key insight: Base tiles are those where (tile_index & tile_incr) == 0.
 * Each base tile (base_row, base_col) forms a group with:
 * - (base_row, base_col | tile_incr)
 * - (base_row | tile_incr, base_col)
 * - (base_row | tile_incr, base_col | tile_incr)
 */
static void h_cross_tile(state_blocked_t *state, dim_t incr, dim_t dim) {
    const dim_t n_blocks = state->n_blocks;
    const dim_t tile_incr = incr / BLOCK_DIM;

    /*
     * Count valid base tiles.
     * A base tile (tb_row, tb_col) must have:
     * 1. (tb_row & tile_incr) == 0  (bit clear for row)
     * 2. (tb_col & tile_incr) == 0  (bit clear for col)
     * 3. (tb_row | tile_incr) < n_blocks  (partner row exists)
     * 4. (tb_col | tile_incr) < n_blocks  (partner col exists)
     * 5. tb_row >= tb_col  (lower triangle of base tiles)
     *
     * For each dimension: tiles with bit clear are 0, 1, ..., tile_incr-1,
     * then 2*tile_incr, 2*tile_incr+1, ..., 3*tile_incr-1, etc.
     *
     * Count base tiles per dimension that have valid partners:
     */
    dim_t n_base_rows = 0;
    for (dim_t tb = 0; tb < n_blocks; ++tb) {
        if ((tb & tile_incr) == 0 && (tb | tile_incr) < n_blocks) {
            n_base_rows++;
        }
    }

    /* Quick exit if no valid groups */
    if (n_base_rows == 0) return;

    /*
     * Build array of valid base tile indices for efficient parallel iteration.
     * This allows us to use a simple parallel for without complex index math.
     */
    dim_t *base_indices = malloc(n_base_rows * sizeof(dim_t));
    if (!base_indices) {
        /* Fallback to sequential if allocation fails */
        for (dim_t base_row = 0; base_row < n_blocks; ++base_row) {
            if (base_row & tile_incr) continue;
            if ((base_row | tile_incr) >= n_blocks) continue;

            for (dim_t base_col = 0; base_col <= base_row; ++base_col) {
                if (base_col & tile_incr) continue;
                if ((base_col | tile_incr) >= n_blocks) continue;

                h_process_tile_group(state, base_row, base_col, tile_incr, dim);
            }
        }
        return;
    }

    /* Fill base indices array */
    dim_t idx = 0;
    for (dim_t tb = 0; tb < n_blocks; ++tb) {
        if ((tb & tile_incr) == 0 && (tb | tile_incr) < n_blocks) {
            base_indices[idx++] = tb;
        }
    }

    /* Total number of groups in lower triangle */
    const dim_t n_groups = n_base_rows * (n_base_rows + 1) / 2;
    const int do_parallel = (dim >= OMP_DIM_THRESHOLD) && (n_groups >= OMP_TILE_THRESHOLD);

    /*
     * Parallelize over independent tile groups.
     * Each group processes 4 tiles that don't overlap with any other group.
     */
    #pragma omp parallel for schedule(guided) if(do_parallel)
    for (dim_t g = 0; g < n_groups; ++g) {
        /* Convert linear index to (br, bc) in lower triangular base grid */
        dim_t br, bc;
        group_coords_from_index(g, n_base_rows, &br, &bc);

        /* Map to actual tile coordinates */
        const dim_t base_row = base_indices[br];
        const dim_t base_col = base_indices[bc];

        /* Process this tile group */
        h_process_tile_group(state, base_row, base_col, tile_incr, dim);
    }

    free(base_indices);
}

void h_blocked(state_blocked_t *state, qubit_t target) {
    const dim_t dim = (dim_t)1 << state->qubits;
    const dim_t n_blocks = state->n_blocks;
    const dim_t n_tiles = n_blocks * (n_blocks + 1) / 2;
    const dim_t incr = (dim_t)1 << target;

    const int do_parallel = (dim >= OMP_DIM_THRESHOLD) && (n_tiles >= OMP_TILE_THRESHOLD);

    if (incr < BLOCK_DIM) {
        /*
         * Within-tile case: all 4 elements of each butterfly are in the same tile.
         * Process each tile independently.
         */
        #pragma omp parallel for schedule(guided) if(do_parallel)
        for (dim_t t = 0; t < n_tiles; ++t) {
            dim_t tb_row, tb_col;
            tile_coords_from_index(t, n_blocks, &tb_row, &tb_col);

            cplx_t *tile = get_tile(state, tb_row, tb_col);
            h_tile_internal(tile, tb_row, tb_col, incr, dim);
        }
    } else {
        /*
         * Cross-tile case: butterflies span multiple tiles.
         * Need to coordinate access across tiles.
         *
         * For simplicity, process sequentially. OpenMP parallelization would
         * require careful partitioning to avoid data races.
         */
        h_cross_tile(state, incr, dim);
    }
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
