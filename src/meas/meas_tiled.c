/**
 * @file meas_tiled.c
 * @brief Expectation value of a diagonal observable for tiled mixed states
 *
 * Computes  <H> = sum_x obs[x] rho_{xx}
 *
 * Diagonal elements live on diagonal tiles (k, k).  Within each tile the
 * diagonal is strided by TILE_DIM + 1 (row-major intra-tile layout).
 *
 * Diagonal tile positions (in TILE_SIZE units): 0, 2, 5, 9, 14, ...
 * Stride between consecutive diagonal tiles: 2, 3, 4, 5, ...
 * (starts at 2, increments by 1 each iteration).
 *
 * Each thread processes a contiguous chunk of diagonal tiles, computing its
 * starting tile offset once and then using the stride-increment pattern.
 * Within each tile, diagonal elements are gathered to a contiguous stack array
 * to enable auto-vectorisation of the dot product.
 */

#include "meas.h"
#include "index.h"

#include <complex.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef OMP_THRESHOLD
#define OMP_THRESHOLD 512
#endif

double mean_tiled(const state_t *state, const double *obs) {
    const idx_t dim = POW2(state->qubits, idx_t);
    const idx_t dim_tile = (dim + TILE_DIM - 1) >> LOG_TILE_DIM;
    const idx_t len_tile = (dim < TILE_DIM) ? dim : TILE_DIM;
    double sum = 0.0;

    #pragma omp parallel reduction(+:sum) if(dim >= OMP_THRESHOLD)
    {
#ifdef _OPENMP
        int tid = omp_get_thread_num();
        int nthreads = omp_get_num_threads();
#else
        int tid = 0;
        int nthreads = 1;
#endif
        idx_t chunk = dim_tile / nthreads;
        idx_t k_start = tid * chunk;
        idx_t k_end = (tid == nthreads - 1) ? dim_tile : k_start + chunk;

        idx_t tile_base = tile_off(k_start, k_start);
        idx_t tile_stride = k_start + 2;

        for (idx_t k = k_start; k < k_end; ++k) {
            /* Gather diagonal elements from tile to contiguous stack array */
            double diag[TILE_DIM];
            for (idx_t j = 0; j < len_tile; ++j) {
                diag[j] = creal(state->data[tile_base + j * (TILE_DIM + 1)]);
            }

            /* Vectorisable contiguous multiply-accumulate */
            const idx_t x_start = k * len_tile;
            for (idx_t j = 0; j < len_tile; ++j) {
                sum += obs[x_start + j] * diag[j];
            }

            tile_base += tile_stride * TILE_SIZE;
            ++tile_stride;
        }
    }
    return sum;
}
