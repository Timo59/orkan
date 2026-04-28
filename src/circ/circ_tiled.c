/**
 * @file circ_tiled.c
 * @brief Circuit operations for tiled mixed states
 *
 * exp_diag:  rho[r][c] *= exp(-i (diag[r] - diag[c]) t)
 *
 * Same precomputation strategy as circ_packed.c.  Iteration follows the tiled
 * lower-triangular layout: tile pairs (tr, tc) with tr >= tc.  Diagonal tiles
 * only touch the strict lower triangle; off-diagonal tiles touch all elements.
 */

#include "circ.h"
#include "index.h"
#include "utils.h"

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef _OPENMP
#include <omp.h>
#endif

void exp_diag_tiled(state_t *state, const double *restrict diag, double t) {
    const idx_t dim      = POW2(state->qubits, idx_t);
    const idx_t dim_tile = (dim + TILE_DIM - 1) >> LOG_TILE_DIM;
    const idx_t len_tile = (dim < TILE_DIM) ? dim : TILE_DIM;
    /* dim is always a power of two, so tiles are never ragged.
     * When dim < TILE_DIM: single tile, len_tile == dim, dim_tile == 1.
     * When dim >= TILE_DIM: len_tile == TILE_DIM, dim divides evenly. */
    cplx_t *restrict data = state->data;

    /* Precompute per-index phases */
    double *restrict cos_d = malloc(dim * sizeof(double));
    double *restrict sin_d = malloc(dim * sizeof(double));
    if (!cos_d || !sin_d) {
        free(cos_d);
        free(sin_d);
        fprintf(stderr, "orkan: circ: exp_diag_tiled: allocation failed\n");
        exit(EXIT_FAILURE);
    }

    #pragma omp parallel for schedule(static) if(dim >= OMP_THRESHOLD_MIXED)
    for (idx_t x = 0; x < dim; ++x) {
        const double angle = diag[x] * t;
        orkan_sincos(angle, &sin_d[x], &cos_d[x]);
    }

    #pragma omp parallel for schedule(dynamic, 1) if(dim >= OMP_THRESHOLD_MIXED)
    for (idx_t tr = 0; tr < dim_tile; ++tr) {
        const idx_t r_base = tr * len_tile;

        for (idx_t tc = 0; tc <= tr; ++tc) {
            const idx_t c_base = tc * len_tile;
            cplx_t *restrict tile = data + tile_off(tr, tc);

            if (tr == tc) {
                /* Diagonal tile: strict lower triangle (skip diagonal).
                 * Mirror each write to the upper triangle so that
                 * within-tile gate paths (which read both triangles)
                 * see consistent post-evolution values. */
                for (idx_t lr = 1; lr < len_tile; ++lr) {
                    const idx_t r = r_base + lr;
                    for (idx_t lc = 0; lc < lr; ++lc) {
                        const idx_t c   = c_base + lc;
                        const idx_t off = elem_off(lr, lc);

                        const double f_re =  cos_d[r] * cos_d[c] + sin_d[r] * sin_d[c];
                        const double f_im = -sin_d[r] * cos_d[c] + cos_d[r] * sin_d[c];

                        const double re = creal(tile[off]);
                        const double im = cimag(tile[off]);
                        tile[off] = CMPLX(f_re * re - f_im * im,
                                          f_re * im + f_im * re);

                        /* Mirror to upper triangle: conj of new lower value */
                        tile[elem_off(lc, lr)] = conj(tile[off]);
                    }
                }
            } else {
                /* Off-diagonal tile: all elements */
                for (idx_t lr = 0; lr < len_tile; ++lr) {
                    const idx_t r = r_base + lr;
                    #pragma omp simd
                    for (idx_t lc = 0; lc < len_tile; ++lc) {
                        const idx_t c   = c_base + lc;
                        const idx_t off = elem_off(lr, lc);

                        const double f_re =  cos_d[r] * cos_d[c] + sin_d[r] * sin_d[c];
                        const double f_im = -sin_d[r] * cos_d[c] + cos_d[r] * sin_d[c];

                        const double re = creal(tile[off]);
                        const double im = cimag(tile[off]);
                        tile[off] = CMPLX(f_re * re - f_im * im,
                                          f_re * im + f_im * re);
                    }
                }
            }
        }
    }

    free(cos_d);
    free(sin_d);
}
