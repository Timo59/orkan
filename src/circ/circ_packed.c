/**
 * @file circ_packed.c
 * @brief Circuit operations for packed mixed states
 *
 * exp_diag:  rho[r][c] *= exp(-i (diag[r] - diag[c]) t)
 *
 * Diagonal elements are unchanged.  Off-diagonal elements within each column
 * are contiguous in packed lower-triangular column-major storage.
 *
 * Trig calls are reduced from O(dim^2) to O(dim) by precomputing per-index
 * cos/sin values and reconstructing the pairwise phase via the product identity
 *   exp(-i(d[r]-d[c])t) = (cos_d[r] - i sin_d[r])(cos_d[c] + i sin_d[c]).
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

void exp_diag_packed(state_t *state, const double *restrict diag, double t) {
    const idx_t dim = POW2(state->qubits, idx_t);
    cplx_t *restrict data = state->data;

    /* Precompute per-index phases: phase[x] = e^{-i diag[x] t} */
    double *restrict cos_d = malloc(dim * sizeof(double));
    double *restrict sin_d = malloc(dim * sizeof(double));
    if (!cos_d || !sin_d) {
        free(cos_d);
        free(sin_d);
        fprintf(stderr, "orkan: circ: exp_diag_packed: allocation failed\n");
        exit(EXIT_FAILURE);
    }

    #pragma omp parallel for schedule(static) if(dim >= OMP_THRESHOLD_MIXED)
    for (idx_t x = 0; x < dim; ++x) {
        const double angle = diag[x] * t;
        orkan_sincos(angle, &sin_d[x], &cos_d[x]);
    }

    /* Iterate columns; elements within a column are contiguous */
    #pragma omp parallel for schedule(guided) if(dim >= OMP_THRESHOLD_MIXED)
    for (idx_t c = 0; c < dim; ++c) {
        const double cc = cos_d[c];
        const double sc = sin_d[c];
        idx_t idx = col_off(dim, c) + 1;        /* first off-diagonal in column */

        #pragma omp simd
        for (idx_t r = c + 1; r < dim; ++r) {
            /* e^{-i(d[r]-d[c])t} = (cr - i sr)(cc + i sc) */
            const double f_re =  cos_d[r] * cc + sin_d[r] * sc;
            const double f_im = -sin_d[r] * cc + cos_d[r] * sc;

            const double re = creal(data[idx]);
            const double im = cimag(data[idx]);
            data[idx] = CMPLX(f_re * re - f_im * im,
                              f_re * im + f_im * re);
            ++idx;
        }
    }

    free(cos_d);
    free(sin_d);
}
