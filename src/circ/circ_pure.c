/**
 * @file circ_pure.c
 * @brief Circuit operations for pure states
 *
 * exp_diag:  psi[x] *= exp(-i diag[x] t)
 *
 * Each iteration computes one cos + sin pair and one complex multiply.
 * The loop body is trig-bound; #pragma omp simd permits the compiler to
 * use vector math (libmvec / SVML) when available.
 */

#include "circ.h"

#include <complex.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef OMP_THRESHOLD
#define OMP_THRESHOLD 4096
#endif

void exp_diag_pure(state_t *state, const double *restrict diag, double t) {
    const idx_t dim = POW2(state->qubits, idx_t);
    cplx_t *restrict data = state->data;

    #pragma omp parallel for schedule(static) if(dim >= OMP_THRESHOLD)
    for (idx_t x = 0; x < dim; ++x) {
        const double angle = diag[x] * t;
        const double c = cos(angle);
        const double s = sin(angle);
        const double re = creal(data[x]);
        const double im = cimag(data[x]);
        data[x] = CMPLX(c * re + s * im, c * im - s * re);
    }
}
