/**
 * @file meas_pure.c
 * @brief Measurement functions for pure states
 *
 * - mean_pure:        <H> = sum_x obs[x] |psi(x)|^2
 * - matel_diag_pure:  <bra|H|ket> = sum_x obs[x] conj(bra[x]) ket[x]
 *
 * Both are length-N reductions, parallelised with OpenMP gated by
 * OMP_THRESHOLD_PURE.
 */

#include "meas.h"

#include <complex.h>

#ifdef _OPENMP
#include <omp.h>
#endif

double mean_pure(const state_t *state, const double *obs) {
    const idx_t dim = POW2(state->qubits, idx_t);
    double sum = 0.0;

    #pragma omp parallel for reduction(+:sum) schedule(static) if(dim >= OMP_THRESHOLD_PURE)
    for (idx_t x = 0; x < dim; ++x) {
        const cplx_t a = state->data[x];
        sum += obs[x] * (creal(a) * creal(a) + cimag(a) * cimag(a));
    }
    return sum;
}

cplx_t matel_diag_pure(const state_t *bra, const state_t *ket, const double *obs) {
    const idx_t dim = POW2(bra->qubits, idx_t);
    double sum_re = 0.0;
    double sum_im = 0.0;

    #pragma omp parallel for reduction(+:sum_re, sum_im) schedule(static) \
            if(dim >= OMP_THRESHOLD_PURE)
    for (idx_t x = 0; x < dim; ++x) {
        const double ar = creal(bra->data[x]), ai = cimag(bra->data[x]);
        const double br = creal(ket->data[x]), bi = cimag(ket->data[x]);
        sum_re += obs[x] * (ar * br + ai * bi);
        sum_im += obs[x] * (ar * bi - ai * br);
    }
    return sum_re + sum_im * I;
}
