/**
 * @file meas_pure.c
 * @brief Expectation value of a diagonal observable for pure states
 *
 * Computes  <H> = sum_x obs[x] |psi(x)|^2
 *
 * The loop body is a contiguous multiply-accumulate, auto-vectorised at -O2.
 * OpenMP parallelises over the state vector with a reduction.
 */

#include "meas.h"

#include <complex.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef OMP_THRESHOLD
#define OMP_THRESHOLD 4096
#endif

double mean_pure(const state_t *state, const double *obs) {
    const idx_t dim = POW2(state->qubits, idx_t);
    double sum = 0.0;

    #pragma omp parallel for reduction(+:sum) schedule(static) if(dim >= OMP_THRESHOLD)
    for (idx_t x = 0; x < dim; ++x) {
        const cplx_t a = state->data[x];
        sum += obs[x] * (creal(a) * creal(a) + cimag(a) * cimag(a));
    }
    return sum;
}
