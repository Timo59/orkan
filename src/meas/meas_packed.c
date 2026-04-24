/**
 * @file meas_packed.c
 * @brief Expectation value of a diagonal observable for packed mixed states
 *
 * Computes  <H> = sum_x obs[x] rho_{xx}
 *
 * Diagonal elements of packed lower-triangular column-major storage sit at
 * indices  0, dim, 2*dim-1, 3*dim-3, ...  with strides  dim, dim-1, dim-2, ...
 *
 * Each thread computes its starting packed index once and then uses the
 * stride-decrement pattern sequentially within its chunk.
 */

#include "meas.h"
#include "index.h"

#include <complex.h>

#ifdef _OPENMP
#include <omp.h>
#endif

double mean_packed(const state_t *state, const double *obs) {
    const idx_t dim = POW2(state->qubits, idx_t);
    double sum = 0.0;

    #pragma omp parallel reduction(+:sum) if(dim >= OMP_THRESHOLD_MIXED)
    {
#ifdef _OPENMP
        int tid = omp_get_thread_num();
        int nthreads = omp_get_num_threads();
#else
        int tid = 0;
        int nthreads = 1;
#endif
        idx_t chunk = dim / nthreads;
        idx_t x_start = tid * chunk;
        idx_t x_end = (tid == nthreads - 1) ? dim : x_start + chunk;

        /* Starting packed index for diagonal element (x_start, x_start) */
        idx_t idx = pack_idx(dim, x_start, x_start);
        idx_t stride = dim - x_start;

        for (idx_t x = x_start; x < x_end; ++x) {
            sum += obs[x] * creal(state->data[idx]);
            idx += stride;
            --stride;
        }
    }
    return sum;
}
