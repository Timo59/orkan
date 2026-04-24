/**
 * @file meas_pure.c
 * @brief Measurement functions for pure states
 *
 * - mean_pure:            <H> = sum_x obs[x] |psi(x)|^2
 * - matel_diag_pure:      <bra|H|ket> = sum_x obs[x] conj(bra[x]) ket[x]
 * - sample_matrix_pure:   S_{i,j} = <iota|V_i^+ H V_j|iota>  (full matrix, column-major)
 *
 * All reduction loops are OpenMP-parallelised, gated by OMP_THRESHOLD_PURE.
 * The sample_matrix outer loops are serial (operations are already parallel).
 */

#include "meas.h"

#include <assert.h>
#include <complex.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/* ====================== mean ====================== */

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

/* ====================== matel_diag ====================== */

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

/* ====================== apply_composite ====================== */

/*
 * @brief   Decode flat index and apply composite operation V_k in-place
 *
 * Decomposes k into multi-index (k_1, ..., k_L) with k_1 least significant
 * (fastest-varying) and applies U_{k_1}^{(1)}, ..., U_{k_L}^{(L)}
 * sequentially.
 */
static inline void apply_composite(state_t *state,
                                   qop_t *const *ops, const idx_t *sizes,
                                   idx_t n_layers, idx_t k)
{
    idx_t remainder = k;
    for (idx_t l = 0; l < n_layers; ++l) {
        const idx_t k_l = remainder % sizes[l];
        remainder /= sizes[l];
        ops[l][k_l](state);
    }
    assert(remainder == 0 && "apply_composite: k >= K (out of range)");
}

/* ====================== sample_matrix ====================== */

void sample_matrix_pure(const state_t *state, const double *obs,
                        qop_t *const *ops, const idx_t *sizes, idx_t n_layers,
                        cplx_t *out)
{
    idx_t K = 1;
    for (idx_t l = 0; l < n_layers; ++l)
        K *= sizes[l];

    for (idx_t j = 0; j < K; ++j) {
        /* evolve ket: phi_j = V_j |iota> */
        state_t phi_j = state_cp(state);
        if (!phi_j.data) goto oom;
        apply_composite(&phi_j, ops, sizes, n_layers, j);

        /* diagonal entry (real-valued): bra == ket */
        out[j * K + j] = matel_diag_pure(&phi_j, &phi_j, obs);

        /* off-diagonal: compute lower, mirror to upper */
        for (idx_t i = j + 1; i < K; ++i) {
            state_t phi_i = state_cp(state);
            if (!phi_i.data) { state_free(&phi_j); goto oom; }
            apply_composite(&phi_i, ops, sizes, n_layers, i);

            const cplx_t val = matel_diag_pure(&phi_i, &phi_j, obs);
            out[i + j * K] = val;        /* lower */
            out[j + i * K] = conj(val);  /* upper */
            state_free(&phi_i);
        }

        state_free(&phi_j);
    }
    return;

oom:
    fprintf(stderr, "qlib: meas: sample_matrix: allocation failed\n");
    exit(EXIT_FAILURE);
}
