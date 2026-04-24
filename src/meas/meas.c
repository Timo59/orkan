/**
 * @file meas.c
 * @brief Dispatch layer for measurement operations
 *
 * Dispatches to type-specific implementations based on state->type.
 */

#include "meas.h"

#include <stdio.h>
#include <stdlib.h>

/*
 * =====================================================================================================================
 * Error handling
 * =====================================================================================================================
 *
 * MEAS_VALIDATE calls exit() on failure.  Passing a null pointer or invalid
 * state type indicates a programming error, not a recoverable runtime
 * condition.  This matches the gate and channel modules.
 */

#define MEAS_VALIDATE(cond, msg) do {                           \
    if (!(cond)) {                                              \
        fprintf(stderr, "qlib: meas: %s\n", (msg));            \
        exit(EXIT_FAILURE);                                     \
    }                                                           \
} while(0)

/*
 * =====================================================================================================================
 * Backend declarations
 * =====================================================================================================================
 */

extern double mean_pure(const state_t *state, const double *obs);
extern double mean_packed(const state_t *state, const double *obs);
extern double mean_tiled(const state_t *state, const double *obs);

extern cplx_t matel_diag_pure(const state_t *bra, const state_t *ket, const double *obs);
extern void   sample_matrix_pure(const state_t *state, const double *obs,
                                 qop_t *const *ops, const idx_t *sizes, idx_t n_layers,
                                 cplx_t *out);

/*
 * =====================================================================================================================
 * Public API
 * =====================================================================================================================
 */

double mean(const state_t *state, const double *obs) {
    MEAS_VALIDATE(state && state->data, "mean: null state or data pointer");
    MEAS_VALIDATE(obs, "mean: null observable pointer");

    switch (state->type) {
        case PURE:
            return mean_pure(state, obs);
        case MIXED_PACKED:
            return mean_packed(state, obs);
        case MIXED_TILED:
            return mean_tiled(state, obs);
        default:
            MEAS_VALIDATE(0, "mean: invalid state type");
            return 0.0;
    }
}

cplx_t matel_diag(const state_t *bra, const state_t *ket, const double *obs) {
    MEAS_VALIDATE(bra && bra->data, "matel_diag: null bra state or data pointer");
    MEAS_VALIDATE(ket && ket->data, "matel_diag: null ket state or data pointer");
    MEAS_VALIDATE(obs, "matel_diag: null observable pointer");
    MEAS_VALIDATE(bra->qubits == ket->qubits, "matel_diag: qubit count mismatch");
    MEAS_VALIDATE(bra->type == PURE, "matel_diag: bra must be PURE (mixed not supported)");
    MEAS_VALIDATE(ket->type == PURE, "matel_diag: ket must be PURE (mixed not supported)");

    return matel_diag_pure(bra, ket, obs);
}

void sample_matrix(const state_t *state, const double *obs,
                   qop_t *const *ops, const idx_t *sizes, idx_t n_layers,
                   cplx_t *out) {
    MEAS_VALIDATE(state && state->data, "sample_matrix: null state or data pointer");
    MEAS_VALIDATE(obs, "sample_matrix: null observable pointer");
    MEAS_VALIDATE(ops, "sample_matrix: null ops pointer");
    MEAS_VALIDATE(sizes, "sample_matrix: null sizes pointer");
    MEAS_VALIDATE(n_layers > 0, "sample_matrix: n_layers must be positive");
    MEAS_VALIDATE(out, "sample_matrix: null output pointer");
    MEAS_VALIDATE(state->type == PURE,
                  "sample_matrix: requires PURE state (mixed not supported)");

    for (idx_t l = 0; l < n_layers; ++l) {
        MEAS_VALIDATE(sizes[l] > 0, "sample_matrix: sizes[l] must be positive");
        MEAS_VALIDATE(ops[l], "sample_matrix: null ops[l] pointer");
    }

    sample_matrix_pure(state, obs, ops, sizes, n_layers, out);
}
