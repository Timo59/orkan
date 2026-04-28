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
        fprintf(stderr, "orkan: meas: %s\n", (msg));           \
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

