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
