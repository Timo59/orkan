/**
 * @file meas.c
 * @brief Dispatch layer for measurement operations
 *
 * Dispatches to type-specific implementations based on state->type.
 */

#include "meas.h"

#include <assert.h>

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
    assert(state && state->data && obs);

    switch (state->type) {
        case PURE:
            return mean_pure(state, obs);
        case MIXED_PACKED:
            return mean_packed(state, obs);
        case MIXED_TILED:
            return mean_tiled(state, obs);
        default:
            assert(0 && "mean(): invalid state type");
            return 0.0;
    }
}
