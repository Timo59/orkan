/*
 * =====================================================================================================================
 *                                              includes
 * =====================================================================================================================
 */

#include "diag.h"

/*
 * =====================================================================================================================
 *                                              apply diagonal operator
 * =====================================================================================================================
 */

void applyOpDiag(state_t* state, const cplx_t op[]) {
    for (dim_t i = 0; i < state->dim; ++i) {
        state->vec[i] *= op[i];
    }
}

/*
 * =====================================================================================================================
 *                                              apply diagonal observable
 * =====================================================================================================================
 */

void applyObsDiag(state_t* state, const double obs[]) {
    for (dim_t i = 0; i < state->dim; ++i) {
        state->vec[i] *= obs[i];
    }
}

void applyObsDiagOmp(state_t* state, const double obs[]) {
#pragma omp parallel for default(none) shared(state, obs)
    for (dim_t i = 0; i < state->dim; ++i) {
        state->vec[i] *= obs[i];
    }
}

/*
 * =====================================================================================================================
 *                                              evolve with diagonal observable
 * =====================================================================================================================
 */

void evolveObsDiag(state_t* state, const double obs[], double angle) {
    for (dim_t i = 0; i < state->dim; ++i) {
        state->vec[i] *= cexp(-I * obs[i] * angle);
    }
}

void evolveObsDiagOmp(state_t* state, const double obs[], double angle) {
#pragma omp parallel for default(none) shared(angle, state, obs)
    for (dim_t i = 0; i < state->dim; ++i) {
        state->vec[i] *= cexp(-I * obs[i] * angle);
    }
}
