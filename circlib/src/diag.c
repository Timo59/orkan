/*
 * =================================================================================================
 *                                              includes
 * =================================================================================================
 */

#include "diag.h"

/*
 * =================================================================================================
 *                                              apply diagonal operator
 * =================================================================================================
 */

void applyOperatorDiag(state_t* state, const cplx_t op[]) {
    for (dim_t i = 0; i < state->dimension; ++i) {
        state->vector[i] *= op[i];
    }
}

/*
 * =================================================================================================
 *                                              apply diagonal observable
 * =================================================================================================
 */

void applyObservableDiag(state_t* state, const double observable[]) {
    for (dim_t i = 0; i < state->dimension; ++i) {
        state->vector[i] *= observable[i];
    }
}

/*
 * =================================================================================================
 *                                              evolve with diagonal observable
 * =================================================================================================
 */

void evolveWithObservableDiag(state_t* state, const double observable[], double angle) {
    for (dim_t i = 0; i < state->dimension; ++i) {
        state->vector[i] *= cexp(-I * observable[i] * angle);
    }
}
