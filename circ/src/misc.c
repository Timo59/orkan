/*
 * =====================================================================================================================
 *                                              includes
 * =====================================================================================================================
 */

#ifndef DIAG_H
#include "misc.h"
#endif

/*
 * =====================================================================================================================
 *                                              apply miscellaneous operator
 * =====================================================================================================================
 */
void applyOpMisc(state_t* state, unitary U) {
    return U(state);
}

/*
 * =====================================================================================================================
 *                                              apply miscellaneous observable
 * =====================================================================================================================
 */
void applyObsMisc(state_t* state, unitary U) {
    return U(state);
}

/*
 * =====================================================================================================================
 *                                              evolve with diagonal observable
 * =====================================================================================================================
 */
void evolveObsMisc(state_t* state, unitary U, double angle) {
    cplx_t* tmp = malloc(sizeof(cplx_t) * state->dim);              // Temporary vector holding state vector multiplied
                                                                    // with cos(angle)
    for (dim_t entry = 0; entry < state->dim; ++entry) {            // Iterate entries of the input's state vector
        tmp[entry] = cos(angle) * state->vec[entry];                // Store it times cos(angle) in tmp
        state->vec[entry] *= -I * sin(angle);                       // Multiply it with -i * sin(angle)
    }
    U(state);                                                       // Apply the unitary to the input state
    for (dim_t entry = 0; entry < state->dim; ++entry) {            // Iterate the entries of the input's state vector
        state->vec[entry] += tmp[entry];                            // Add the tmp's entry to it
    }
    free(tmp);
}
