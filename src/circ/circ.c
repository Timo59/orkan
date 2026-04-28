/**
 * @file circ.c
 * @brief Dispatch layer for circuit operations
 *
 * Dispatches to type-specific implementations based on state->type.
 */

#include "circ.h"

#include <stdio.h>
#include <stdlib.h>

/*
 * =====================================================================================================================
 * Error handling
 * =====================================================================================================================
 *
 * CIRC_VALIDATE calls exit() on failure.  Passing a null pointer or invalid
 * state type indicates a programming error, not a recoverable runtime
 * condition.  This matches the gate, channel, and meas modules.
 */

#define CIRC_VALIDATE(cond, msg) do {                           \
    if (!(cond)) {                                              \
        fprintf(stderr, "orkan: circ: %s\n", (msg));           \
        exit(EXIT_FAILURE);                                     \
    }                                                           \
} while(0)

/*
 * =====================================================================================================================
 * Backend declarations
 * =====================================================================================================================
 */

extern void exp_diag_pure(state_t *state, const double *diag, double t);
extern void exp_diag_packed(state_t *state, const double *diag, double t);
extern void exp_diag_tiled(state_t *state, const double *diag, double t);

/*
 * =====================================================================================================================
 * Public API
 * =====================================================================================================================
 */

void exp_diag(state_t *state, const double *diag, double t) {
    CIRC_VALIDATE(state && state->data, "exp_diag: null state or data pointer");
    CIRC_VALIDATE(diag, "exp_diag: null diagonal pointer");

    if (t == 0.0) return;

    switch (state->type) {
        case PURE:
            exp_diag_pure(state, diag, t);
            break;
        case MIXED_PACKED:
            exp_diag_packed(state, diag, t);
            break;
        case MIXED_TILED:
            exp_diag_tiled(state, diag, t);
            break;
        default:
            CIRC_VALIDATE(0, "exp_diag: invalid state type");
    }
}

