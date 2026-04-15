/**
 * @file test_circ.h
 * @brief Shared helpers for circuit module tests
 */

#ifndef TEST_CIRC_H
#define TEST_CIRC_H

#include "q_types.h"

static inline void fill_ramp_diag(double *diag, idx_t dim) {
    for (idx_t x = 0; x < dim; ++x)
        diag[x] = (double)x;
}

static inline void fill_zero_diag(double *diag, idx_t dim) {
    for (idx_t x = 0; x < dim; ++x)
        diag[x] = 0.0;
}

static inline void fill_uniform_diag(double *diag, idx_t dim) {
    for (idx_t x = 0; x < dim; ++x)
        diag[x] = 1.0;
}

#endif /* TEST_CIRC_H */
