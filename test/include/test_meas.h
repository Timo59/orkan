/**
 * @file test_meas.h
 * @brief Shared helpers for measurement module tests
 */

#ifndef TEST_MEAS_H
#define TEST_MEAS_H

#include "q_types.h"

static inline void fill_identity_obs(double *obs, idx_t dim) {
    for (idx_t x = 0; x < dim; ++x)
        obs[x] = (double)x;
}

static inline void fill_uniform_obs(double *obs, idx_t dim) {
    for (idx_t x = 0; x < dim; ++x)
        obs[x] = 1.0;
}

#endif /* TEST_MEAS_H */
