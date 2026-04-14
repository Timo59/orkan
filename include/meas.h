#ifndef MEAS_H
#define MEAS_H

/*
 * =====================================================================================================================
 * includes
 * =====================================================================================================================
 */

#include "state.h"

/*
 * =====================================================================================================================
 * C++ check
 * =====================================================================================================================
 */

#ifdef __cplusplus
extern "C" {
#endif

/*
 * =====================================================================================================================
 * Public API
 * =====================================================================================================================
 */

/*
 * @brief   Expectation value of a diagonal observable
 *
 * Computes  <H> = sum_x obs[x] |psi(x)|^2          (PURE)
 *           <H> = sum_x obs[x] rho_{xx}             (MIXED_PACKED / MIXED_TILED)
 *
 * @param[in]   state   Quantum state (any representation)
 * @param[in]   obs     Diagonal coefficients, length 2^(state->qubits), caller-owned
 * @return      Real expectation value Tr(rho H)
 */
double mean(const state_t *state, const double *obs);

#ifdef __cplusplus
}
#endif

#endif  /* MEAS_H */
