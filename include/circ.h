#ifndef CIRC_H
#define CIRC_H

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
 * @brief   Diagonal unitary evolution: U = exp(-i D t)
 *
 * Applies a state-dependent phase to each basis state:
 *   PURE:   psi[x]    *= exp(-i diag[x] t)
 *   MIXED:  rho[r][c] *= exp(-i (diag[r] - diag[c]) t)
 *
 * For mixed states, diagonal density-matrix elements are unchanged.
 *
 * @param[in,out]  state  Quantum state (any representation)
 * @param[in]      diag   Diagonal coefficients, length 2^(state->qubits), caller-owned
 * @param[in]      t      Evolution parameter (no implicit factor of 1/2)
 */
void exp_diag(state_t *state, const double *diag, double t);

#ifdef __cplusplus
}
#endif

#endif /* CIRC_H */
