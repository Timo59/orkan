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

/*
 * @brief   Matrix element of a diagonal observable between two pure states
 *
 * Computes  <bra| H |ket> = sum_x obs[x] conj(bra[x]) ket[x]
 * for a diagonal observable H = diag(obs).
 *
 * When bra == ket this reduces to mean(), but returns a complex value
 * (imaginary part zero to machine precision).
 *
 * Only PURE states are supported.  Mixed states cause program termination
 * with a diagnostic message.
 *
 * @param[in]   bra     Pure quantum state (bra vector)
 * @param[in]   ket     Pure quantum state (ket vector)
 * @param[in]   obs     Diagonal coefficients, length 2^(bra->qubits), caller-owned
 * @return      Complex matrix element <bra| H |ket>
 */
cplx_t matel_diag(const state_t *bra, const state_t *ket, const double *obs);

#ifdef __cplusplus
}
#endif

#endif  /* MEAS_H */
