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
 * @brief   Function pointer type for quantum operations
 *
 * A quantum operation modifies a state in-place.  Operations with additional
 * parameters (target qubit, rotation angle, ...) must be wrapped in a
 * function matching this signature.
 */
typedef void (*qop_t)(state_t *state);

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

/*
 * @brief   Sample matrix of a diagonal observable under layered operations
 *
 * Computes the Hermitian sample matrix S with entries
 *
 *      S_{i,j} = <iota| V_i^+ H V_j |iota>
 *
 * where V_k = U_{k_L}^{(L)} ... U_{k_1}^{(1)} is the composite operation
 * for multi-index k, and H = diag(obs) is a diagonal observable.
 *
 * The flat index k maps to multi-index (k_1, ..., k_L) via mixed-radix
 * decomposition with k_1 least significant (fastest-varying):
 *
 *      k = k_L * (n_{L-1} * ... * n_1) + ... + k_2 * n_1 + k_1
 *
 * Output is the packed lower triangle of S in column-major order,
 * K(K+1)/2 elements where K = sizes[0] * ... * sizes[n_layers-1].
 *
 * Memory: holds at most two state-vector copies simultaneously.
 *
 * Only PURE states are supported.  Mixed states cause program termination
 * with a diagnostic message.
 *
 * @param[in]   state       Pure quantum state |iota>
 * @param[in]   obs         Diagonal coefficients, length 2^(state->qubits)
 * @param[in]   ops         Layered operations: ops[l][k] is the k-th operation
 *                          in layer l (l = 0..n_layers-1, k = 0..sizes[l]-1)
 * @param[in]   sizes       Number of operations per layer, length n_layers
 * @param[in]   n_layers    Number of operation layers (L)
 * @param[out]  out         Full matrix (column-major), caller-allocated,
 *                          length K*K
 */
void sample_matrix(const state_t *state, const double *obs,
                   qop_t *const *ops, const idx_t *sizes, idx_t n_layers,
                   cplx_t *out);

#ifdef __cplusplus
}
#endif

#endif  /* MEAS_H */
