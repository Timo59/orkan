/**
 * @file gate_tiled_p.c
 * @brief Phase gate P(theta) for tiled density matrices
 *
 * P(theta) = diag(1, e^(i*theta)).  Conjugation: rho' = P rho P-dagger.
 *
 * The off-diagonal blocks transform as:
 *   rho'_10 = e^(+i*theta) * rho_10
 *   rho'_01 = e^(-i*theta) * rho_01
 * while diagonal blocks (rho_00, rho_11) are unchanged.
 *
 * This is identical to the Rz(theta) conjugation, which multiplies
 * the same e^(+/-i*theta) phases on the (1,0) and (0,1) sub-blocks.
 * The global phase difference between P and Rz (P has no phase on |0>,
 * while Rz has e^(-i*theta/2)) cancels in the conjugation P rho P-dagger
 * vs Rz rho Rz-dagger.  Therefore p_tiled simply delegates to rz_tiled.
 */

#include "gate_tiled.h"

void p_tiled(state_t *state, const qubit_t target, const double theta) {
    rz_tiled(state, target, theta);
}
