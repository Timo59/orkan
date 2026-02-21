/**
 * @file gate_tiled.h
 * @brief Shared inline functions for tiled density matrix gate implementations
 *
 * Tiled storage partitions the density matrix into TILE_DIM x TILE_DIM tiles.
 * Only lower-triangular tiles are stored (at tile level). Within each tile,
 * elements are stored row-major.
 *
 * This header provides common inline helpers shared across gate_tiled*.c files:
 *   - tile_off(): tile base offset in the data array
 *   - elem_off(): element offset within a tile
 *   - dtile_read/dtile_write(): canonical access for diagonal tiles
 *
 * insertBit0() and insertBits2_0() are defined in gate.h.
 */

#ifndef GATE_TILED_H
#define GATE_TILED_H

#include "gate.h"
#include <complex.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/* Minimum dimension to enable OpenMP parallelization.
 * At dim < 512 (n < 9), OpenMP thread-launch overhead (~32µs) dominates
 * the actual computation time. See profile/REPORT_x_tiled.md §5.1. */
#ifndef OMP_THRESHOLD
#define OMP_THRESHOLD 512
#endif

/*
 * =====================================================================================================================
 * Tiled index helpers
 * =====================================================================================================================
 */

static inline gate_idx_t tile_off(gate_idx_t tr, gate_idx_t tc) {
    return (tr * (tr + 1) / 2 + tc) * TILE_SIZE;
}

static inline gate_idx_t elem_off(gate_idx_t lr, gate_idx_t lc) {
    return lr * TILE_DIM + lc;
}

/*
 * =====================================================================================================================
 * Diagonal tile helpers
 * =====================================================================================================================
 *
 * Helpers for reading/writing from diagonal tiles through the canonical
 * (lower-triangle) position. In diagonal tiles (tr == tc), elements at
 * (lr, lc) and (lc, lr) are conjugates. During swap operations, one
 * element of a conjugate pair may be modified before the other is read,
 * causing stale data. These helpers always access the lower-triangle
 * position, avoiding stale upper-triangle reads.
 */

static inline cplx_t dtile_read(cplx_t *tile, gate_idx_t lr, gate_idx_t lc) {
    return (lr >= lc) ? tile[elem_off(lr, lc)] : conj(tile[elem_off(lc, lr)]);
}

static inline void dtile_write(cplx_t *tile, gate_idx_t lr, gate_idx_t lc, cplx_t val) {
    if (lr >= lc) {
        tile[elem_off(lr, lc)] = val;
    } else {
        tile[elem_off(lc, lr)] = conj(val);
    }
}

/*
 * =====================================================================================================================
 * Tiled gate function declarations
 * =====================================================================================================================
 */

/* Arbitrary 2x2 unitary gate applied as ρ' = U ρ U† (col-major mat) */
void single_from_mat(state_t *state, qubit_t target, const cplx_t *mat);

/* Single-qubit gates (no parameter) */
void h_tiled(state_t *state, const qubit_t target);
void hy_tiled(state_t *state, const qubit_t target);
void s_tiled(state_t *state, const qubit_t target);
void sdg_tiled(state_t *state, const qubit_t target);
void t_tiled(state_t *state, const qubit_t target);
void tdg_tiled(state_t *state, const qubit_t target);
void x_tiled(state_t *state, const qubit_t target);
void y_tiled(state_t *state, const qubit_t target);
void z_tiled(state_t *state, const qubit_t target);

/* Single-qubit rotation gates (with angle parameter) */
void p_tiled(state_t *state, const qubit_t target, const double theta);
void rx_tiled(state_t *state, const qubit_t target, const double theta);
void ry_tiled(state_t *state, const qubit_t target, const double theta);
void rz_tiled(state_t *state, const qubit_t target, const double theta);

/* Two-qubit gates */
void cx_tiled(state_t *state, const qubit_t control, const qubit_t target);
void cy_tiled(state_t *state, const qubit_t control, const qubit_t target);
void cz_tiled(state_t *state, const qubit_t control, const qubit_t target);
void swap_tiled(state_t *state, const qubit_t q1, const qubit_t q2);

/* Three-qubit gates */
void ccx_tiled(state_t *state, const qubit_t control1, const qubit_t control2, const qubit_t target);

#endif /* GATE_TILED_H */
