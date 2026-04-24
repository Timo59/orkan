/**
 * @file gate_tiled.h
 * @brief Shared inline functions for tiled density matrix gate implementations
 *
 * Tiled storage partitions the density matrix into TILE_DIM x TILE_DIM tiles.
 * Only lower-triangular tiles are stored (at tile level). Within each tile,
 * elements are stored row-major.
 *
 * Shared declarations and OpenMP configuration for gate_tiled_*.c files.
 * Index helpers (insertBit0, tile_off, elem_off, etc.) are in index.h.
 */

#ifndef GATE_TILED_H
#define GATE_TILED_H

#include "gate.h"
#include "index.h"
#include <complex.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/* Minimum dimension to enable OpenMP parallelization.
 * At dim < 512 (n < 9), OpenMP thread-launch overhead (~32µs) dominates
 * the actual computation time. See profile/REPORT_x_tiled.md §5.1. */


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
