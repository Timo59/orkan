/**
 * @file state_internal.h
 * @brief Type-specific state function declarations (internal to libq)
 *
 * These functions are the per-backend implementations that state.c dispatches
 * to based on state->type. They are not part of the public API — users should
 * call the generic state_get/state_set/state_len/etc. from state.h instead.
 */

#ifndef STATE_INTERNAL_H
#define STATE_INTERNAL_H

#include "state.h"

/*
 * =====================================================================================================================
 * PURE
 * =====================================================================================================================
 */

idx_t  state_pure_len(qubit_t qubits);
void   state_pure_plus(state_t *state, qubit_t qubits);
cplx_t state_pure_get(const state_t *state, idx_t row);
void   state_pure_set(state_t *state, idx_t row, cplx_t val);
void   state_pure_print(const state_t *state);

/*
 * =====================================================================================================================
 * MIXED_PACKED
 * =====================================================================================================================
 */

idx_t  state_packed_len(qubit_t qubits);
void   state_packed_plus(state_t *state, qubit_t qubits);
cplx_t state_packed_get(const state_t *state, idx_t row, idx_t col);
void   state_packed_set(state_t *state, idx_t row, idx_t col, cplx_t val);
void   state_packed_print(const state_t *state);

/*
 * =====================================================================================================================
 * MIXED_TILED
 *
 * WARNING: Intra-tile layout is row-major, NOT column-major like LAPACK packed.
 * Do NOT pass tile pointers to LAPACK routines expecting column-major storage.
 * =====================================================================================================================
 */

idx_t  state_tiled_len(qubit_t qubits);
void   state_tiled_plus(state_t *state, qubit_t qubits);
cplx_t state_tiled_get(const state_t *state, idx_t row, idx_t col);
void   state_tiled_set(state_t *state, idx_t row, idx_t col, cplx_t val);
void   state_tiled_print(const state_t *state);

#endif /* STATE_INTERNAL_H */
