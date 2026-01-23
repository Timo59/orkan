#ifndef MHIPSTER_BLOCK_H
#define MHIPSTER_BLOCK_H

/**
 * @file mhipster_block.h
 * @brief Quantum gates for blocked density matrices
 *
 * Blocked gate implementations optimized for cache locality and SIMD.
 * Each gate transforms ρ → UρU† using tile-based processing.
 */

#include "state_blocked.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * =====================================================================================================================
 * Single-qubit gates (blocked format)
 * =====================================================================================================================
 */

/**
 * @brief Pauli-X gate: ρ → XρX
 */
void x_blocked(state_blocked_t *state, qubit_t target);

/**
 * @brief Pauli-Y gate: ρ → YρY†
 */
void y_blocked(state_blocked_t *state, qubit_t target);

/**
 * @brief Pauli-Z gate: ρ → ZρZ
 */
void z_blocked(state_blocked_t *state, qubit_t target);

/**
 * @brief Hadamard gate: ρ → HρH
 */
void h_blocked(state_blocked_t *state, qubit_t target);

/**
 * @brief S gate (Phase gate): ρ → SρS†
 */
void s_blocked(state_blocked_t *state, qubit_t target);

/**
 * @brief S-dagger gate: ρ → S†ρS
 */
void sdg_blocked(state_blocked_t *state, qubit_t target);

/**
 * @brief T gate (π/8 gate): ρ → TρT†
 */
void t_blocked(state_blocked_t *state, qubit_t target);

/**
 * @brief T-dagger gate: ρ → T†ρT
 */
void tdg_blocked(state_blocked_t *state, qubit_t target);

#ifdef __cplusplus
}
#endif

#endif /* MHIPSTER_BLOCK_H */
