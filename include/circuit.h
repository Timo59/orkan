// circuit.h - Circuit-level abstractions for quantum simulation
//
// Provides Pauli string types and operations for efficient multi-qubit Pauli
// application. The pauli_t structure uses bitmask representations (x_mask,
// z_mask, y_mask) to enable O(2^n) in-place application without per-qubit
// gate dispatch overhead.

#ifndef CIRCUIT_H
#define CIRCUIT_H

/*
 * =====================================================================================================================
 * Includes
 * =====================================================================================================================
 */
#include "state.h"

#include <stdint.h>

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
 * Type definitions
 * =====================================================================================================================
 */

/**
 * @brief Single-qubit Pauli label
 *
 * Encoding chosen so that bit 0 flags X-component and bit 1 flags Z-component:
 *   I = 0b00, X = 0b01, Y = 0b11 (both X and Z), Z = 0b10
 * This makes mask construction branchless.
 */
typedef enum {
    PAULI_I = 0,
    PAULI_X = 1,
    PAULI_Y = 2,
    PAULI_Z = 3
} pauli_label_t;

/**
 * @brief Multi-qubit Pauli string P = P_{n-1} otimes ... otimes P_0
 *
 * Stores both a dense label array (for introspection) and precomputed bitmasks
 * (for fast application). The bitmask convention is:
 *   x_mask: bit j set iff labels[j] in {X, Y}  (Pauli flips qubit j)
 *   z_mask: bit j set iff labels[j] in {Y, Z}  (Pauli phases qubit j)
 *   y_mask: bit j set iff labels[j] = Y
 */
typedef struct {
    pauli_label_t *labels;      /* heap-allocated, length n_qubits */
    qubit_t        n_qubits;
    uint64_t       x_mask;      /* bit j set iff labels[j] in {X, Y} */
    uint64_t       z_mask;      /* bit j set iff labels[j] in {Y, Z} */
    uint64_t       y_mask;      /* bit j set iff labels[j] = Y        */
} pauli_t;

/**
 * @brief Phase factor returned by pauli_phase
 *
 * Represents a phase from {+1, -1, +i, -i} as (real, imag) with integer
 * components in {-1, 0, +1}.
 */
typedef struct {
    int real;
    int imag;
} pauli_phase_t;

/*
 * =====================================================================================================================
 * Pauli string API
 * =====================================================================================================================
 */

/**
 * @brief Allocate a Pauli string from a dense label array
 *
 * @param[in] labels  Array of n_qubits Pauli labels (copied internally)
 * @param[in] n_qubits  Number of qubits (must be <= 64)
 * @return Heap-allocated pauli_t, or NULL on OOM or n_qubits > 64
 */
pauli_t *pauli_new(const pauli_label_t *labels, qubit_t n_qubits);

/**
 * @brief Parse a Pauli string from an "IXYZ..." character string
 *
 * The string is read left-to-right as qubit n-1 down to qubit 0 (big-endian
 * qubit ordering, matching standard physics notation).
 *
 * @param[in] s         Null-terminated string of 'I','X','Y','Z' characters
 * @param[in] n_qubits  Expected number of qubits (must match strlen(s))
 * @return Heap-allocated pauli_t, or NULL on invalid char, length mismatch, or OOM
 */
pauli_t *pauli_from_str(const char *s, qubit_t n_qubits);

/**
 * @brief Free a Pauli string and its internal label array
 *
 * Safe to call with NULL.
 *
 * @param[in,out] p  Pauli string to free
 */
void pauli_free(pauli_t *p);

/**
 * @brief Compute the phase factor P|x> = eps(x) * |x XOR x_mask>
 *
 * eps(x) = i^{popcount(y_mask)} * (-1)^{popcount(x & z_mask)}
 *
 * @param[in] P  Pauli string
 * @param[in] x  Computational basis index
 * @return Phase factor as (real, imag) with integer components
 */
pauli_phase_t pauli_phase(const pauli_t *P, uint64_t x);

/**
 * @brief Apply a Pauli string to a pure state in-place: |psi> -> P|psi>
 *
 * Uses bitmask representation for efficient traversal. For diagonal Paulis
 * (x_mask == 0), applies phase flips in O(2^n). For non-diagonal Paulis,
 * iterates over pairs (j, j XOR x_mask) to swap and phase amplitudes.
 *
 * @param[in,out] state  Pure quantum state (must have type == PURE)
 * @param[in]     P      Pauli string (masks must fit within state dimension)
 */
void pauli_apply_pure(state_t *state, const pauli_t *P);

/**
 * @brief Apply the Pauli gadget exp(-i*theta/2 * P) to a pure state in-place
 *
 * Implements the Euler decomposition:
 *   U(theta) = exp(-i*theta/2 * P) = cos(theta/2) * I  -  i * sin(theta/2) * P
 *
 * Three cases dispatched by mask structure:
 *
 * - Identity (x_mask == 0 && z_mask == 0): global phase e^{-i*theta/2} applied
 *   to every amplitude. Not a physical no-op; differs from the density-matrix case.
 *
 * - Diagonal (x_mask == 0): per-amplitude scalar multiply by e^{-i*theta/2}
 *   if epsilon(i) = +1, or e^{+i*theta/2} if epsilon(i) = -1.
 *
 * - Non-diagonal (x_mask != 0): paired update over (j, j XOR x_mask).
 *   Both amplitudes are read before either is written (no aliasing hazard).
 *
 * OpenMP-parallelised when dim >= OMP_THRESHOLD (default 4096).
 *
 * @param[in,out] state  Pure quantum state (must have type == PURE)
 * @param[in]     P      Pauli string (masks must fit within state dimension)
 * @param[in]     theta  Rotation angle in radians; U = exp(-i*theta/2 * P)
 */
void pauli_exp_pure(state_t *state, const pauli_t *P, double theta);

#ifdef __cplusplus
}
#endif

#endif /* CIRCUIT_H */
