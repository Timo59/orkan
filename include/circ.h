// circ.h - Circuit-level abstractions for quantum simulation
//
// Provides Pauli string types and operations for efficient multi-qubit Pauli
// application. The pauli_t structure uses bitmask representations (x_mask,
// z_mask, y_mask) to enable O(2^n) in-place application without per-qubit
// gate dispatch overhead.

#ifndef CIRC_H
#define CIRC_H

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
 * @brief Multi-qubit Pauli string P = P_{n-1} otimes ... otimes P_0
 *
 * Stores precomputed bitmasks for O(2^n) application:
 *   x_mask: bit j set iff qubit j in {X, Y}  (Pauli flips qubit j)
 *   z_mask: bit j set iff qubit j in {Y, Z}  (Pauli phases qubit j)
 *   y_mask: bit j set iff qubit j = Y
 * Invariants: y_mask ⊆ x_mask, y_mask ⊆ z_mask. n_qubits <= 64.
 */
typedef struct {
    qubit_t  n_qubits;
    uint64_t x_mask;      /* bit j set iff qubit j in {X, Y} */
    uint64_t z_mask;      /* bit j set iff qubit j in {Y, Z} */
    uint64_t y_mask;      /* bit j set iff qubit j = Y        */
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
 * @brief Parse a Pauli string from an "IXYZ..." character string
 *
 * The string is read left-to-right as qubit n-1 down to qubit 0 (big-endian
 * qubit ordering, matching standard physics notation).
 *
 * Aborts via GATE_VALIDATE on NULL s, length mismatch, or invalid character.
 *
 * @param[in] s         Null-terminated string of 'I','X','Y','Z' characters
 * @param[in] n_qubits  Expected number of qubits (must match strlen(s), <= 64)
 * @return Heap-allocated pauli_t
 */
pauli_t *pauli_from_str(const char *s, qubit_t n_qubits);

/**
 * @brief Free a Pauli string. Safe to call with NULL.
 *
 * @param[in,out] p  Pauli string to free
 */
void pauli_free(pauli_t *p);

/**
 * @brief Reconstruct a Pauli string as a human-readable character string
 *
 * Returns a heap-allocated null-terminated string of length n_qubits written
 * big-endian: s[0] = qubit n-1, s[n-1] = qubit 0. The caller must free() it.
 *
 * @param[in] P  Pauli string to serialise
 * @return Heap-allocated string of length n_qubits
 */
char *pauli_to_str(const pauli_t *P);

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
 * @brief Apply a Pauli string to a quantum state in-place: |psi> -> P|psi>
 *
 * Dispatches to the correct backend based on state->type.
 * Aborts via GATE_VALIDATE on null pointers or mask overflow.
 *
 * @param[in,out] state  Quantum state
 * @param[in]     P      Pauli string (masks must fit within state dimension)
 */
void pauli_apply(state_t *state, const pauli_t *P);

/**
 * @brief Apply the Pauli gadget exp(-i*theta/2 * P) to a quantum state in-place
 *
 * U(theta) = exp(-i*theta/2 * P) = cos(theta/2) * I  -  i * sin(theta/2) * P
 *
 * Dispatches to the correct backend based on state->type.
 * Aborts via GATE_VALIDATE on null pointers or mask overflow.
 *
 * @param[in,out] state  Quantum state
 * @param[in]     P      Pauli string (masks must fit within state dimension)
 * @param[in]     theta  Rotation angle in radians; U = exp(-i*theta/2 * P)
 */
void pauli_exp(state_t *state, const pauli_t *P, double theta);

#ifdef __cplusplus
}
#endif

#endif /* CIRC_H */
