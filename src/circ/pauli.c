// pauli.c - Pauli string construction and state-type dispatchers
//
// A Pauli string P = P_{n-1} (x) ... (x) P_0 is stored as bitmasks for
// efficient O(2^n) in-place application. The key identity is:
//
//   P|x> = i^{popcount(y_mask)} * (-1)^{popcount(x & z_mask)} * |x XOR x_mask>
//
// This avoids per-qubit gate dispatch and enables vectorized traversal.

#include "circ.h"
#include "gate.h"       /* GATE_VALIDATE */

#include <string.h>     /* strlen */

/*
 * =====================================================================================================================
 * Pure state implementations (pauli_pure.c)
 * =====================================================================================================================
 */

extern void pauli_apply_pure(state_t *state, const pauli_t *P);
extern void pauli_exp_pure(state_t *state, const pauli_t *P, double theta);

/*
 * =====================================================================================================================
 * Construction and destruction
 * =====================================================================================================================
 */


pauli_t *pauli_from_str(const char *s, qubit_t n_qubits) {
    GATE_VALIDATE(s, "pauli_from_str: null string");

    const size_t len = strlen(s);
    GATE_VALIDATE(len == (size_t)n_qubits,
                  "pauli_from_str: string length does not match n_qubits");
    GATE_VALIDATE(n_qubits <= 64,
                  "pauli_from_str: n_qubits exceeds maximum of 64");

    pauli_t *p = malloc(sizeof(pauli_t));
    GATE_VALIDATE(p, "pauli_from_str: out of memory");

    p->n_qubits = n_qubits;
    p->x_mask = 0;
    p->z_mask = 0;
    p->y_mask = 0;

    /* String is big-endian: s[0] = qubit n-1, s[n-1] = qubit 0 */
    for (size_t i = 0; i < len; ++i) {
        const qubit_t j = (qubit_t)(n_qubits - 1 - i);
        const uint64_t bit = (uint64_t)1 << j;
        switch (s[i]) {
            case 'I': break;
            case 'X': p->x_mask |= bit; break;
            case 'Y': p->x_mask |= bit; p->z_mask |= bit; p->y_mask |= bit; break;
            case 'Z': p->z_mask |= bit; break;
            default:
                GATE_VALIDATE(0, "pauli_from_str: invalid character (expected 'I','X','Y','Z')");
        }
    }

    return p;
}

char *pauli_to_str(const pauli_t *P) {
    GATE_VALIDATE(P, "pauli_to_str: null Pauli pointer");

    char *s = malloc((size_t)P->n_qubits + 1);
    GATE_VALIDATE(s, "pauli_to_str: out of memory");

    /* Big-endian: s[0] = qubit n-1, s[n-1] = qubit 0 */
    for (qubit_t j = 0; j < P->n_qubits; ++j) {
        const int has_x = (int)((P->x_mask >> j) & 1u);
        const int has_z = (int)((P->z_mask >> j) & 1u);
        char c;
        if      (has_x && has_z) c = 'Y';
        else if (has_x)          c = 'X';
        else if (has_z)          c = 'Z';
        else                     c = 'I';
        s[P->n_qubits - 1 - j] = c;
    }
    s[P->n_qubits] = '\0';
    return s;
}

void pauli_free(pauli_t *p) {
    if (!p) return;
    free(p);
}

/*
 * =====================================================================================================================
 * Phase computation
 * =====================================================================================================================
 */

pauli_phase_t pauli_phase(const pauli_t *P, uint64_t x) {
    /* i^{popcount(y_mask)} * (-1)^{popcount(x & z_mask)} */
    const unsigned y_pop = (unsigned)(__builtin_popcountll(P->y_mask) & 3u);
    const unsigned z_parity = (unsigned)(__builtin_popcountll(x & P->z_mask) & 1u);

    /* i^y_pop gives: 0->+1, 1->+i, 2->-1, 3->-i */
    static const int i_pow_real[] = { 1, 0, -1,  0};
    static const int i_pow_imag[] = { 0, 1,  0, -1};

    const int sign = 1 - 2 * (int)z_parity;   /* (-1)^z_parity */

    pauli_phase_t result;
    result.real = sign * i_pow_real[y_pop];
    result.imag = sign * i_pow_imag[y_pop];
    return result;
}

/*
 * =====================================================================================================================
 * Dispatchers
 * =====================================================================================================================
 */

void pauli_apply(state_t *state, const pauli_t *P) {
    GATE_VALIDATE(state && state->data,
                  "pauli_apply: null state or data pointer");
    GATE_VALIDATE(P,
                  "pauli_apply: null Pauli pointer");
    GATE_VALIDATE(state->qubits >= 64 ||
                  (P->x_mask | P->z_mask) >> state->qubits == 0,
                  "pauli_apply: Pauli mask exceeds state dimension");
    switch (state->type) {
        case PURE:         pauli_apply_pure(state, P);   break;
        case MIXED_PACKED: GATE_VALIDATE(0, "pauli_apply: not implemented for MIXED_PACKED"); break;
        case MIXED_TILED:  GATE_VALIDATE(0, "pauli_apply: not implemented for MIXED_TILED");  break;
        default:           GATE_VALIDATE(0, "pauli_apply: unknown state type");
    }
}

void pauli_exp(state_t *state, const pauli_t *P, double theta) {
    GATE_VALIDATE(state && state->data,
                  "pauli_exp: null state or data pointer");
    GATE_VALIDATE(P,
                  "pauli_exp: null Pauli pointer");
    GATE_VALIDATE(state->qubits >= 64 ||
                  (P->x_mask | P->z_mask) >> state->qubits == 0,
                  "pauli_exp: Pauli mask exceeds state dimension");
    switch (state->type) {
        case PURE:         pauli_exp_pure(state, P, theta);   break;
        case MIXED_PACKED: GATE_VALIDATE(0, "pauli_exp: not implemented for MIXED_PACKED"); break;
        case MIXED_TILED:  GATE_VALIDATE(0, "pauli_exp: not implemented for MIXED_TILED");  break;
        default:           GATE_VALIDATE(0, "pauli_exp: unknown state type");
    }
}
