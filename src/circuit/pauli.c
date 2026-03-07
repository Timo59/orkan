// pauli.c - Pauli string construction and application to pure states
//
// A Pauli string P = P_{n-1} (x) ... (x) P_0 is stored as bitmasks for
// efficient O(2^n) in-place application. The key identity is:
//
//   P|x> = i^{popcount(y_mask)} * (-1)^{popcount(x & z_mask)} * |x XOR x_mask>
//
// This avoids per-qubit gate dispatch and enables vectorized traversal.

#include "circuit.h"
#include "gate.h"       /* insertBit0, gate_idx_t, GATE_VALIDATE */

#include <complex.h>    /* creal, cimag, CMPLX */
#include <math.h>       /* cos, sin */
#include <string.h>     /* strlen, memcpy */

#ifdef _OPENMP
#include <omp.h>
#endif

/* Minimum dimension to enable OpenMP parallelization.
 * Pauli application is O(dim) with low per-element cost, so we need a
 * reasonable threshold to offset thread-spawn overhead.
 */
#ifndef OMP_THRESHOLD
#define OMP_THRESHOLD 4096
#endif

/*
 * =====================================================================================================================
 * Construction and destruction
 * =====================================================================================================================
 */

pauli_t *pauli_new(const pauli_label_t *labels, qubit_t n_qubits) {
    GATE_VALIDATE(n_qubits <= 64,
                  "pauli_new: n_qubits exceeds maximum of 64");
    GATE_VALIDATE(n_qubits == 0 || labels != NULL,
                  "pauli_new: null labels with n_qubits > 0");

    pauli_t *p = malloc(sizeof(pauli_t));
    GATE_VALIDATE(p, "pauli_new: out of memory");

    p->n_qubits = n_qubits;
    p->x_mask = 0;
    p->z_mask = 0;
    p->y_mask = 0;

    if (n_qubits == 0) return p;

    for (qubit_t j = 0; j < n_qubits; ++j) {
        const uint64_t bit = (uint64_t)1 << j;
        switch (labels[j]) {
            case PAULI_I: break;
            case PAULI_X: p->x_mask |= bit; break;
            case PAULI_Y: p->x_mask |= bit; p->z_mask |= bit; p->y_mask |= bit; break;
            case PAULI_Z: p->z_mask |= bit; break;
        }
    }

    return p;
}

pauli_t *pauli_from_str(const char *s, qubit_t n_qubits) {
    GATE_VALIDATE(s, "pauli_from_str: null string");

    const size_t len = strlen(s);
    GATE_VALIDATE(len == (size_t)n_qubits,
                  "pauli_from_str: string length does not match n_qubits");

    if (n_qubits == 0) return pauli_new(NULL, 0);

    pauli_label_t *labels = malloc((size_t)n_qubits * sizeof(pauli_label_t));
    GATE_VALIDATE(labels, "pauli_from_str: out of memory");

    /* String is big-endian: s[0] = qubit n-1, s[n-1] = qubit 0 */
    for (size_t i = 0; i < len; ++i) {
        const qubit_t qubit_idx = (qubit_t)(n_qubits - 1 - i);
        switch (s[i]) {
            case 'I': labels[qubit_idx] = PAULI_I; break;
            case 'X': labels[qubit_idx] = PAULI_X; break;
            case 'Y': labels[qubit_idx] = PAULI_Y; break;
            case 'Z': labels[qubit_idx] = PAULI_Z; break;
            default:
                GATE_VALIDATE(0, "pauli_from_str: invalid character (expected 'I','X','Y','Z')");
        }
    }

    pauli_t *p = pauli_new(labels, n_qubits);
    free(labels);
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
 * Pure state application
 * =====================================================================================================================
 */

/**
 * @brief Multiply a complex number by i^exp (exp mod 4)
 *
 * i^0 = 1, i^1 = i, i^2 = -1, i^3 = -i
 * Each case is a simple coordinate rotation -- no multiplications needed.
 */
static inline cplx_t apply_ipow(const cplx_t z, const unsigned exp) {
    const double re = creal(z);
    const double im = cimag(z);
    switch (exp & 3u) {
        case 0: return z;
        case 1: return CMPLX(-im,  re);
        case 2: return CMPLX(-re, -im);
        case 3: return CMPLX( im, -re);
    }
    __builtin_unreachable();
}

/**
 * @brief Apply Pauli string P to pure state in-place: |psi> -> P|psi>
 *
 * Two cases dispatched by mask structure:
 *
 * 1. Diagonal (x_mask == 0): identity returns immediately; pure Z-string
 *    applies real sign flips (-1)^{popcount(i & z_mask)}.
 *    (y_mask != 0 with x_mask == 0 is structurally impossible: PAULI_Y always
 *    sets both x_mask and y_mask, so y_mask ⊆ x_mask by construction.)
 *
 * 2. Non-diagonal (x_mask != 0): swap + phase on pairs (j, j XOR x_mask).
 *    Enumerate j by inserting 0 at the LOWEST set bit of x_mask.
 *    This gives exactly dim/2 unique pairs, each (j, jp) differing at all
 *    x_mask bit positions. XOR with the full x_mask then applies all X-flips.
 *
 *    Key insight: n_pairs = dim/2, not dim/2^{popcount(x_mask)}.
 *    Only ONE insertBit0 call is needed — at the lowest x_mask bit.
 */
void pauli_apply_pure(state_t *state, const pauli_t *P) {
    GATE_VALIDATE(state && state->data,
                  "pauli_apply_pure: null state or data pointer");
    GATE_VALIDATE(P,
                  "pauli_apply_pure: null Pauli pointer");
    GATE_VALIDATE(state->type == PURE,
                  "pauli_apply_pure: expected PURE state");
    GATE_VALIDATE(state->qubits >= 64 ||
                  (P->x_mask | P->z_mask) >> state->qubits == 0,
                  "pauli_apply_pure: Pauli mask exceeds state dimension");

    const gate_idx_t dim = (gate_idx_t)1 << state->qubits;
    cplx_t *restrict data = state->data;

    const uint64_t x_mask = P->x_mask;
    const uint64_t z_mask = P->z_mask;
    const uint64_t y_mask = P->y_mask;

    if (x_mask == 0) {
        if (z_mask == 0) return;  /* identity: no-op */
        /* Pure Z-string (y_mask == 0 guaranteed): real sign flips only */
        #pragma omp parallel for if(dim >= OMP_THRESHOLD)
        for (gate_idx_t i = 0; i < dim; ++i) {
            const unsigned parity = (unsigned)(__builtin_popcountll(i & z_mask) & 1u);
            const double sign = 1.0 - 2.0 * (double)parity;
            data[i] = data[i] * sign;
        }
    } else {
        /* Non-diagonal: swap-and-phase on pairs (j, j XOR x_mask).
         *
         * Insert 0 at the lowest set bit of x_mask to enumerate exactly dim/2
         * indices j where that bit is 0; the partner jp = j XOR x_mask flips
         * all X-bits at once (including the split bit), ensuring jp != j and
         * each pair is visited exactly once.
         */
        const unsigned eta_exp   = (unsigned)(__builtin_popcountll(y_mask) & 3u);
        const unsigned y_parity  = (unsigned)(__builtin_popcountll(y_mask) & 1u);
        const qubit_t  split_bit = (qubit_t)__builtin_ctzll(x_mask);
        const gate_idx_t n_pairs = dim >> 1;

        #pragma omp parallel for if(dim >= OMP_THRESHOLD)
        for (gate_idx_t k = 0; k < n_pairs; ++k) {
            const gate_idx_t j  = insertBit0(k, split_bit);
            const gate_idx_t jp = j ^ x_mask;

            const unsigned parity_j = (unsigned)(__builtin_popcountll(j & z_mask) & 1u);
            const unsigned exp_j  = (eta_exp + 2u * parity_j) & 3u;
            const unsigned exp_jp = (exp_j + 2u * y_parity) & 3u;

            const cplx_t tmp = data[j];
            data[j]  = apply_ipow(data[jp], exp_jp);
            data[jp] = apply_ipow(tmp,      exp_j);
        }
    }
}

/**
 * @brief Apply Pauli gadget U(theta) = exp(-i*theta/2 * P) to a pure state in-place
 *
 * Uses the Euler decomposition:
 *   U(theta) = cos(theta/2) * I  -  i * sin(theta/2) * P
 *
 * Two cases dispatched by mask structure:
 *
 * 1. Diagonal (x_mask == 0):
 *    - Identity (z_mask == 0): global phase e^{-i*theta/2}. All amplitudes
 *      scaled by CMPLX(cos(theta/2), -sin(theta/2)). Not a no-op.
 *    - Pure Z-string (z_mask != 0): two precomputed scalars e^{±i*theta/2}
 *      selected per amplitude by (-1)^{popcount(i & z_mask)}.
 *      (y_mask != 0 with x_mask == 0 is structurally impossible.)
 *
 * 2. Non-diagonal (x_mask != 0): paired update, reading both old values first.
 *    new_j  = c*old_j  + s*apply_ipow(old_jp, (exp_jp+3)&3)
 *    new_jp = c*old_jp + s*apply_ipow(old_j,  (exp_j +3)&3)
 *    Derivation: -i*i^k = i^3 * i^k = i^{(k+3)&3}. This avoids all complex
 *    multiplications — only real scaling and coordinate rotations are used.
 *
 *    Pair enumeration identical to pauli_apply_pure (insertBit0 at split_bit).
 */
void pauli_exp_pure(state_t *state, const pauli_t *P, double theta) {
    GATE_VALIDATE(state && state->data,
                  "pauli_exp_pure: null state or data pointer");
    GATE_VALIDATE(P,
                  "pauli_exp_pure: null Pauli pointer");
    GATE_VALIDATE(state->type == PURE,
                  "pauli_exp_pure: expected PURE state");
    GATE_VALIDATE(state->qubits >= 64 ||
                  (P->x_mask | P->z_mask) >> state->qubits == 0,
                  "pauli_exp_pure: Pauli mask exceeds state dimension");

    const gate_idx_t dim = (gate_idx_t)1 << state->qubits;
    cplx_t *restrict data = state->data;

    const uint64_t x_mask = P->x_mask;
    const uint64_t z_mask = P->z_mask;
    const uint64_t y_mask = P->y_mask;

    const double half = theta * 0.5;
    const double c    = cos(half);
    const double s    = sin(half);

    if (x_mask == 0) {
        if (z_mask == 0) {
            /* Global phase: all amplitudes multiplied by e^{-i*theta/2} = c - i*s */
            const cplx_t phase = CMPLX(c, -s);
            #pragma omp parallel for if(dim >= OMP_THRESHOLD)
            for (gate_idx_t i = 0; i < dim; ++i) {
                data[i] = data[i] * phase;
            }
            return;
        }
        /* Pure Z-string (y_mask == 0 guaranteed): eta_exp = 0, parity_i in {0, 1}.
         * parity=0: e^{-i*theta/2}*data  (phase +1 eigenstate)
         * parity=1: e^{+i*theta/2}*data  (phase -1 eigenstate)
         * Precompute both scalars to avoid per-iteration apply_ipow dispatch.
         */
        const cplx_t fac_pos = CMPLX(c, -s);   /* epsilon = +1: e^{-i*theta/2} */
        const cplx_t fac_neg = CMPLX(c,  s);   /* epsilon = -1: e^{+i*theta/2} */
        #pragma omp parallel for if(dim >= OMP_THRESHOLD)
        for (gate_idx_t i = 0; i < dim; ++i) {
            const unsigned parity = (unsigned)(__builtin_popcountll(i & z_mask) & 1u);
            data[i] = data[i] * (parity ? fac_neg : fac_pos);
        }
    } else {
        /* Non-diagonal: paired update on (j, jp = j XOR x_mask).
         *
         * (P|psi>)_j  = eps(jp) * psi[jp]   where eps(x) = i^{exp_x}
         * new_j  = c*old_j  - i*s*eps(jp)*old_jp = c*old_j  + s*apply_ipow(old_jp, (exp_jp+3)&3)
         * new_jp = c*old_jp - i*s*eps(j) *old_j  = c*old_jp + s*apply_ipow(old_j,  (exp_j +3)&3)
         *
         * Both old values must be read before either write.
         * Phase exponents computed identically to pauli_apply_pure.
         */
        const unsigned eta_exp   = (unsigned)(__builtin_popcountll(y_mask) & 3u);
        const unsigned y_parity  = (unsigned)(__builtin_popcountll(y_mask) & 1u);
        const qubit_t  split_bit = (qubit_t)__builtin_ctzll(x_mask);
        const gate_idx_t n_pairs = dim >> 1;

        #pragma omp parallel for if(dim >= OMP_THRESHOLD)
        for (gate_idx_t k = 0; k < n_pairs; ++k) {
            const gate_idx_t j  = insertBit0(k, split_bit);
            const gate_idx_t jp = j ^ x_mask;

            const unsigned parity_j = (unsigned)(__builtin_popcountll(j & z_mask) & 1u);
            const unsigned exp_j    = (eta_exp + 2u * parity_j) & 3u;
            const unsigned exp_jp   = (exp_j + 2u * y_parity) & 3u;

            const cplx_t old_j  = data[j];
            const cplx_t old_jp = data[jp];

            data[j]  = c * old_j  + s * apply_ipow(old_jp, (exp_jp + 3u) & 3u);
            data[jp] = c * old_jp + s * apply_ipow(old_j,  (exp_j  + 3u) & 3u);
        }
    }
}
