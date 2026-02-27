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
    if (n_qubits > 64 || (n_qubits > 0 && labels == NULL)) {
        return NULL;
    }

    pauli_t *p = malloc(sizeof(pauli_t));
    if (!p) return NULL;

    p->n_qubits = n_qubits;
    p->x_mask = 0;
    p->z_mask = 0;
    p->y_mask = 0;

    if (n_qubits == 0) {
        p->labels = NULL;
        return p;
    }

    p->labels = malloc((size_t)n_qubits * sizeof(pauli_label_t));
    if (!p->labels) {
        free(p);
        return NULL;
    }

    memcpy(p->labels, labels, (size_t)n_qubits * sizeof(pauli_label_t));

    for (qubit_t j = 0; j < n_qubits; ++j) {
        const uint64_t bit = (uint64_t)1 << j;
        switch (labels[j]) {
            case PAULI_I:
                break;
            case PAULI_X:
                p->x_mask |= bit;
                break;
            case PAULI_Y:
                p->x_mask |= bit;
                p->z_mask |= bit;
                p->y_mask |= bit;
                break;
            case PAULI_Z:
                p->z_mask |= bit;
                break;
        }
    }

    return p;
}

pauli_t *pauli_from_str(const char *s, qubit_t n_qubits) {
    if (!s) return NULL;

    const size_t len = strlen(s);
    if (len != (size_t)n_qubits) return NULL;

    pauli_label_t *labels = malloc((size_t)n_qubits * sizeof(pauli_label_t));
    if (!labels) return NULL;

    /* String is big-endian: s[0] = qubit n-1, s[n-1] = qubit 0 */
    for (size_t i = 0; i < len; ++i) {
        const qubit_t qubit_idx = (qubit_t)(n_qubits - 1 - i);
        switch (s[i]) {
            case 'I': labels[qubit_idx] = PAULI_I; break;
            case 'X': labels[qubit_idx] = PAULI_X; break;
            case 'Y': labels[qubit_idx] = PAULI_Y; break;
            case 'Z': labels[qubit_idx] = PAULI_Z; break;
            default:
                free(labels);
                return NULL;
        }
    }

    pauli_t *p = pauli_new(labels, n_qubits);
    free(labels);
    return p;
}

void pauli_free(pauli_t *p) {
    if (!p) return;
    free(p->labels);
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
 * Three cases, dispatched by mask structure:
 *
 * 1. Identity (x_mask == 0 && z_mask == 0): no-op.
 *
 * 2. Diagonal (x_mask == 0): pure phase application.
 *    - Z-only (y_mask == 0): real sign flips (-1)^{popcount(i & z_mask)}.
 *    - With Y factors (unreachable: y_mask != 0 implies x_mask != 0, so
 *      this branch is dead, but kept for completeness).
 *
 * 3. Non-diagonal (x_mask != 0): swap + phase on pairs (j, j XOR x_mask).
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

    if (x_mask == 0 && z_mask == 0) {
        return;
    }

    const unsigned eta_exp = (unsigned)(__builtin_popcountll(y_mask) & 3u);

    if (x_mask == 0) {
        if (y_mask == 0) {
            /* Pure Z-string: real sign flips only */
            #pragma omp parallel for if(dim >= OMP_THRESHOLD)
            for (gate_idx_t i = 0; i < dim; ++i) {
                const unsigned parity = (unsigned)(__builtin_popcountll(i & z_mask) & 1u);
                const double sign = 1.0 - 2.0 * (double)parity;
                data[i] = data[i] * sign;
            }
        } else {
            /* Diagonal with Y factors: multiply by i^{eta_exp + 2*parity} */
            #pragma omp parallel for if(dim >= OMP_THRESHOLD)
            for (gate_idx_t i = 0; i < dim; ++i) {
                const unsigned parity = (unsigned)(__builtin_popcountll(i & z_mask) & 1u);
                const unsigned exp = (eta_exp + 2u * parity) & 3u;
                data[i] = apply_ipow(data[i], exp);
            }
        }
    } else {
        /* Non-diagonal: swap-and-phase on pairs (j, j XOR x_mask).
         *
         * Insert 0 at the lowest set bit of x_mask to enumerate exactly dim/2
         * indices j where that bit is 0; the partner jp = j XOR x_mask flips
         * all X-bits at once (including the split bit), ensuring jp != j and
         * each pair is visited exactly once.
         */
        const unsigned y_parity = (unsigned)(__builtin_popcountll(y_mask) & 1u);
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
 * Three cases, matching the structure of pauli_apply_pure:
 *
 * 1. Identity (x_mask == 0 && z_mask == 0): global phase e^{-i*theta/2}.
 *    All amplitudes scaled by CMPLX(cos(theta/2), -sin(theta/2)).
 *
 * 2. Diagonal (x_mask == 0): per-amplitude scalar multiply.
 *    data[i] = c*data[i] + s*apply_ipow(data[i], (exp_i+3)&3)
 *    where exp_i = (eta_exp + 2*parity_i) & 3 and -i corresponds to i^3.
 *
 *    For pure Z (y_mask == 0), only parity in {0,1} and eta_exp=0, so
 *    two precomputed scalars e^{-i*theta/2} and e^{+i*theta/2} suffice.
 *
 * 3. Non-diagonal (x_mask != 0): paired update, reading both old values first.
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

    if (x_mask == 0 && z_mask == 0) {
        /* Global phase: all amplitudes multiplied by e^{-i*theta/2} = c - i*s */
        const cplx_t phase = CMPLX(c, -s);
        #pragma omp parallel for if(dim >= OMP_THRESHOLD)
        for (gate_idx_t i = 0; i < dim; ++i) {
            data[i] = data[i] * phase;
        }
        return;
    }

    const unsigned eta_exp = (unsigned)(__builtin_popcountll(y_mask) & 3u);

    if (x_mask == 0) {
        if (y_mask == 0) {
            /* Pure Z-string: eta_exp = 0, parity_i in {0, 1} only.
             * parity=0 => exp_i=0 => -i * i^0 = i^3 => factor = c*data + s*i^3*data
             *          = (c - is)*data = e^{-i*theta/2}*data   (phase +1 eigenstate)
             * parity=1 => exp_i=2 => -i * i^2 = i^1 => factor = c*data + s*i^1*data
             *          = (c + is)*data = e^{+i*theta/2}*data   (phase -1 eigenstate)
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
            /* Diagonal with Y factors (structurally dead: y_mask != 0 implies x_mask != 0,
             * kept for completeness matching pauli_apply_pure).
             * data[i] = c*data[i] + s*apply_ipow(data[i], (exp_i+3)&3)
             */
            #pragma omp parallel for if(dim >= OMP_THRESHOLD)
            for (gate_idx_t i = 0; i < dim; ++i) {
                const unsigned parity = (unsigned)(__builtin_popcountll(i & z_mask) & 1u);
                const unsigned exp_i  = (eta_exp + 2u * parity) & 3u;
                data[i] = c * data[i] + s * apply_ipow(data[i], (exp_i + 3u) & 3u);
            }
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
