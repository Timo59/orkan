// pauli_pure.c - Pauli string operations on PURE states
//
// Implementations for the PURE backend. Input validation is performed by the
// dispatchers in pauli.c before these functions are called.

#include "circ.h"
#include "gate.h"       /* insertBit0, idx_t */
#include "index.h"

#include <complex.h>    /* creal, cimag, CMPLX */
#include <math.h>       /* cos, sin */

#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef OMP_THRESHOLD
#define OMP_THRESHOLD 4096
#endif

/*
 * =====================================================================================================================
 * Internal helpers
 * =====================================================================================================================
 */

/**
 * @brief Multiply a complex number by i^exp (exp mod 4)
 *
 * i^0 = 1, i^1 = i, i^2 = -1, i^3 = -i
 * Each case is a coordinate rotation — no floating-point multiplications.
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

/*
 * =====================================================================================================================
 * Pure-state Pauli application
 * =====================================================================================================================
 */

/**
 * @brief Apply Pauli string P to pure state in-place: |psi> -> P|psi>
 *
 * Two cases dispatched by mask structure:
 *
 * 1. Diagonal (x_mask == 0): identity returns immediately; pure Z-string
 *    applies real sign flips (-1)^{popcount(i & z_mask)}.
 *
 * 2. Non-diagonal (x_mask != 0): swap + phase on pairs (j, j XOR x_mask).
 *    Enumerate j by inserting 0 at the LOWEST set bit of x_mask.
 *    This gives exactly dim/2 unique pairs, each (j, jp) differing at all
 *    x_mask bit positions.
 */
void pauli_apply_pure(state_t *state, const pauli_t *P) {
    const idx_t dim = (idx_t)1 << state->qubits;
    cplx_t *restrict data = state->data;

    const uint64_t x_mask = P->x_mask;
    const uint64_t z_mask = P->z_mask;
    const uint64_t y_mask = P->y_mask;

    if (x_mask == 0) {
        if (z_mask == 0) return;  /* identity: no-op */
        /* Pure Z-string (y_mask == 0 guaranteed): real sign flips only */
        #pragma omp parallel for if(dim >= OMP_THRESHOLD)
        for (idx_t i = 0; i < dim; ++i) {
            const unsigned parity = (unsigned)(__builtin_popcountll(i & z_mask) & 1u);
            const double sign = 1.0 - 2.0 * (double)parity;
            data[i] = data[i] * sign;
        }
    } else {
        const unsigned eta_exp   = (unsigned)(__builtin_popcountll(y_mask) & 3u);
        const unsigned y_parity  = (unsigned)(__builtin_popcountll(y_mask) & 1u);
        const qubit_t  split_bit = (qubit_t)__builtin_ctzll(x_mask);
        const idx_t n_pairs = dim >> 1;

        #pragma omp parallel for if(dim >= OMP_THRESHOLD)
        for (idx_t k = 0; k < n_pairs; ++k) {
            const idx_t j  = insertBit0(k, split_bit);
            const idx_t jp = j ^ x_mask;

            const unsigned parity_j = (unsigned)(__builtin_popcountll(j & z_mask) & 1u);
            const unsigned exp_j  = (eta_exp + 2u * parity_j) & 3u;
            const unsigned exp_jp = (exp_j + 2u * y_parity) & 3u;

            const cplx_t tmp = data[j];
            data[j]  = apply_ipow(data[jp], exp_jp);
            data[jp] = apply_ipow(tmp,      exp_j);
        }
    }
}

/*
 * =====================================================================================================================
 * Pure-state Pauli exponential
 * =====================================================================================================================
 */

/**
 * @brief Apply U(theta) = exp(-i*theta/2 * P) to a pure state in-place
 *
 * Uses the Euler decomposition:
 *   U(theta) = cos(theta/2) * I  -  i * sin(theta/2) * P
 *
 * Three cases dispatched by mask structure:
 *
 * 1. Identity (x_mask == 0 && z_mask == 0): global phase e^{-i*theta/2}.
 * 2. Diagonal (x_mask == 0): per-amplitude phase by eigenvalue parity.
 * 3. Non-diagonal (x_mask != 0): paired update; both old values read before write.
 */
void pauli_exp_pure(state_t *state, const pauli_t *P, double theta) {
    const idx_t dim = (idx_t)1 << state->qubits;
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
            for (idx_t i = 0; i < dim; ++i) {
                data[i] = data[i] * phase;
            }
            return;
        }
        /* Pure Z-string: two precomputed scalars e^{±i*theta/2} selected by parity */
        const cplx_t fac_pos = CMPLX(c, -s);   /* epsilon = +1: e^{-i*theta/2} */
        const cplx_t fac_neg = CMPLX(c,  s);   /* epsilon = -1: e^{+i*theta/2} */
        #pragma omp parallel for if(dim >= OMP_THRESHOLD)
        for (idx_t i = 0; i < dim; ++i) {
            const unsigned parity = (unsigned)(__builtin_popcountll(i & z_mask) & 1u);
            data[i] = data[i] * (parity ? fac_neg : fac_pos);
        }
    } else {
        /* Non-diagonal: paired update on (j, jp = j XOR x_mask).
         * Both old values must be read before either write.
         * new_j  = c*old_j  + s*apply_ipow(old_jp, (exp_jp+3)&3)
         * new_jp = c*old_jp + s*apply_ipow(old_j,  (exp_j +3)&3)
         */
        const unsigned eta_exp   = (unsigned)(__builtin_popcountll(y_mask) & 3u);
        const unsigned y_parity  = (unsigned)(__builtin_popcountll(y_mask) & 1u);
        const qubit_t  split_bit = (qubit_t)__builtin_ctzll(x_mask);
        const idx_t n_pairs = dim >> 1;

        #pragma omp parallel for if(dim >= OMP_THRESHOLD)
        for (idx_t k = 0; k < n_pairs; ++k) {
            const idx_t j  = insertBit0(k, split_bit);
            const idx_t jp = j ^ x_mask;

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
