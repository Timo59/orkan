/**
 * @file gate_packed_cx.c
 * @brief CX (CNOT) gate for mixed states in packed lower-triangular storage
 *
 * CX is a zero-flop permutation gate (self-adjoint): ρ' = CX·ρ·CX.
 * perm(x) = x ^ (1<<target) when bit control is set, else x.
 * ρ'[r,c] = ρ[perm(r), perm(c)].  Involution: each element is fixed or in a 2-cycle.
 *
 * Per-pair operations (6 swap pairs):
 *   A (r10,c00 <-> r11,c00): plain swap    E (r10,c10 <-> r11,c11): plain swap
 *   B (r11,c01 <-> r10,c01): plain swap    F (r11,c10 <-> r10,c11): plain swap
 *   C,D (r00/r01,c10 <-> c11): plain swap
 *
 * For each block (br, bc) the 6 swap pairs split into:
 *   - Fast path (r00 > c11): all lower-tri, plain swaps
 *   - Slow path: up to 4 pairs may cross the diagonal, requiring conjugation
 *   - Diagonal blocks (bc == br): skip pairs C,D (Hermitian pairs share storage)
 */

#include "gate.h"
#include "index.h"
#include <complex.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef OMP_THRESHOLD
#define OMP_THRESHOLD 512
#endif

void cx_packed(state_t * restrict state, const qubit_t control, const qubit_t target) {
    const idx_t dim = (idx_t)1 << state->qubits;
    const qubit_t hi = (control > target) ? control : target;
    const qubit_t lo = (control < target) ? control : target;
    const idx_t incr_ctrl = (idx_t)1 << control;
    const idx_t incr_tgt = (idx_t)1 << target;
    cplx_t * restrict data = state->data;
    const idx_t quarter_dim = dim >> 2;

    #pragma omp parallel for schedule(static, 1) if(dim >= OMP_THRESHOLD)
    for (idx_t bc = 0; bc < quarter_dim; ++bc) {
        const idx_t c00 = insertBits2_0(bc, lo, hi);
        const idx_t c01 = c00 | incr_tgt, c10 = c00 | incr_ctrl, c11 = c10 | incr_tgt;

        /* Precompute offsets of columns */
        idx_t offset_c00 = c00 * (2 * dim - c00 + 1) / 2;
        idx_t offset_c01 = c01 * (2 * dim - c01 + 1) / 2;
        idx_t offset_c10 = c10 * (2 * dim - c10 + 1) / 2;
        idx_t offset_c11 = c11 * (2 * dim - c11 + 1) / 2;

        for (idx_t br = bc; br < quarter_dim; ++br) {
            const idx_t r00 = insertBits2_0(br, lo, hi);
            const idx_t r01 = r00 | incr_tgt, r10 = r00 | incr_ctrl, r11 = r10 | incr_tgt;

            /* Pair A: (r10,c00) <-> (r11,c00) — always lower tri; plain swap */
            cplx_t tmp = data[(r10 - c00) + offset_c00];
            data[(r10 - c00) + offset_c00] = data[(r11 - c00) + offset_c00];
            data[(r11 - c00) + offset_c00] = tmp;

            /* Pair E: (r10,c10) <-> (r11,c11) — always lower tri; plain swap */
            tmp = data[(r10 - c10) + offset_c10];
            data[(r10 - c10) + offset_c10] = data[(r11 - c11) + offset_c11];
            data[(r11 - c11) + offset_c11] = tmp;

            /* Fast path: all remaining pairs in lower triangle */
            if (r00 > c11) {
                /* Pair B: (r11,c01) <-> (r10,c01); plain swap */
                tmp = data[(r11 - c01) + offset_c01];
                data[(r11 - c01) + offset_c01] = data[(r10 - c01) + offset_c01];
                data[(r10 - c01) + offset_c01] = tmp;

                /* Pair F: (r11,c10) <-> (r10,c11); plain swap */
                tmp = data[(r11 - c10) + offset_c10];
                data[(r11 - c10) + offset_c10] = data[(r10 - c11) + offset_c11];
                data[(r10 - c11) + offset_c11] = tmp;

                /* Pair D: (r01,c10) <-> (r01,c11); plain swap */
                tmp = data[(r01 - c10) + offset_c10];
                data[(r01 - c10) + offset_c10] = data[(r01 - c11) + offset_c11];
                data[(r01 - c11) + offset_c11] = tmp;

                /* Pair C: (r00,c10) <-> (r00,c11); plain swap */
                tmp = data[(r00 - c10) + offset_c10];
                data[(r00 - c10) + offset_c10] = data[(r00 - c11) + offset_c11];
                data[(r00 - c11) + offset_c11] = tmp;
            }
            else {
                if (r10 > c11) {
                    /* Pair F: plain swap, both lower tri */
                    tmp = data[(r11 - c10) + offset_c10];
                    data[(r11 - c10) + offset_c10] = data[(r10 - c11) + offset_c11];
                    data[(r10 - c11) + offset_c11] = tmp;

                    /* Pair B: plain swap, both lower tri */
                    tmp = data[(r11 - c01) + offset_c01];
                    data[(r11 - c01) + offset_c01] = data[(r10 - c01) + offset_c01];
                    data[(r10 - c01) + offset_c01] = tmp;
                }
                else {
                    /* Pair F: (r11,c10) lower, (r10,c11) upper — conj swap */
                    tmp = data[(r11 - c10) + offset_c10];
                    data[(r11 - c10) + offset_c10] = conj(data[pack_idx(dim, c11, r10)]);
                    data[pack_idx(dim, c11, r10)] = conj(tmp);

                    /* Pair B: (r11,c01) lower; (r10,c01) may cross diagonal */
                    tmp = data[(r11 - c01) + offset_c01];
                    if (r10 > c01) {
                        data[(r11 - c01) + offset_c01] = data[(r10 - c01) + offset_c01];
                        data[(r10 - c01) + offset_c01] = tmp;
                    }
                    else {
                        data[(r11 - c01) + offset_c01] = conj(data[pack_idx(dim, c01, r10)]);
                        data[pack_idx(dim, c01, r10)] = conj(tmp);
                    }
                }

                if (bc != br) {
                    if (r00 > c10) {
                        /* Pair C: (r00,c10) lower, (r00,c11) upper */
                        /* r00 < c11, otherwise fast path would have been taken */
                        tmp = data[(r00 - c10) + offset_c10];
                        data[(r00 - c10) + offset_c10] = conj(data[pack_idx(dim, c11, r00)]);
                        data[pack_idx(dim, c11, r00)] = conj(tmp);

                        /* Pair D: both lower tri (r00 > c10 implies r01 > c11) */
                        tmp = data[(r01 - c10) + offset_c10];
                        data[(r01 - c10) + offset_c10] = data[(r01 - c11) + offset_c11];
                        data[(r01 - c11) + offset_c11] = tmp;
                    }
                    else {
                        /* Pair C: both upper tri */
                        tmp = data[pack_idx(dim, c10, r00)];
                        data[pack_idx(dim, c10, r00)] = data[pack_idx(dim, c11, r00)];
                        data[pack_idx(dim, c11, r00)] = tmp;

                        /* Pair D: (r01,c11) always upper; (r01,c10) may cross */
                        /* r00 < c10 implies r01 < c11 */
                        tmp = data[pack_idx(dim, c11, r01)];
                        if (r01 > c10) {
                            data[pack_idx(dim, c11, r01)] = conj(data[(r01 - c10) + offset_c10]);
                            data[(r01 - c10) + offset_c10] = conj(tmp);
                        }
                        else {
                            data[pack_idx(dim, c11, r01)] = data[pack_idx(dim, c10, r01)];
                            data[pack_idx(dim, c10, r01)] = tmp;
                        }
                    }
                }
            }
        }
    }
}
