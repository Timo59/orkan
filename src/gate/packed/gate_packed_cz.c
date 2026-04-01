/**
 * @file gate_packed_cz.c
 * @brief CZ gate for mixed states in packed lower-triangular storage
 *
 * CZ: diag(1,1,1,-1), zero-flop diagonal gate — only sign flips.
 *
 * Negate elements where exactly one of (row, col) has both ctrl+tgt bits set:
 *   Group A (row=11): (r11,c00), (r11,c01), (r11,c10) — always lower-tri
 *   Group B (col=11): (r00,c11), (r01,c11), (r10,c11) — may cross diagonal
 *
 * Diagonal blocks (bc == br): skip group B (same storage as A, avoids double-negate).
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

void cz_packed(state_t *state, const qubit_t control, const qubit_t target) {
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

        idx_t offset_c00 = c00 * (2 * dim - c00 + 1) / 2;
        idx_t offset_c01 = c01 * (2 * dim - c01 + 1) / 2;
        idx_t offset_c10 = c10 * (2 * dim - c10 + 1) / 2;
        idx_t offset_c11 = c11 * (2 * dim - c11 + 1) / 2;

        for (idx_t br = bc; br < quarter_dim; ++br) {
            const idx_t r00 = insertBits2_0(br, lo, hi);
            const idx_t r01 = r00 | incr_tgt, r10 = r00 | incr_ctrl, r11 = r10 | incr_tgt;

            /* Group A: row = 11 — always in lower triangle */
            data[(r11 - c00) + offset_c00] = -data[(r11 - c00) + offset_c00];
            data[(r11 - c01) + offset_c01] = -data[(r11 - c01) + offset_c01];
            data[(r11 - c10) + offset_c10] = -data[(r11 - c10) + offset_c10];

            /* Group B: col = 11 — skip on diagonal block (Hermitian pairs already negated) */
            if (bc != br) {
                if (r00 > c11) {
                    /* Fast path: all in lower triangle */
                    data[(r00 - c11) + offset_c11] = -data[(r00 - c11) + offset_c11];
                    data[(r01 - c11) + offset_c11] = -data[(r01 - c11) + offset_c11];
                    data[(r10 - c11) + offset_c11] = -data[(r10 - c11) + offset_c11];
                }
                else {
                    /* (r00, c11): always upper tri here — negate stored at (c11, r00) */
                    data[pack_idx(dim, c11, r00)] = -data[pack_idx(dim, c11, r00)];

                    /* (r01, c11) */
                    if (r01 > c11) {
                        data[(r01 - c11) + offset_c11] = -data[(r01 - c11) + offset_c11];
                    } else {
                        data[pack_idx(dim, c11, r01)] = -data[pack_idx(dim, c11, r01)];
                    }

                    /* (r10, c11) */
                    if (r10 > c11) {
                        data[(r10 - c11) + offset_c11] = -data[(r10 - c11) + offset_c11];
                    } else {
                        data[pack_idx(dim, c11, r10)] = -data[pack_idx(dim, c11, r10)];
                    }
                }
            }
        }
    }
}
