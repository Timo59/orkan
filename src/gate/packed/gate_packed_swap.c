/**
 * @file gate_packed_swap.c
 * @brief SWAP gate for mixed states in packed lower-triangular storage
 *
 * SWAP: zero-flop permutation gate.  perm(x) = x with bits lo and hi exchanged.
 * rho'[r,c] = rho[perm(r), perm(c)].
 *
 * Per-pair operations (6 swap pairs):
 *   A (r01,c00 <-> r10,c00): plain swap    D (r10,c01 <-> r01,c10): plain swap
 *   B (r01,c01 <-> r10,c10): plain swap    E (r00,c01 <-> r00,c10): plain swap
 *   C (r11,c01 <-> r11,c10): plain swap    F (r10,c11 <-> r01,c11): plain swap
 *
 * For each block (br, bc) the 6 swap pairs split into:
 *   - Fast path (r00 > c10): all lower-tri, plain swaps
 *   - Slow path: up to 3 pairs may cross the diagonal, requiring conjugation
 *   - Diagonal blocks (bc == br): skip pairs E,F when both upper tri
 */

#include "gate.h"
#include <complex.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#define OMP_THRESHOLD 512

void swap_packed(state_t *state, const qubit_t q1, const qubit_t q2) {
    const gate_idx_t dim = (gate_idx_t)1 << state->qubits;
    const qubit_t hi = q1 > q2 ? q1 : q2;
    const qubit_t lo = q1 < q2 ? q1 : q2;
    const gate_idx_t incr_hi = (gate_idx_t)1 << hi;
    const gate_idx_t incr_lo = (gate_idx_t)1 << lo;
    cplx_t * restrict data = state->data;
    const gate_idx_t quarter_dim = dim >> 2;

    #pragma omp parallel for schedule(static, 1) if(dim >= OMP_THRESHOLD)
    for (gate_idx_t bc = 0; bc < quarter_dim; ++bc) {
        const gate_idx_t c00 = insertBits2_0(bc, lo, hi);
        const gate_idx_t c01 = c00 | incr_lo, c10 = c00 | incr_hi, c11 = c10 | incr_lo;

        /* Precompute offsets of columns */
        gate_idx_t offset_c00 = c00 * (2 * dim - c00 + 1) / 2;
        gate_idx_t offset_c01 = c01 * (2 * dim - c01 + 1) / 2;
        gate_idx_t offset_c10 = c10 * (2 * dim - c10 + 1) / 2;
        gate_idx_t offset_c11 = c11 * (2 * dim - c11 + 1) / 2;


        for (gate_idx_t br = bc; br < quarter_dim; ++br) {
            const gate_idx_t r00 = insertBits2_0(br, lo, hi);
            const gate_idx_t r01 = r00 | incr_lo, r10 = r00 | incr_hi, r11 = r10 | incr_lo;

            /* Pair A: (r01,c00) <-> (r10,c00) — always lower tri; plain swap */
            cplx_t tmp = data[(r01 - c00) + offset_c00];
            data[(r01 - c00) + offset_c00] = data[(r10 - c00) + offset_c00];
            data[(r10 - c00) + offset_c00] = tmp;

            /* Pair B: (r01,c01) <-> (r10,c10) — always lower tri; plain swap */
            tmp = data[(r01 - c01) + offset_c01];
            data[(r01 - c01) + offset_c01] = data[(r10 - c10) + offset_c10];
            data[(r10 - c10) + offset_c10] = tmp;

            /* Pair C: (r11,c01) <-> (r11,c10) — always lower tri; plain swap */
            tmp = data[(r11 - c01) + offset_c01];
            data[(r11 - c01) + offset_c01] = data[(r11 - c10) + offset_c10];
            data[(r11 - c10) + offset_c10] = tmp;

            /* Fast path: all remaining pairs in lower triangle */
            if (r00 > c10) {
                /* Pair E: (r00,c01) <-> (r00,c10); plain swap */
                tmp = data[(r00 - c01) + offset_c01];
                data[(r00 - c01) + offset_c01] = data[(r00 - c10) + offset_c10];
                data[(r00 - c10) + offset_c10] = tmp;

                /* Pair F: (r01,c11) <-> (r10,c11); plain swap */
                tmp = data[(r01 - c11) + offset_c11];
                data[(r01 - c11) + offset_c11] = data[(r10 - c11) + offset_c11];
                data[(r10 - c11) + offset_c11] = tmp;

                /* Pair D: (r10,c01) <-> (r01,c10); plain swap */
                tmp = data[(r01 - c10) + offset_c10];
                data[(r01 - c10) + offset_c10] = data[(r10 - c01) + offset_c01];
                data[(r10 - c01) + offset_c01] = tmp;
            }
            else {
                /* Pair D: (r10,c01) lower; (r01,c10) may cross diagonal */
                tmp = data[(r10 - c01) + offset_c01];
                if (r01 > c10) {
                    data[(r10 - c01) + offset_c01] = data[(r01 - c10) + offset_c10];
                    data[(r01 - c10) + offset_c10] = tmp;
                }
                else {
                    data[(r10 - c01) + offset_c01] = conj(data[pack_idx(dim, c10, r01)]);
                    data[pack_idx(dim, c10, r01)] = conj(tmp);
                }

                if (r00 > c01) {
                    /* Pair E: (r00,c01) lower, (r00,c10) upper — conj swap */
                    tmp = data[(r00 - c01) + offset_c01];
                    data[(r00 - c01) + offset_c01] = conj(data[pack_idx(dim, c10, r00)]);
                    data[pack_idx(dim, c10, r00)] = conj(tmp);

                    /* Pair F: (r10,c11) lower, (r01,c11) upper — conj swap */
                    tmp = data[(r10 - c11) + offset_c11];
                    data[(r10 - c11) + offset_c11] = conj(data[pack_idx(dim, c11, r01)]);
                    data[pack_idx(dim, c11, r01)] = conj(tmp);
                }
                else if (bc != br) {
                    /* Pair E: (r00,c01) <-> (r00,c10) — both upper tri */
                    tmp = data[pack_idx(dim, c01, r00)];
                    data[pack_idx(dim, c01, r00)] = data[pack_idx(dim, c10, r00)];
                    data[pack_idx(dim, c10, r00)] = tmp;

                    /* Pair F: (r01,c11) <-> (r10,c11) — both upper tri */
                    tmp = data[pack_idx(dim, c11, r10)];
                    data[pack_idx(dim, c11, r10)] = data[pack_idx(dim, c11, r01)];
                    data[pack_idx(dim, c11, r01)] = tmp;
                }
            }
        }
    }
}
