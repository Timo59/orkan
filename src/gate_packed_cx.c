/**
 * @file gate_packed_cx.c
 * @brief CX (CNOT) gate for mixed states in packed lower-triangular storage
 *
 * CX is a zero-flop permutation gate (self-adjoint): ρ' = CX·ρ·CX.
 * perm(x) = x ^ (1<<target) when bit control is set, else x.
 * ρ'[r,c] = ρ[perm(r), perm(c)].  Involution: each element is fixed or in a 2-cycle.
 *
 * For each block (br, bc) the 6 swap pairs split into:
 *   - Fast path (r00 > c11): all lower-tri, plain swaps
 *   - Slow path: up to 4 pairs may cross the diagonal, requiring conjugation
 *   - Diagonal blocks (bc == br): skip group B (Hermitian pairs share storage)
 */

#include "gate.h"
#include <complex.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#define OMP_THRESHOLD 64

void cx_packed(state_t * restrict state, const qubit_t control, const qubit_t target) {
    const gate_idx_t dim = (gate_idx_t)1 << state->qubits;
    const qubit_t hi = (control > target) ? control : target;
    const qubit_t lo = (control < target) ? control : target;
    const gate_idx_t incr_ctrl = (gate_idx_t)1 << control;
    const gate_idx_t incr_tgt = (gate_idx_t)1 << target;
    cplx_t * restrict data = state->data;
    const gate_idx_t quarter_dim = dim >> 2;

    #pragma omp parallel for schedule(static, 1) if(dim >= OMP_THRESHOLD)
    for (gate_idx_t bc = 0; bc < quarter_dim; ++bc) {
        const gate_idx_t c00 = insertBits2_0(bc, lo, hi);
        const gate_idx_t c01 = c00 | incr_tgt, c10 = c00 | incr_ctrl, c11 = c10 | incr_tgt;

        /* Precompute offsets of columns */
        gate_idx_t offset_c00 = c00 * (2 * dim - c00 + 1) / 2;
        gate_idx_t offset_c01 = c01 * (2 * dim - c01 + 1) / 2;
        gate_idx_t offset_c10 = c10 * (2 * dim - c10 + 1) / 2;
        gate_idx_t offset_c11 = c11 * (2 * dim - c11 + 1) / 2;

        for (gate_idx_t br = bc; br < quarter_dim; ++br) {
            const gate_idx_t r00 = insertBits2_0(br, lo, hi);
            const gate_idx_t r01 = r00 | incr_tgt, r10 = r00 | incr_ctrl, r11 = r10 | incr_tgt;

            /* The first two pairs are always in the lower triangle */
            /* |c=1,t=0><c=0,t=0| <--> |c=1,t=1><c=0,t=0| */
            cplx_t tmp = data[(r10 - c00) + offset_c00];
            data[(r10 - c00) + offset_c00] = data[(r11 - c00) + offset_c00];
            data[(r11 - c00) + offset_c00] = tmp;

            /* |c=1,t=0><c=1,t=0| <--> |c=1,t=1><c=1,t=1| */
            tmp = data[(r10 - c10) + offset_c10];
            data[(r10 - c10) + offset_c10] = data[(r11 - c11) + offset_c11];
            data[(r11 - c11) + offset_c11] = tmp;

            /* Fast-path: Remaining 4 pairs are in the lower triangle. */
            if (r00 > c11) {
                /* |c=1,t=1><c=0,t=1| <--> |c=1,t=0><c=0,t=1| */
                tmp = data[(r11 - c01) + offset_c01];
                data[(r11 - c01) + offset_c01] = data[(r10 - c01) + offset_c01];
                data[(r10 - c01) + offset_c01] = tmp;

                /* |c=1,t=1><c=1,t=0| <--> |c=1,t=0><c=1,t=1| */
                tmp = data[(r11 - c10) + offset_c10];
                data[(r11 - c10) + offset_c10] = data[(r10 - c11) + offset_c11];
                data[(r10 - c11) + offset_c11] = tmp;

                /* |c=0,t=1><c=1,t=0| <--> |c=0,t=1><c=1,t=1| */
                tmp = data[(r01 - c10) + offset_c10];
                data[(r01 - c10) + offset_c10] = data[(r01 - c11) + offset_c11];
                data[(r01 - c11) + offset_c11] = tmp;

                /* |c=0,t=0><c=1,t=0| <--> |c=0,t=0><c=1,t=1| */
                tmp = data[(r00 - c10) + offset_c10];
                data[(r00 - c10) + offset_c10] = data[(r00 - c11) + offset_c11];
                data[(r00 - c11) + offset_c11] = tmp;
            }
            else {
                /* r10 > c11 implies r10 > c01 */
                if (r10 > c11) {
                    /* |c=1,t=1><c=1,t=0| <--> |c=1,t=0><c=1,t=1| */
                    tmp = data[(r11 - c10) + offset_c10];
                    data[(r11 - c10) + offset_c10] = data[(r10 - c11) + offset_c11];
                    data[(r10 - c11) + offset_c11] = tmp;

                    /* |c=1,t=1><c=0,t=1| <--> |c=1,t=0><c=0,t=1| */
                    tmp = data[(r11 - c01) + offset_c01];
                    data[(r11 - c01) + offset_c01] = data[(r10 - c01) + offset_c01];
                    data[(r10 - c01) + offset_c01] = tmp;
                }
                else {
                    /* |c=1,t=1><c=1,t=0| <--> |c=1,t=0><c=1,t=1| */
                    tmp = data[(r11 - c10) + offset_c10];
                    data[(r11 - c10) + offset_c10] = conj(data[pack_idx(dim, c11, r10)]);
                    data[pack_idx(dim, c11, r10)] = conj(tmp);

                    /* |c=1,t=1><c=0,t=1| <--> |c=1,t=0><c=0,t=1| */
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
                    /* r00 > c10 implies r01 > c10 */
                    if (r00 > c10) {
                        /* |c=0,t=0><c=1,t=0| <--> |c=0,t=0><c=1,t=1| */
                        /* r00 < c11, otherwise fast path would have been taken */
                        /* (r00, c11) is in the upper triangle => conj(c11, r00) */
                        tmp = data[(r00 - c10) + offset_c10];
                        data[(r00 - c10) + offset_c10] = conj(data[pack_idx(dim, c11, r00)]);
                        data[pack_idx(dim, c11, r00)] = conj(tmp);

                        /* |c=0,t=1><c=1,t=0| <--> |c=0,t=1><c=1,t=1| */
                        /* r00 > c10 also implies r01 > c11 */
                        tmp = data[(r01 - c10) + offset_c10];
                        data[(r01 - c10) + offset_c10] = data[(r01 - c11) + offset_c11];
                        data[(r01 - c11) + offset_c11] = tmp;
                    }
                    else {
                        /* |c=0,t=0><c=1,t=0| <--> |c=0,t=0><c=1,t=1| */
                        /* r00 < c10 implies both elements are in the upper triangle of a lower triangle block */
                        tmp = data[pack_idx(dim, c10, r00)];
                        data[pack_idx(dim, c10, r00)] = data[pack_idx(dim, c11, r00)];
                        data[pack_idx(dim, c11, r00)] = tmp;

                        /* |c=0,t=1><c=1,t=0| <--> |c=0,t=1><c=1,t=1| */
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
