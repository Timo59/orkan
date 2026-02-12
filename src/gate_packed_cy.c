/**
 * @file gate_packed_cy.c
 * @brief CY (controlled-Y) gate for mixed states in packed lower-triangular storage
 *
 * Same permutation as CX with per-element phase factors.
 * CY|10> = i|11>, CY|11> = -i|10>.  rho'[a,b] = phi(a)*phi*(b)*rho[perm(a), perm(b)]
 * where phi(00)=phi(01)=1, phi(10)=-i, phi(11)=i.
 *
 * Per-pair phases (same 6 pairs as CX):
 *   A (r10,c00 <-> r11,c00): x(-i, i)     E (r10,c10 <-> r11,c11): plain swap
 *   B (r11,c01 <-> r10,c01): x(i, -i)     F (r11,c10 <-> r10,c11): negate-swap
 *   C,D (r00/r01,c10 <-> c11): x(i, -i)
 *
 * Pitfall: when a pair crosses the diagonal, phase and conj() compose
 * differently.  E.g. xi on lower-tri value v: CMPLX(-Im(v), Re(v)); but
 * xi on upper-tri stored s (where rho = conj(s)): result stored as
 * CMPLX(Im(s), Re(s)), not CMPLX(-Im(s), Re(s)).
 */

#include "gate.h"
#include <complex.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#define OMP_THRESHOLD 64

void cy_packed(state_t * restrict state, const qubit_t control, const qubit_t target) {
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

        gate_idx_t offset_c00 = c00 * (2 * dim - c00 + 1) / 2;
        gate_idx_t offset_c01 = c01 * (2 * dim - c01 + 1) / 2;
        gate_idx_t offset_c10 = c10 * (2 * dim - c10 + 1) / 2;
        gate_idx_t offset_c11 = c11 * (2 * dim - c11 + 1) / 2;

        for (gate_idx_t br = bc; br < quarter_dim; ++br) {
            const gate_idx_t r00 = insertBits2_0(br, lo, hi);
            const gate_idx_t r01 = r00 | incr_tgt, r10 = r00 | incr_ctrl, r11 = r10 | incr_tgt;

            /* Pair A: (r10,c00) <-> (r11,c00) — always lower tri; phases (-i, i) */
            cplx_t tmp = data[(r10 - c00) + offset_c00];
            data[(r10 - c00) + offset_c00] = CMPLX(cimag(data[(r11 - c00) + offset_c00]), -creal(data[(r11 - c00) + offset_c00]));
            data[(r11 - c00) + offset_c00] = CMPLX(-cimag(tmp), creal(tmp));

            /* Pair E: (r10,c10) <-> (r11,c11) — always lower tri; plain swap */
            tmp = data[(r10 - c10) + offset_c10];
            data[(r10 - c10) + offset_c10] = data[(r11 - c11) + offset_c11];
            data[(r11 - c11) + offset_c11] = tmp;

            /* Fast path: all remaining pairs in lower triangle */
            if (r00 > c11) {
                /* Pair B: (r11,c01) <-> (r10,c01); phases (i, -i) */
                tmp = data[(r11 - c01) + offset_c01];
                data[(r11 - c01) + offset_c01] = CMPLX(-cimag(data[(r10 - c01) + offset_c01]), creal(data[(r10 - c01) + offset_c01]));
                data[(r10 - c01) + offset_c01] = CMPLX(cimag(tmp), -creal(tmp));

                /* Pair F: (r11,c10) <-> (r10,c11); negate-swap */
                tmp = data[(r11 - c10) + offset_c10];
                data[(r11 - c10) + offset_c10] = -data[(r10 - c11) + offset_c11];
                data[(r10 - c11) + offset_c11] = -tmp;

                /* Pair D: (r01,c10) <-> (r01,c11); phases (i, -i) */
                tmp = data[(r01 - c10) + offset_c10];
                data[(r01 - c10) + offset_c10] = CMPLX(-cimag(data[(r01 - c11) + offset_c11]), creal(data[(r01 - c11) + offset_c11]));
                data[(r01 - c11) + offset_c11] = CMPLX(cimag(tmp), -creal(tmp));

                /* Pair C: (r00,c10) <-> (r00,c11); phases (i, -i) */
                tmp = data[(r00 - c10) + offset_c10];
                data[(r00 - c10) + offset_c10] = CMPLX(-cimag(data[(r00 - c11) + offset_c11]), creal(data[(r00 - c11) + offset_c11]));
                data[(r00 - c11) + offset_c11] = CMPLX(cimag(tmp), -creal(tmp));
            }
            else {
                if (r10 > c11) {
                    /* Pair F: negate-swap, both lower tri */
                    tmp = data[(r11 - c10) + offset_c10];
                    data[(r11 - c10) + offset_c10] = -data[(r10 - c11) + offset_c11];
                    data[(r10 - c11) + offset_c11] = -tmp;

                    /* Pair B: phases (i, -i), both lower tri */
                    tmp = data[(r11 - c01) + offset_c01];
                    data[(r11 - c01) + offset_c01] = CMPLX(-cimag(data[(r10 - c01) + offset_c01]), creal(data[(r10 - c01) + offset_c01]));
                    data[(r10 - c01) + offset_c01] = CMPLX(cimag(tmp), -creal(tmp));
                }
                else {
                    /* Pair F: (r11,c10) lower, (r10,c11) upper — negate + conj */
                    tmp = data[(r11 - c10) + offset_c10];
                    data[(r11 - c10) + offset_c10] = -conj(data[pack_idx(dim, c11, r10)]);
                    data[pack_idx(dim, c11, r10)] = -conj(tmp);

                    /* Pair B: (r11,c01) lower; (r10,c01) may cross diagonal */
                    tmp = data[(r11 - c01) + offset_c01];
                    if (r10 > c01) {
                        data[(r11 - c01) + offset_c01] = CMPLX(-cimag(data[(r10 - c01) + offset_c01]), creal(data[(r10 - c01) + offset_c01]));
                        data[(r10 - c01) + offset_c01] = CMPLX(cimag(tmp), -creal(tmp));
                    }
                    else {
                        data[(r11 - c01) + offset_c01] = CMPLX(cimag(data[pack_idx(dim, c01, r10)]), creal(data[pack_idx(dim, c01, r10)]));
                        data[pack_idx(dim, c01, r10)] = CMPLX(cimag(tmp), creal(tmp));
                    }
                }

                if (bc != br) {
                    if (r00 > c10) {
                        /* Pair C: (r00,c10) lower, (r00,c11) upper */
                        /* r00 < c11, otherwise fast path would have been taken */
                        tmp = data[(r00 - c10) + offset_c10];
                        data[(r00 - c10) + offset_c10] = CMPLX(cimag(data[pack_idx(dim, c11, r00)]), creal(data[pack_idx(dim, c11, r00)]));
                        data[pack_idx(dim, c11, r00)] = CMPLX(cimag(tmp), creal(tmp));

                        /* Pair D: both lower tri (r00 > c10 implies r01 > c11) */
                        tmp = data[(r01 - c10) + offset_c10];
                        data[(r01 - c10) + offset_c10] = CMPLX(-cimag(data[(r01 - c11) + offset_c11]), creal(data[(r01 - c11) + offset_c11]));
                        data[(r01 - c11) + offset_c11] = CMPLX(cimag(tmp), -creal(tmp));
                    }
                    else {
                        /* Pair C: both upper tri */
                        tmp = data[pack_idx(dim, c10, r00)];
                        data[pack_idx(dim, c10, r00)] = CMPLX(cimag(data[pack_idx(dim, c11, r00)]), -creal(data[pack_idx(dim, c11, r00)]));
                        data[pack_idx(dim, c11, r00)] = CMPLX(-cimag(tmp), creal(tmp));

                        /* Pair D: (r01,c11) always upper; (r01,c10) may cross */
                        /* r00 < c10 implies r01 < c11 */
                        tmp = data[pack_idx(dim, c11, r01)];
                        if (r01 > c10) {
                            data[pack_idx(dim, c11, r01)] = CMPLX(cimag(data[(r01 - c10) + offset_c10]), creal(data[(r01 - c10) + offset_c10]));
                            data[(r01 - c10) + offset_c10] = CMPLX(cimag(tmp), creal(tmp));
                        }
                        else {
                            data[pack_idx(dim, c11, r01)] = CMPLX(-cimag(data[pack_idx(dim, c10, r01)]), creal(data[pack_idx(dim, c10, r01)]));
                            data[pack_idx(dim, c10, r01)] = CMPLX(cimag(tmp), -creal(tmp));
                        }
                    }
                }
            }
        }
    }
}
