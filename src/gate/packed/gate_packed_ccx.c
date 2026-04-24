/**
 * @file gate_packed_ccx.c
 * @brief CCX (Toffoli) gate for mixed states in packed lower-triangular storage
 *
 * CCX is a zero-flop permutation gate (self-adjoint): ρ' = CCX·ρ·CCX.
 * perm(x) = x ^ (1<<target) when bits ctrl1 AND ctrl2 are both set, else x.
 * ρ'[r,c] = ρ[perm(r), perm(c)].  Involution: each element is fixed or in a 2-cycle.
 *
 * Per-pair operations (14 swap pairs):
 *   A (r110,c000 <-> r111,c000): plain swap   E (r111,c001 <-> r110,c001): plain swap
 *   B (r110,c010 <-> r111,c010): plain swap   F (r111,c011 <-> r110,c011): plain swap
 *   C (r110,c100 <-> r111,c100): plain swap   G (r111,c101 <-> r110,c101): plain swap
 *   D (r110,c110 <-> r111,c111): plain swap   H (r111,c110 <-> r110,c111): plain swap
 *   I-N (rXXX,c110 <-> rXXX,c111): plain swap, rXXX in {r000..r101}
 *
 * For each block (br, bc) the 14 swap pairs split into:
 *   - Fast path (r000 > c111): all lower-tri, plain swaps
 *   - Slow path: pairs E-N may cross the diagonal, requiring conjugation
 *   - Diagonal blocks (bc == br): skip pairs I-N (Hermitian pairs share storage)
 */

#include "gate.h"
#include "index.h"
#include <complex.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/* Minimum dimension to enable OpenMP parallelization (avoid thread overhead for small systems) */


void ccx_packed(state_t *state, const qubit_t ctrl1, const qubit_t ctrl2, const qubit_t target) {
    const idx_t dim = (idx_t)1 << state->qubits;
    cplx_t * restrict data = state->data;
    const idx_t n_base = dim >> 3;

    /* Mutable copies for potential swap */
    qubit_t c1 = ctrl1, c2 = ctrl2;

    /* Enforce c1 > c2 */
    if (c1 < c2) {
        const qubit_t tmp = c1;
        c1 = c2;
        c2 = tmp;
    }
    idx_t incr_ctrl1 = (idx_t)1 << c1;
    idx_t incr_ctrl2 = (idx_t)1 << c2;
    idx_t incr_target = (idx_t)1 << target;

    qubit_t hi, med, lo;
    if (c2 > target) {
        hi = c1;
        med = c2;
        lo = target;
    }
    else if (c1 > target) {
        hi = c1;
        med = target;
        lo = c2;
    }
    else {
        hi = target;
        med = c1;
        lo = c2;
    }

    #pragma omp parallel for schedule(static, 1) if(dim >= OMP_THRESHOLD_MIXED)
    for (idx_t bc = 0; bc < n_base; ++bc) {
        idx_t c000 = insertBit0(insertBit0(insertBit0(bc, lo), med), hi);
        idx_t c001 = c000 | incr_target, c010 = c000 | incr_ctrl2, c100 = c000 | incr_ctrl1;
        idx_t c011 = c001 | incr_ctrl2, c101 = c001 | incr_ctrl1, c110 = c010 | incr_ctrl1;
        idx_t c111 = c011 | incr_ctrl1;

        const idx_t offset_c000 = c000 * (2 * dim - c000 + 1) / 2;
        const idx_t offset_c001 = c001 * (2 * dim - c001 + 1) / 2;
        const idx_t offset_c010 = c010 * (2 * dim - c010 + 1) / 2;
        const idx_t offset_c011 = c011 * (2 * dim - c011 + 1) / 2;
        const idx_t offset_c100 = c100 * (2 * dim - c100 + 1) / 2;
        const idx_t offset_c101 = c101 * (2 * dim - c101 + 1) / 2;
        const idx_t offset_c110 = c110 * (2 * dim - c110 + 1) / 2;
        const idx_t offset_c111 = c111 * (2 * dim - c111 + 1) / 2;

        for (idx_t br = bc ; br < n_base; ++br) {
            idx_t r000 = insertBit0(insertBit0(insertBit0(br, lo), med), hi);
            idx_t r001 = r000 | incr_target, r010 = r000 | incr_ctrl2, r100 = r000 | incr_ctrl1;
            idx_t r011 = r001 | incr_ctrl2, r101 = r001 | incr_ctrl1, r110 = r010 | incr_ctrl1;
            idx_t r111 = r011 | incr_ctrl1;

            /* Pair A: (r110,c000) <-> (r111,c000) — always lower tri; plain swap */
            cplx_t tmp = data[(r110 - c000) + offset_c000];
            data[(r110 - c000) + offset_c000] = data[(r111 - c000) + offset_c000];
            data[(r111 - c000) + offset_c000] = tmp;

            /* Pair B: (r110,c010) <-> (r111,c010) — always lower tri; plain swap */
            tmp = data[(r110 - c010) + offset_c010];
            data[(r110 - c010) + offset_c010] = data[(r111 - c010) + offset_c010];
            data[(r111 - c010) + offset_c010] = tmp;

            /* Pair C: (r110,c100) <-> (r111,c100) — always lower tri; plain swap */
            tmp = data[(r110 - c100) + offset_c100];
            data[(r110 - c100) + offset_c100] = data[(r111 - c100) + offset_c100];
            data[(r111 - c100) + offset_c100] = tmp;

            /* Pair D: (r110,c110) <-> (r111,c111) — always lower tri; plain swap */
            tmp = data[(r110 - c110) + offset_c110];
            data[(r110 - c110) + offset_c110] = data[(r111 - c111) + offset_c111];
            data[(r111 - c111) + offset_c111] = tmp;

            /* Fast path: all remaining pairs in lower triangle */
            if (r000 > c111) {
                /* Pair E: (r111,c001) <-> (r110,c001); plain swap */
                tmp = data[(r111 - c001) + offset_c001];
                data[(r111 - c001) + offset_c001] = data[(r110 - c001) + offset_c001];
                data[(r110 - c001) + offset_c001] = tmp;

                /* Pair F: (r111,c011) <-> (r110,c011); plain swap */
                tmp = data[(r111 - c011) + offset_c011];
                data[(r111 - c011) + offset_c011] = data[(r110 - c011) + offset_c011];
                data[(r110 - c011) + offset_c011] = tmp;

                /* Pair G: (r111,c101) <-> (r110,c101); plain swap */
                tmp = data[(r111 - c101) + offset_c101];
                data[(r111 - c101) + offset_c101] = data[(r110 - c101) + offset_c101];
                data[(r110 - c101) + offset_c101] = tmp;

                /* Pair H: (r111,c110) <-> (r110,c111); plain swap */
                tmp = data[(r111 - c110) + offset_c110];
                data[(r111 - c110) + offset_c110] = data[(r110 - c111) + offset_c111];
                data[(r110 - c111) + offset_c111] = tmp;

                /* Pair I: (r000,c110) <-> (r000,c111); plain swap */
                tmp = data[(r000 - c110) + offset_c110];
                data[(r000 - c110) + offset_c110] = data[(r000 - c111) + offset_c111];
                data[(r000 - c111) + offset_c111] = tmp;

                /* Pair J: (r001,c110) <-> (r001,c111); plain swap */
                tmp = data[(r001 - c110) + offset_c110];
                data[(r001 - c110) + offset_c110] = data[(r001 - c111) + offset_c111];
                data[(r001 - c111) + offset_c111] = tmp;

                /* Pair K: (r010,c110) <-> (r010,c111); plain swap */
                tmp = data[(r010 - c110) + offset_c110];
                data[(r010 - c110) + offset_c110] = data[(r010 - c111) + offset_c111];
                data[(r010 - c111) + offset_c111] = tmp;

                /* Pair L: (r011,c110) <-> (r011,c111); plain swap */
                tmp = data[(r011 - c110) + offset_c110];
                data[(r011 - c110) + offset_c110] = data[(r011 - c111) + offset_c111];
                data[(r011 - c111) + offset_c111] = tmp;

                /* Pair M: (r100,c110) <-> (r100,c111); plain swap */
                tmp = data[(r100 - c110) + offset_c110];
                data[(r100 - c110) + offset_c110] = data[(r100 - c111) + offset_c111];
                data[(r100 - c111) + offset_c111] = tmp;

                /* Pair N: (r101,c110) <-> (r101,c111); plain swap */
                tmp = data[(r101 - c110) + offset_c110];
                data[(r101 - c110) + offset_c110] = data[(r101 - c111) + offset_c111];
                data[(r101 - c111) + offset_c111] = tmp;
            }
            else {
                /* Pair E: (r111,c001) lower; (r110,c001) may cross diagonal */
                tmp = data[(r111 - c001) + offset_c001];
                if (r110 > c001) {
                    data[(r111 - c001) + offset_c001] = data[(r110 - c001) + offset_c001];
                    data[(r110 - c001) + offset_c001] = tmp;
                }
                else {
                    data[(r111 - c001) + offset_c001] = conj(data[pack_idx(dim, c001, r110)]);
                    data[pack_idx(dim, c001, r110)] = conj(tmp);
                }

                /* Pair F: (r111,c011) lower; (r110,c011) may cross diagonal */
                tmp = data[(r111 - c011) + offset_c011];
                if (r110 > c011) {
                    data[(r111 - c011) + offset_c011] = data[(r110 - c011) + offset_c011];
                    data[(r110 - c011) + offset_c011] = tmp;
                }
                else {
                    data[(r111 - c011) + offset_c011] = conj(data[pack_idx(dim, c011, r110)]);
                    data[pack_idx(dim, c011, r110)] = conj(tmp);
                }

                /* Pair G: (r111,c101) lower; (r110,c101) may cross diagonal */
                tmp = data[(r111 - c101) + offset_c101];
                if (r110 > c101) {
                    data[(r111 - c101) + offset_c101] = data[(r110 - c101) + offset_c101];
                    data[(r110 - c101) + offset_c101] = tmp;
                }
                else {
                    data[(r111 - c101) + offset_c101] = conj(data[pack_idx(dim, c101, r110)]);
                    data[pack_idx(dim, c101, r110)] = conj(tmp);
                }

                /* Pair H: (r111,c110) lower; (r110,c111) may cross diagonal */
                tmp = data[(r111 - c110) + offset_c110];
                if (r110 > c111) {
                    data[(r111 - c110) + offset_c110] = data[(r110 - c111) + offset_c111];
                    data[(r110 - c111) + offset_c111] = tmp;
                }
                else {
                    data[(r111 - c110) + offset_c110] = conj(data[pack_idx(dim, c111, r110)]);
                    data[pack_idx(dim, c111, r110)] = conj(tmp);
                }

                if (bc != br) {
                    /* Pair I: (r000,c110) <-> (r000,c111); may cross diagonal */
                    if (r000 > c110) {
                        tmp = data[(r000 - c110) + offset_c110];
                        data[(r000 - c110) + offset_c110] = conj(data[pack_idx(dim, c111, r000)]);
                        data[pack_idx(dim, c111, r000)] = conj(tmp);
                    }
                    else {
                        tmp = data[pack_idx(dim, c110, r000)];
                        data[pack_idx(dim, c110, r000)] = data[pack_idx(dim, c111, r000)];
                        data[pack_idx(dim, c111, r000)] = tmp;
                    }

                    /* Pair J: (r001,c110) <-> (r001,c111); may cross diagonal */
                    if (r001 > c111) {
                        tmp = data[(r001 - c110) + offset_c110];
                        data[(r001 - c110) + offset_c110] = data[(r001 - c111) + offset_c111];
                        data[(r001 - c111) + offset_c111] = tmp;
                    }
                    else if (r001 > c110) {
                        tmp = data[(r001 - c110) + offset_c110];
                        data[(r001 - c110) + offset_c110] = conj(data[pack_idx(dim, c111, r001)]);
                        data[pack_idx(dim, c111, r001)] = conj(tmp);
                    }
                    else {
                        tmp = data[pack_idx(dim, c110, r001)];
                        data[pack_idx(dim, c110, r001)] = data[pack_idx(dim, c111, r001)];
                        data[pack_idx(dim, c111, r001)] = tmp;
                    }

                    /* Pair K: (r010,c110) <-> (r010,c111); may cross diagonal */
                    if (r010 > c111) {
                        tmp = data[(r010 - c110) + offset_c110];
                        data[(r010 - c110) + offset_c110] = data[(r010 - c111) + offset_c111];
                        data[(r010 - c111) + offset_c111] = tmp;
                    }
                    else if (r010 > c110) {
                        tmp = data[(r010 - c110) + offset_c110];
                        data[(r010 - c110) + offset_c110] = conj(data[pack_idx(dim, c111, r010)]);
                        data[pack_idx(dim, c111, r010)] = conj(tmp);
                    }
                    else {
                        tmp = data[pack_idx(dim, c110, r010)];
                        data[pack_idx(dim, c110, r010)] = data[pack_idx(dim, c111, r010)];
                        data[pack_idx(dim, c111, r010)] = tmp;
                    }

                    /* Pair L: (r011,c110) <-> (r011,c111); may cross diagonal */
                    if (r011 > c111) {
                        tmp = data[(r011 - c110) + offset_c110];
                        data[(r011 - c110) + offset_c110] = data[(r011 - c111) + offset_c111];
                        data[(r011 - c111) + offset_c111] = tmp;
                    }
                    else if (r011 > c110) {
                        tmp = data[(r011 - c110) + offset_c110];
                        data[(r011 - c110) + offset_c110] = conj(data[pack_idx(dim, c111, r011)]);
                        data[pack_idx(dim, c111, r011)] = conj(tmp);
                    }
                    else {
                        tmp = data[pack_idx(dim, c110, r011)];
                        data[pack_idx(dim, c110, r011)] = data[pack_idx(dim, c111, r011)];
                        data[pack_idx(dim, c111, r011)] = tmp;
                    }

                    /* Pair M: (r100,c110) <-> (r100,c111); may cross diagonal */
                    if (r100 > c111) {
                        tmp = data[(r100 - c110) + offset_c110];
                        data[(r100 - c110) + offset_c110] = data[(r100 - c111) + offset_c111];
                        data[(r100 - c111) + offset_c111] = tmp;
                    }
                    else if (r100 > c110) {
                        tmp = data[(r100 - c110) + offset_c110];
                        data[(r100 - c110) + offset_c110] = conj(data[pack_idx(dim, c111, r100)]);
                        data[pack_idx(dim, c111, r100)] = conj(tmp);
                    }
                    else {
                        tmp = data[pack_idx(dim, c110, r100)];
                        data[pack_idx(dim, c110, r100)] = data[pack_idx(dim, c111, r100)];
                        data[pack_idx(dim, c111, r100)] = tmp;
                    }

                    /* Pair N: (r101,c110) <-> (r101,c111); may cross diagonal */
                    if (r101 > c111) {
                        tmp = data[(r101 - c110) + offset_c110];
                        data[(r101 - c110) + offset_c110] = data[(r101 - c111) + offset_c111];
                        data[(r101 - c111) + offset_c111] = tmp;
                    }
                    else if (r101 > c110) {
                        tmp = data[(r101 - c110) + offset_c110];
                        data[(r101 - c110) + offset_c110] = conj(data[pack_idx(dim, c111, r101)]);
                        data[pack_idx(dim, c111, r101)] = conj(tmp);
                    }
                    else {
                        tmp = data[pack_idx(dim, c110, r101)];
                        data[pack_idx(dim, c110, r101)] = data[pack_idx(dim, c111, r101)];
                        data[pack_idx(dim, c111, r101)] = tmp;
                    }
                }
            }
        }
    }
}
