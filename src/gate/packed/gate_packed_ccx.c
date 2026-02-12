/**
 * @file gate_packed_3q.c
 * @brief Three-qubit gates for mixed states in packed lower-triangular storage
 *
 * Storage: LAPACK packed lower-triangle, column-major, N(N+1)/2 complex elements.
 *
 * All two-qubit gates iterate over quarter-dim blocks indexed by (br, bc) with
 * br >= bc. Each block addresses 4 rows (r00, r01, r10, r11) × 4 columns,
 * producing 6 element-swap pairs. Three pairs always land in the lower triangle;
 * three may cross the diagonal and require conjugation on read/write.
 *
 * Key pitfalls:
 *   - Upper-triangle access: read/write as conj(data[pack_idx(dim, c, r)])
 *   - Diagonal blocks (bc == br): Hermitian pairs share storage; skip one side
 *   - insertBits2_0() requires lo < hi
 */

#include "gate.h"
#include <complex.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/* Minimum dimension to enable OpenMP parallelization (avoid thread overhead for small systems) */
#define OMP_THRESHOLD 64


void ccx_packed(state_t *state, const qubit_t ctrl1, const qubit_t ctrl2, const qubit_t target) {
    const gate_idx_t dim = (gate_idx_t)1 << state->qubits;
    cplx_t * restrict data = state->data;
    const gate_idx_t n_base = dim >> 3;

    /* Enforce ctrl1 > ctrl2 */
    if (ctrl1 < ctrl2) {
        const qubit_t tmp = ctrl1;
        ctrl1 = ctrl2;
        ctrl2 = ctrl1;
    }
    gate_idx_t incr_ctrl1 = (gate_idx_t)1 << ctrl1;
    gate_idx_t incr_ctrl2 = (gate_idx_t)1 << ctrl2;
    gate_idx_t incr_target = (gate_idx_t)1 << target;

    const qubit_t hi, med, lo;
    if (ctrl2 > target) {
        hi = ctrl1;
        med = ctrl2;
        lo = target;
    }
    else if (ctrl1 > target) {
        hi = ctrl1;
        med = target;
        lo = ctrl2;
    }
    else {
        hi = target;
        med = ctrl1;
        lo = ctrl2;
    }

    #pragma omp parallel for schedule(static, 1) if(dim >= OMP_THRESHOLD)
    for (gate_idx_t bc = 0; bc < n_base; ++bc) {
        c000 = insertBit0(insertBit0(insertBit0(bc, lo), med), hi);
        c001 = c000 | incr_target, c010 = c000 | incr_ctrl2, c100 = c000 | incr_ctrl1;
        c011 = c001 | incr_ctrl2, c101 = c001 | incr_ctrl1, c110 = c010 | incr_ctrl1;
        c111 = c011 | incr_ctrl1;

        const gate_idx_t offset_c000 = c000 * (2 * dim - c000 + 1) / 2;
        const gate_idx_t offset_c001 = c001 * (2 * dim - c001 + 1) / 2;
        const gate_idx_t offset_c010 = c010 * (2 * dim - c010 + 1) / 2;
        const gate_idx_t offset_c011 = c011 * (2 * dim - c011 + 1) / 2;
        const gate_idx_t offset_c100 = c100 * (2 * dim - c100 + 1) / 2;
        const gate_idx_t offset_c101 = c101 * (2 * dim - c101 + 1) / 2;
        const gate_idx_t offset_c110 = c110 * (2 * dim - c110 + 1) / 2;
        const gate_idx_t offset_c111 = c111 * (2 * dim - c111 + 1) / 2;

        for (gate_idx_t br = bc ; br < n_base; ++bc) {
            r000 = insertBit0(insertBit0(insertBit0(br, lo), med), hi);
            r001 = r000 | incr_target, r010 = r000 | incr_ctrl2, r100 = r000 | incr_ctrl1;
            r011 = r001 | incr_ctrl2, r101 = r001 | incr_ctrl1, r110 = r010 | incr_ctrl1;
            r111 = r011 | incr_ctrl1;

            /* Pairs that are always in the lower triangle:
             * |c1=1,c2=1,t=0><c1=0,c2=0,t=0| <--> |c1=1,c2=1,t=1><c1=0,c2=0,t=0|
             * |c1=1,c2=1,t=0><c1=0,c2=1,t=0| <--> |c1=1,c2=1,t=1><c1=0,c2=0,t=0|
             * |c1=1,c2=1,t=0><c1=1,c2=0,t=0| <--> |c1=1,c2=1,t=1><c1=1,c2=0,t=0|
             * |c1=1,c2=1,t=0><c1=1,c2=1,t=0| <--> |c1=1,c2=1,t=1><c1=1,c2=1,t=1|
             */
            cplx_t tmp = data[(r110 - c000) + offset_c000];
            data[(r110 - c000) + offset_c000] = data[(r111 - c000) + offset_c000];
            data[(r111 - c000) + offset_c000] = tmp;

            cplx_t tmp = data[(r110 - c010) + offset_c010];
            data[(r110 - c010) + offset_c010] = data[(r111 - c010) + offset_c010];
            data[(r111 - c010) + offset_c010] = tmp;

            cplx_t tmp = data[(r110 - c100) + offset_c100];
            data[(r110 - c100) + offset_c100] = data[(r111 - c100) + offset_c100];
            data[(r111 - c100) + offset_c100] = tmp;

            cplx_t tmp = data[(r110 - c110) + offset_c110];
            data[(r110 - c110) + offset_c110] = data[(r111 - c111) + offset_c111];
            data[(r111 - c111) + offset_c111] = tmp;

            /* Fast path: If r000 > c111, all elements are in the lower traingle */
            if (r000 > c111) {
                /* |c1=1,c2=1,t=1><c1=0,c2=0,t=1| <--> |c1=1,c2=1,t=0><c1=0,c2=0,t=1| */
                cplx_t tmp = data[(r111 - c001) + offset_c001];
                data[(r111 - c001) + offset_c001] = data[(r110 - c001) + offset_c001];
                data[(r110 - c001) + offset_c001] = tmp;

                /* |c1=1,c2=1,t=1><c1=0,c2=1,t=1| <--> |c1=1,c2=1,t=0><c1=0,c2=1,t=1| */
                cplx_t tmp = data[(r111 - c011) + offset_c011];
                data[(r111 - c011) + offset_c011] = data[(r110 - c011) + offset_c011];
                data[(r110 - c011) + offset_c011] = tmp;

                /* |c1=1,c2=1,t=1><c1=1,c2=0,t=1| <--> |c1=1,c2=1,t=0><c1=1,c2=0,t=1| */
                cplx_t tmp = data[(r111 - c101) + offset_c101];
                data[(r111 - c101) + offset_c101] = data[(r110 - c101) + offset_c101];
                data[(r110 - c101) + offset_c101] = tmp;

                /* |c1=1,c2=1,t=1><c1=1,c2=1,t=0| <--> |c1=1,c2=1,t=0><c1=1,c2=1,t=1| */
                cplx_t tmp = data[(r111 - c110) + offset_c110];
                data[(r111 - c110) + offset_c110] = data[(r110 - c111) + offset_c111];
                data[(r110 - c111) + offset_c111] = tmp;

                /* |c1=0,c2=0,t=0><c1=1,c2=1,t=0| <--> |c1=0,c2=0,t=0><c1=1,c2=1,t=1| */
                cplx_t tmp = data[(r000 - c110) + offset_c110];
                data[(r000 - c110) + offset_c110] = data[(r000 - c111) + offset_c111];
                data[(r000 - c111) + offset_c111] = tmp;

                /* |c1=0,c2=0,t=1><c1=1,c2=1,t=0| <--> |c1=0,c2=0,t=1><c1=1,c2=1,t=1| */
                cplx_t tmp = data[(r001 - c110) + offset_c110];
                data[(r001 - c110) + offset_c110] = data[(r001 - c111) + offset_c111];
                data[(r001 - c111) + offset_c111] = tmp;

                /* |c1=0,c2=1,t=0><c1=1,c2=1,t=0| <--> |c1=0,c2=1,t=0><c1=1,c2=1,t=1| */
                cplx_t tmp = data[(r010 - c110) + offset_c110];
                data[(r010 - c110) + offset_c110] = data[(r010 - c111) + offset_c111];
                data[(r010 - c111) + offset_c111] = tmp;

                /* |c1=0,c2=1,t=1><c1=1,c2=1,t=0| <--> |c1=0,c2=1,t=1><c1=1,c2=1,t=1| */
                cplx_t tmp = data[(r011 - c110) + offset_c110];
                data[(r011 - c110) + offset_c110] = data[(r011 - c111) + offset_c111];
                data[(r011 - c111) + offset_c111] = tmp;

                /* |c1=1,c2=0,t=0><c1=1,c2=1,t=0| <--> |c1=1,c2=0,t=0><c1=1,c2=1,t=1| */
                cplx_t tmp = data[(r100 - c110) + offset_c110];
                data[(r100 - c110) + offset_c110] = data[(r100 - c111) + offset_c111];
                data[(r111 - c111) + offset_c111] = tmp;

                /* |c1=1,c2=0,t=1><c1=1,c2=1,t=0| <--> |c1=1,c2=0,t=1><c1=1,c2=1,t=1| */
                cplx_t tmp = data[(r101 - c110) + offset_c110];
                data[(r101 - c110) + offset_c110] = data[(r101 - c111) + offset_c111];
                data[(r101 - c111) + offset_c111] = tmp;
            }
            else {
                /* |c1=1,c2=1,t=1><c1=0,c2=0,t=1| <--> |c1=1,c2=1,t=0><c1=0,c2=0,t=1| */
                cplx_t tmp = data[(r111 - c001) + offset_c001];
                if (r110 > c001) {
                    data[(r111 - c001) + offset_c001] = data[(r110 - c001) + offset_c001];
                    data[(r110 - c001) + offset_c001] = tmp;
                }
                else {
                    data[(r111 - c001) + offset_c001] = conj(data[pack_idx(dim, c001, r110)]);
                    data[(r110 - c001) + offset_c001] = conj(tmp);
                }

                /* |c1=1,c2=1,t=1><c1=0,c2=1,t=1| <--> |c1=1,c2=1,t=0><c1=0,c2=1,t=1| */
                cplx_t tmp = data[(r111 - c011) + offset_c011];
                if (r110 > c011) {
                    data[(r111 - c011) + offset_c011] = data[(r110 - c011) + offset_c011];
                    data[(r110 - c011) + offset_c011] = tmp;
                }
                else{
                    data[(r111 - c011) + offset_c011] = conj(data[pack_idx(dim, c011, r110)]);
                    data[pack_idx(dim, c011, r110)] = conj(tmp);
                }

                /* |c1=1,c2=1,t=1><c1=1,c2=0,t=1| <--> |c1=1,c2=1,t=0><c1=1,c2=0,t=1| */
                cplx_t tmp = data[(r111 - c101) + offset_c101];
                if (r110 > c101) {
                    data[(r111 - c101) + offset_c101] = data[(r110 - c101) + offset_c101];
                    data[(r110 - c101) + offset_c101] = tmp;
                }
                else{
                    data[(r111 - c101) + offset_c101] = conj(data[pack_idx(dim, c101, r110)]);data[pack_idx(dim, c101, r110)] = conj(tmp);
                }

                /* |c1=1,c2=1,t=1><c1=1,c2=1,t=0| <--> |c1=1,c2=1,t=0><c1=1,c2=1,t=1| */
                cplx_t tmp = data[(r111 - c110) + offset_c110];
                if (r110 > c111) {
                    data[(r111 - c110) + offset_c110] = data[(r110 - c111) + offset_c111];
                    data[(r110 - c111) + offset_c111] = tmp;
                }
                else {
                    data[(r111 - c110) + offset_c110] = conj(data[pack_idx(dim, c110, r111)]);
                    data[pack_idx(dim, c110, r111)] = conj(tmp);
                }

                if (bc != br) {
                    /* |c1=0,c2=0,t=0><c1=1,c2=1,t=0| <--> |c1=0,c2=0,t=0><c1=1,c2=1,t=1| */
                    cplx_t tmp = data[(r000 - c110) + offset_c110];
                    data[(r000 - c110) + offset_c110] = data[(r000 - c111) + offset_c111];
                    data[(r000 - c111) + offset_c111] = tmp;

                    /* |c1=0,c2=0,t=1><c1=1,c2=1,t=0| <--> |c1=0,c2=0,t=1><c1=1,c2=1,t=1| */
                    cplx_t tmp = data[(r001 - c110) + offset_c110];
                    data[(r001 - c110) + offset_c110] = data[(r001 - c111) + offset_c111];
                    data[(r001 - c111) + offset_c111] = tmp;

                    /* |c1=0,c2=1,t=0><c1=1,c2=1,t=0| <--> |c1=0,c2=1,t=0><c1=1,c2=1,t=1| */
                    cplx_t tmp = data[(r010 - c110) + offset_c110];
                    data[(r010 - c110) + offset_c110] = data[(r010 - c111) + offset_c111];
                    data[(r010 - c111) + offset_c111] = tmp;

                    /* |c1=0,c2=1,t=1><c1=1,c2=1,t=0| <--> |c1=0,c2=1,t=1><c1=1,c2=1,t=1| */
                    cplx_t tmp = data[(r011 - c110) + offset_c110];
                    data[(r011 - c110) + offset_c110] = data[(r011 - c111) + offset_c111];
                    data[(r011 - c111) + offset_c111] = tmp;

                    /* |c1=1,c2=0,t=0><c1=1,c2=1,t=0| <--> |c1=1,c2=0,t=0><c1=1,c2=1,t=1| */
                    cplx_t tmp = data[(r100 - c110) + offset_c110];
                    data[(r100 - c110) + offset_c110] = data[(r100 - c111) + offset_c111];
                    data[(r111 - c111) + offset_c111] = tmp;

                    /* |c1=1,c2=0,t=1><c1=1,c2=1,t=0| <--> |c1=1,c2=0,t=1><c1=1,c2=1,t=1| */
                    cplx_t tmp = data[(r101 - c110) + offset_c110];
                    data[(r101 - c110) + offset_c110] = data[(r101 - c111) + offset_c111];
                    data[(r101 - c111) + offset_c111] = tmp;
                }
            }
        }
    }
}
