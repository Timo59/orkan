/**
 * @file gate_tiled_ccx.c
 * @brief CCX (Toffoli) gate for tiled density matrices
 */

#include "gate_tiled.h"

void ccx_tiled(state_t *state, const qubit_t control1, const qubit_t control2, const qubit_t target) {
    const qubit_t ctrl1 = (control1 > control2) ? control1 : control2;
    const qubit_t ctrl2 = (control1 < control2) ? control1 : control2;
    const idx_t dim = (idx_t)1 << state->qubits;
    const idx_t dim_tile = (dim + TILE_DIM - 1) >> LOG_TILE_DIM;
    cplx_t * restrict data = state->data;

    /* Sort all three qubit positions: lo < mid < hi */
    qubit_t lo, mid, hi;
    if (ctrl2 > target) {
        lo = target; mid = ctrl2; hi = ctrl1;
    }
    else if (ctrl1 > target) {
        lo = ctrl2; mid = target; hi = ctrl1;
    }
    else {
        lo = ctrl2; mid = ctrl1; hi = target;
    }

    if (ctrl1 < LOG_TILE_DIM) {
        /*
         * Case 1:
         * All three qubit indices are smaller than LOG_TILE_DIM.
         *
         * All swap pairs are contained in the same tile.
         */
        if(target < LOG_TILE_DIM) {
            const idx_t n_tiles = dim_tile * (dim_tile + 1) / 2;
            const idx_t n_base = (dim < TILE_DIM) ? dim >> 3 : TILE_DIM >> 3;
            const idx_t incr_ctrl1 = (idx_t)1 << ctrl1;
            const idx_t incr_ctrl2 = (idx_t)1 << ctrl2;
            const idx_t incr_tgt   = (idx_t)1 << target;

            #pragma omp parallel for schedule(static) if(dim >= OMP_THRESHOLD)
            for (idx_t i = 0; i < n_tiles; ++i) {
                cplx_t * restrict tile = data + i * TILE_SIZE;

                for (idx_t blr = 0; blr < n_base; ++blr) {
                    /* insertBits3_0: insert 0-bits at positions lo, mid, hi */
                    const idx_t r000 = insertBit0(insertBits2_0(blr, lo, mid), hi);
                    const idx_t r_c1   = r000 | incr_ctrl1;                   /* ctrl1=1, ctrl2=0, tgt=0 */
                    const idx_t r_c2   = r000 | incr_ctrl2;                   /* ctrl1=0, ctrl2=1, tgt=0 */
                    const idx_t r_t    = r000 | incr_tgt;                     /* ctrl1=0, ctrl2=0, tgt=1 */
                    const idx_t r_c1c2 = r000 | incr_ctrl1 | incr_ctrl2;     /* ctrl1=1, ctrl2=1, tgt=0 */
                    const idx_t r_c1t  = r000 | incr_ctrl1 | incr_tgt;       /* ctrl1=1, ctrl2=0, tgt=1 */
                    const idx_t r_c2t  = r000 | incr_ctrl2 | incr_tgt;       /* ctrl1=0, ctrl2=1, tgt=1 */
                    const idx_t r_all  = r000 | incr_ctrl1 | incr_ctrl2 | incr_tgt; /* ctrl1=1, ctrl2=1, tgt=1 */

                    for (idx_t blc = 0; blc < n_base; ++blc) {
                        const idx_t c000 = insertBit0(insertBits2_0(blc, lo, mid), hi);
                        const idx_t c_c1   = c000 | incr_ctrl1;
                        const idx_t c_c2   = c000 | incr_ctrl2;
                        const idx_t c_t    = c000 | incr_tgt;
                        const idx_t c_c1c2 = c000 | incr_ctrl1 | incr_ctrl2;
                        const idx_t c_c1t  = c000 | incr_ctrl1 | incr_tgt;
                        const idx_t c_c2t  = c000 | incr_ctrl2 | incr_tgt;
                        const idx_t c_all  = c000 | incr_ctrl1 | incr_ctrl2 | incr_tgt;

                        cplx_t tmp;

                        /*
                         * 14 swap pairs from ρ'[r,c] = ρ[perm(r), perm(c)]:
                         *   perm flips target bit when both ctrl1=1 and ctrl2=1.
                         *   Row swaps: r_c1c2 <-> r_all (ctrl1=1,ctrl2=1,tgt=0 <-> ctrl1=1,ctrl2=1,tgt=1)
                         *   Column swaps: c_c1c2 <-> c_all
                         */

                        /* ===== Pairs A-C, E-G: Row changes (r_c1c2 <-> r_all), column fixed ===== */

                        /* Pair A: |c1=1,c2=1,t=0><c1=0,c2=0,t=0| <--> |c1=1,c2=1,t=1><c1=0,c2=0,t=0| */
                        tmp = tile[elem_off(r_c1c2, c000)];
                        tile[elem_off(r_c1c2, c000)] = tile[elem_off(r_all, c000)];
                        tile[elem_off(r_all, c000)] = tmp;

                        /* Pair B: |c1=1,c2=1,t=0><c1=0,c2=1,t=0| <--> |c1=1,c2=1,t=1><c1=0,c2=1,t=0| */
                        tmp = tile[elem_off(r_c1c2, c_c2)];
                        tile[elem_off(r_c1c2, c_c2)] = tile[elem_off(r_all, c_c2)];
                        tile[elem_off(r_all, c_c2)] = tmp;

                        /* Pair C: |c1=1,c2=1,t=0><c1=1,c2=0,t=0| <--> |c1=1,c2=1,t=1><c1=1,c2=0,t=0| */
                        tmp = tile[elem_off(r_c1c2, c_c1)];
                        tile[elem_off(r_c1c2, c_c1)] = tile[elem_off(r_all, c_c1)];
                        tile[elem_off(r_all, c_c1)] = tmp;

                        /* Pair E: |c1=1,c2=1,t=1><c1=0,c2=0,t=1| <--> |c1=1,c2=1,t=0><c1=0,c2=0,t=1| */
                        tmp = tile[elem_off(r_all, c_t)];
                        tile[elem_off(r_all, c_t)] = tile[elem_off(r_c1c2, c_t)];
                        tile[elem_off(r_c1c2, c_t)] = tmp;

                        /* Pair F: |c1=1,c2=1,t=1><c1=0,c2=1,t=1| <--> |c1=1,c2=1,t=0><c1=0,c2=1,t=1| */
                        tmp = tile[elem_off(r_all, c_c2t)];
                        tile[elem_off(r_all, c_c2t)] = tile[elem_off(r_c1c2, c_c2t)];
                        tile[elem_off(r_c1c2, c_c2t)] = tmp;

                        /* Pair G: |c1=1,c2=1,t=1><c1=1,c2=0,t=1| <--> |c1=1,c2=1,t=0><c1=1,c2=0,t=1| */
                        tmp = tile[elem_off(r_all, c_c1t)];
                        tile[elem_off(r_all, c_c1t)] = tile[elem_off(r_c1c2, c_c1t)];
                        tile[elem_off(r_c1c2, c_c1t)] = tmp;

                        /* ===== Pairs D, H: Both row and column change ===== */

                        /* Pair D: |c1=1,c2=1,t=0><c1=1,c2=1,t=0| <--> |c1=1,c2=1,t=1><c1=1,c2=1,t=1| */
                        tmp = tile[elem_off(r_c1c2, c_c1c2)];
                        tile[elem_off(r_c1c2, c_c1c2)] = tile[elem_off(r_all, c_all)];
                        tile[elem_off(r_all, c_all)] = tmp;

                        /* Pair H: |c1=1,c2=1,t=1><c1=1,c2=1,t=0| <--> |c1=1,c2=1,t=0><c1=1,c2=1,t=1| */
                        tmp = tile[elem_off(r_all, c_c1c2)];
                        tile[elem_off(r_all, c_c1c2)] = tile[elem_off(r_c1c2, c_all)];
                        tile[elem_off(r_c1c2, c_all)] = tmp;

                        /* ===== Pairs I-N: Column changes (c_c1c2 <-> c_all), row fixed ===== */

                        /* Pair I: |c1=0,c2=0,t=0><c1=1,c2=1,t=0| <--> |c1=0,c2=0,t=0><c1=1,c2=1,t=1| */
                        tmp = tile[elem_off(r000, c_c1c2)];
                        tile[elem_off(r000, c_c1c2)] = tile[elem_off(r000, c_all)];
                        tile[elem_off(r000, c_all)] = tmp;

                        /* Pair J: |c1=0,c2=0,t=1><c1=1,c2=1,t=0| <--> |c1=0,c2=0,t=1><c1=1,c2=1,t=1| */
                        tmp = tile[elem_off(r_t, c_c1c2)];
                        tile[elem_off(r_t, c_c1c2)] = tile[elem_off(r_t, c_all)];
                        tile[elem_off(r_t, c_all)] = tmp;

                        /* Pair K: |c1=0,c2=1,t=0><c1=1,c2=1,t=0| <--> |c1=0,c2=1,t=0><c1=1,c2=1,t=1| */
                        tmp = tile[elem_off(r_c2, c_c1c2)];
                        tile[elem_off(r_c2, c_c1c2)] = tile[elem_off(r_c2, c_all)];
                        tile[elem_off(r_c2, c_all)] = tmp;

                        /* Pair L: |c1=0,c2=1,t=1><c1=1,c2=1,t=0| <--> |c1=0,c2=1,t=1><c1=1,c2=1,t=1| */
                        tmp = tile[elem_off(r_c2t, c_c1c2)];
                        tile[elem_off(r_c2t, c_c1c2)] = tile[elem_off(r_c2t, c_all)];
                        tile[elem_off(r_c2t, c_all)] = tmp;

                        /* Pair M: |c1=1,c2=0,t=0><c1=1,c2=1,t=0| <--> |c1=1,c2=0,t=0><c1=1,c2=1,t=1| */
                        tmp = tile[elem_off(r_c1, c_c1c2)];
                        tile[elem_off(r_c1, c_c1c2)] = tile[elem_off(r_c1, c_all)];
                        tile[elem_off(r_c1, c_all)] = tmp;

                        /* Pair N: |c1=1,c2=0,t=1><c1=1,c2=1,t=0| <--> |c1=1,c2=0,t=1><c1=1,c2=1,t=1| */
                        tmp = tile[elem_off(r_c1t, c_c1c2)];
                        tile[elem_off(r_c1t, c_c1c2)] = tile[elem_off(r_c1t, c_all)];
                        tile[elem_off(r_c1t, c_all)] = tmp;
                    }
                }
            }
        }

        /*
         * Case 2:
         * Control indices are smaller than LOG_TILE_DIM & target index is larger than LOG_TILE_DIM.
         *
         * The target qubit selects tile rows/columns (tr0=tgt0, tr1=tgt1).
         * Both controls are local bits within each tile.
         *
         * CCX flips the target bit only when ctrl1=ctrl2=1.  Since the target
         * is a tile-level bit, elements whose row/col have ctrl1=ctrl2=1 move
         * between tiles tr0<->tr1 or tc0<->tc1.
         *
         * The 14 swap pairs (from Case 1) split into 5 groups by tile pair:
         *
         *   Group 1 (A,B,C):  t(tr0,tc0) <-> t(tr1,tc0)  -- row flips, col in tc0
         *   Group 2 (E,F,G):  t(tr1,tc1) <-> t(tr0,tc1)  -- row flips, col in tc1
         *   Group 3 (D,H):    t(tr0,tc0) <-> t(tr1,tc1)  and  t(tr1,tc0) <-> t(tr0,tc1)
         *   Group 4 (I,K,M):  t(tr0,tc0) <-> t(tr0,tc1)  -- col flips, row in tr0
         *   Group 5 (J,L,N):  t(tr1,tc0) <-> t(tr1,tc1)  -- col flips, row in tr1
         *
         * Groups 1, 5, and pair D are always in the lower triangle.
         * Groups 2, 4, and pair H involve tile t(tr0,tc1) which may be upper triangle.
         */
        else {
            const idx_t n_bt = dim_tile >> 1;
            const qubit_t tile_bit = target - LOG_TILE_DIM;
            const idx_t incr_tile = (idx_t)1 << tile_bit;
            const idx_t n_bl = TILE_DIM >> 2;
            const qubit_t lo_local = ctrl2;  /* ctrl2 < ctrl1, both < LOG_TILE_DIM */
            const qubit_t hi_local = ctrl1;
            const idx_t incr_ctrl1 = (idx_t)1 << ctrl1;
            const idx_t incr_ctrl2 = (idx_t)1 << ctrl2;

            #pragma omp parallel for schedule(static, 1) if(dim >= OMP_THRESHOLD)
            for (idx_t btr = 0; btr < n_bt; ++btr) {
                const idx_t tr0 = insertBit0(btr, tile_bit);
                const idx_t tr1 = tr0 | incr_tile;

                for (idx_t btc = 0; btc <= btr; ++btc) {
                    const idx_t tc0 = insertBit0(btc, tile_bit);
                    const idx_t tc1 = tc0 | incr_tile;

                    cplx_t * restrict t00 = data + tile_off(tr0, tc0);
                    cplx_t * restrict t10 = data + tile_off(tr1, tc0);
                    cplx_t * restrict t11 = data + tile_off(tr1, tc1);

                    /* Fast path: tr0 > tc1 means all four tiles strictly in lower triangle */
                    if (tr0 > tc1) {
                        cplx_t * restrict t01 = data + tile_off(tr0, tc1);

                        for (idx_t blr = 0; blr < n_bl; ++blr) {
                            const idx_t r000   = insertBits2_0(blr, lo_local, hi_local);
                            const idx_t r_c1   = r000 | incr_ctrl1;
                            const idx_t r_c2   = r000 | incr_ctrl2;
                            const idx_t r_c1c2 = r000 | incr_ctrl1 | incr_ctrl2;

                            for (idx_t blc = 0; blc < n_bl; ++blc) {
                                const idx_t c000   = insertBits2_0(blc, lo_local, hi_local);
                                const idx_t c_c1   = c000 | incr_ctrl1;
                                const idx_t c_c2   = c000 | incr_ctrl2;
                                const idx_t c_c1c2 = c000 | incr_ctrl1 | incr_ctrl2;

                                cplx_t tmp;

                                /* ===== Group 1: t00 <-> t10, row (c1=1,c2=1) flips tile, col fixed in tc0 ===== */

                                /* Pair A: row c1c2 tgt=0 <-> tgt=1, col (0,0) */
                                tmp = t00[elem_off(r_c1c2, c000)];
                                t00[elem_off(r_c1c2, c000)] = t10[elem_off(r_c1c2, c000)];
                                t10[elem_off(r_c1c2, c000)] = tmp;

                                /* Pair B: row c1c2 tgt=0 <-> tgt=1, col (0,1) */
                                tmp = t00[elem_off(r_c1c2, c_c2)];
                                t00[elem_off(r_c1c2, c_c2)] = t10[elem_off(r_c1c2, c_c2)];
                                t10[elem_off(r_c1c2, c_c2)] = tmp;

                                /* Pair C: row c1c2 tgt=0 <-> tgt=1, col (1,0) */
                                tmp = t00[elem_off(r_c1c2, c_c1)];
                                t00[elem_off(r_c1c2, c_c1)] = t10[elem_off(r_c1c2, c_c1)];
                                t10[elem_off(r_c1c2, c_c1)] = tmp;

                                /* ===== Group 5: t10 <-> t11, col (c1=1,c2=1) flips tile, row fixed in tr1 ===== */

                                /* Pair J: row (0,0) in tr1, col c1c2 tc0 <-> tc1 */
                                tmp = t10[elem_off(r000, c_c1c2)];
                                t10[elem_off(r000, c_c1c2)] = t11[elem_off(r000, c_c1c2)];
                                t11[elem_off(r000, c_c1c2)] = tmp;

                                /* Pair L: row (0,1) in tr1, col c1c2 tc0 <-> tc1 */
                                tmp = t10[elem_off(r_c2, c_c1c2)];
                                t10[elem_off(r_c2, c_c1c2)] = t11[elem_off(r_c2, c_c1c2)];
                                t11[elem_off(r_c2, c_c1c2)] = tmp;

                                /* Pair N: row (1,0) in tr1, col c1c2 tc0 <-> tc1 */
                                tmp = t10[elem_off(r_c1, c_c1c2)];
                                t10[elem_off(r_c1, c_c1c2)] = t11[elem_off(r_c1, c_c1c2)];
                                t11[elem_off(r_c1, c_c1c2)] = tmp;

                                /* ===== Group 3, Pair D: both flip, t00 <-> t11 ===== */

                                tmp = t00[elem_off(r_c1c2, c_c1c2)];
                                t00[elem_off(r_c1c2, c_c1c2)] = t11[elem_off(r_c1c2, c_c1c2)];
                                t11[elem_off(r_c1c2, c_c1c2)] = tmp;

                                /* ===== Group 2: t11 <-> t01, row (c1=1,c2=1) flips tile, col fixed in tc1 ===== */

                                /* Pair E: row c1c2 tgt=1 <-> tgt=0, col (0,0) in tc1 */
                                tmp = t11[elem_off(r_c1c2, c000)];
                                t11[elem_off(r_c1c2, c000)] = t01[elem_off(r_c1c2, c000)];
                                t01[elem_off(r_c1c2, c000)] = tmp;

                                /* Pair F: row c1c2 tgt=1 <-> tgt=0, col (0,1) in tc1 */
                                tmp = t11[elem_off(r_c1c2, c_c2)];
                                t11[elem_off(r_c1c2, c_c2)] = t01[elem_off(r_c1c2, c_c2)];
                                t01[elem_off(r_c1c2, c_c2)] = tmp;

                                /* Pair G: row c1c2 tgt=1 <-> tgt=0, col (1,0) in tc1 */
                                tmp = t11[elem_off(r_c1c2, c_c1)];
                                t11[elem_off(r_c1c2, c_c1)] = t01[elem_off(r_c1c2, c_c1)];
                                t01[elem_off(r_c1c2, c_c1)] = tmp;

                                /* ===== Group 3, Pair H: both flip, t10 <-> t01 ===== */

                                tmp = t10[elem_off(r_c1c2, c_c1c2)];
                                t10[elem_off(r_c1c2, c_c1c2)] = t01[elem_off(r_c1c2, c_c1c2)];
                                t01[elem_off(r_c1c2, c_c1c2)] = tmp;

                                /* ===== Group 4: t00 <-> t01, col (c1=1,c2=1) flips tile, row fixed in tr0 ===== */

                                /* Pair I: row (0,0) in tr0, col c1c2 tc0 <-> tc1 */
                                tmp = t00[elem_off(r000, c_c1c2)];
                                t00[elem_off(r000, c_c1c2)] = t01[elem_off(r000, c_c1c2)];
                                t01[elem_off(r000, c_c1c2)] = tmp;

                                /* Pair K: row (0,1) in tr0, col c1c2 tc0 <-> tc1 */
                                tmp = t00[elem_off(r_c2, c_c1c2)];
                                t00[elem_off(r_c2, c_c1c2)] = t01[elem_off(r_c2, c_c1c2)];
                                t01[elem_off(r_c2, c_c1c2)] = tmp;

                                /* Pair M: row (1,0) in tr0, col c1c2 tc0 <-> tc1 */
                                tmp = t00[elem_off(r_c1, c_c1c2)];
                                t00[elem_off(r_c1, c_c1c2)] = t01[elem_off(r_c1, c_c1c2)];
                                t01[elem_off(r_c1, c_c1c2)] = tmp;
                            }
                        }
                    }

                    /* t01 = t(tr0,tc1) is upper triangle; use its adjoint t(tc1,tr0) */
                    else if (btr != btc) {
                        cplx_t * restrict t01 = data + tile_off(tc1, tr0);

                        for (idx_t blr = 0; blr < n_bl; ++blr) {
                            const idx_t r000   = insertBits2_0(blr, lo_local, hi_local);
                            const idx_t r_c1   = r000 | incr_ctrl1;
                            const idx_t r_c2   = r000 | incr_ctrl2;
                            const idx_t r_c1c2 = r000 | incr_ctrl1 | incr_ctrl2;

                            for (idx_t blc = 0; blc < n_bl; ++blc) {
                                const idx_t c000   = insertBits2_0(blc, lo_local, hi_local);
                                const idx_t c_c1   = c000 | incr_ctrl1;
                                const idx_t c_c2   = c000 | incr_ctrl2;
                                const idx_t c_c1c2 = c000 | incr_ctrl1 | incr_ctrl2;

                                cplx_t tmp;

                                /* ===== Group 1: t00 <-> t10, always lower tri ===== */

                                /* Pair A */
                                tmp = t00[elem_off(r_c1c2, c000)];
                                t00[elem_off(r_c1c2, c000)] = t10[elem_off(r_c1c2, c000)];
                                t10[elem_off(r_c1c2, c000)] = tmp;

                                /* Pair B */
                                tmp = t00[elem_off(r_c1c2, c_c2)];
                                t00[elem_off(r_c1c2, c_c2)] = t10[elem_off(r_c1c2, c_c2)];
                                t10[elem_off(r_c1c2, c_c2)] = tmp;

                                /* Pair C */
                                tmp = t00[elem_off(r_c1c2, c_c1)];
                                t00[elem_off(r_c1c2, c_c1)] = t10[elem_off(r_c1c2, c_c1)];
                                t10[elem_off(r_c1c2, c_c1)] = tmp;

                                /* ===== Group 5: t10 <-> t11, always lower tri ===== */

                                /* Pair J */
                                tmp = t10[elem_off(r000, c_c1c2)];
                                t10[elem_off(r000, c_c1c2)] = t11[elem_off(r000, c_c1c2)];
                                t11[elem_off(r000, c_c1c2)] = tmp;

                                /* Pair L */
                                tmp = t10[elem_off(r_c2, c_c1c2)];
                                t10[elem_off(r_c2, c_c1c2)] = t11[elem_off(r_c2, c_c1c2)];
                                t11[elem_off(r_c2, c_c1c2)] = tmp;

                                /* Pair N */
                                tmp = t10[elem_off(r_c1, c_c1c2)];
                                t10[elem_off(r_c1, c_c1c2)] = t11[elem_off(r_c1, c_c1c2)];
                                t11[elem_off(r_c1, c_c1c2)] = tmp;

                                /* ===== Group 3, Pair D: t00 <-> t11, always lower tri ===== */

                                tmp = t00[elem_off(r_c1c2, c_c1c2)];
                                t00[elem_off(r_c1c2, c_c1c2)] = t11[elem_off(r_c1c2, c_c1c2)];
                                t11[elem_off(r_c1c2, c_c1c2)] = tmp;

                                /* ===== Group 2: t11[..] <-> conj(t01[col,row]), t01 is adjoint ===== */

                                /* Pair E: row c1c2, col (0,0) */
                                tmp = t11[elem_off(r_c1c2, c000)];
                                t11[elem_off(r_c1c2, c000)] = conj(t01[elem_off(c000, r_c1c2)]);
                                t01[elem_off(c000, r_c1c2)] = conj(tmp);

                                /* Pair F: row c1c2, col (0,1) */
                                tmp = t11[elem_off(r_c1c2, c_c2)];
                                t11[elem_off(r_c1c2, c_c2)] = conj(t01[elem_off(c_c2, r_c1c2)]);
                                t01[elem_off(c_c2, r_c1c2)] = conj(tmp);

                                /* Pair G: row c1c2, col (1,0) */
                                tmp = t11[elem_off(r_c1c2, c_c1)];
                                t11[elem_off(r_c1c2, c_c1)] = conj(t01[elem_off(c_c1, r_c1c2)]);
                                t01[elem_off(c_c1, r_c1c2)] = conj(tmp);

                                /* ===== Group 3, Pair H: t10[..] <-> conj(t01[col,row]) ===== */

                                tmp = t10[elem_off(r_c1c2, c_c1c2)];
                                t10[elem_off(r_c1c2, c_c1c2)] = conj(t01[elem_off(c_c1c2, r_c1c2)]);
                                t01[elem_off(c_c1c2, r_c1c2)] = conj(tmp);

                                /* ===== Group 4: t00[..] <-> conj(t01[col,row]) ===== */

                                /* Pair I: row (0,0), col c1c2 */
                                tmp = t00[elem_off(r000, c_c1c2)];
                                t00[elem_off(r000, c_c1c2)] = conj(t01[elem_off(c_c1c2, r000)]);
                                t01[elem_off(c_c1c2, r000)] = conj(tmp);

                                /* Pair K: row (0,1), col c1c2 */
                                tmp = t00[elem_off(r_c2, c_c1c2)];
                                t00[elem_off(r_c2, c_c1c2)] = conj(t01[elem_off(c_c1c2, r_c2)]);
                                t01[elem_off(c_c1c2, r_c2)] = conj(tmp);

                                /* Pair M: row (1,0), col c1c2 */
                                tmp = t00[elem_off(r_c1, c_c1c2)];
                                t00[elem_off(r_c1, c_c1c2)] = conj(t01[elem_off(c_c1c2, r_c1)]);
                                t01[elem_off(c_c1c2, r_c1)] = conj(tmp);
                            }
                        }
                    }

                    /* Diagonal block: btr == btc, so t(tr0,tc1) = t(tr0,tr1) is the adjoint of t(tr1,tr0) = t10 */
                    else {
                        /*
                         * Phase 1: one-way writes for pairs involving t01.
                         * Since t01 is the adjoint of t10, reading t01[lr,lc] = conj(t10[lc,lr]).
                         * But t10 is about to be modified by Groups 1/5/H, so we must read from
                         * t10 first (as the adjoint source) before t10 is overwritten.
                         *
                         * Group 2 (E,F,G): t11[r_c1c2, c_xx] = conj(t10[c_xx, r_c1c2])
                         * Group 4 (I,K,M): t00[r_xx, c_c1c2] = conj(t10[c_c1c2, r_xx])
                         */
                        for (idx_t blr = 0; blr < n_bl; ++blr) {
                            const idx_t r000   = insertBits2_0(blr, lo_local, hi_local);
                            const idx_t r_c1   = r000 | incr_ctrl1;
                            const idx_t r_c2   = r000 | incr_ctrl2;
                            const idx_t r_c1c2 = r000 | incr_ctrl1 | incr_ctrl2;

                            for (idx_t blc = 0; blc < n_bl; ++blc) {
                                const idx_t c000   = insertBits2_0(blc, lo_local, hi_local);
                                const idx_t c_c1   = c000 | incr_ctrl1;
                                const idx_t c_c2   = c000 | incr_ctrl2;
                                const idx_t c_c1c2 = c000 | incr_ctrl1 | incr_ctrl2;

                                /* Group 2 (E,F,G): t11 receives from adjoint of t10 */
                                /* Pair E */
                                t11[elem_off(r_c1c2, c000)] = conj(t10[elem_off(c000, r_c1c2)]);
                                /* Pair F */
                                t11[elem_off(r_c1c2, c_c2)] = conj(t10[elem_off(c_c2, r_c1c2)]);
                                /* Pair G */
                                t11[elem_off(r_c1c2, c_c1)] = conj(t10[elem_off(c_c1, r_c1c2)]);

                                /* Group 4 (I,K,M): t00 receives from adjoint of t10 */
                                /* Pair I */
                                t00[elem_off(r000, c_c1c2)] = conj(t10[elem_off(c_c1c2, r000)]);
                                /* Pair K */
                                t00[elem_off(r_c2, c_c1c2)] = conj(t10[elem_off(c_c1c2, r_c2)]);
                                /* Pair M */
                                t00[elem_off(r_c1, c_c1c2)] = conj(t10[elem_off(c_c1c2, r_c1)]);
                            }
                        }

                        /*
                         * Phase 2: two-way swaps on tiles that are not aliased with t01.
                         * Groups 1 (A,B,C), 5 (J,L,N), pair D always in lower tri.
                         * Pair H: t10 <-> t01 = conj(t10^T), self-adjoint swap within t10.
                         */
                        for (idx_t blr = 0; blr < n_bl; ++blr) {
                            const idx_t r000   = insertBits2_0(blr, lo_local, hi_local);
                            const idx_t r_c1   = r000 | incr_ctrl1;
                            const idx_t r_c2   = r000 | incr_ctrl2;
                            const idx_t r_c1c2 = r000 | incr_ctrl1 | incr_ctrl2;

                            for (idx_t blc = 0; blc < n_bl; ++blc) {
                                const idx_t c000   = insertBits2_0(blc, lo_local, hi_local);
                                const idx_t c_c1   = c000 | incr_ctrl1;
                                const idx_t c_c2   = c000 | incr_ctrl2;
                                const idx_t c_c1c2 = c000 | incr_ctrl1 | incr_ctrl2;

                                cplx_t tmp;

                                /* ===== Group 1: t00 <-> t10 ===== */

                                /* Pair A */
                                tmp = t00[elem_off(r_c1c2, c000)];
                                t00[elem_off(r_c1c2, c000)] = t10[elem_off(r_c1c2, c000)];
                                t10[elem_off(r_c1c2, c000)] = tmp;

                                /* Pair B */
                                tmp = t00[elem_off(r_c1c2, c_c2)];
                                t00[elem_off(r_c1c2, c_c2)] = t10[elem_off(r_c1c2, c_c2)];
                                t10[elem_off(r_c1c2, c_c2)] = tmp;

                                /* Pair C */
                                tmp = t00[elem_off(r_c1c2, c_c1)];
                                t00[elem_off(r_c1c2, c_c1)] = t10[elem_off(r_c1c2, c_c1)];
                                t10[elem_off(r_c1c2, c_c1)] = tmp;

                                /* ===== Group 5: t10 <-> t11 ===== */

                                /* Pair J */
                                tmp = t10[elem_off(r000, c_c1c2)];
                                t10[elem_off(r000, c_c1c2)] = t11[elem_off(r000, c_c1c2)];
                                t11[elem_off(r000, c_c1c2)] = tmp;

                                /* Pair L */
                                tmp = t10[elem_off(r_c2, c_c1c2)];
                                t10[elem_off(r_c2, c_c1c2)] = t11[elem_off(r_c2, c_c1c2)];
                                t11[elem_off(r_c2, c_c1c2)] = tmp;

                                /* Pair N */
                                tmp = t10[elem_off(r_c1, c_c1c2)];
                                t10[elem_off(r_c1, c_c1c2)] = t11[elem_off(r_c1, c_c1c2)];
                                t11[elem_off(r_c1, c_c1c2)] = tmp;

                                /* ===== Group 3, Pair D: t00 <-> t11 ===== */

                                tmp = t00[elem_off(r_c1c2, c_c1c2)];
                                t00[elem_off(r_c1c2, c_c1c2)] = t11[elem_off(r_c1c2, c_c1c2)];
                                t11[elem_off(r_c1c2, c_c1c2)] = tmp;

                                /* ===== Group 3, Pair H: t10[r_c1c2,c_c1c2] <-> conj(t10[c_c1c2,r_c1c2]) ===== */
                                /* Only process lower triangle to avoid double-swap */
                                if (blc <= blr) {
                                    tmp = t10[elem_off(r_c1c2, c_c1c2)];
                                    t10[elem_off(r_c1c2, c_c1c2)] = conj(t10[elem_off(c_c1c2, r_c1c2)]);
                                    t10[elem_off(c_c1c2, r_c1c2)] = conj(tmp);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    else if (ctrl2 < LOG_TILE_DIM) {
        /*
         * Case 3:
         * One control index is larger than LOG_TILE_DIM & the other control and the target index are smaller than
         * LOG_TILE_DIM.
         *
         * Swap pairs are spread accross three tiles t10=(ctrl1=1|ctrl1=0), t11=(ctrl1=1|ctrl1=1) and
         * t01=(ctrl1=0|ctrl1=1), but each pair is contained in a single tile. Tile t01 may be in the upper traingle.
         */
        if (target < LOG_TILE_DIM) {
            /*
             * Case 3:
             * ctrl1 >= LOG_TILE_DIM (tile-level bit), ctrl2 < LOG_TILE_DIM, target < LOG_TILE_DIM.
             *
             * The 14 swap pairs distribute across three tiles:
             *   t10 = (tr1, tc0): pairs A, B, E, F — row swaps r_c2 <-> r_c2t
             *   t11 = (tr1, tc1): pairs C, D, G, H, M, N — mixed row/col/cross swaps
             *   t01 = (tr0, tc1): pairs I, J, K, L — col swaps c_c2 <-> c_c2t
             *
             * Each swap pair is self-contained within a single tile (no inter-tile swaps).
             * Tile t01 may be in the upper triangle when tr0 < tc1.
             */
            const idx_t n_bt = dim_tile >> 1;
            const qubit_t tile_bit = ctrl1 - LOG_TILE_DIM;
            const idx_t incr_tile = (idx_t)1 << tile_bit;
            const idx_t n_bl = TILE_DIM >> 2;
            const idx_t incr_ctrl2 = (idx_t)1 << ctrl2;
            const idx_t incr_tgt   = (idx_t)1 << target;
            const qubit_t lo_local = (ctrl2 < target) ? ctrl2 : target;
            const qubit_t hi_local = (ctrl2 > target) ? ctrl2 : target;

            #pragma omp parallel for schedule(static, 1) if(dim >= OMP_THRESHOLD)
            for (idx_t btr = 0; btr < n_bt; ++btr) {
                const idx_t tr0 = insertBit0(btr, tile_bit);
                const idx_t tr1 = tr0 | incr_tile;

                for (idx_t btc = 0; btc <= btr; ++btc) {
                    const idx_t tc0 = insertBit0(btc, tile_bit);
                    const idx_t tc1 = tc0 | incr_tile;

                    cplx_t * restrict t10 = data + tile_off(tr1, tc0);
                    cplx_t * restrict t11 = data + tile_off(tr1, tc1);

                    /* Fast path: tr0 > tc1 means t01 is in the lower triangle */
                    if (tr0 > tc1) {
                        cplx_t * restrict t01 = data + tile_off(tr0, tc1);

                        for (idx_t blr = 0; blr < n_bl; ++blr) {
                            const idx_t r000  = insertBits2_0(blr, lo_local, hi_local);
                            const idx_t r_c2  = r000 | incr_ctrl2;
                            const idx_t r_t   = r000 | incr_tgt;
                            const idx_t r_c2t = r000 | incr_ctrl2 | incr_tgt;

                            for (idx_t blc = 0; blc < n_bl; ++blc) {
                                const idx_t c000  = insertBits2_0(blc, lo_local, hi_local);
                                const idx_t c_c2  = c000 | incr_ctrl2;
                                const idx_t c_t   = c000 | incr_tgt;
                                const idx_t c_c2t = c000 | incr_ctrl2 | incr_tgt;

                                cplx_t tmp;

                                /* ===== t10 = (tr1, tc0): Pairs A, B, E, F — row swap r_c2 <-> r_c2t ===== */

                                /* Pair A: row c2=1,tgt=0 <-> c2=1,tgt=1, col c2=0,tgt=0 */
                                tmp = t10[elem_off(r_c2, c000)];
                                t10[elem_off(r_c2, c000)] = t10[elem_off(r_c2t, c000)];
                                t10[elem_off(r_c2t, c000)] = tmp;

                                /* Pair B: row c2=1,tgt=0 <-> c2=1,tgt=1, col c2=1,tgt=0 */
                                tmp = t10[elem_off(r_c2, c_c2)];
                                t10[elem_off(r_c2, c_c2)] = t10[elem_off(r_c2t, c_c2)];
                                t10[elem_off(r_c2t, c_c2)] = tmp;

                                /* Pair E: row c2=1,tgt=0 <-> c2=1,tgt=1, col c2=0,tgt=1 */
                                tmp = t10[elem_off(r_c2, c_t)];
                                t10[elem_off(r_c2, c_t)] = t10[elem_off(r_c2t, c_t)];
                                t10[elem_off(r_c2t, c_t)] = tmp;

                                /* Pair F: row c2=1,tgt=0 <-> c2=1,tgt=1, col c2=1,tgt=1 */
                                tmp = t10[elem_off(r_c2, c_c2t)];
                                t10[elem_off(r_c2, c_c2t)] = t10[elem_off(r_c2t, c_c2t)];
                                t10[elem_off(r_c2t, c_c2t)] = tmp;

                                /* ===== t11 = (tr1, tc1): Pairs C, G — row swap r_c2 <-> r_c2t ===== */

                                /* Pair C: row c2=1,tgt=0 <-> c2=1,tgt=1, col c2=0,tgt=0 */
                                tmp = t11[elem_off(r_c2, c000)];
                                t11[elem_off(r_c2, c000)] = t11[elem_off(r_c2t, c000)];
                                t11[elem_off(r_c2t, c000)] = tmp;

                                /* Pair G: row c2=1,tgt=0 <-> c2=1,tgt=1, col c2=0,tgt=1 */
                                tmp = t11[elem_off(r_c2, c_t)];
                                t11[elem_off(r_c2, c_t)] = t11[elem_off(r_c2t, c_t)];
                                t11[elem_off(r_c2t, c_t)] = tmp;

                                /* ===== t11: Pairs M, N — col swap c_c2 <-> c_c2t ===== */

                                /* Pair M: col c2=1,tgt=0 <-> c2=1,tgt=1, row c2=0,tgt=0 */
                                tmp = t11[elem_off(r000, c_c2)];
                                t11[elem_off(r000, c_c2)] = t11[elem_off(r000, c_c2t)];
                                t11[elem_off(r000, c_c2t)] = tmp;

                                /* Pair N: col c2=1,tgt=0 <-> c2=1,tgt=1, row c2=0,tgt=1 */
                                tmp = t11[elem_off(r_t, c_c2)];
                                t11[elem_off(r_t, c_c2)] = t11[elem_off(r_t, c_c2t)];
                                t11[elem_off(r_t, c_c2t)] = tmp;

                                /* ===== t11: Pairs D, H — both row and col change ===== */

                                /* Pair D: (r_c2, c_c2) <-> (r_c2t, c_c2t) */
                                tmp = t11[elem_off(r_c2, c_c2)];
                                t11[elem_off(r_c2, c_c2)] = t11[elem_off(r_c2t, c_c2t)];
                                t11[elem_off(r_c2t, c_c2t)] = tmp;

                                /* Pair H: (r_c2t, c_c2) <-> (r_c2, c_c2t) */
                                tmp = t11[elem_off(r_c2t, c_c2)];
                                t11[elem_off(r_c2t, c_c2)] = t11[elem_off(r_c2, c_c2t)];
                                t11[elem_off(r_c2, c_c2t)] = tmp;

                                /* ===== t01 = (tr0, tc1): Pairs I, J, K, L — col swap c_c2 <-> c_c2t ===== */

                                /* Pair I: col c2=1,tgt=0 <-> c2=1,tgt=1, row c2=0,tgt=0 */
                                tmp = t01[elem_off(r000, c_c2)];
                                t01[elem_off(r000, c_c2)] = t01[elem_off(r000, c_c2t)];
                                t01[elem_off(r000, c_c2t)] = tmp;

                                /* Pair J: col c2=1,tgt=0 <-> c2=1,tgt=1, row c2=0,tgt=1 */
                                tmp = t01[elem_off(r_t, c_c2)];
                                t01[elem_off(r_t, c_c2)] = t01[elem_off(r_t, c_c2t)];
                                t01[elem_off(r_t, c_c2t)] = tmp;

                                /* Pair K: col c2=1,tgt=0 <-> c2=1,tgt=1, row c2=1,tgt=0 */
                                tmp = t01[elem_off(r_c2, c_c2)];
                                t01[elem_off(r_c2, c_c2)] = t01[elem_off(r_c2, c_c2t)];
                                t01[elem_off(r_c2, c_c2t)] = tmp;

                                /* Pair L: col c2=1,tgt=0 <-> c2=1,tgt=1, row c2=1,tgt=1 */
                                tmp = t01[elem_off(r_c2t, c_c2)];
                                t01[elem_off(r_c2t, c_c2)] = t01[elem_off(r_c2t, c_c2t)];
                                t01[elem_off(r_c2t, c_c2t)] = tmp;
                            }
                        }
                    }

                    /* t01 in upper triangle, use adjoint tile(tc1, tr0) */
                    else if (btr != btc) {
                        cplx_t * restrict t01 = data + tile_off(tc1, tr0);

                        for (idx_t blr = 0; blr < n_bl; ++blr) {
                            const idx_t r000  = insertBits2_0(blr, lo_local, hi_local);
                            const idx_t r_c2  = r000 | incr_ctrl2;
                            const idx_t r_t   = r000 | incr_tgt;
                            const idx_t r_c2t = r000 | incr_ctrl2 | incr_tgt;

                            for (idx_t blc = 0; blc < n_bl; ++blc) {
                                const idx_t c000  = insertBits2_0(blc, lo_local, hi_local);
                                const idx_t c_c2  = c000 | incr_ctrl2;
                                const idx_t c_t   = c000 | incr_tgt;
                                const idx_t c_c2t = c000 | incr_ctrl2 | incr_tgt;

                                cplx_t tmp;

                                /* ===== t10: Pairs A, B, E, F — row swap (normal, always lower tri) ===== */

                                /* Pair A */
                                tmp = t10[elem_off(r_c2, c000)];
                                t10[elem_off(r_c2, c000)] = t10[elem_off(r_c2t, c000)];
                                t10[elem_off(r_c2t, c000)] = tmp;

                                /* Pair B */
                                tmp = t10[elem_off(r_c2, c_c2)];
                                t10[elem_off(r_c2, c_c2)] = t10[elem_off(r_c2t, c_c2)];
                                t10[elem_off(r_c2t, c_c2)] = tmp;

                                /* Pair E */
                                tmp = t10[elem_off(r_c2, c_t)];
                                t10[elem_off(r_c2, c_t)] = t10[elem_off(r_c2t, c_t)];
                                t10[elem_off(r_c2t, c_t)] = tmp;

                                /* Pair F */
                                tmp = t10[elem_off(r_c2, c_c2t)];
                                t10[elem_off(r_c2, c_c2t)] = t10[elem_off(r_c2t, c_c2t)];
                                t10[elem_off(r_c2t, c_c2t)] = tmp;

                                /* ===== t11: Pairs C, G, M, N, D, H — normal swaps (always lower tri) ===== */

                                /* Pair C */
                                tmp = t11[elem_off(r_c2, c000)];
                                t11[elem_off(r_c2, c000)] = t11[elem_off(r_c2t, c000)];
                                t11[elem_off(r_c2t, c000)] = tmp;

                                /* Pair G */
                                tmp = t11[elem_off(r_c2, c_t)];
                                t11[elem_off(r_c2, c_t)] = t11[elem_off(r_c2t, c_t)];
                                t11[elem_off(r_c2t, c_t)] = tmp;

                                /* Pair M */
                                tmp = t11[elem_off(r000, c_c2)];
                                t11[elem_off(r000, c_c2)] = t11[elem_off(r000, c_c2t)];
                                t11[elem_off(r000, c_c2t)] = tmp;

                                /* Pair N */
                                tmp = t11[elem_off(r_t, c_c2)];
                                t11[elem_off(r_t, c_c2)] = t11[elem_off(r_t, c_c2t)];
                                t11[elem_off(r_t, c_c2t)] = tmp;

                                /* Pair D */
                                tmp = t11[elem_off(r_c2, c_c2)];
                                t11[elem_off(r_c2, c_c2)] = t11[elem_off(r_c2t, c_c2t)];
                                t11[elem_off(r_c2t, c_c2t)] = tmp;

                                /* Pair H */
                                tmp = t11[elem_off(r_c2t, c_c2)];
                                t11[elem_off(r_c2t, c_c2)] = t11[elem_off(r_c2, c_c2t)];
                                t11[elem_off(r_c2, c_c2t)] = tmp;

                                /* ===== t01 (adjoint): Pairs I, J, K, L — col swap via conj-transpose ===== */
                                /* t01[lr,lc] = conj(adjt01[lc,lr]) where adjt01 = tile(tc1, tr0) */

                                /* Pair I: col c_c2 <-> c_c2t, row r000 */
                                tmp = t01[elem_off(c_c2, r000)];
                                t01[elem_off(c_c2, r000)] = t01[elem_off(c_c2t, r000)];
                                t01[elem_off(c_c2t, r000)] = tmp;

                                /* Pair J: col c_c2 <-> c_c2t, row r_t */
                                tmp = t01[elem_off(c_c2, r_t)];
                                t01[elem_off(c_c2, r_t)] = t01[elem_off(c_c2t, r_t)];
                                t01[elem_off(c_c2t, r_t)] = tmp;

                                /* Pair K: col c_c2 <-> c_c2t, row r_c2 */
                                tmp = t01[elem_off(c_c2, r_c2)];
                                t01[elem_off(c_c2, r_c2)] = t01[elem_off(c_c2t, r_c2)];
                                t01[elem_off(c_c2t, r_c2)] = tmp;

                                /* Pair L: col c_c2 <-> c_c2t, row r_c2t */
                                tmp = t01[elem_off(c_c2, r_c2t)];
                                t01[elem_off(c_c2, r_c2t)] = t01[elem_off(c_c2t, r_c2t)];
                                t01[elem_off(c_c2t, r_c2t)] = tmp;
                            }
                        }
                    }

                    /*
                     * Diagonal: btr == btc, so t01 = (tr0, tc1) is upper triangle
                     * and its adjoint is t(tc1, tr0) = t(tr1, tc0) = t10.
                     *
                     * Pairs I,J,K,L (col swaps in t01) translate to row swaps
                     * in t10 with transposed loop indices — identical to pairs
                     * A,B,E,F. Doing both would double-swap (cancel out), so we
                     * skip I,J,K,L entirely.
                     *
                     * t11 = (tr1, tc1) is a diagonal tile (tr1 == tc1 since
                     * btr == btc), but all 6 pairs (C,D,G,H,M,N) write to
                     * distinct memory locations and maintain Hermiticity
                     * together, so no special guards are needed.
                     */
                    else {
                        for (idx_t blr = 0; blr < n_bl; ++blr) {
                            const idx_t r000  = insertBits2_0(blr, lo_local, hi_local);
                            const idx_t r_c2  = r000 | incr_ctrl2;
                            const idx_t r_t   = r000 | incr_tgt;
                            const idx_t r_c2t = r000 | incr_ctrl2 | incr_tgt;

                            for (idx_t blc = 0; blc < n_bl; ++blc) {
                                const idx_t c000  = insertBits2_0(blc, lo_local, hi_local);
                                const idx_t c_c2  = c000 | incr_ctrl2;
                                const idx_t c_t   = c000 | incr_tgt;
                                const idx_t c_c2t = c000 | incr_ctrl2 | incr_tgt;

                                cplx_t tmp;

                                /* ===== t10: Pairs A, B, E, F — row swap r_c2 <-> r_c2t ===== */

                                /* Pair A */
                                tmp = t10[elem_off(r_c2, c000)];
                                t10[elem_off(r_c2, c000)] = t10[elem_off(r_c2t, c000)];
                                t10[elem_off(r_c2t, c000)] = tmp;

                                /* Pair B */
                                tmp = t10[elem_off(r_c2, c_c2)];
                                t10[elem_off(r_c2, c_c2)] = t10[elem_off(r_c2t, c_c2)];
                                t10[elem_off(r_c2t, c_c2)] = tmp;

                                /* Pair E */
                                tmp = t10[elem_off(r_c2, c_t)];
                                t10[elem_off(r_c2, c_t)] = t10[elem_off(r_c2t, c_t)];
                                t10[elem_off(r_c2t, c_t)] = tmp;

                                /* Pair F */
                                tmp = t10[elem_off(r_c2, c_c2t)];
                                t10[elem_off(r_c2, c_c2t)] = t10[elem_off(r_c2t, c_c2t)];
                                t10[elem_off(r_c2t, c_c2t)] = tmp;

                                /* ===== t11: Pairs C, G, M, N, D, H ===== */

                                /* Pair C */
                                tmp = t11[elem_off(r_c2, c000)];
                                t11[elem_off(r_c2, c000)] = t11[elem_off(r_c2t, c000)];
                                t11[elem_off(r_c2t, c000)] = tmp;

                                /* Pair G */
                                tmp = t11[elem_off(r_c2, c_t)];
                                t11[elem_off(r_c2, c_t)] = t11[elem_off(r_c2t, c_t)];
                                t11[elem_off(r_c2t, c_t)] = tmp;

                                /* Pair M */
                                tmp = t11[elem_off(r000, c_c2)];
                                t11[elem_off(r000, c_c2)] = t11[elem_off(r000, c_c2t)];
                                t11[elem_off(r000, c_c2t)] = tmp;

                                /* Pair N */
                                tmp = t11[elem_off(r_t, c_c2)];
                                t11[elem_off(r_t, c_c2)] = t11[elem_off(r_t, c_c2t)];
                                t11[elem_off(r_t, c_c2t)] = tmp;

                                /* Pair D */
                                tmp = t11[elem_off(r_c2, c_c2)];
                                t11[elem_off(r_c2, c_c2)] = t11[elem_off(r_c2t, c_c2t)];
                                t11[elem_off(r_c2t, c_c2t)] = tmp;

                                /* Pair H */
                                tmp = t11[elem_off(r_c2t, c_c2)];
                                t11[elem_off(r_c2t, c_c2)] = t11[elem_off(r_c2, c_c2t)];
                                t11[elem_off(r_c2, c_c2t)] = tmp;
                            }
                        }
                    }
                }
            }
        }

        /*
         * Case 4:
         * ctrl2 < LOG_TILE_DIM (local bit), ctrl1 >= LOG_TILE_DIM and target >= LOG_TILE_DIM (tile bits).
         *
         * Two tile-level bits (ctrl1, target) create a 4x4 grid of tile pairs:
         *   trXY = tile row with ctrl1=X, target=Y
         *   tcXY = tile col with ctrl1=X, target=Y
         *
         * The permutation flips the target tile-bit when ctrl1 tile-bit=1 AND
         * ctrl2 local-bit=1.  This distributes the 14 swap pairs across 10
         * distinct tile pairs:
         *
         *   Always lower triangle:
         *     (tr10,tc00) <-> (tr11,tc00): Pairs A, B
         *     (tr10,tc10) <-> (tr11,tc10): Pair C
         *     (tr10,tc10) <-> (tr11,tc11): Pair D
         *     (tr11,tc10) <-> (tr11,tc11): Pair N
         *
         *   May cross the diagonal:
         *     (tr11,tc01) <-> (tr10,tc01): Pairs E, F
         *     (tr11,tc11) <-> (tr10,tc11): Pair G
         *     (tr11,tc10) <-> (tr10,tc11): Pair H
         *     (tr00,tc10) <-> (tr00,tc11): Pairs I, K
         *     (tr01,tc10) <-> (tr01,tc11): Pairs J, L
         *     (tr10,tc10) <-> (tr10,tc11): Pair M
         */
        else {
            const idx_t n_bt = dim_tile >> 2;
            const qubit_t ctrl1_tile = ctrl1 - LOG_TILE_DIM;
            const qubit_t tgt_tile   = target - LOG_TILE_DIM;
            const qubit_t lo_tile = (ctrl1_tile < tgt_tile) ? ctrl1_tile : tgt_tile;
            const qubit_t hi_tile = (ctrl1_tile > tgt_tile) ? ctrl1_tile : tgt_tile;
            const idx_t incr_ctrl1_tile = (idx_t)1 << ctrl1_tile;
            const idx_t incr_tgt_tile   = (idx_t)1 << tgt_tile;
            const idx_t n_bl = TILE_DIM >> 1;
            const idx_t incr_ctrl2 = (idx_t)1 << ctrl2;

            #pragma omp parallel for schedule(static, 1) if(dim >= OMP_THRESHOLD)
            for (idx_t btr = 0; btr < n_bt; ++btr) {
                const idx_t tr00 = insertBits2_0(btr, lo_tile, hi_tile);
                const idx_t tr01 = tr00 | incr_tgt_tile;
                const idx_t tr10 = tr00 | incr_ctrl1_tile;
                const idx_t tr11 = tr10 | incr_tgt_tile;

                for (idx_t btc = 0; btc <= btr; ++btc) {
                    const idx_t tc00 = insertBits2_0(btc, lo_tile, hi_tile);
                    const idx_t tc01 = tc00 | incr_tgt_tile;
                    const idx_t tc10 = tc00 | incr_ctrl1_tile;
                    const idx_t tc11 = tc10 | incr_tgt_tile;

                    /* Pointers to always-lower-triangle tiles */
                    cplx_t * restrict t10_00 = data + tile_off(tr10, tc00);
                    cplx_t * restrict t11_00 = data + tile_off(tr11, tc00);
                    cplx_t * restrict t10_10 = data + tile_off(tr10, tc10);
                    cplx_t * restrict t11_10 = data + tile_off(tr11, tc10);
                    cplx_t * restrict t11_11 = data + tile_off(tr11, tc11);

                    /* Fast path: tr00 > tc11 means all tiles are in the lower triangle */
                    if (tr00 > tc11) {
                        cplx_t * restrict t11_01 = data + tile_off(tr11, tc01);
                        cplx_t * restrict t10_01 = data + tile_off(tr10, tc01);
                        cplx_t * restrict t10_11 = data + tile_off(tr10, tc11);
                        cplx_t * restrict t00_10 = data + tile_off(tr00, tc10);
                        cplx_t * restrict t00_11 = data + tile_off(tr00, tc11);
                        cplx_t * restrict t01_10 = data + tile_off(tr01, tc10);
                        cplx_t * restrict t01_11 = data + tile_off(tr01, tc11);
                        for (idx_t blr = 0; blr < n_bl; ++blr) {
                            const idx_t r000 = insertBit0(blr, ctrl2);
                            const idx_t r_c2 = r000 | incr_ctrl2;

                            for (idx_t blc = 0; blc < n_bl; ++blc) {
                                const idx_t c000 = insertBit0(blc, ctrl2);
                                const idx_t c_c2 = c000 | incr_ctrl2;

                                cplx_t tmp;

                                /* ===== Pairs A, B: (tr10,tc00) <-> (tr11,tc00) ===== */

                                /* Pair A: [r_c2, c000] */
                                tmp = t10_00[elem_off(r_c2, c000)];
                                t10_00[elem_off(r_c2, c000)] = t11_00[elem_off(r_c2, c000)];
                                t11_00[elem_off(r_c2, c000)] = tmp;

                                /* Pair B: [r_c2, c_c2] */
                                tmp = t10_00[elem_off(r_c2, c_c2)];
                                t10_00[elem_off(r_c2, c_c2)] = t11_00[elem_off(r_c2, c_c2)];
                                t11_00[elem_off(r_c2, c_c2)] = tmp;

                                /* ===== Pair C: (tr10,tc10) <-> (tr11,tc10) [r_c2, c000] ===== */

                                tmp = t10_10[elem_off(r_c2, c000)];
                                t10_10[elem_off(r_c2, c000)] = t11_10[elem_off(r_c2, c000)];
                                t11_10[elem_off(r_c2, c000)] = tmp;

                                /* ===== Pair D: (tr10,tc10) <-> (tr11,tc11) [r_c2, c_c2] ===== */

                                tmp = t10_10[elem_off(r_c2, c_c2)];
                                t10_10[elem_off(r_c2, c_c2)] = t11_11[elem_off(r_c2, c_c2)];
                                t11_11[elem_off(r_c2, c_c2)] = tmp;

                                /* ===== Pairs E, F: (tr11,tc01) <-> (tr10,tc01) ===== */

                                /* Pair E: [r_c2, c000] */
                                tmp = t11_01[elem_off(r_c2, c000)];
                                t11_01[elem_off(r_c2, c000)] = t10_01[elem_off(r_c2, c000)];
                                t10_01[elem_off(r_c2, c000)] = tmp;

                                /* Pair F: [r_c2, c_c2] */
                                tmp = t11_01[elem_off(r_c2, c_c2)];
                                t11_01[elem_off(r_c2, c_c2)] = t10_01[elem_off(r_c2, c_c2)];
                                t10_01[elem_off(r_c2, c_c2)] = tmp;

                                /* ===== Pair G: (tr11,tc11) <-> (tr10,tc11) [r_c2, c000] ===== */

                                tmp = t11_11[elem_off(r_c2, c000)];
                                t11_11[elem_off(r_c2, c000)] = t10_11[elem_off(r_c2, c000)];
                                t10_11[elem_off(r_c2, c000)] = tmp;

                                /* ===== Pair H: (tr11,tc10) <-> (tr10,tc11) [r_c2, c_c2] ===== */

                                tmp = t11_10[elem_off(r_c2, c_c2)];
                                t11_10[elem_off(r_c2, c_c2)] = t10_11[elem_off(r_c2, c_c2)];
                                t10_11[elem_off(r_c2, c_c2)] = tmp;

                                /* ===== Pairs I, K: (tr00,tc10) <-> (tr00,tc11) ===== */

                                /* Pair I: [r000, c_c2] */
                                tmp = t00_10[elem_off(r000, c_c2)];
                                t00_10[elem_off(r000, c_c2)] = t00_11[elem_off(r000, c_c2)];
                                t00_11[elem_off(r000, c_c2)] = tmp;

                                /* Pair K: [r_c2, c_c2] */
                                tmp = t00_10[elem_off(r_c2, c_c2)];
                                t00_10[elem_off(r_c2, c_c2)] = t00_11[elem_off(r_c2, c_c2)];
                                t00_11[elem_off(r_c2, c_c2)] = tmp;

                                /* ===== Pairs J, L: (tr01,tc10) <-> (tr01,tc11) ===== */

                                /* Pair J: [r000, c_c2] */
                                tmp = t01_10[elem_off(r000, c_c2)];
                                t01_10[elem_off(r000, c_c2)] = t01_11[elem_off(r000, c_c2)];
                                t01_11[elem_off(r000, c_c2)] = tmp;

                                /* Pair L: [r_c2, c_c2] */
                                tmp = t01_10[elem_off(r_c2, c_c2)];
                                t01_10[elem_off(r_c2, c_c2)] = t01_11[elem_off(r_c2, c_c2)];
                                t01_11[elem_off(r_c2, c_c2)] = tmp;

                                /* ===== Pair M: (tr10,tc10) <-> (tr10,tc11) [r000, c_c2] ===== */

                                tmp = t10_10[elem_off(r000, c_c2)];
                                t10_10[elem_off(r000, c_c2)] = t10_11[elem_off(r000, c_c2)];
                                t10_11[elem_off(r000, c_c2)] = tmp;

                                /* ===== Pair N: (tr11,tc10) <-> (tr11,tc11) [r000, c_c2] ===== */

                                tmp = t11_10[elem_off(r000, c_c2)];
                                t11_10[elem_off(r000, c_c2)] = t11_11[elem_off(r000, c_c2)];
                                t11_11[elem_off(r000, c_c2)] = tmp;
                            }
                        }
                    }

                    /*
                     * Off-diagonal: btr != btc but some tiles may be in the upper triangle.
                     *
                     * Upper-triangle tiles (tc > tr) are accessed via their adjoint:
                     *   tile(tr,tc)[lr,lc] = conj(tile(tc,tr)[lc,lr])
                     *
                     * Tiles that may be upper tri: (tr10,tc11), (tr10,tc01), (tr00,tc10),
                     * (tr00,tc11), (tr01,tc10), (tr01,tc11).
                     *
                     * Note: (tr10,tc01) is upper iff incr_tgt_tile > incr_ctrl1_tile
                     *   (i.e., target > ctrl1 in tile-bit space).
                     *   (tr01,tc10) is upper iff incr_ctrl1_tile > incr_tgt_tile.
                     *   Exactly one of these two is upper when they differ.
                     */
                    else if (btr != btc) {
                        for (idx_t blr = 0; blr < n_bl; ++blr) {
                            const idx_t r000 = insertBit0(blr, ctrl2);
                            const idx_t r_c2 = r000 | incr_ctrl2;

                            for (idx_t blc = 0; blc < n_bl; ++blc) {
                                const idx_t c000 = insertBit0(blc, ctrl2);
                                const idx_t c_c2 = c000 | incr_ctrl2;

                                cplx_t tmp;

                                /* ===== Always-lower-tri pairs ===== */

                                /* Pairs A, B: (tr10,tc00) <-> (tr11,tc00) */
                                tmp = t10_00[elem_off(r_c2, c000)];
                                t10_00[elem_off(r_c2, c000)] = t11_00[elem_off(r_c2, c000)];
                                t11_00[elem_off(r_c2, c000)] = tmp;

                                tmp = t10_00[elem_off(r_c2, c_c2)];
                                t10_00[elem_off(r_c2, c_c2)] = t11_00[elem_off(r_c2, c_c2)];
                                t11_00[elem_off(r_c2, c_c2)] = tmp;

                                /* Pair C: (tr10,tc10) <-> (tr11,tc10) [r_c2, c000] */
                                tmp = t10_10[elem_off(r_c2, c000)];
                                t10_10[elem_off(r_c2, c000)] = t11_10[elem_off(r_c2, c000)];
                                t11_10[elem_off(r_c2, c000)] = tmp;

                                /* Pair D: (tr10,tc10) <-> (tr11,tc11) [r_c2, c_c2] */
                                tmp = t10_10[elem_off(r_c2, c_c2)];
                                t10_10[elem_off(r_c2, c_c2)] = t11_11[elem_off(r_c2, c_c2)];
                                t11_11[elem_off(r_c2, c_c2)] = tmp;

                                /* Pair N: (tr11,tc10) <-> (tr11,tc11) [r000, c_c2] */
                                tmp = t11_10[elem_off(r000, c_c2)];
                                t11_10[elem_off(r000, c_c2)] = t11_11[elem_off(r000, c_c2)];
                                t11_11[elem_off(r000, c_c2)] = tmp;
                            }
                        }

                        /*
                         * Pairs E, F: (tr11,tc01) <-> (tr10,tc01)
                         * tr11 > tc01 always. tr10 vs tc01 depends on bit ordering.
                         */
                        if (tr10 > tc01) {
                            cplx_t * restrict t11_01 = data + tile_off(tr11, tc01);
                            cplx_t * restrict t10_01 = data + tile_off(tr10, tc01);
                            for (idx_t blr = 0; blr < n_bl; ++blr) {
                                const idx_t r_c2 = insertBit0(blr, ctrl2) | incr_ctrl2;
                                for (idx_t blc = 0; blc < n_bl; ++blc) {
                                    const idx_t c000 = insertBit0(blc, ctrl2);
                                    const idx_t c_c2 = c000 | incr_ctrl2;
                                    cplx_t tmp;
                                    /* Pair E */
                                    tmp = t11_01[elem_off(r_c2, c000)];
                                    t11_01[elem_off(r_c2, c000)] = t10_01[elem_off(r_c2, c000)];
                                    t10_01[elem_off(r_c2, c000)] = tmp;
                                    /* Pair F */
                                    tmp = t11_01[elem_off(r_c2, c_c2)];
                                    t11_01[elem_off(r_c2, c_c2)] = t10_01[elem_off(r_c2, c_c2)];
                                    t10_01[elem_off(r_c2, c_c2)] = tmp;
                                }
                            }
                        }
                        else {
                            /* (tr10,tc01) is upper tri; adjoint is (tc01,tr10) */
                            cplx_t * restrict t11_01 = data + tile_off(tr11, tc01);
                            cplx_t * restrict adj10_01 = data + tile_off(tc01, tr10);
                            for (idx_t blr = 0; blr < n_bl; ++blr) {
                                const idx_t r_c2 = insertBit0(blr, ctrl2) | incr_ctrl2;
                                for (idx_t blc = 0; blc < n_bl; ++blc) {
                                    const idx_t c000 = insertBit0(blc, ctrl2);
                                    const idx_t c_c2 = c000 | incr_ctrl2;
                                    cplx_t tmp;
                                    /* Pair E: t11_01[r_c2,c000] <-> conj(adj10_01[c000,r_c2]) */
                                    tmp = t11_01[elem_off(r_c2, c000)];
                                    t11_01[elem_off(r_c2, c000)] = conj(adj10_01[elem_off(c000, r_c2)]);
                                    adj10_01[elem_off(c000, r_c2)] = conj(tmp);
                                    /* Pair F: t11_01[r_c2,c_c2] <-> conj(adj10_01[c_c2,r_c2]) */
                                    tmp = t11_01[elem_off(r_c2, c_c2)];
                                    t11_01[elem_off(r_c2, c_c2)] = conj(adj10_01[elem_off(c_c2, r_c2)]);
                                    adj10_01[elem_off(c_c2, r_c2)] = conj(tmp);
                                }
                            }
                        }

                        /*
                         * Pairs G, H, M: involve (tr10,tc11) which is always upper tri
                         * when btr != btc (since tc11 > tr10 requires tc11 = tc10|incr_tgt > tr10,
                         * but with btr > btc the base difference may make tr10 > tc11).
                         *
                         * Check tr10 > tc11 to determine.
                         */
                        if (tr10 > tc11) {
                            cplx_t * restrict t10_11 = data + tile_off(tr10, tc11);
                            for (idx_t blr = 0; blr < n_bl; ++blr) {
                                const idx_t r000 = insertBit0(blr, ctrl2);
                                const idx_t r_c2 = r000 | incr_ctrl2;
                                for (idx_t blc = 0; blc < n_bl; ++blc) {
                                    const idx_t c000 = insertBit0(blc, ctrl2);
                                    const idx_t c_c2 = c000 | incr_ctrl2;
                                    cplx_t tmp;
                                    /* Pair G: (tr11,tc11)[r_c2,c000] <-> (tr10,tc11)[r_c2,c000] */
                                    tmp = t11_11[elem_off(r_c2, c000)];
                                    t11_11[elem_off(r_c2, c000)] = t10_11[elem_off(r_c2, c000)];
                                    t10_11[elem_off(r_c2, c000)] = tmp;
                                    /* Pair H: (tr11,tc10)[r_c2,c_c2] <-> (tr10,tc11)[r_c2,c_c2] */
                                    tmp = t11_10[elem_off(r_c2, c_c2)];
                                    t11_10[elem_off(r_c2, c_c2)] = t10_11[elem_off(r_c2, c_c2)];
                                    t10_11[elem_off(r_c2, c_c2)] = tmp;
                                    /* Pair M: (tr10,tc10)[r000,c_c2] <-> (tr10,tc11)[r000,c_c2] */
                                    tmp = t10_10[elem_off(r000, c_c2)];
                                    t10_10[elem_off(r000, c_c2)] = t10_11[elem_off(r000, c_c2)];
                                    t10_11[elem_off(r000, c_c2)] = tmp;
                                }
                            }
                        }
                        else {
                            /* (tr10,tc11) is upper tri; adjoint is (tc11,tr10) */
                            cplx_t * restrict adj10_11 = data + tile_off(tc11, tr10);
                            for (idx_t blr = 0; blr < n_bl; ++blr) {
                                const idx_t r000 = insertBit0(blr, ctrl2);
                                const idx_t r_c2 = r000 | incr_ctrl2;
                                for (idx_t blc = 0; blc < n_bl; ++blc) {
                                    const idx_t c000 = insertBit0(blc, ctrl2);
                                    const idx_t c_c2 = c000 | incr_ctrl2;
                                    cplx_t tmp;
                                    /* Pair G: t11_11[r_c2,c000] <-> conj(adj10_11[c000,r_c2]) */
                                    tmp = t11_11[elem_off(r_c2, c000)];
                                    t11_11[elem_off(r_c2, c000)] = conj(adj10_11[elem_off(c000, r_c2)]);
                                    adj10_11[elem_off(c000, r_c2)] = conj(tmp);
                                    /* Pair H: t11_10[r_c2,c_c2] <-> conj(adj10_11[c_c2,r_c2]) */
                                    tmp = t11_10[elem_off(r_c2, c_c2)];
                                    t11_10[elem_off(r_c2, c_c2)] = conj(adj10_11[elem_off(c_c2, r_c2)]);
                                    adj10_11[elem_off(c_c2, r_c2)] = conj(tmp);
                                    /* Pair M: t10_10[r000,c_c2] <-> conj(adj10_11[c_c2,r000]) */
                                    tmp = t10_10[elem_off(r000, c_c2)];
                                    t10_10[elem_off(r000, c_c2)] = conj(adj10_11[elem_off(c_c2, r000)]);
                                    adj10_11[elem_off(c_c2, r000)] = conj(tmp);
                                }
                            }
                        }

                        /*
                         * Pairs I, K: (tr00,tc10) <-> (tr00,tc11)
                         * Pairs J, L: (tr01,tc10) <-> (tr01,tc11)
                         *
                         * tc10 and tc11 have ctrl1 bit set; tr00 and tr01 do not.
                         * These tiles may be in the upper triangle depending on
                         * btr vs btc magnitudes.
                         *
                         * Key relationships (off-diagonal, btr > btc):
                         *   tr00 > tc10 => tr01 > tc10, tr01 > tc11
                         *     (adding tgt bit to both sides preserves ordering)
                         *   tr00 <= tc10 => (tr00,tc10) upper. Also (tr00,tc11) upper.
                         *     tr01 vs tc10 and tr01 vs tc11 require runtime checks.
                         */
                        if (tr00 > tc10) {
                            /* (tr00,tc10): lower. (tr00,tc11): upper (since tr00 <= tc11).
                             * (tr01,tc10): lower (tr01 > tc10 since tr00 > tc10).
                             * (tr01,tc11): lower (tr01 > tc11 by same argument). */
                            cplx_t * restrict t00_10 = data + tile_off(tr00, tc10);
                            cplx_t * restrict adj00_11 = data + tile_off(tc11, tr00);
                            cplx_t * restrict t01_10 = data + tile_off(tr01, tc10);
                            cplx_t * restrict t01_11 = data + tile_off(tr01, tc11);

                            for (idx_t blr = 0; blr < n_bl; ++blr) {
                                const idx_t r000 = insertBit0(blr, ctrl2);
                                const idx_t r_c2 = r000 | incr_ctrl2;
                                for (idx_t blc = 0; blc < n_bl; ++blc) {
                                    const idx_t c_c2 = insertBit0(blc, ctrl2) | incr_ctrl2;
                                    cplx_t tmp;

                                    /* Pair I: t00_10[r000,c_c2] <-> conj(adj00_11[c_c2,r000]) */
                                    tmp = t00_10[elem_off(r000, c_c2)];
                                    t00_10[elem_off(r000, c_c2)] = conj(adj00_11[elem_off(c_c2, r000)]);
                                    adj00_11[elem_off(c_c2, r000)] = conj(tmp);

                                    /* Pair K: t00_10[r_c2,c_c2] <-> conj(adj00_11[c_c2,r_c2]) */
                                    tmp = t00_10[elem_off(r_c2, c_c2)];
                                    t00_10[elem_off(r_c2, c_c2)] = conj(adj00_11[elem_off(c_c2, r_c2)]);
                                    adj00_11[elem_off(c_c2, r_c2)] = conj(tmp);

                                    /* Pair J: t01_10[r000,c_c2] <-> t01_11[r000,c_c2] */
                                    tmp = t01_10[elem_off(r000, c_c2)];
                                    t01_10[elem_off(r000, c_c2)] = t01_11[elem_off(r000, c_c2)];
                                    t01_11[elem_off(r000, c_c2)] = tmp;

                                    /* Pair L: t01_10[r_c2,c_c2] <-> t01_11[r_c2,c_c2] */
                                    tmp = t01_10[elem_off(r_c2, c_c2)];
                                    t01_10[elem_off(r_c2, c_c2)] = t01_11[elem_off(r_c2, c_c2)];
                                    t01_11[elem_off(r_c2, c_c2)] = tmp;
                                }
                            }
                        }
                        else {
                            /* tr00 <= tc10: (tr00,tc10) and (tr00,tc11) are upper tri.
                             * adj(tr00,tc10) = (tc10,tr00), adj(tr00,tc11) = (tc11,tr00). */
                            cplx_t * restrict adj00_10 = data + tile_off(tc10, tr00);
                            cplx_t * restrict adj00_11 = data + tile_off(tc11, tr00);

                            /* Pairs I, K: swap between adjoints directly */
                            for (idx_t blr = 0; blr < n_bl; ++blr) {
                                const idx_t r000 = insertBit0(blr, ctrl2);
                                const idx_t r_c2 = r000 | incr_ctrl2;
                                for (idx_t blc = 0; blc < n_bl; ++blc) {
                                    const idx_t c_c2 = insertBit0(blc, ctrl2) | incr_ctrl2;
                                    cplx_t tmp;
                                    /* Pair I: adj00_10[c_c2,r000] <-> adj00_11[c_c2,r000] */
                                    tmp = adj00_10[elem_off(c_c2, r000)];
                                    adj00_10[elem_off(c_c2, r000)] = adj00_11[elem_off(c_c2, r000)];
                                    adj00_11[elem_off(c_c2, r000)] = tmp;
                                    /* Pair K: adj00_10[c_c2,r_c2] <-> adj00_11[c_c2,r_c2] */
                                    tmp = adj00_10[elem_off(c_c2, r_c2)];
                                    adj00_10[elem_off(c_c2, r_c2)] = adj00_11[elem_off(c_c2, r_c2)];
                                    adj00_11[elem_off(c_c2, r_c2)] = tmp;
                                }
                            }

                            /* Pairs J, L: (tr01,tc10) <-> (tr01,tc11)
                             * tr01 vs tc10: runtime check needed. */
                            if (tr01 > tc11) {
                                /* Both (tr01,tc10) and (tr01,tc11) are lower tri */
                                cplx_t * restrict t01_10 = data + tile_off(tr01, tc10);
                                cplx_t * restrict t01_11 = data + tile_off(tr01, tc11);
                                for (idx_t blr = 0; blr < n_bl; ++blr) {
                                    const idx_t r000 = insertBit0(blr, ctrl2);
                                    const idx_t r_c2 = r000 | incr_ctrl2;
                                    for (idx_t blc = 0; blc < n_bl; ++blc) {
                                        const idx_t c_c2 = insertBit0(blc, ctrl2) | incr_ctrl2;
                                        cplx_t tmp;
                                        tmp = t01_10[elem_off(r000, c_c2)];
                                        t01_10[elem_off(r000, c_c2)] = t01_11[elem_off(r000, c_c2)];
                                        t01_11[elem_off(r000, c_c2)] = tmp;
                                        tmp = t01_10[elem_off(r_c2, c_c2)];
                                        t01_10[elem_off(r_c2, c_c2)] = t01_11[elem_off(r_c2, c_c2)];
                                        t01_11[elem_off(r_c2, c_c2)] = tmp;
                                    }
                                }
                            }
                            else if (tr01 > tc10) {
                                /* (tr01,tc10) lower, (tr01,tc11) upper */
                                cplx_t * restrict t01_10 = data + tile_off(tr01, tc10);
                                cplx_t * restrict adj01_11 = data + tile_off(tc11, tr01);
                                for (idx_t blr = 0; blr < n_bl; ++blr) {
                                    const idx_t r000 = insertBit0(blr, ctrl2);
                                    const idx_t r_c2 = r000 | incr_ctrl2;
                                    for (idx_t blc = 0; blc < n_bl; ++blc) {
                                        const idx_t c_c2 = insertBit0(blc, ctrl2) | incr_ctrl2;
                                        cplx_t tmp;
                                        /* Pair J */
                                        tmp = t01_10[elem_off(r000, c_c2)];
                                        t01_10[elem_off(r000, c_c2)] = conj(adj01_11[elem_off(c_c2, r000)]);
                                        adj01_11[elem_off(c_c2, r000)] = conj(tmp);
                                        /* Pair L */
                                        tmp = t01_10[elem_off(r_c2, c_c2)];
                                        t01_10[elem_off(r_c2, c_c2)] = conj(adj01_11[elem_off(c_c2, r_c2)]);
                                        adj01_11[elem_off(c_c2, r_c2)] = conj(tmp);
                                    }
                                }
                            }
                            else {
                                /* Both (tr01,tc10) and (tr01,tc11) are upper tri */
                                cplx_t * restrict adj01_10 = data + tile_off(tc10, tr01);
                                cplx_t * restrict adj01_11 = data + tile_off(tc11, tr01);
                                for (idx_t blr = 0; blr < n_bl; ++blr) {
                                    const idx_t r000 = insertBit0(blr, ctrl2);
                                    const idx_t r_c2 = r000 | incr_ctrl2;
                                    for (idx_t blc = 0; blc < n_bl; ++blc) {
                                        const idx_t c_c2 = insertBit0(blc, ctrl2) | incr_ctrl2;
                                        cplx_t tmp;
                                        tmp = adj01_10[elem_off(c_c2, r000)];
                                        adj01_10[elem_off(c_c2, r000)] = adj01_11[elem_off(c_c2, r000)];
                                        adj01_11[elem_off(c_c2, r000)] = tmp;
                                        tmp = adj01_10[elem_off(c_c2, r_c2)];
                                        adj01_10[elem_off(c_c2, r_c2)] = adj01_11[elem_off(c_c2, r_c2)];
                                        adj01_11[elem_off(c_c2, r_c2)] = tmp;
                                    }
                                }
                            }
                        }
                    }

                    /*
                     * Diagonal: btr == btc.
                     *
                     * When btr == btc: tr00=tc00, tr01=tc01, tr10=tc10, tr11=tc11.
                     * Upper-triangle tiles alias lower-triangle tiles in the grid:
                     *   adj(tr10,tc11) = (tr11,tc10) = t11_10
                     *   adj(tr00,tc10) = (tr10,tc00) = t10_00
                     *   adj(tr00,tc11) = (tr11,tc00) = t11_00
                     *   adj(tr01,tc11) = (tr11,tc01)
                     *   One of (tr01,tc10)/(tr10,tc01) is upper; its adj is the other.
                     *
                     * Pairs I,K are the Hermitian conjugates of Pairs A,B (they
                     * access the same tiles at transposed indices and coincide when
                     * blr==blc). Skip I,K on diagonal.
                     *
                     * Pairs J,L are similarly redundant with Pairs E,F. Skip J,L.
                     *
                     * Pairs G,M read from t11_10 as adjoint of (tr10,tc11), but
                     * t11_10 is modified by Pairs C,N,H. Phase 1 handles G,M
                     * as one-way writes before t11_10 is modified.
                     *
                     * Pair H is a self-adjoint swap within t11_10 (both sides map
                     * to the same tile with transposed indices).
                     *
                     * Pairs E,F access tiles t11_01 and t10_01 (or its adjoint
                     * t01_10), which are not touched by any other pair. Done as
                     * independent two-way swaps.
                     */
                    else {
                        /* Phase 1: one-way writes for G and M (read from t11_10 before it's modified) */
                        for (idx_t blr = 0; blr < n_bl; ++blr) {
                            const idx_t r000 = insertBit0(blr, ctrl2);
                            const idx_t r_c2 = r000 | incr_ctrl2;

                            for (idx_t blc = 0; blc < n_bl; ++blc) {
                                const idx_t c000 = insertBit0(blc, ctrl2);
                                const idx_t c_c2 = c000 | incr_ctrl2;

                                /* Pair G: t11_11[r_c2,c000] = conj(t11_10[c000,r_c2]) */
                                t11_11[elem_off(r_c2, c000)] = conj(t11_10[elem_off(c000, r_c2)]);

                                /* Pair M: t10_10[r000,c_c2] = conj(t11_10[c_c2,r000]) */
                                t10_10[elem_off(r000, c_c2)] = conj(t11_10[elem_off(c_c2, r000)]);
                            }
                        }

                        /*
                         * Pairs E, F: (tr11,tc01) <-> (tr10,tc01).
                         * Tiles t11_01 and t10_01 (or t01_10 adjoint) are NOT modified
                         * by any other pair in either phase, so E/F can be done as
                         * independent two-way swaps at any point.
                         *
                         * If ctrl1 > target in tile space: both (tr11,tc01) and (tr10,tc01)
                         * are lower tri (ctrl1 bit > tgt bit). Direct swap.
                         *
                         * If target > ctrl1 in tile space: (tr10,tc01) is upper tri,
                         * adjoint = (tr01,tc10) = t01_10. Two-way swap with conj-transpose.
                         *
                         * J/L are skipped on diagonal (redundant with E/F).
                         */
                        {
                            cplx_t * restrict t11_01 = data + tile_off(tr11, tc01);
                            if (incr_ctrl1_tile > incr_tgt_tile) {
                                /* ctrl1 > target: both tiles lower tri, direct swap */
                                cplx_t * restrict t10_01 = data + tile_off(tr10, tc01);
                                for (idx_t blr = 0; blr < n_bl; ++blr) {
                                    const idx_t r_c2 = insertBit0(blr, ctrl2) | incr_ctrl2;
                                    for (idx_t blc = 0; blc < n_bl; ++blc) {
                                        const idx_t c000 = insertBit0(blc, ctrl2);
                                        const idx_t c_c2 = c000 | incr_ctrl2;
                                        cplx_t tmp;
                                        /* Pair E */
                                        tmp = t11_01[elem_off(r_c2, c000)];
                                        t11_01[elem_off(r_c2, c000)] = t10_01[elem_off(r_c2, c000)];
                                        t10_01[elem_off(r_c2, c000)] = tmp;
                                        /* Pair F */
                                        tmp = t11_01[elem_off(r_c2, c_c2)];
                                        t11_01[elem_off(r_c2, c_c2)] = t10_01[elem_off(r_c2, c_c2)];
                                        t10_01[elem_off(r_c2, c_c2)] = tmp;
                                    }
                                }
                            }
                            else {
                                /* target > ctrl1: (tr10,tc01) upper, adj = (tr01,tc10) */
                                cplx_t * restrict t01_10 = data + tile_off(tr01, tc10);
                                for (idx_t blr = 0; blr < n_bl; ++blr) {
                                    const idx_t r_c2 = insertBit0(blr, ctrl2) | incr_ctrl2;
                                    for (idx_t blc = 0; blc < n_bl; ++blc) {
                                        const idx_t c000 = insertBit0(blc, ctrl2);
                                        const idx_t c_c2 = c000 | incr_ctrl2;
                                        cplx_t tmp;
                                        /* Pair E: t11_01[r_c2,c000] <-> conj(t01_10[c000,r_c2]) */
                                        tmp = t11_01[elem_off(r_c2, c000)];
                                        t11_01[elem_off(r_c2, c000)] = conj(t01_10[elem_off(c000, r_c2)]);
                                        t01_10[elem_off(c000, r_c2)] = conj(tmp);
                                        /* Pair F: t11_01[r_c2,c_c2] <-> conj(t01_10[c_c2,r_c2]) */
                                        tmp = t11_01[elem_off(r_c2, c_c2)];
                                        t11_01[elem_off(r_c2, c_c2)] = conj(t01_10[elem_off(c_c2, r_c2)]);
                                        t01_10[elem_off(c_c2, r_c2)] = conj(tmp);
                                    }
                                }
                            }
                        }

                        /* Phase 2: two-way swaps */
                        for (idx_t blr = 0; blr < n_bl; ++blr) {
                            const idx_t r000 = insertBit0(blr, ctrl2);
                            const idx_t r_c2 = r000 | incr_ctrl2;

                            for (idx_t blc = 0; blc < n_bl; ++blc) {
                                const idx_t c000 = insertBit0(blc, ctrl2);
                                const idx_t c_c2 = c000 | incr_ctrl2;

                                cplx_t tmp;

                                /* Pairs A, B: (tr10,tc00) <-> (tr11,tc00) — always lower */
                                tmp = t10_00[elem_off(r_c2, c000)];
                                t10_00[elem_off(r_c2, c000)] = t11_00[elem_off(r_c2, c000)];
                                t11_00[elem_off(r_c2, c000)] = tmp;

                                tmp = t10_00[elem_off(r_c2, c_c2)];
                                t10_00[elem_off(r_c2, c_c2)] = t11_00[elem_off(r_c2, c_c2)];
                                t11_00[elem_off(r_c2, c_c2)] = tmp;

                                /* Pair C: (tr10,tc10) <-> (tr11,tc10) [r_c2, c000] */
                                tmp = t10_10[elem_off(r_c2, c000)];
                                t10_10[elem_off(r_c2, c000)] = t11_10[elem_off(r_c2, c000)];
                                t11_10[elem_off(r_c2, c000)] = tmp;

                                /* Pair D: (tr10,tc10) <-> (tr11,tc11) [r_c2, c_c2] */
                                tmp = t10_10[elem_off(r_c2, c_c2)];
                                t10_10[elem_off(r_c2, c_c2)] = t11_11[elem_off(r_c2, c_c2)];
                                t11_11[elem_off(r_c2, c_c2)] = tmp;

                                /* Pair N: (tr11,tc10) <-> (tr11,tc11) [r000, c_c2] */
                                tmp = t11_10[elem_off(r000, c_c2)];
                                t11_10[elem_off(r000, c_c2)] = t11_11[elem_off(r000, c_c2)];
                                t11_11[elem_off(r000, c_c2)] = tmp;

                                /* Pair H: self-adjoint swap within t11_10
                                 * t11_10[r_c2,c_c2] <-> conj(t11_10[c_c2,r_c2])
                                 * Only process lower triangle (blc <= blr) to avoid double-swap */
                                if (blc <= blr) {
                                    tmp = t11_10[elem_off(r_c2, c_c2)];
                                    t11_10[elem_off(r_c2, c_c2)] = conj(t11_10[elem_off(c_c2, r_c2)]);
                                    t11_10[elem_off(c_c2, r_c2)] = conj(tmp);
                                }
                            }
                        }

                        /* Note: Pairs I, K, J, L are SKIPPED on diagonal — they are the
                         * Hermitian conjugates of Pairs A, B, E, F respectively, and performing
                         * both would double-swap (cancel). */
                    }
                }
            }
        }
    }
    else {
        /* Case 5:
         * Both control indices are at least LOG_TILE_DIM & target index is smaller than LOG_TILE_DIM.
         *
         * Swap pairs are spread across 7 tiles t1100=(ctrl1=1,ctrl2=1|ctrl1=0,ctrl2=0), t1101, t1110, t1111, t1011,
         * t0111, t0011, but each pair is contained in a single tile. Tiles my be in the upper triangle.
         */
        if (target < LOG_TILE_DIM) {
            /*
             * Case 5:
             * ctrl1 > ctrl2 >= LOG_TILE_DIM (both tile-level bits), target < LOG_TILE_DIM (local bit).
             *
             * The two tile-level control bits create a 4x4 grid of tile rows/cols:
             *   trXY = tile row with ctrl1=X, ctrl2=Y (and similarly tcXY).
             *
             * The permutation flips the target LOCAL bit when the tile row/col
             * has both ctrl1=1 and ctrl2=1 (i.e., in tr11 or tc11).  Since the
             * target bit is local, the tile address is unchanged by the permutation
             * -- all 14 swap pairs are intra-tile swaps.
             *
             * The 14 pairs distribute across 7 tile locations:
             *
             *   Always lower triangle (tr11 > tcXX for XX != 11):
             *     (tr11, tc00): Pairs A, E — row swap lr0 <-> lr1, col fixed
             *     (tr11, tc01): Pairs B, F — row swap lr0 <-> lr1, col fixed
             *     (tr11, tc10): Pairs C, G — row swap lr0 <-> lr1, col fixed
             *
             *   Lower triangle when btr >= btc:
             *     (tr11, tc11): Pairs D, H — both row and col swap
             *
             *   May be upper triangle:
             *     (tr00, tc11): Pairs I, J — col swap lc0 <-> lc1, row fixed
             *     (tr01, tc11): Pairs K, L — col swap lc0 <-> lc1, row fixed
             *     (tr10, tc11): Pairs M, N — col swap lc0 <-> lc1, row fixed
             *
             * Row swap in tile: elem(lr0, lc) <-> elem(lr1, lc)
             *   where lr0 = insertBit0(blr, target), lr1 = lr0 | incr_tgt
             * Col swap in tile: elem(lr, lc0) <-> elem(lr, lc1)
             *   where lc0 = insertBit0(blc, target), lc1 = lc0 | incr_tgt
             * Both swap in diagonal tile (D,H): lr0,lc0 <-> lr1,lc1 and lr1,lc0 <-> lr0,lc1
             *
             * For upper-triangle tiles (trXX, tc11) where trXX < tc11:
             *   Access via adjoint tile(tc11, trXX) with transposed, conjugated elements.
             *   A col swap lc0<->lc1 at fixed row lr in the logical tile becomes a row
             *   swap in the adjoint: adj[lc0,lr] <-> adj[lc1,lr] (no conjugation needed
             *   since both values are just exchanged).
             *
             * On the diagonal (btr == btc), tiles (tr00,tc11), (tr01,tc11), (tr10,tc11)
             * are upper triangle with adjoints (tr11,tc00), (tr11,tc01), (tr11,tc10)
             * respectively. Pairs I,J / K,L / M,N (col swaps in the adjoint) coincide
             * with Pairs A,E / B,F / C,G (row swaps) — they are Hermitian conjugate
             * pairs and performing both would double-swap. Skip I,J,K,L,M,N on diagonal.
             */
            const idx_t n_bt = dim_tile >> 2;
            const qubit_t ctrl1_tile = ctrl1 - LOG_TILE_DIM;
            const qubit_t ctrl2_tile = ctrl2 - LOG_TILE_DIM;
            const qubit_t lo_tile = ctrl2_tile;   /* ctrl2_tile < ctrl1_tile */
            const qubit_t hi_tile = ctrl1_tile;
            const idx_t incr_ctrl1_tile = (idx_t)1 << ctrl1_tile;
            const idx_t incr_ctrl2_tile = (idx_t)1 << ctrl2_tile;
            const idx_t n_bl = TILE_DIM >> 1;
            const idx_t incr_tgt = (idx_t)1 << target;

            #pragma omp parallel for schedule(static, 1) if(dim >= OMP_THRESHOLD)
            for (idx_t btr = 0; btr < n_bt; ++btr) {
                const idx_t tr00 = insertBits2_0(btr, lo_tile, hi_tile);
                const idx_t tr01 = tr00 | incr_ctrl2_tile;
                const idx_t tr10 = tr00 | incr_ctrl1_tile;
                const idx_t tr11 = tr00 | incr_ctrl1_tile | incr_ctrl2_tile;

                for (idx_t btc = 0; btc <= btr; ++btc) {
                    const idx_t tc00 = insertBits2_0(btc, lo_tile, hi_tile);
                    const idx_t tc01 = tc00 | incr_ctrl2_tile;
                    const idx_t tc10 = tc00 | incr_ctrl1_tile;
                    const idx_t tc11 = tc00 | incr_ctrl1_tile | incr_ctrl2_tile;

                    /* Pointers to always-lower-triangle tiles */
                    cplx_t * restrict t11_00 = data + tile_off(tr11, tc00);
                    cplx_t * restrict t11_01 = data + tile_off(tr11, tc01);
                    cplx_t * restrict t11_10 = data + tile_off(tr11, tc10);
                    cplx_t * restrict t11_11 = data + tile_off(tr11, tc11);

                    /* Fast path: tr00 > tc11 means all 7 tiles are in the lower triangle */
                    if (tr00 > tc11) {
                        cplx_t * restrict t00_11 = data + tile_off(tr00, tc11);
                        cplx_t * restrict t01_11 = data + tile_off(tr01, tc11);
                        cplx_t * restrict t10_11 = data + tile_off(tr10, tc11);

                        for (idx_t blr = 0; blr < n_bl; ++blr) {
                            const idx_t lr0 = insertBit0(blr, target);
                            const idx_t lr1 = lr0 | incr_tgt;

                            for (idx_t blc = 0; blc < n_bl; ++blc) {
                                const idx_t lc0 = insertBit0(blc, target);
                                const idx_t lc1 = lc0 | incr_tgt;

                                cplx_t tmp;

                                /* ===== (tr11, tc00): Pairs A, E — row swap lr0 <-> lr1 ===== */

                                /* Pair A: col lc0 */
                                tmp = t11_00[elem_off(lr0, lc0)];
                                t11_00[elem_off(lr0, lc0)] = t11_00[elem_off(lr1, lc0)];
                                t11_00[elem_off(lr1, lc0)] = tmp;

                                /* Pair E: col lc1 */
                                tmp = t11_00[elem_off(lr0, lc1)];
                                t11_00[elem_off(lr0, lc1)] = t11_00[elem_off(lr1, lc1)];
                                t11_00[elem_off(lr1, lc1)] = tmp;

                                /* ===== (tr11, tc01): Pairs B, F — row swap lr0 <-> lr1 ===== */

                                /* Pair B: col lc0 */
                                tmp = t11_01[elem_off(lr0, lc0)];
                                t11_01[elem_off(lr0, lc0)] = t11_01[elem_off(lr1, lc0)];
                                t11_01[elem_off(lr1, lc0)] = tmp;

                                /* Pair F: col lc1 */
                                tmp = t11_01[elem_off(lr0, lc1)];
                                t11_01[elem_off(lr0, lc1)] = t11_01[elem_off(lr1, lc1)];
                                t11_01[elem_off(lr1, lc1)] = tmp;

                                /* ===== (tr11, tc10): Pairs C, G — row swap lr0 <-> lr1 ===== */

                                /* Pair C: col lc0 */
                                tmp = t11_10[elem_off(lr0, lc0)];
                                t11_10[elem_off(lr0, lc0)] = t11_10[elem_off(lr1, lc0)];
                                t11_10[elem_off(lr1, lc0)] = tmp;

                                /* Pair G: col lc1 */
                                tmp = t11_10[elem_off(lr0, lc1)];
                                t11_10[elem_off(lr0, lc1)] = t11_10[elem_off(lr1, lc1)];
                                t11_10[elem_off(lr1, lc1)] = tmp;

                                /* ===== (tr11, tc11): Pairs D, H — both row and col swap ===== */

                                /* Pair D: (lr0, lc0) <-> (lr1, lc1) */
                                tmp = t11_11[elem_off(lr0, lc0)];
                                t11_11[elem_off(lr0, lc0)] = t11_11[elem_off(lr1, lc1)];
                                t11_11[elem_off(lr1, lc1)] = tmp;

                                /* Pair H: (lr1, lc0) <-> (lr0, lc1) */
                                tmp = t11_11[elem_off(lr1, lc0)];
                                t11_11[elem_off(lr1, lc0)] = t11_11[elem_off(lr0, lc1)];
                                t11_11[elem_off(lr0, lc1)] = tmp;

                                /* ===== (tr00, tc11): Pairs I, J — col swap lc0 <-> lc1 ===== */

                                /* Pair I: row lr0 */
                                tmp = t00_11[elem_off(lr0, lc0)];
                                t00_11[elem_off(lr0, lc0)] = t00_11[elem_off(lr0, lc1)];
                                t00_11[elem_off(lr0, lc1)] = tmp;

                                /* Pair J: row lr1 */
                                tmp = t00_11[elem_off(lr1, lc0)];
                                t00_11[elem_off(lr1, lc0)] = t00_11[elem_off(lr1, lc1)];
                                t00_11[elem_off(lr1, lc1)] = tmp;

                                /* ===== (tr01, tc11): Pairs K, L — col swap lc0 <-> lc1 ===== */

                                /* Pair K: row lr0 */
                                tmp = t01_11[elem_off(lr0, lc0)];
                                t01_11[elem_off(lr0, lc0)] = t01_11[elem_off(lr0, lc1)];
                                t01_11[elem_off(lr0, lc1)] = tmp;

                                /* Pair L: row lr1 */
                                tmp = t01_11[elem_off(lr1, lc0)];
                                t01_11[elem_off(lr1, lc0)] = t01_11[elem_off(lr1, lc1)];
                                t01_11[elem_off(lr1, lc1)] = tmp;

                                /* ===== (tr10, tc11): Pairs M, N — col swap lc0 <-> lc1 ===== */

                                /* Pair M: row lr0 */
                                tmp = t10_11[elem_off(lr0, lc0)];
                                t10_11[elem_off(lr0, lc0)] = t10_11[elem_off(lr0, lc1)];
                                t10_11[elem_off(lr0, lc1)] = tmp;

                                /* Pair N: row lr1 */
                                tmp = t10_11[elem_off(lr1, lc0)];
                                t10_11[elem_off(lr1, lc0)] = t10_11[elem_off(lr1, lc1)];
                                t10_11[elem_off(lr1, lc1)] = tmp;
                            }
                        }
                    }

                    /*
                     * Off-diagonal: btr != btc, but some tiles may be in the upper triangle.
                     *
                     * Tiles (tr11, tcXX) for XX in {00,01,10} are always lower tri.
                     * Tile (tr11, tc11) is lower tri since btr > btc => tr11 > tc11.
                     *
                     * The three col-swap tiles may be upper tri:
                     *   (tr00, tc11): always upper in this branch (tr00 <= tc11)
                     *   (tr01, tc11): upper iff tr01 <= tc11, equivalently tr00 <= tc10
                     *   (tr10, tc11): upper iff tr10 <= tc11, equivalently tr00 <= tc01
                     *
                     * When (trXX, tc11) is upper, access its adjoint (tc11, trXX)
                     * which is a stored lower-triangle tile.  A col swap lc0<->lc1
                     * at row lr in the logical tile becomes a row swap in the adjoint:
                     *   adj[lc0, lr] <-> adj[lc1, lr]
                     */
                    else if (btr != btc) {
                        /* Always-lower-triangle row swaps and diagonal-tile swaps */
                        for (idx_t blr = 0; blr < n_bl; ++blr) {
                            const idx_t lr0 = insertBit0(blr, target);
                            const idx_t lr1 = lr0 | incr_tgt;

                            for (idx_t blc = 0; blc < n_bl; ++blc) {
                                const idx_t lc0 = insertBit0(blc, target);
                                const idx_t lc1 = lc0 | incr_tgt;

                                cplx_t tmp;

                                /* ===== (tr11, tc00): Pairs A, E ===== */
                                tmp = t11_00[elem_off(lr0, lc0)];
                                t11_00[elem_off(lr0, lc0)] = t11_00[elem_off(lr1, lc0)];
                                t11_00[elem_off(lr1, lc0)] = tmp;

                                tmp = t11_00[elem_off(lr0, lc1)];
                                t11_00[elem_off(lr0, lc1)] = t11_00[elem_off(lr1, lc1)];
                                t11_00[elem_off(lr1, lc1)] = tmp;

                                /* ===== (tr11, tc01): Pairs B, F ===== */
                                tmp = t11_01[elem_off(lr0, lc0)];
                                t11_01[elem_off(lr0, lc0)] = t11_01[elem_off(lr1, lc0)];
                                t11_01[elem_off(lr1, lc0)] = tmp;

                                tmp = t11_01[elem_off(lr0, lc1)];
                                t11_01[elem_off(lr0, lc1)] = t11_01[elem_off(lr1, lc1)];
                                t11_01[elem_off(lr1, lc1)] = tmp;

                                /* ===== (tr11, tc10): Pairs C, G ===== */
                                tmp = t11_10[elem_off(lr0, lc0)];
                                t11_10[elem_off(lr0, lc0)] = t11_10[elem_off(lr1, lc0)];
                                t11_10[elem_off(lr1, lc0)] = tmp;

                                tmp = t11_10[elem_off(lr0, lc1)];
                                t11_10[elem_off(lr0, lc1)] = t11_10[elem_off(lr1, lc1)];
                                t11_10[elem_off(lr1, lc1)] = tmp;

                                /* ===== (tr11, tc11): Pairs D, H ===== */
                                tmp = t11_11[elem_off(lr0, lc0)];
                                t11_11[elem_off(lr0, lc0)] = t11_11[elem_off(lr1, lc1)];
                                t11_11[elem_off(lr1, lc1)] = tmp;

                                tmp = t11_11[elem_off(lr1, lc0)];
                                t11_11[elem_off(lr1, lc0)] = t11_11[elem_off(lr0, lc1)];
                                t11_11[elem_off(lr0, lc1)] = tmp;
                            }
                        }

                        /*
                         * Pairs I, J: col swap in (tr00, tc11) — always upper in this branch.
                         * Adjoint is (tc11, tr00).  Col swap lc0<->lc1 at row lr becomes
                         * row swap adj[lc0,lr] <-> adj[lc1,lr].
                         */
                        {
                            cplx_t * restrict adj00_11 = data + tile_off(tc11, tr00);
                            for (idx_t blr = 0; blr < n_bl; ++blr) {
                                const idx_t lr0 = insertBit0(blr, target);
                                const idx_t lr1 = lr0 | incr_tgt;
                                for (idx_t blc = 0; blc < n_bl; ++blc) {
                                    const idx_t lc0 = insertBit0(blc, target);
                                    const idx_t lc1 = lc0 | incr_tgt;
                                    cplx_t tmp;
                                    /* Pair I: adj[lc0, lr0] <-> adj[lc1, lr0] */
                                    tmp = adj00_11[elem_off(lc0, lr0)];
                                    adj00_11[elem_off(lc0, lr0)] = adj00_11[elem_off(lc1, lr0)];
                                    adj00_11[elem_off(lc1, lr0)] = tmp;
                                    /* Pair J: adj[lc0, lr1] <-> adj[lc1, lr1] */
                                    tmp = adj00_11[elem_off(lc0, lr1)];
                                    adj00_11[elem_off(lc0, lr1)] = adj00_11[elem_off(lc1, lr1)];
                                    adj00_11[elem_off(lc1, lr1)] = tmp;
                                }
                            }
                        }

                        /*
                         * Pairs K, L: col swap in (tr01, tc11).
                         * Upper iff tr01 <= tc11 (equivalently tr00 <= tc10).
                         */
                        if (tr01 > tc11) {
                            /* (tr01, tc11) is lower tri — direct col swap */
                            cplx_t * restrict t01_11 = data + tile_off(tr01, tc11);
                            for (idx_t blr = 0; blr < n_bl; ++blr) {
                                const idx_t lr0 = insertBit0(blr, target);
                                const idx_t lr1 = lr0 | incr_tgt;
                                for (idx_t blc = 0; blc < n_bl; ++blc) {
                                    const idx_t lc0 = insertBit0(blc, target);
                                    const idx_t lc1 = lc0 | incr_tgt;
                                    cplx_t tmp;
                                    /* Pair K */
                                    tmp = t01_11[elem_off(lr0, lc0)];
                                    t01_11[elem_off(lr0, lc0)] = t01_11[elem_off(lr0, lc1)];
                                    t01_11[elem_off(lr0, lc1)] = tmp;
                                    /* Pair L */
                                    tmp = t01_11[elem_off(lr1, lc0)];
                                    t01_11[elem_off(lr1, lc0)] = t01_11[elem_off(lr1, lc1)];
                                    t01_11[elem_off(lr1, lc1)] = tmp;
                                }
                            }
                        }
                        else {
                            /* (tr01, tc11) is upper tri — row swap in adjoint (tc11, tr01) */
                            cplx_t * restrict adj01_11 = data + tile_off(tc11, tr01);
                            for (idx_t blr = 0; blr < n_bl; ++blr) {
                                const idx_t lr0 = insertBit0(blr, target);
                                const idx_t lr1 = lr0 | incr_tgt;
                                for (idx_t blc = 0; blc < n_bl; ++blc) {
                                    const idx_t lc0 = insertBit0(blc, target);
                                    const idx_t lc1 = lc0 | incr_tgt;
                                    cplx_t tmp;
                                    /* Pair K: adj[lc0, lr0] <-> adj[lc1, lr0] */
                                    tmp = adj01_11[elem_off(lc0, lr0)];
                                    adj01_11[elem_off(lc0, lr0)] = adj01_11[elem_off(lc1, lr0)];
                                    adj01_11[elem_off(lc1, lr0)] = tmp;
                                    /* Pair L: adj[lc0, lr1] <-> adj[lc1, lr1] */
                                    tmp = adj01_11[elem_off(lc0, lr1)];
                                    adj01_11[elem_off(lc0, lr1)] = adj01_11[elem_off(lc1, lr1)];
                                    adj01_11[elem_off(lc1, lr1)] = tmp;
                                }
                            }
                        }

                        /*
                         * Pairs M, N: col swap in (tr10, tc11).
                         * Upper iff tr10 <= tc11 (equivalently tr00 <= tc01).
                         */
                        if (tr10 > tc11) {
                            /* (tr10, tc11) is lower tri — direct col swap */
                            cplx_t * restrict t10_11 = data + tile_off(tr10, tc11);
                            for (idx_t blr = 0; blr < n_bl; ++blr) {
                                const idx_t lr0 = insertBit0(blr, target);
                                const idx_t lr1 = lr0 | incr_tgt;
                                for (idx_t blc = 0; blc < n_bl; ++blc) {
                                    const idx_t lc0 = insertBit0(blc, target);
                                    const idx_t lc1 = lc0 | incr_tgt;
                                    cplx_t tmp;
                                    /* Pair M */
                                    tmp = t10_11[elem_off(lr0, lc0)];
                                    t10_11[elem_off(lr0, lc0)] = t10_11[elem_off(lr0, lc1)];
                                    t10_11[elem_off(lr0, lc1)] = tmp;
                                    /* Pair N */
                                    tmp = t10_11[elem_off(lr1, lc0)];
                                    t10_11[elem_off(lr1, lc0)] = t10_11[elem_off(lr1, lc1)];
                                    t10_11[elem_off(lr1, lc1)] = tmp;
                                }
                            }
                        }
                        else {
                            /* (tr10, tc11) is upper tri — row swap in adjoint (tc11, tr10) */
                            cplx_t * restrict adj10_11 = data + tile_off(tc11, tr10);
                            for (idx_t blr = 0; blr < n_bl; ++blr) {
                                const idx_t lr0 = insertBit0(blr, target);
                                const idx_t lr1 = lr0 | incr_tgt;
                                for (idx_t blc = 0; blc < n_bl; ++blc) {
                                    const idx_t lc0 = insertBit0(blc, target);
                                    const idx_t lc1 = lc0 | incr_tgt;
                                    cplx_t tmp;
                                    /* Pair M: adj[lc0, lr0] <-> adj[lc1, lr0] */
                                    tmp = adj10_11[elem_off(lc0, lr0)];
                                    adj10_11[elem_off(lc0, lr0)] = adj10_11[elem_off(lc1, lr0)];
                                    adj10_11[elem_off(lc1, lr0)] = tmp;
                                    /* Pair N: adj[lc0, lr1] <-> adj[lc1, lr1] */
                                    tmp = adj10_11[elem_off(lc0, lr1)];
                                    adj10_11[elem_off(lc0, lr1)] = adj10_11[elem_off(lc1, lr1)];
                                    adj10_11[elem_off(lc1, lr1)] = tmp;
                                }
                            }
                        }
                    }

                    /*
                     * Diagonal: btr == btc.
                     *
                     * When btr == btc: tr00=tc00, tr01=tc01, tr10=tc10, tr11=tc11.
                     *
                     * Upper-triangle tiles alias lower-triangle tiles:
                     *   adj(tr00, tc11) = (tc11, tr00) = (tr11, tc00) = t11_00
                     *   adj(tr01, tc11) = (tc11, tr01) = (tr11, tc01) = t11_01
                     *   adj(tr10, tc11) = (tc11, tr10) = (tr11, tc10) = t11_10
                     *
                     * Col swaps I,J in adj(tr00,tc11) = row swaps in t11_00 = Pairs A,E.
                     * Col swaps K,L in adj(tr01,tc11) = row swaps in t11_01 = Pairs B,F.
                     * Col swaps M,N in adj(tr10,tc11) = row swaps in t11_10 = Pairs C,G.
                     *
                     * These are Hermitian conjugate pairs — performing both would
                     * double-swap (cancel). Skip I,J,K,L,M,N on diagonal.
                     *
                     * (tr11, tc11) is a diagonal tile (tr11 == tc11). Pairs D and H
                     * swap elements that are conjugate-related:
                     *   D: (lr0,lc0) <-> (lr1,lc1) — symmetric pair (lr0<lc0 iff lr1<lc1)
                     *   H: (lr1,lc0) <-> (lr0,lc1) — anti-diagonal pair
                     * Both pairs write to distinct memory locations within the diagonal
                     * tile and jointly preserve Hermiticity, so no special guards needed.
                     */
                    else {
                        for (idx_t blr = 0; blr < n_bl; ++blr) {
                            const idx_t lr0 = insertBit0(blr, target);
                            const idx_t lr1 = lr0 | incr_tgt;

                            for (idx_t blc = 0; blc < n_bl; ++blc) {
                                const idx_t lc0 = insertBit0(blc, target);
                                const idx_t lc1 = lc0 | incr_tgt;

                                cplx_t tmp;

                                /* ===== (tr11, tc00): Pairs A, E — row swap ===== */

                                /* Pair A */
                                tmp = t11_00[elem_off(lr0, lc0)];
                                t11_00[elem_off(lr0, lc0)] = t11_00[elem_off(lr1, lc0)];
                                t11_00[elem_off(lr1, lc0)] = tmp;

                                /* Pair E */
                                tmp = t11_00[elem_off(lr0, lc1)];
                                t11_00[elem_off(lr0, lc1)] = t11_00[elem_off(lr1, lc1)];
                                t11_00[elem_off(lr1, lc1)] = tmp;

                                /* ===== (tr11, tc01): Pairs B, F — row swap ===== */

                                /* Pair B */
                                tmp = t11_01[elem_off(lr0, lc0)];
                                t11_01[elem_off(lr0, lc0)] = t11_01[elem_off(lr1, lc0)];
                                t11_01[elem_off(lr1, lc0)] = tmp;

                                /* Pair F */
                                tmp = t11_01[elem_off(lr0, lc1)];
                                t11_01[elem_off(lr0, lc1)] = t11_01[elem_off(lr1, lc1)];
                                t11_01[elem_off(lr1, lc1)] = tmp;

                                /* ===== (tr11, tc10): Pairs C, G — row swap ===== */

                                /* Pair C */
                                tmp = t11_10[elem_off(lr0, lc0)];
                                t11_10[elem_off(lr0, lc0)] = t11_10[elem_off(lr1, lc0)];
                                t11_10[elem_off(lr1, lc0)] = tmp;

                                /* Pair G */
                                tmp = t11_10[elem_off(lr0, lc1)];
                                t11_10[elem_off(lr0, lc1)] = t11_10[elem_off(lr1, lc1)];
                                t11_10[elem_off(lr1, lc1)] = tmp;

                                /* ===== (tr11, tc11): Pairs D, H — both row and col swap ===== */

                                /* Pair D: (lr0, lc0) <-> (lr1, lc1) */
                                tmp = t11_11[elem_off(lr0, lc0)];
                                t11_11[elem_off(lr0, lc0)] = t11_11[elem_off(lr1, lc1)];
                                t11_11[elem_off(lr1, lc1)] = tmp;

                                /* Pair H: (lr1, lc0) <-> (lr0, lc1) */
                                tmp = t11_11[elem_off(lr1, lc0)];
                                t11_11[elem_off(lr1, lc0)] = t11_11[elem_off(lr0, lc1)];
                                t11_11[elem_off(lr0, lc1)] = tmp;
                            }
                        }

                        /* Note: Pairs I,J,K,L,M,N are SKIPPED on diagonal — they are the
                         * Hermitian conjugates of Pairs A,E,B,F,C,G respectively, and
                         * performing both would double-swap (cancel). */
                    }
                }
            }
        }

        /*
         * Case 6:
         * All three qubit indices are at least LOG_TILE_DIM (all tile-level bits).
         *
         * Three tile bits create an 8x8 grid of tile row/col pairs. The base
         * loop variable btr runs up to dim_tile >> 3.  Each swap pair operates
         * on entire tiles (TILE_SIZE elements).
         *
         * The 14 swap pairs at tile level:
         *
         *   Row changes (r_c1c2 <-> r_all), column fixed:
         *     A: tile(tr_c1c2, tc000)  <-> tile(tr_all, tc000)   -- always lower tri
         *     B: tile(tr_c1c2, tc_c2)  <-> tile(tr_all, tc_c2)   -- always lower tri
         *     C: tile(tr_c1c2, tc_c1)  <-> tile(tr_all, tc_c1)   -- always lower tri
         *     E: tile(tr_c1c2, tc_t)   <-> tile(tr_all, tc_t)    -- may cross diagonal
         *     F: tile(tr_c1c2, tc_c2t) <-> tile(tr_all, tc_c2t)  -- may cross diagonal
         *     G: tile(tr_c1c2, tc_c1t) <-> tile(tr_all, tc_c1t)  -- may cross diagonal
         *
         *   Both row and column change:
         *     D: tile(tr_c1c2, tc_c1c2) <-> tile(tr_all, tc_all)   -- always lower tri
         *     H: tile(tr_all, tc_c1c2)  <-> tile(tr_c1c2, tc_all)  -- may cross diagonal
         *
         *   Column changes (c_c1c2 <-> c_all), row fixed:
         *     I: tile(tr000, tc_c1c2) <-> tile(tr000, tc_all)   -- may cross diagonal
         *     J: tile(tr_t,   tc_c1c2) <-> tile(tr_t,   tc_all) -- may cross diagonal
         *     K: tile(tr_c2,  tc_c1c2) <-> tile(tr_c2,  tc_all) -- may cross diagonal
         *     L: tile(tr_c2t, tc_c1c2) <-> tile(tr_c2t, tc_all) -- may cross diagonal
         *     M: tile(tr_c1,  tc_c1c2) <-> tile(tr_c1,  tc_all) -- may cross diagonal
         *     N: tile(tr_c1t, tc_c1c2) <-> tile(tr_c1t, tc_all) -- may cross diagonal
         *
         * Triangle analysis:
         *
         * "Always lower tri" pairs (A-C, D):
         *   tr_c1c2 has both ctrl1 and ctrl2 tile bits set; tr_all has all three.
         *   For cols tc000, tc_c2, tc_c1 (missing at least one ctrl bit), the
         *   row always dominates:
         *     - tr_c1c2 > tc000 (both ctrl bits set vs none)
         *     - tr_c1c2 > tc_c2 (ctrl1 bit difference)
         *     - tr_c1c2 > tc_c1 (ctrl2 bit difference)
         *   Similarly tr_all > tc_XX for all these columns.
         *   For D: tr_c1c2 >= tc_c1c2 and tr_all >= tc_all (same bits, btr >= btc).
         *
         * "May cross diagonal" pairs (E-G):
         *   tr_c1c2 vs tc_t/tc_c2t/tc_c1t: depends on relative magnitudes of
         *   the tile bits. When tgt_tile > ctrl1_tile, incr_t > incr_c1 > incr_c2,
         *   so tc_t > tr_c1c2 is possible. Runtime check needed.
         *   However, tr_all (all 3 bits set) always > tc_XX (< 3 bits) from btr >= btc.
         *
         * "May cross diagonal" pairs (H, I-N):
         *   The column is tc_all or tc_c1c2 which may exceed the row.
         *
         * Key structural property: on the btr==btc diagonal, tr_XX != tc_YY
         * for any XX != YY because insertBits3_0 forces the three marker bits
         * to distinct patterns. This eliminates diagonal density-matrix tiles
         * in cross-label pairs and simplifies the diagonal case.
         */
        else {
            const idx_t ctrl1_tile = ctrl1 - LOG_TILE_DIM;
            const idx_t ctrl2_tile = ctrl2 - LOG_TILE_DIM;
            const idx_t tgt_tile   = target - LOG_TILE_DIM;

            /* Sort tile bits: lo_tile < mid_tile < hi_tile */
            qubit_t lo_tile, mid_tile, hi_tile;
            if (ctrl2_tile > tgt_tile) {
                lo_tile = tgt_tile; mid_tile = ctrl2_tile; hi_tile = ctrl1_tile;
            }
            else if (ctrl1_tile > tgt_tile) {
                lo_tile = ctrl2_tile; mid_tile = tgt_tile; hi_tile = ctrl1_tile;
            }
            else {
                lo_tile = ctrl2_tile; mid_tile = ctrl1_tile; hi_tile = tgt_tile;
            }

            const idx_t incr_c1 = (idx_t)1 << ctrl1_tile;
            const idx_t incr_c2 = (idx_t)1 << ctrl2_tile;
            const idx_t incr_t  = (idx_t)1 << tgt_tile;

            const idx_t n_bt = dim_tile >> 3;

            #pragma omp parallel for schedule(static, 1) if(dim >= OMP_THRESHOLD)
            for (idx_t btr = 0; btr < n_bt; ++btr) {
                /* insertBits3_0: insert 0-bits at positions lo_tile, mid_tile, hi_tile */
                const idx_t tr000  = insertBit0(insertBits2_0(btr, lo_tile, mid_tile), hi_tile);
                const idx_t tr_c1  = tr000 | incr_c1;
                const idx_t tr_c2  = tr000 | incr_c2;
                const idx_t tr_t   = tr000 | incr_t;
                const idx_t tr_c1c2 = tr000 | incr_c1 | incr_c2;
                const idx_t tr_c1t  = tr000 | incr_c1 | incr_t;
                const idx_t tr_c2t  = tr000 | incr_c2 | incr_t;
                const idx_t tr_all  = tr000 | incr_c1 | incr_c2 | incr_t;

                for (idx_t btc = 0; btc <= btr; ++btc) {
                    const idx_t tc000  = insertBit0(insertBits2_0(btc, lo_tile, mid_tile), hi_tile);
                    const idx_t tc_c1  = tc000 | incr_c1;
                    const idx_t tc_c2  = tc000 | incr_c2;
                    const idx_t tc_t   = tc000 | incr_t;
                    const idx_t tc_c1c2 = tc000 | incr_c1 | incr_c2;
                    const idx_t tc_c1t  = tc000 | incr_c1 | incr_t;
                    const idx_t tc_c2t  = tc000 | incr_c2 | incr_t;
                    const idx_t tc_all  = tc000 | incr_c1 | incr_c2 | incr_t;

                    /*
                     * =====================================================================
                     * Pairs A-C: Always lower triangle (row swap, whole-tile memswap).
                     *
                     * For these pairs, both tiles have tr_c1c2 or tr_all as tile
                     * row and tc000/tc_c2/tc_c1 as tile column. Since:
                     *   tr_c1c2 = tr000 | incr_c1 | incr_c2  (both ctrl bits set)
                     *   tc000/tc_c2/tc_c1 lack at least one ctrl bit
                     * the row always exceeds the column (ctrl1 is the highest bit,
                     * so having both ctrl bits set guarantees dominance over any
                     * column missing at least one ctrl bit). tr_all > tr_c1c2
                     * trivially.
                     * =====================================================================
                     */

                    /* Pair A: tile(tr_c1c2, tc000) <-> tile(tr_all, tc000) */
                    {
                        cplx_t * restrict tA = data + tile_off(tr_c1c2, tc000);
                        cplx_t * restrict tB = data + tile_off(tr_all, tc000);
                        for (idx_t i = 0; i < TILE_SIZE; ++i) {
                            cplx_t tmp = tA[i]; tA[i] = tB[i]; tB[i] = tmp;
                        }
                    }

                    /* Pair B: tile(tr_c1c2, tc_c2) <-> tile(tr_all, tc_c2) */
                    {
                        cplx_t * restrict tA = data + tile_off(tr_c1c2, tc_c2);
                        cplx_t * restrict tB = data + tile_off(tr_all, tc_c2);
                        for (idx_t i = 0; i < TILE_SIZE; ++i) {
                            cplx_t tmp = tA[i]; tA[i] = tB[i]; tB[i] = tmp;
                        }
                    }

                    /* Pair C: tile(tr_c1c2, tc_c1) <-> tile(tr_all, tc_c1) */
                    {
                        cplx_t * restrict tA = data + tile_off(tr_c1c2, tc_c1);
                        cplx_t * restrict tB = data + tile_off(tr_all, tc_c1);
                        for (idx_t i = 0; i < TILE_SIZE; ++i) {
                            cplx_t tmp = tA[i]; tA[i] = tB[i]; tB[i] = tmp;
                        }
                    }

                    /*
                     * =====================================================================
                     * Pair D: tile(tr_c1c2, tc_c1c2) <-> tile(tr_all, tc_all)
                     *
                     * Always lower tri: tr_c1c2 >= tc_c1c2 (from btr >= btc,
                     * same bits set) and tr_all >= tc_all (same argument).
                     * On diagonal (btr == btc): both are diagonal tiles, so
                     * we use the lower-triangle-only swap pattern.
                     * =====================================================================
                     */
                    if (btr != btc) {
                        cplx_t * restrict tA = data + tile_off(tr_c1c2, tc_c1c2);
                        cplx_t * restrict tB = data + tile_off(tr_all, tc_all);
                        for (idx_t i = 0; i < TILE_SIZE; ++i) {
                            cplx_t tmp = tA[i]; tA[i] = tB[i]; tB[i] = tmp;
                        }
                    }
                    else {
                        /* Diagonal: tr_c1c2 == tc_c1c2 and tr_all == tc_all,
                         * both tiles sit on the matrix diagonal. Swap only
                         * lower-triangle elements to preserve Hermitian storage. */
                        cplx_t * restrict tA = data + tile_off(tr_c1c2, tc_c1c2);
                        cplx_t * restrict tB = data + tile_off(tr_all, tc_all);
                        for (idx_t lr = 0; lr < TILE_DIM; ++lr) {
                            for (idx_t lc = 0; lc <= lr; ++lc) {
                                cplx_t tmp = tA[lr * TILE_DIM + lc];
                                tA[lr * TILE_DIM + lc] = tB[lr * TILE_DIM + lc];
                                tB[lr * TILE_DIM + lc] = tmp;
                            }
                        }
                    }

                    /*
                     * =====================================================================
                     * Fast path: tr000 > tc_all means ALL tiles are in the lower
                     * triangle.  The row base (tr000, with no bits set) exceeds
                     * the largest column (tc_all, with all bits set), so every
                     * (tr_XX, tc_YY) satisfies tr_XX >= tr000 > tc_all >= tc_YY.
                     *
                     * This covers pairs E-G, H, I-N as simple whole-tile memswaps.
                     * =====================================================================
                     */
                    if (tr000 > tc_all) {

                        /* Pair E: tile(tr_c1c2, tc_t) <-> tile(tr_all, tc_t) */
                        {
                            cplx_t * restrict tA = data + tile_off(tr_c1c2, tc_t);
                            cplx_t * restrict tB = data + tile_off(tr_all, tc_t);
                            for (idx_t i = 0; i < TILE_SIZE; ++i) {
                                cplx_t tmp = tA[i]; tA[i] = tB[i]; tB[i] = tmp;
                            }
                        }

                        /* Pair F: tile(tr_c1c2, tc_c2t) <-> tile(tr_all, tc_c2t) */
                        {
                            cplx_t * restrict tA = data + tile_off(tr_c1c2, tc_c2t);
                            cplx_t * restrict tB = data + tile_off(tr_all, tc_c2t);
                            for (idx_t i = 0; i < TILE_SIZE; ++i) {
                                cplx_t tmp = tA[i]; tA[i] = tB[i]; tB[i] = tmp;
                            }
                        }

                        /* Pair G: tile(tr_c1c2, tc_c1t) <-> tile(tr_all, tc_c1t) */
                        {
                            cplx_t * restrict tA = data + tile_off(tr_c1c2, tc_c1t);
                            cplx_t * restrict tB = data + tile_off(tr_all, tc_c1t);
                            for (idx_t i = 0; i < TILE_SIZE; ++i) {
                                cplx_t tmp = tA[i]; tA[i] = tB[i]; tB[i] = tmp;
                            }
                        }

                        /* Pair H: tile(tr_all, tc_c1c2) <-> tile(tr_c1c2, tc_all) */
                        {
                            cplx_t * restrict tA = data + tile_off(tr_all, tc_c1c2);
                            cplx_t * restrict tB = data + tile_off(tr_c1c2, tc_all);
                            for (idx_t i = 0; i < TILE_SIZE; ++i) {
                                cplx_t tmp = tA[i]; tA[i] = tB[i]; tB[i] = tmp;
                            }
                        }

                        /* Pair I: tile(tr000, tc_c1c2) <-> tile(tr000, tc_all) */
                        {
                            cplx_t * restrict tA = data + tile_off(tr000, tc_c1c2);
                            cplx_t * restrict tB = data + tile_off(tr000, tc_all);
                            for (idx_t i = 0; i < TILE_SIZE; ++i) {
                                cplx_t tmp = tA[i]; tA[i] = tB[i]; tB[i] = tmp;
                            }
                        }

                        /* Pair J: tile(tr_t, tc_c1c2) <-> tile(tr_t, tc_all) */
                        {
                            cplx_t * restrict tA = data + tile_off(tr_t, tc_c1c2);
                            cplx_t * restrict tB = data + tile_off(tr_t, tc_all);
                            for (idx_t i = 0; i < TILE_SIZE; ++i) {
                                cplx_t tmp = tA[i]; tA[i] = tB[i]; tB[i] = tmp;
                            }
                        }

                        /* Pair K: tile(tr_c2, tc_c1c2) <-> tile(tr_c2, tc_all) */
                        {
                            cplx_t * restrict tA = data + tile_off(tr_c2, tc_c1c2);
                            cplx_t * restrict tB = data + tile_off(tr_c2, tc_all);
                            for (idx_t i = 0; i < TILE_SIZE; ++i) {
                                cplx_t tmp = tA[i]; tA[i] = tB[i]; tB[i] = tmp;
                            }
                        }

                        /* Pair L: tile(tr_c2t, tc_c1c2) <-> tile(tr_c2t, tc_all) */
                        {
                            cplx_t * restrict tA = data + tile_off(tr_c2t, tc_c1c2);
                            cplx_t * restrict tB = data + tile_off(tr_c2t, tc_all);
                            for (idx_t i = 0; i < TILE_SIZE; ++i) {
                                cplx_t tmp = tA[i]; tA[i] = tB[i]; tB[i] = tmp;
                            }
                        }

                        /* Pair M: tile(tr_c1, tc_c1c2) <-> tile(tr_c1, tc_all) */
                        {
                            cplx_t * restrict tA = data + tile_off(tr_c1, tc_c1c2);
                            cplx_t * restrict tB = data + tile_off(tr_c1, tc_all);
                            for (idx_t i = 0; i < TILE_SIZE; ++i) {
                                cplx_t tmp = tA[i]; tA[i] = tB[i]; tB[i] = tmp;
                            }
                        }

                        /* Pair N: tile(tr_c1t, tc_c1c2) <-> tile(tr_c1t, tc_all) */
                        {
                            cplx_t * restrict tA = data + tile_off(tr_c1t, tc_c1c2);
                            cplx_t * restrict tB = data + tile_off(tr_c1t, tc_all);
                            for (idx_t i = 0; i < TILE_SIZE; ++i) {
                                cplx_t tmp = tA[i]; tA[i] = tB[i]; tB[i] = tmp;
                            }
                        }
                    }

                    /*
                     * =====================================================================
                     * General off-diagonal case (btr != btc, not all tiles lower tri).
                     *
                     * Pairs E, F, G: swap tile(tr_c1c2, tc_XX) <-> tile(tr_all, tc_XX).
                     *   tile(tr_all, tc_XX) is always lower tri (tr_all has all 3 bits
                     *   set and at least c1+c2 > any single/double bit lacking both ctrls).
                     *   tile(tr_c1c2, tc_XX) may be upper tri when tgt_tile is large
                     *   enough that tc_XX > tr_c1c2.
                     *
                     * Pair H: tile(tr_all, tc_c1c2) <-> tile(tr_c1c2, tc_all).
                     *   tile(tr_all, tc_c1c2) always lower tri. tile(tr_c1c2, tc_all)
                     *   may be upper tri.
                     *
                     * Pairs I-N: tile(tr_XX, tc_c1c2) <-> tile(tr_XX, tc_all).
                     *   Both may be upper tri.
                     * =====================================================================
                     */
                    else if (btr != btc) {

                        /*
                         * Pairs E, F, G: The tile(tr_all, tc_XX) side is always
                         * lower tri. The tile(tr_c1c2, tc_XX) side needs a check.
                         * If tr_c1c2 > tc_XX: both lower tri, direct swap.
                         * If tr_c1c2 <= tc_XX: tile(tr_c1c2, tc_XX) is upper tri;
                         *   use its adjoint tile(tc_XX, tr_c1c2) for the swap.
                         */

                        /* Pair E: tile(tr_c1c2, tc_t) <-> tile(tr_all, tc_t) */
                        if (tr_c1c2 > tc_t) {
                            cplx_t * restrict tA = data + tile_off(tr_c1c2, tc_t);
                            cplx_t * restrict tB = data + tile_off(tr_all, tc_t);
                            for (idx_t i = 0; i < TILE_SIZE; ++i) {
                                cplx_t tmp = tA[i]; tA[i] = tB[i]; tB[i] = tmp;
                            }
                        }
                        else {
                            cplx_t * restrict tA = data + tile_off(tc_t, tr_c1c2);
                            cplx_t * restrict tB = data + tile_off(tr_all, tc_t);
                            for (idx_t lr = 0; lr < TILE_DIM; ++lr) {
                                for (idx_t lc = 0; lc < TILE_DIM; ++lc) {
                                    cplx_t tmp = tB[lr * TILE_DIM + lc];
                                    tB[lr * TILE_DIM + lc] = conj(tA[lc * TILE_DIM + lr]);
                                    tA[lc * TILE_DIM + lr] = conj(tmp);
                                }
                            }
                        }

                        /* Pair F: tile(tr_c1c2, tc_c2t) <-> tile(tr_all, tc_c2t) */
                        if (tr_c1c2 > tc_c2t) {
                            cplx_t * restrict tA = data + tile_off(tr_c1c2, tc_c2t);
                            cplx_t * restrict tB = data + tile_off(tr_all, tc_c2t);
                            for (idx_t i = 0; i < TILE_SIZE; ++i) {
                                cplx_t tmp = tA[i]; tA[i] = tB[i]; tB[i] = tmp;
                            }
                        }
                        else {
                            cplx_t * restrict tA = data + tile_off(tc_c2t, tr_c1c2);
                            cplx_t * restrict tB = data + tile_off(tr_all, tc_c2t);
                            for (idx_t lr = 0; lr < TILE_DIM; ++lr) {
                                for (idx_t lc = 0; lc < TILE_DIM; ++lc) {
                                    cplx_t tmp = tB[lr * TILE_DIM + lc];
                                    tB[lr * TILE_DIM + lc] = conj(tA[lc * TILE_DIM + lr]);
                                    tA[lc * TILE_DIM + lr] = conj(tmp);
                                }
                            }
                        }

                        /* Pair G: tile(tr_c1c2, tc_c1t) <-> tile(tr_all, tc_c1t) */
                        if (tr_c1c2 > tc_c1t) {
                            cplx_t * restrict tA = data + tile_off(tr_c1c2, tc_c1t);
                            cplx_t * restrict tB = data + tile_off(tr_all, tc_c1t);
                            for (idx_t i = 0; i < TILE_SIZE; ++i) {
                                cplx_t tmp = tA[i]; tA[i] = tB[i]; tB[i] = tmp;
                            }
                        }
                        else {
                            cplx_t * restrict tA = data + tile_off(tc_c1t, tr_c1c2);
                            cplx_t * restrict tB = data + tile_off(tr_all, tc_c1t);
                            for (idx_t lr = 0; lr < TILE_DIM; ++lr) {
                                for (idx_t lc = 0; lc < TILE_DIM; ++lc) {
                                    cplx_t tmp = tB[lr * TILE_DIM + lc];
                                    tB[lr * TILE_DIM + lc] = conj(tA[lc * TILE_DIM + lr]);
                                    tA[lc * TILE_DIM + lr] = conj(tmp);
                                }
                            }
                        }

                        /*
                         * Pair H: tile(tr_all, tc_c1c2) <-> tile(tr_c1c2, tc_all)
                         *
                         * tile(tr_all, tc_c1c2) is always lower tri (btr > btc and
                         * tr_all has tgt bit extra over tc_c1c2).
                         * tile(tr_c1c2, tc_all): check tr_c1c2 > tc_all.
                         */
                        if (tr_c1c2 > tc_all) {
                            cplx_t * restrict tA = data + tile_off(tr_all, tc_c1c2);
                            cplx_t * restrict tB = data + tile_off(tr_c1c2, tc_all);
                            for (idx_t i = 0; i < TILE_SIZE; ++i) {
                                cplx_t tmp = tA[i]; tA[i] = tB[i]; tB[i] = tmp;
                            }
                        }
                        else {
                            /* tile(tr_c1c2, tc_all) upper tri; adjoint = tile(tc_all, tr_c1c2) */
                            cplx_t * restrict tA = data + tile_off(tr_all, tc_c1c2);
                            cplx_t * restrict tB = data + tile_off(tc_all, tr_c1c2);
                            for (idx_t lr = 0; lr < TILE_DIM; ++lr) {
                                for (idx_t lc = 0; lc < TILE_DIM; ++lc) {
                                    cplx_t tmp = tA[lr * TILE_DIM + lc];
                                    tA[lr * TILE_DIM + lc] = conj(tB[lc * TILE_DIM + lr]);
                                    tB[lc * TILE_DIM + lr] = conj(tmp);
                                }
                            }
                        }

                        /*
                         * Pairs I-N: column swaps (tc_c1c2 <-> tc_all), row fixed.
                         *
                         * For each pair, both tiles share the same row index tr_XX.
                         * Three cases:
                         *   (a) tr_XX > tc_all: both lower tri -> direct memswap
                         *   (b) tr_XX > tc_c1c2 but tr_XX <= tc_all: first lower,
                         *       second upper -> cross-diagonal swap
                         *   (c) tr_XX <= tc_c1c2: both upper tri -> swap adjoints
                         *
                         * tc_all > tc_c1c2 always (tc_all = tc_c1c2 | incr_t).
                         */

                        /* Pair I: tile(tr000, tc_c1c2) <-> tile(tr000, tc_all) */
                        if (tr000 > tc_all) {
                            cplx_t * restrict tA = data + tile_off(tr000, tc_c1c2);
                            cplx_t * restrict tB = data + tile_off(tr000, tc_all);
                            for (idx_t i = 0; i < TILE_SIZE; ++i) {
                                cplx_t tmp = tA[i]; tA[i] = tB[i]; tB[i] = tmp;
                            }
                        }
                        else if (tr000 > tc_c1c2) {
                            cplx_t * restrict tA = data + tile_off(tr000, tc_c1c2);
                            cplx_t * restrict tB = data + tile_off(tc_all, tr000);
                            for (idx_t lr = 0; lr < TILE_DIM; ++lr) {
                                for (idx_t lc = 0; lc < TILE_DIM; ++lc) {
                                    cplx_t tmp = tA[lr * TILE_DIM + lc];
                                    tA[lr * TILE_DIM + lc] = conj(tB[lc * TILE_DIM + lr]);
                                    tB[lc * TILE_DIM + lr] = conj(tmp);
                                }
                            }
                        }
                        else {
                            cplx_t * restrict tA = data + tile_off(tc_c1c2, tr000);
                            cplx_t * restrict tB = data + tile_off(tc_all, tr000);
                            for (idx_t i = 0; i < TILE_SIZE; ++i) {
                                cplx_t tmp = tA[i]; tA[i] = tB[i]; tB[i] = tmp;
                            }
                        }

                        /* Pair J: tile(tr_t, tc_c1c2) <-> tile(tr_t, tc_all) */
                        if (tr_t > tc_all) {
                            cplx_t * restrict tA = data + tile_off(tr_t, tc_c1c2);
                            cplx_t * restrict tB = data + tile_off(tr_t, tc_all);
                            for (idx_t i = 0; i < TILE_SIZE; ++i) {
                                cplx_t tmp = tA[i]; tA[i] = tB[i]; tB[i] = tmp;
                            }
                        }
                        else if (tr_t > tc_c1c2) {
                            cplx_t * restrict tA = data + tile_off(tr_t, tc_c1c2);
                            cplx_t * restrict tB = data + tile_off(tc_all, tr_t);
                            for (idx_t lr = 0; lr < TILE_DIM; ++lr) {
                                for (idx_t lc = 0; lc < TILE_DIM; ++lc) {
                                    cplx_t tmp = tA[lr * TILE_DIM + lc];
                                    tA[lr * TILE_DIM + lc] = conj(tB[lc * TILE_DIM + lr]);
                                    tB[lc * TILE_DIM + lr] = conj(tmp);
                                }
                            }
                        }
                        else {
                            cplx_t * restrict tA = data + tile_off(tc_c1c2, tr_t);
                            cplx_t * restrict tB = data + tile_off(tc_all, tr_t);
                            for (idx_t i = 0; i < TILE_SIZE; ++i) {
                                cplx_t tmp = tA[i]; tA[i] = tB[i]; tB[i] = tmp;
                            }
                        }

                        /* Pair K: tile(tr_c2, tc_c1c2) <-> tile(tr_c2, tc_all) */
                        if (tr_c2 > tc_all) {
                            cplx_t * restrict tA = data + tile_off(tr_c2, tc_c1c2);
                            cplx_t * restrict tB = data + tile_off(tr_c2, tc_all);
                            for (idx_t i = 0; i < TILE_SIZE; ++i) {
                                cplx_t tmp = tA[i]; tA[i] = tB[i]; tB[i] = tmp;
                            }
                        }
                        else if (tr_c2 > tc_c1c2) {
                            cplx_t * restrict tA = data + tile_off(tr_c2, tc_c1c2);
                            cplx_t * restrict tB = data + tile_off(tc_all, tr_c2);
                            for (idx_t lr = 0; lr < TILE_DIM; ++lr) {
                                for (idx_t lc = 0; lc < TILE_DIM; ++lc) {
                                    cplx_t tmp = tA[lr * TILE_DIM + lc];
                                    tA[lr * TILE_DIM + lc] = conj(tB[lc * TILE_DIM + lr]);
                                    tB[lc * TILE_DIM + lr] = conj(tmp);
                                }
                            }
                        }
                        else {
                            cplx_t * restrict tA = data + tile_off(tc_c1c2, tr_c2);
                            cplx_t * restrict tB = data + tile_off(tc_all, tr_c2);
                            for (idx_t i = 0; i < TILE_SIZE; ++i) {
                                cplx_t tmp = tA[i]; tA[i] = tB[i]; tB[i] = tmp;
                            }
                        }

                        /* Pair L: tile(tr_c2t, tc_c1c2) <-> tile(tr_c2t, tc_all) */
                        if (tr_c2t > tc_all) {
                            cplx_t * restrict tA = data + tile_off(tr_c2t, tc_c1c2);
                            cplx_t * restrict tB = data + tile_off(tr_c2t, tc_all);
                            for (idx_t i = 0; i < TILE_SIZE; ++i) {
                                cplx_t tmp = tA[i]; tA[i] = tB[i]; tB[i] = tmp;
                            }
                        }
                        else if (tr_c2t > tc_c1c2) {
                            cplx_t * restrict tA = data + tile_off(tr_c2t, tc_c1c2);
                            cplx_t * restrict tB = data + tile_off(tc_all, tr_c2t);
                            for (idx_t lr = 0; lr < TILE_DIM; ++lr) {
                                for (idx_t lc = 0; lc < TILE_DIM; ++lc) {
                                    cplx_t tmp = tA[lr * TILE_DIM + lc];
                                    tA[lr * TILE_DIM + lc] = conj(tB[lc * TILE_DIM + lr]);
                                    tB[lc * TILE_DIM + lr] = conj(tmp);
                                }
                            }
                        }
                        else {
                            cplx_t * restrict tA = data + tile_off(tc_c1c2, tr_c2t);
                            cplx_t * restrict tB = data + tile_off(tc_all, tr_c2t);
                            for (idx_t i = 0; i < TILE_SIZE; ++i) {
                                cplx_t tmp = tA[i]; tA[i] = tB[i]; tB[i] = tmp;
                            }
                        }

                        /* Pair M: tile(tr_c1, tc_c1c2) <-> tile(tr_c1, tc_all) */
                        if (tr_c1 > tc_all) {
                            cplx_t * restrict tA = data + tile_off(tr_c1, tc_c1c2);
                            cplx_t * restrict tB = data + tile_off(tr_c1, tc_all);
                            for (idx_t i = 0; i < TILE_SIZE; ++i) {
                                cplx_t tmp = tA[i]; tA[i] = tB[i]; tB[i] = tmp;
                            }
                        }
                        else if (tr_c1 > tc_c1c2) {
                            cplx_t * restrict tA = data + tile_off(tr_c1, tc_c1c2);
                            cplx_t * restrict tB = data + tile_off(tc_all, tr_c1);
                            for (idx_t lr = 0; lr < TILE_DIM; ++lr) {
                                for (idx_t lc = 0; lc < TILE_DIM; ++lc) {
                                    cplx_t tmp = tA[lr * TILE_DIM + lc];
                                    tA[lr * TILE_DIM + lc] = conj(tB[lc * TILE_DIM + lr]);
                                    tB[lc * TILE_DIM + lr] = conj(tmp);
                                }
                            }
                        }
                        else {
                            cplx_t * restrict tA = data + tile_off(tc_c1c2, tr_c1);
                            cplx_t * restrict tB = data + tile_off(tc_all, tr_c1);
                            for (idx_t i = 0; i < TILE_SIZE; ++i) {
                                cplx_t tmp = tA[i]; tA[i] = tB[i]; tB[i] = tmp;
                            }
                        }

                        /* Pair N: tile(tr_c1t, tc_c1c2) <-> tile(tr_c1t, tc_all) */
                        if (tr_c1t > tc_all) {
                            cplx_t * restrict tA = data + tile_off(tr_c1t, tc_c1c2);
                            cplx_t * restrict tB = data + tile_off(tr_c1t, tc_all);
                            for (idx_t i = 0; i < TILE_SIZE; ++i) {
                                cplx_t tmp = tA[i]; tA[i] = tB[i]; tB[i] = tmp;
                            }
                        }
                        else if (tr_c1t > tc_c1c2) {
                            cplx_t * restrict tA = data + tile_off(tr_c1t, tc_c1c2);
                            cplx_t * restrict tB = data + tile_off(tc_all, tr_c1t);
                            for (idx_t lr = 0; lr < TILE_DIM; ++lr) {
                                for (idx_t lc = 0; lc < TILE_DIM; ++lc) {
                                    cplx_t tmp = tA[lr * TILE_DIM + lc];
                                    tA[lr * TILE_DIM + lc] = conj(tB[lc * TILE_DIM + lr]);
                                    tB[lc * TILE_DIM + lr] = conj(tmp);
                                }
                            }
                        }
                        else {
                            cplx_t * restrict tA = data + tile_off(tc_c1c2, tr_c1t);
                            cplx_t * restrict tB = data + tile_off(tc_all, tr_c1t);
                            for (idx_t i = 0; i < TILE_SIZE; ++i) {
                                cplx_t tmp = tA[i]; tA[i] = tB[i]; tB[i] = tmp;
                            }
                        }
                    }

                    /*
                     * =====================================================================
                     * Diagonal case: btr == btc.
                     *
                     * When btr == btc: tr000 == tc000, so tr_XX == tc_XX for all XX.
                     *
                     * Pairs A-G (row swaps) and I-N (column swaps) are conjugate pairs:
                     *   A <-> I, B <-> K, C <-> M, E <-> J, F <-> L, G <-> N
                     * Performing both would double-swap (cancel). Pairs A-C are already
                     * handled unconditionally above and are always lower tri. On the
                     * diagonal, I is the adjoint of A, K of B, M of C -- these three
                     * are automatically correct. But E-G may involve upper-tri tiles
                     * on diagonal, and their adjoints J-L-N would handle those same
                     * positions. We process whichever side has the lower-tri tile.
                     *
                     * Since the permutation is self-inverse, and on diagonal the tile
                     * at (c1c2, XX) and (XX, c1c2) are adjoints, we must pick ONE
                     * of {E,J}, {F,L}, {G,N} to avoid double-swap.
                     *
                     * Pairs E, F, G each swap tile(tr_c1c2, tc_XX) <-> tile(tr_all, tc_XX).
                     * On diagonal: tr_all = tc_all, tr_c1c2 = tc_c1c2.
                     *   tile(tr_all, tc_XX) is always lower tri (tr_all has all 3 bits).
                     *   tile(tr_c1c2, tc_XX): lower iff tr_c1c2 > tc_XX.
                     *     If lower: direct whole-tile swap.
                     *     If upper: tile(tr_c1c2, tc_XX) = adj(tile(tc_XX, tc_c1c2)).
                     *       But tile(tc_XX, tc_c1c2) is the I/J/K/L/M/N counterpart
                     *       which we would skip. So we use the adjoint access pattern.
                     *
                     * Pair D is already handled above.
                     *
                     * Pair H: tile(tr_all, tc_c1c2) <-> tile(tr_c1c2, tc_all).
                     * On diagonal: tile(tr_all, tr_c1c2) <-> tile(tr_c1c2, tr_all).
                     * These are adjoint tiles. Self-adjoint swap in single tile.
                     * =====================================================================
                     */
                    else {

                        /*
                         * Pairs E, F, G on diagonal.
                         * tile(tr_all, tc_XX) is always lower tri.
                         * tile(tr_c1c2, tc_XX): lower iff tr_c1c2 > tc_XX.
                         *
                         * Note: tr_c1c2 != tc_XX is guaranteed on diagonal because
                         * the insertBits3_0 structure forces distinct 3-bit patterns
                         * at the marker positions (c1c2 = {1,1,0} vs t = {0,0,1},
                         * c2t = {0,1,1}, c1t = {1,0,1}). So strict > or < always.
                         *
                         *   If tr_c1c2 > tc_XX: both lower tri, direct swap.
                         *   If tr_c1c2 < tc_XX: tile(tr_c1c2, tc_XX) is upper tri.
                         *     Use adjoint tile(tc_XX, tr_c1c2).
                         */

                        /* Pair E: tile(tr_c1c2, tc_t) <-> tile(tr_all, tc_t) */
                        if (tr_c1c2 > tc_t) {
                            cplx_t * restrict tA = data + tile_off(tr_c1c2, tc_t);
                            cplx_t * restrict tB = data + tile_off(tr_all, tc_t);
                            for (idx_t i = 0; i < TILE_SIZE; ++i) {
                                cplx_t tmp = tA[i]; tA[i] = tB[i]; tB[i] = tmp;
                            }
                        }
                        else {
                            /* Upper tri: use adjoint tile(tc_t, tr_c1c2) */
                            cplx_t * restrict tA = data + tile_off(tc_t, tr_c1c2);
                            cplx_t * restrict tB = data + tile_off(tr_all, tc_t);
                            for (idx_t lr = 0; lr < TILE_DIM; ++lr) {
                                for (idx_t lc = 0; lc < TILE_DIM; ++lc) {
                                    cplx_t tmp = tB[lr * TILE_DIM + lc];
                                    tB[lr * TILE_DIM + lc] = conj(tA[lc * TILE_DIM + lr]);
                                    tA[lc * TILE_DIM + lr] = conj(tmp);
                                }
                            }
                        }

                        /* Pair F: tile(tr_c1c2, tc_c2t) <-> tile(tr_all, tc_c2t) */
                        if (tr_c1c2 > tc_c2t) {
                            cplx_t * restrict tA = data + tile_off(tr_c1c2, tc_c2t);
                            cplx_t * restrict tB = data + tile_off(tr_all, tc_c2t);
                            for (idx_t i = 0; i < TILE_SIZE; ++i) {
                                cplx_t tmp = tA[i]; tA[i] = tB[i]; tB[i] = tmp;
                            }
                        }
                        else {
                            cplx_t * restrict tA = data + tile_off(tc_c2t, tr_c1c2);
                            cplx_t * restrict tB = data + tile_off(tr_all, tc_c2t);
                            for (idx_t lr = 0; lr < TILE_DIM; ++lr) {
                                for (idx_t lc = 0; lc < TILE_DIM; ++lc) {
                                    cplx_t tmp = tB[lr * TILE_DIM + lc];
                                    tB[lr * TILE_DIM + lc] = conj(tA[lc * TILE_DIM + lr]);
                                    tA[lc * TILE_DIM + lr] = conj(tmp);
                                }
                            }
                        }

                        /* Pair G: tile(tr_c1c2, tc_c1t) <-> tile(tr_all, tc_c1t) */
                        if (tr_c1c2 > tc_c1t) {
                            cplx_t * restrict tA = data + tile_off(tr_c1c2, tc_c1t);
                            cplx_t * restrict tB = data + tile_off(tr_all, tc_c1t);
                            for (idx_t i = 0; i < TILE_SIZE; ++i) {
                                cplx_t tmp = tA[i]; tA[i] = tB[i]; tB[i] = tmp;
                            }
                        }
                        else {
                            cplx_t * restrict tA = data + tile_off(tc_c1t, tr_c1c2);
                            cplx_t * restrict tB = data + tile_off(tr_all, tc_c1t);
                            for (idx_t lr = 0; lr < TILE_DIM; ++lr) {
                                for (idx_t lc = 0; lc < TILE_DIM; ++lc) {
                                    cplx_t tmp = tB[lr * TILE_DIM + lc];
                                    tB[lr * TILE_DIM + lc] = conj(tA[lc * TILE_DIM + lr]);
                                    tA[lc * TILE_DIM + lr] = conj(tmp);
                                }
                            }
                        }

                        /* Pair H: self-adjoint swap within tile(tr_all, tr_c1c2).
                         * Since tr_all = tr_c1c2 | incr_t, we have tr_all > tr_c1c2,
                         * so tile(tr_all, tr_c1c2) is the stored lower-tri tile.
                         * tile(tr_c1c2, tr_all) is its adjoint: element [lr,lc] of
                         * the adjoint = conj(tile(tr_all, tr_c1c2)[lc,lr]).
                         * The swap becomes a self-conjugation within a single tile:
                         *   T[lr,lc] <-> conj(T[lc,lr])
                         * Only process lr >= lc to avoid double-swapping. */
                        cplx_t * restrict tH = data + tile_off(tr_all, tr_c1c2);
                        for (idx_t lr = 0; lr < TILE_DIM; ++lr) {
                            for (idx_t lc = 0; lc <= lr; ++lc) {
                                cplx_t tmp = tH[lr * TILE_DIM + lc];
                                tH[lr * TILE_DIM + lc] = conj(tH[lc * TILE_DIM + lr]);
                                tH[lc * TILE_DIM + lr] = conj(tmp);
                            }
                        }

                        /* Pairs I-N: SKIPPED on diagonal -- they are the Hermitian
                         * conjugates of Pairs A,E,B,F,C,G respectively, and
                         * performing both would double-swap (cancel). */
                    }
                }
            }
        }
    }

}
