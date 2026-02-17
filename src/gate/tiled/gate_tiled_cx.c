/**
 * @file gate_tiled_cx.c
 * @brief CX (CNOT) gate for tiled density matrices
 *
 * CX is a zero-flop permutation gate (self-adjoint): ρ' = CX·ρ·CX.
 * perm(x) = x ^ (1<<target) when bit control is set, else x.
 * ρ'[r,c] = ρ[perm(r), perm(c)].  Involution: each element is fixed or in a 2-cycle.
 *
 * For each block (br, bc), row indices r00,r01,r10,r11 and col indices c00,c01,c10,c11
 * (where bit pattern ab means control=a, target=b), the 6 swap pairs are:
 *   Pair A: (r10,c00) <-> (r11,c00)   -- always lower-tri
 *   Pair E: (r10,c10) <-> (r11,c11)   -- always lower-tri
 *   Pair B: (r11,c01) <-> (r10,c01)   -- may cross diagonal
 *   Pair F: (r11,c10) <-> (r10,c11)   -- may cross diagonal
 *   Pair C: (r00,c10) <-> (r00,c11)   -- may cross diagonal
 *   Pair D: (r01,c10) <-> (r01,c11)   -- may cross diagonal
 *
 * Three cases based on where lo=min(ctrl,tgt) and hi=max(ctrl,tgt) fall
 * relative to LOG_TILE_DIM, following the same structure as swap_tiled().
 */

#include "gate_tiled.h"

void cx_tiled(state_t *state, const qubit_t control, const qubit_t target) {
    const gate_idx_t dim = (gate_idx_t)1 << state->qubits;
    const gate_idx_t dim_tile = (dim + TILE_DIM - 1) >> LOG_TILE_DIM;
    cplx_t * restrict data = state->data;

    if (control < LOG_TILE_DIM) {
        /* Case 1:
         * Control and target smaller than LOG_TILE_DIM
         *
         * All swap pairs are contained within the same tile:
         *  |ctrl=1,tgt=0><ctrl=0,tgt=0| <--> |ctrl=1,tgt=1><ctrl=0,tgt=0|
         *  |ctrl=1,tgt=0><ctrl=1,tgt=0| <--> |ctrl=1,tgt=1><ctrl=1,tgt=1|
         *  |ctrl=1,tgt=0><ctrl=0,tgt=1| <--> |ctrl=1,tgt=1><ctrl=0,tgt=1|
         *  |ctrl=1,tgt=0><ctrl=1,tgt=1| <--> |ctrl=1,tgt=1><ctrl=1,tgt=0|
         *  |ctrl=0,tgt=1><ctrl=1,tgt=0| <--> |ctrl=0,tgt=1><ctrl=1,tgt=1|
         *  |ctrl=0,tgt=0><ctrl=1,tgt=0| <--> |ctrl=0,tgt=0><ctrl=1,tgt=1|
         */
        if (target < LOG_TILE_DIM) {
            const gate_idx_t n_tiles = dim_tile * (dim_tile + 1) / 2;
            const gate_idx_t n_base = (dim < TILE_DIM) ? dim >> 2 : TILE_DIM >> 2;
            const qubit_t lo = (control < target) ? control : target;
            const qubit_t hi = (control > target) ? control : target;
            const gate_idx_t incr_ctrl = (gate_idx_t)1 << control;
            const gate_idx_t incr_tgt = (gate_idx_t)1 << target;

            #pragma omp parallel for schedule(static) if(dim >= OMP_THRESHOLD)
            for (gate_idx_t i = 0; i < n_tiles; ++i) {
                cplx_t * restrict tile = data + i * TILE_SIZE;

                for (gate_idx_t blr = 0; blr < n_base; ++blr) {
                    const gate_idx_t r00 = insertBits2_0(blr, lo, hi);
                    const gate_idx_t r01 = r00 | incr_tgt, r10 = r00 | incr_ctrl, r11 = r10 | incr_tgt;

                    for (gate_idx_t blc = 0; blc < n_base; ++blc) {
                        const gate_idx_t c00 = insertBits2_0(blc, lo, hi);
                        const gate_idx_t c01 = c00 | incr_tgt, c10 = c00 | incr_ctrl, c11 = c10 | incr_tgt;

                        cplx_t tmp = tile[elem_off(r10, c00)];
                        tile[elem_off(r10, c00)] = tile[elem_off(r11, c00)];
                        tile[elem_off(r11, c00)] = tmp;

                        tmp = tile[elem_off(r10, c10)];
                        tile[elem_off(r10, c10)] = tile[elem_off(r11, c11)];
                        tile[elem_off(r11, c11)] = tmp;

                        tmp = tile[elem_off(r11, c01)];
                        tile[elem_off(r11, c01)] = tile[elem_off(r10, c01)];
                        tile[elem_off(r10, c01)] = tmp;

                        tmp = tile[elem_off(r11, c10)];
                        tile[elem_off(r11, c10)] = tile[elem_off(r10, c11)];
                        tile[elem_off(r10, c11)] = tmp;

                        tmp = tile[elem_off(r01, c10)];
                        tile[elem_off(r01, c10)] = tile[elem_off(r01, c11)];
                        tile[elem_off(r01, c11)] = tmp;

                        tmp = tile[elem_off(r00, c10)];
                        tile[elem_off(r00, c10)] = tile[elem_off(r00, c11)];
                        tile[elem_off(r00, c11)] = tmp;
                    }
                }
            }
        }

        /* Case 2:
         * Target is at least LOG_TILE_DIM & control is smaller than LOG_TILE_DIM
         *
         * All swap pairs are spread across four tiles t00=(tgt=0,tgt=0), t10=(tgt=1,tgt=0), t11=(tgt=1,tgt=1),
         * t01=(tgt=0,tgt=0) and each swap is an inter-tile swap.
         *  |ctrl=1,tgt=0><ctrl=0,tgt=0| <--> |ctrl=1,tgt=1><ctrl=0,tgt=0|  -- always lower tri
         *  |ctrl=1,tgt=0><ctrl=1,tgt=0| <--> |ctrl=1,tgt=1><ctrl=1,tgt=1|  -- always lower tri
         *  |ctrl=1,tgt=0><ctrl=0,tgt=1| <--> |ctrl=1,tgt=1><ctrl=0,tgt=1|  -- may cross diagonal
         *  |ctrl=1,tgt=0><ctrl=1,tgt=1| <--> |ctrl=1,tgt=1><ctrl=1,tgt=0|  -- may cross diagonal
         *  |ctrl=0,tgt=1><ctrl=1,tgt=0| <--> |ctrl=0,tgt=1><ctrl=1,tgt=1|  -- always lower tri
         *  |ctrl=0,tgt=0><ctrl=1,tgt=0| <--> |ctrl=0,tgt=0><ctrl=1,tgt=1|  -- may cross diagonal
         */
        else {
            const gate_idx_t n_bt = dim_tile >> 1;
            const qubit_t tile_bit = target - LOG_TILE_DIM;
            const gate_idx_t incr_tile = (gate_idx_t)1 << tile_bit;
            const gate_idx_t n_bl = TILE_DIM >> 1;
            const gate_idx_t incr_local = (gate_idx_t)1 << control;

            #pragma omp parallel for schedule(static, 1) if(dim >= OMP_THRESHOLD)
            for (gate_idx_t btr = 0; btr < n_bt; ++btr) {
                const gate_idx_t tr0 = insertBit0(btr, tile_bit);
                const gate_idx_t tr1 = tr0 | incr_tile;

                for (gate_idx_t btc = 0; btc <= btr; ++btc) {
                    const gate_idx_t tc0 = insertBit0(btc, tile_bit);
                    const gate_idx_t tc1 = tc0 | incr_tile;

                    cplx_t * restrict t00 = data + tile_off(tr0, tc0);
                    cplx_t * restrict t10 = data + tile_off(tr1, tc0);
                    cplx_t * restrict t11 = data + tile_off(tr1, tc1);

                    /* Fast path: All tiles are in the lower triangle */
                    if (tr0 > tc1) {
                        cplx_t * restrict t01 = data + tile_off(tr0, tc1);

                        for (gate_idx_t blr = 0; blr < n_bl; ++blr) {
                            const gate_idx_t r0 = insertBit0(blr, control);
                            const gate_idx_t r1 = r0 | incr_local;

                            for (gate_idx_t blc = 0; blc < n_bl; ++blc) {
                                const gate_idx_t c0 = insertBit0(blc, control);
                                const gate_idx_t c1 = c0 | incr_local;

                                /* Pair A: |ctrl=1,tgt=0><ctrl=0,tgt=0| <--> |ctrl=1,tgt=1><ctrl=0,tgt=0| */
                                cplx_t tmp = t00[elem_off(r1, c0)];
                                t00[elem_off(r1, c0)] = t10[elem_off(r1, c0)];
                                t10[elem_off(r1, c0)] = tmp;

                                /* Pair B: |ctrl=1,tgt=0><ctrl=1,tgt=0| <--> |ctrl=1,tgt=1><ctrl=1,tgt=1| */
                                tmp = t00[elem_off(r1, c1)];
                                t00[elem_off(r1, c1)] = t11[elem_off(r1, c1)];
                                t11[elem_off(r1, c1)] = tmp;

                                /* Pair E: |ctrl=0,tgt=1><ctrl=1,tgt=0| <--> |ctrl=0,tgt=1><ctrl=1,tgt=1| */
                                tmp = t10[elem_off(r0, c1)];
                                t10[elem_off(r0, c1)] = t11[elem_off(r0, c1)];
                                t11[elem_off(r0, c1)] = tmp;

                                /* Pair C: |ctrl=1,tgt=0><ctrl=0,tgt=1| <--> |ctrl=1,tgt=1><ctrl=0,tgt=1| */
                                tmp = t11[elem_off(r1, c0)];
                                t11[elem_off(r1, c0)] = t01[elem_off(r1, c0)];
                                t01[elem_off(r1, c0)] = tmp;

                                /* Pair D: |ctrl=1,tgt=0><ctrl=1,tgt=1| <--> |ctrl=1,tgt=1><ctrl=1,tgt=0| */
                                tmp = t10[elem_off(r1, c1)];
                                t10[elem_off(r1, c1)] = t01[elem_off(r1, c1)];
                                t01[elem_off(r1, c1)] = tmp;

                                /* Pair F: |ctrl=0,tgt=0><ctrl=1,tgt=0| <--> |ctrl=0,tgt=0><ctrl=1,tgt=1| */
                                tmp = t00[elem_off(r0, c1)];
                                t00[elem_off(r0, c1)] = t01[elem_off(r0, c1)];
                                t01[elem_off(r0, c1)] = tmp;
                            }
                        }
                    }

                    /* Original tile t01 is in the upper triangle but its adjoint is contained in the lower triangle */
                    else if (btr != btc) {
                        cplx_t * restrict t01 = data + tile_off(tc1, tr0);

                        for (gate_idx_t blr = 0; blr < n_bl; ++blr) {
                            const gate_idx_t r0 = insertBit0(blr, control);
                            const gate_idx_t r1 = r0 | incr_local;

                            for (gate_idx_t blc = 0; blc < n_bl; ++blc) {
                                const gate_idx_t c0 = insertBit0(blc, control);
                                const gate_idx_t c1 = c0 | incr_local;

                                /* Pair A: |ctrl=1,tgt=0><ctrl=0,tgt=0| <--> |ctrl=1,tgt=1><ctrl=0,tgt=0| */
                                cplx_t tmp = t00[elem_off(r1, c0)];
                                t00[elem_off(r1, c0)] = t10[elem_off(r1, c0)];
                                t10[elem_off(r1, c0)] = tmp;

                                /* Pair B: |ctrl=1,tgt=0><ctrl=1,tgt=0| <--> |ctrl=1,tgt=1><ctrl=1,tgt=1| */
                                tmp = t00[elem_off(r1, c1)];
                                t00[elem_off(r1, c1)] = t11[elem_off(r1, c1)];
                                t11[elem_off(r1, c1)] = tmp;

                                /* Pair E: |ctrl=0,tgt=1><ctrl=1,tgt=0| <--> |ctrl=0,tgt=1><ctrl=1,tgt=1| */
                                tmp = t10[elem_off(r0, c1)];
                                t10[elem_off(r0, c1)] = t11[elem_off(r0, c1)];
                                t11[elem_off(r0, c1)] = tmp;

                                /* Pair C: |ctrl=1,tgt=0><ctrl=0,tgt=1| <--> |ctrl=1,tgt=1><ctrl=0,tgt=1| */
                                tmp = t11[elem_off(r1, c0)];
                                t11[elem_off(r1, c0)] = conj(t01[elem_off(c0, r1)]);
                                t01[elem_off(c0, r1)] = conj(tmp);

                                /* Pair D: |ctrl=1,tgt=0><ctrl=1,tgt=1| <--> |ctrl=1,tgt=1><ctrl=1,tgt=0| */
                                tmp = t10[elem_off(r1, c1)];
                                t10[elem_off(r1, c1)] = conj(t01[elem_off(c1, r1)]);
                                t01[elem_off(c1, r1)] = conj(tmp);

                                /* Pair F: |ctrl=0,tgt=0><ctrl=1,tgt=0| <--> |ctrl=0,tgt=0><ctrl=1,tgt=1| */
                                tmp = t00[elem_off(r0, c1)];
                                t00[elem_off(r0, c1)] = conj(t01[elem_off(c1, r0)]);
                                t01[elem_off(c1, r0)] = conj(tmp);
                            }
                        }
                    }

                    /* Original tile t01 is in the upper triangle and its adjoint is t10 */
                    else {
                        /* One-way writes: Elements in the lower triangle receive values t01 (from its adjoint t10), but
                         * do not overwrite to avoid double write
                         */
                        for (gate_idx_t blr = 0; blr < n_bl; ++blr) {
                            const gate_idx_t r0 = insertBit0(blr, control);
                            const gate_idx_t r1 = r0 | incr_local;

                            for (gate_idx_t blc = 0; blc < n_bl; ++blc) {
                                const gate_idx_t c0 = insertBit0(blc, control);
                                const gate_idx_t c1 = c0 | incr_local;

                                /* Pair C: |ctrl=1,tgt=0><ctrl=0,tgt=1| <--> |ctrl=1,tgt=1><ctrl=0,tgt=1| */
                                t11[elem_off(r1, c0)] = conj(t10[elem_off(c0, r1)]);

                                /* Pair F: |ctrl=0,tgt=0><ctrl=1,tgt=0| <--> |ctrl=0,tgt=0><ctrl=1,tgt=1| */
                                t00[elem_off(r0, c1)] = conj(t10[elem_off(c1, r0)]);
                            }
                        }

                        for (gate_idx_t blr = 0; blr < n_bl; ++blr) {
                            const gate_idx_t r0 = insertBit0(blr, control);
                            const gate_idx_t r1 = r0 | incr_local;

                            for (gate_idx_t blc = 0; blc < n_bl; ++blc) {
                                const gate_idx_t c0 = insertBit0(blc, control);
                                const gate_idx_t c1 = c0 | incr_local;

                                /* Pair A: |ctrl=1,tgt=0><ctrl=0,tgt=0| <--> |ctrl=1,tgt=1><ctrl=0,tgt=0| */
                                cplx_t tmp = t00[elem_off(r1, c0)];
                                t00[elem_off(r1, c0)] = t10[elem_off(r1, c0)];
                                t10[elem_off(r1, c0)] = tmp;

                                /* Pair B: |ctrl=1,tgt=0><ctrl=1,tgt=0| <--> |ctrl=1,tgt=1><ctrl=1,tgt=1| */
                                tmp = t00[elem_off(r1, c1)];
                                t00[elem_off(r1, c1)] = t11[elem_off(r1, c1)];
                                t11[elem_off(r1, c1)] = tmp;

                                /* Pair E: |ctrl=0,tgt=1><ctrl=1,tgt=0| <--> |ctrl=0,tgt=1><ctrl=1,tgt=1| */
                                tmp = t10[elem_off(r0, c1)];
                                t10[elem_off(r0, c1)] = t11[elem_off(r0, c1)];
                                t11[elem_off(r0, c1)] = tmp;

                                /* Pair D: |ctrl=1,tgt=0><ctrl=1,tgt=1| <--> |ctrl=1,tgt=1><ctrl=1,tgt=0| */
                                if (blc <= blr) {
                                    tmp = t10[elem_off(r1, c1)];
                                    t10[elem_off(r1, c1)] = conj(t10[elem_off(c1, r1)]);
                                    t10[elem_off(c1, r1)] = conj(tmp);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    else {
        /* Case 3:
         * Control is at least LOG_TILE_DIM & target is smaller than LOG_TILE_DIM
         *
         * All swap pairs are spread across three tiles t10=(ctrl=1,ctrl=0), t11=(ctrl=1,ctrl=1), t01=(ctrl=0,ctrl=1),
         * but each swap pair is contained in one tile.
         *  |ctrl=1,tgt=0><ctrl=0,tgt=0| <--> |ctrl=1,tgt=1><ctrl=0,tgt=0|  -- always lower tri
         *  |ctrl=1,tgt=0><ctrl=1,tgt=0| <--> |ctrl=1,tgt=1><ctrl=1,tgt=1|  -- always lower tri
         *  |ctrl=1,tgt=0><ctrl=0,tgt=1| <--> |ctrl=1,tgt=1><ctrl=0,tgt=1|  -- always lower tri
         *  |ctrl=1,tgt=0><ctrl=1,tgt=1| <--> |ctrl=1,tgt=1><ctrl=1,tgt=0|  -- always lower tri
         *  |ctrl=0,tgt=1><ctrl=1,tgt=0| <--> |ctrl=0,tgt=1><ctrl=1,tgt=1|  -- may be in upper triangle
         *  |ctrl=0,tgt=0><ctrl=1,tgt=0| <--> |ctrl=0,tgt=0><ctrl=1,tgt=1|  -- may be in upper triangle
         */
        if (target < LOG_TILE_DIM) {
            const gate_idx_t n_bt = dim_tile >> 1;
            const gate_idx_t tile_bit = control - LOG_TILE_DIM;
            const gate_idx_t incr_tile = (gate_idx_t)1 << tile_bit;
            const gate_idx_t n_bl = TILE_DIM >> 1;
            const gate_idx_t incr_local = (gate_idx_t)1 << target;

            #pragma omp parallel for schedule(static, 1) if(dim >= OMP_THRESHOLD)
            for (gate_idx_t btr = 0; btr < n_bt; ++btr) {
                const gate_idx_t tr0 = insertBit0(btr, tile_bit);
                const gate_idx_t tr1 = tr0 | incr_tile;

                for (gate_idx_t btc = 0; btc <= btr; ++btc) {
                    const gate_idx_t tc0 = insertBit0(btc, tile_bit);
                    const gate_idx_t tc1 = tc0 | incr_tile;

                    cplx_t * restrict t10 = data + tile_off(tr1, tc0);
                    cplx_t * restrict t11 = data + tile_off(tr1, tc1);

                    if (tr0 > tc1) {
                      cplx_t * restrict t01 = data + tile_off(tr0, tc1);

                        for (gate_idx_t blr = 0; blr < n_bl; ++blr) {
                            const gate_idx_t lr0 = insertBit0(blr, target);
                            const gate_idx_t lr1 = lr0 | incr_local;

                            for (gate_idx_t blc = 0; blc < n_bl; ++blc) {
                                const gate_idx_t lc0 = insertBit0(blc, target);
                                const gate_idx_t lc1 = lc0 | incr_local;

                                /* Pair A: |ctrl=1,tgt=0><ctrl=0,tgt=0| <--> |ctrl=1,tgt=1><ctrl=0,tgt=0| */
                                cplx_t tmp = t10[elem_off(lr0, lc0)];
                                t10[elem_off(lr0, lc0)] = t10[elem_off(lr1, lc0)];
                                t10[elem_off(lr1, lc0)] = tmp;

                                /* Pair C: |ctrl=1,tgt=0><ctrl=0,tgt=1| <--> |ctrl=1,tgt=1><ctrl=0,tgt=1| */
                                tmp = t10[elem_off(lr0, lc1)];
                                t10[elem_off(lr0, lc1)] = t10[elem_off(lr1, lc1)];
                                t10[elem_off(lr1, lc1)] = tmp;

                                /* Pair B: |ctrl=1,tgt=0><ctrl=1,tgt=0| <--> |ctrl=1,tgt=1><ctrl=1,tgt=1| */
                                tmp = t11[elem_off(lr0, lc0)];
                                t11[elem_off(lr0, lc0)] = t11[elem_off(lr1, lc1)];
                                t11[elem_off(lr1, lc1)] = tmp;

                                /* Pair D: |ctrl=1,tgt=0><ctrl=1,tgt=1| <--> |ctrl=1,tgt=1><ctrl=1,tgt=0| */
                                tmp = t11[elem_off(lr0, lc1)];
                                t11[elem_off(lr0, lc1)] = t11[elem_off(lr1, lc0)];
                                t11[elem_off(lr1, lc0)] = tmp;

                                /* Pair E: |ctrl=0,tgt=1><ctrl=1,tgt=0| <--> |ctrl=0,tgt=1><ctrl=1,tgt=1| */
                                tmp = t01[elem_off(lr1, lc0)];
                                t01[elem_off(lr1, lc0)] = t01[elem_off(lr1, lc1)];
                                t01[elem_off(lr1, lc1)] = tmp;

                                /* Pair F: |ctrl=0,tgt=0><ctrl=1,tgt=0| <--> |ctrl=0,tgt=0><ctrl=1,tgt=1| */
                                tmp = t01[elem_off(lr0, lc0)];
                                t01[elem_off(lr0, lc0)] = t01[elem_off(lr0, lc1)];
                                t01[elem_off(lr0, lc1)] = tmp;
                            }
                        }
                    }

                    else if (btr != btc) {
                        cplx_t * restrict t01 = data + tile_off(tc1, tr0);   

                        for (gate_idx_t blr = 0; blr < n_bl; ++blr) {
                            const gate_idx_t lr0 = insertBit0(blr, target);
                            const gate_idx_t lr1 = lr0 | incr_local;

                            for (gate_idx_t blc = 0; blc < n_bl; ++blc) {
                                const gate_idx_t lc0 = insertBit0(blc, target);
                                const gate_idx_t lc1 = lc0 | incr_local;

                                /* Pair A: |ctrl=1,tgt=0><ctrl=0,tgt=0| <--> |ctrl=1,tgt=1><ctrl=0,tgt=0| */
                                cplx_t tmp = t10[elem_off(lr0, lc0)];
                                t10[elem_off(lr0, lc0)] = t10[elem_off(lr1, lc0)];
                                t10[elem_off(lr1, lc0)] = tmp;

                                /* Pair C: |ctrl=1,tgt=0><ctrl=0,tgt=1| <--> |ctrl=1,tgt=1><ctrl=0,tgt=1| */
                                tmp = t10[elem_off(lr0, lc1)];
                                t10[elem_off(lr0, lc1)] = t10[elem_off(lr1, lc1)];
                                t10[elem_off(lr1, lc1)] = tmp;

                                /* Pair B: |ctrl=1,tgt=0><ctrl=1,tgt=0| <--> |ctrl=1,tgt=1><ctrl=1,tgt=1| */
                                tmp = t11[elem_off(lr0, lc0)];
                                t11[elem_off(lr0, lc0)] = t11[elem_off(lr1, lc1)];
                                t11[elem_off(lr1, lc1)] = tmp;

                                /* Pair D: |ctrl=1,tgt=0><ctrl=1,tgt=1| <--> |ctrl=1,tgt=1><ctrl=1,tgt=0| */
                                tmp = t11[elem_off(lr0, lc1)];
                                t11[elem_off(lr0, lc1)] = t11[elem_off(lr1, lc0)];
                                t11[elem_off(lr1, lc0)] = tmp;

                                /* Pair E: |ctrl=0,tgt=1><ctrl=1,tgt=0| <--> |ctrl=0,tgt=1><ctrl=1,tgt=1| */
                                tmp = t01[elem_off(lr1, lc0)];
                                t01[elem_off(lr1, lc0)] = t01[elem_off(lr1, lc1)];
                                t01[elem_off(lr1, lc1)] = tmp;

                                /* Pair F: |ctrl=0,tgt=0><ctrl=1,tgt=0| <--> |ctrl=0,tgt=0><ctrl=1,tgt=1| */
                                tmp = t01[elem_off(lr0, lc0)];
                                t01[elem_off(lr0, lc0)] = t01[elem_off(lr0, lc1)];
                                t01[elem_off(lr0, lc1)] = tmp;
                            }
                        }
                    }

                    else {
                        for (gate_idx_t blr = 0; blr < n_bl; ++blr) {
                            const gate_idx_t lr0 = insertBit0(blr, target);
                            const gate_idx_t lr1 = lr0 | incr_local;

                            for (gate_idx_t blc = 0; blc < n_bl; ++blc) {
                                const gate_idx_t lc0 = insertBit0(blc, target);
                                const gate_idx_t lc1 = lc0 | incr_local;

                                /* Pair A: |ctrl=1,tgt=0><ctrl=0,tgt=0| <--> |ctrl=1,tgt=1><ctrl=0,tgt=0| */
                                cplx_t tmp = t10[elem_off(lr0, lc0)];
                                t10[elem_off(lr0, lc0)] = t10[elem_off(lr1, lc0)];
                                t10[elem_off(lr1, lc0)] = tmp;

                                /* Pair C: |ctrl=1,tgt=0><ctrl=0,tgt=1| <--> |ctrl=1,tgt=1><ctrl=0,tgt=1| */
                                tmp = t10[elem_off(lr0, lc1)];
                                t10[elem_off(lr0, lc1)] = t10[elem_off(lr1, lc1)];
                                t10[elem_off(lr1, lc1)] = tmp;

                                /* Pair B: |ctrl=1,tgt=0><ctrl=1,tgt=0| <--> |ctrl=1,tgt=1><ctrl=1,tgt=1| */
                                tmp = t11[elem_off(lr0, lc0)];
                                t11[elem_off(lr0, lc0)] = t11[elem_off(lr1, lc1)];
                                t11[elem_off(lr1, lc1)] = tmp;

                                /* Pair D: |ctrl=1,tgt=0><ctrl=1,tgt=1| <--> |ctrl=1,tgt=1><ctrl=1,tgt=0| */
                                tmp = t11[elem_off(lr0, lc1)];
                                t11[elem_off(lr0, lc1)] = t11[elem_off(lr1, lc0)];
                                t11[elem_off(lr1, lc0)] = tmp;
                            }
                        }
                    }
                }
            }
        }

        /* Case 4:
         * Control and target are at least LOG_TILE_DIM
         *
         * All swap pairs are spread across 12 tiles and swap the entire tiles.
         *  |ctrl=1,tgt=0><ctrl=0,tgt=0| <--> |ctrl=1,tgt=1><ctrl=0,tgt=0|  -- always lower tri
         *  |ctrl=1,tgt=0><ctrl=1,tgt=0| <--> |ctrl=1,tgt=1><ctrl=1,tgt=1|  -- always lower tri
         *  |ctrl=1,tgt=0><ctrl=0,tgt=1| <--> |ctrl=1,tgt=1><ctrl=0,tgt=1|  -- may cross diagonal
         *  |ctrl=1,tgt=0><ctrl=1,tgt=1| <--> |ctrl=1,tgt=1><ctrl=1,tgt=0|  -- may cross diagonal
         *  |ctrl=0,tgt=1><ctrl=1,tgt=0| <--> |ctrl=0,tgt=1><ctrl=1,tgt=1|  -- may cross diagonal
         *  |ctrl=0,tgt=0><ctrl=1,tgt=0| <--> |ctrl=0,tgt=0><ctrl=1,tgt=1|  -- may cross diagonal
         */
        else {
            const gate_idx_t n_bt = dim_tile >> 2;
            const qubit_t lo = (control < target) ? control - LOG_TILE_DIM : target - LOG_TILE_DIM;
            const qubit_t hi = (control > target) ? control - LOG_TILE_DIM : target - LOG_TILE_DIM;
            const gate_idx_t incr_ctrl = (gate_idx_t)1 << (control - LOG_TILE_DIM);
            const gate_idx_t incr_tgt = (gate_idx_t)1 << (target - LOG_TILE_DIM);

            #pragma omp parallel for schedule(static, 1) if(dim >= OMP_THRESHOLD)
            for (gate_idx_t btr = 0; btr < n_bt; ++btr) {
                const gate_idx_t tr00 = insertBits2_0(btr, lo, hi);
                const gate_idx_t tr01 = tr00 | incr_tgt, tr10 = tr00 | incr_ctrl, tr11 = tr10 | incr_tgt;

                for (gate_idx_t btc = 0; btc <= btr; ++btc) {
                    const gate_idx_t tc00 = insertBits2_0(btc, lo, hi);
                    const gate_idx_t tc01 = tc00 | incr_tgt, tc10 = tc00 | incr_ctrl, tc11 = tc10 | incr_tgt;

                    /* Pair A: |ctrl=1,tgt=0><ctrl=0,tgt=0| <--> |ctrl=1,tgt=1><ctrl=0,tgt=0| */
                    {
                        cplx_t * restrict tA = data + tile_off(tr10, tc00);
                        cplx_t * restrict tB = data + tile_off(tr11, tc00);
                        for (gate_idx_t i = 0; i < TILE_SIZE; ++i) {
                            cplx_t tmp = tA[i]; tA[i] = tB[i]; tB[i] = tmp;
                        }
                    }

                    /* Pair B: |ctrl=1,tgt=0><ctrl=1,tgt=0| <--> |ctrl=1,tgt=1><ctrl=1,tgt=1| */
                    {
                        cplx_t * restrict tA = data + tile_off(tr10, tc10);
                        cplx_t * restrict tB = data + tile_off(tr11, tc11);
                        for (gate_idx_t i = 0; i < TILE_SIZE; ++i) {
                            cplx_t tmp = tA[i]; tA[i] = tB[i]; tB[i] = tmp;
                        }
                    }

                    /* Fast path: If tr00 > tc11, all tiles are in the lower triangle */
                    if (tr00 > tc11) {
                        /* Pair C: |ctrl=1,tgt=0><ctrl=0,tgt=1| <--> |ctrl=1,tgt=1><ctrl=0,tgt=1| */
                        {
                            cplx_t * restrict tA = data + tile_off(tr11, tc01);
                            cplx_t * restrict tB = data + tile_off(tr10, tc01);
                            for (gate_idx_t i = 0; i < TILE_SIZE; ++i) {
                                cplx_t tmp = tA[i]; tA[i] = tB[i]; tB[i] = tmp;
                            }
                        }

                        /* Pair D: |ctrl=1,tgt=0><ctrl=1,tgt=1| <--> |ctrl=1,tgt=1><ctrl=1,tgt=0| */
                        {
                            cplx_t * restrict tA = data + tile_off(tr11, tc10);
                            cplx_t * restrict tB = data + tile_off(tr10, tc11);
                            for (gate_idx_t i = 0; i < TILE_SIZE; ++i) {
                                cplx_t tmp = tA[i]; tA[i] = tB[i]; tB[i] = tmp;
                            }
                        }

                        /* Pair E: |ctrl=0,tgt=1><ctrl=1,tgt=0| <--> |ctrl=0,tgt=1><ctrl=1,tgt=1| */
                        {
                            cplx_t * restrict tA = data + tile_off(tr01, tc10);
                            cplx_t * restrict tB = data + tile_off(tr01, tc11);
                            for (gate_idx_t i = 0; i < TILE_SIZE; ++i) {
                                cplx_t tmp = tA[i]; tA[i] = tB[i]; tB[i] = tmp;
                            }
                        }

                        /* Pair F: |ctrl=0,tgt=0><ctrl=1,tgt=0| <--> |ctrl=0,tgt=0><ctrl=1,tgt=1| */
                        {
                            cplx_t * restrict tA = data + tile_off(tr00, tc10);
                            cplx_t * restrict tB = data + tile_off(tr00, tc11);
                            for (gate_idx_t i = 0; i < TILE_SIZE; ++i) {
                                cplx_t tmp = tA[i]; tA[i] = tB[i]; tB[i] = tmp;
                            }
                        }
                    }

                    else {
                        if (tr10 > tc11) {
                            {
                                /* Pair C: |ctrl=1,tgt=0><ctrl=0,tgt=1| <--> |ctrl=1,tgt=1><ctrl=0,tgt=1| */
                                cplx_t * restrict tA = data + tile_off(tr11, tc01);
                                cplx_t * restrict tB = data + tile_off(tr10, tc01);
                                for (gate_idx_t i = 0; i < TILE_SIZE; ++i) {
                                    cplx_t tmp = tA[i]; tA[i] = tB[i]; tB[i] = tmp;
                                }
                            }

                            {
                                /* Pair D: |ctrl=1,tgt=0><ctrl=1,tgt=1| <--> |ctrl=1,tgt=1><ctrl=1,tgt=0| */
                                cplx_t * restrict tA = data + tile_off(tr11, tc10);
                                cplx_t * restrict tB = data + tile_off(tr10, tc11);
                                for (gate_idx_t i = 0; i < TILE_SIZE; ++i) {
                                    cplx_t tmp = tA[i]; tA[i] = tB[i]; tB[i] = tmp;
                                }
                            }
                        }

                        else {
                            {
                                /* Pair C: |ctrl=1,tgt=0><ctrl=0,tgt=1| <--> |ctrl=1,tgt=1><ctrl=0,tgt=1| */
                                cplx_t * restrict tA = data + tile_off(tr11, tc01);
                                cplx_t * restrict tB = data + tile_off(tc01, tr10);
                                for (gate_idx_t lr = 0; lr < TILE_DIM; ++lr) {
                                    for (gate_idx_t lc = 0; lc < TILE_DIM; ++lc) {
                                        cplx_t tmp = tA[lr * TILE_DIM + lc];
                                        tA[lr * TILE_DIM + lc] = conj(tB[lc * TILE_DIM + lr]);
                                        tB[lc * TILE_DIM + lr] = conj(tmp);
                                    }
                                }
                            }

                            if (btr != btc) {
                                /* Pair D: |ctrl=1,tgt=0><ctrl=1,tgt=1| <--> |ctrl=1,tgt=1><ctrl=1,tgt=0| */
                                cplx_t * restrict tA = data + tile_off(tr11, tc10);
                                cplx_t * restrict tB = data + tile_off(tc11, tr10);
                                for (gate_idx_t lr = 0; lr < TILE_DIM; ++lr) {
                                    for (gate_idx_t lc = 0; lc < TILE_DIM; ++lc) {
                                        cplx_t tmp = tA[lr * TILE_DIM + lc];
                                        tA[lr * TILE_DIM + lc] = conj(tB[lc * TILE_DIM + lr]);
                                        tB[lc * TILE_DIM + lr] = conj(tmp);
                                    }
                                }
                            }

                            else {
                                /* Pair D: |ctrl=1,tgt=0><ctrl=1,tgt=1| <--> |ctrl=1,tgt=1><ctrl=1,tgt=0| */
                                cplx_t * restrict tA = data + tile_off(tr11, tc10);
                                cplx_t * restrict tB = data + tile_off(tc11, tr10);
                                for (gate_idx_t lr = 0; lr < TILE_DIM; ++lr) {
                                    for (gate_idx_t lc = 0; lc <= lr; ++lc) {
                                        cplx_t tmp = tA[lr * TILE_DIM + lc];
                                        tA[lr * TILE_DIM + lc] = conj(tB[lc * TILE_DIM + lr]);
                                        tB[lc * TILE_DIM + lr] = conj(tmp);
                                    }
                                }
                            }
                        }

                        if (btr != btc) {
                            if (tr00 > tc10) {
                                /* Pair E: |ctrl=0,tgt=1><ctrl=1,tgt=0| <--> |ctrl=0,tgt=1><ctrl=1,tgt=1| */
                                {
                                    cplx_t * restrict tA = data + tile_off(tr01, tc10);
                                    cplx_t * restrict tB = data + tile_off(tr01, tc11);
                                    for (gate_idx_t i = 0; i < TILE_SIZE; ++i) {
                                        cplx_t tmp = tA[i]; tA[i] = tB[i]; tB[i] = tmp;
                                    }
                                }

                                /* Pair F: |ctrl=0,tgt=0><ctrl=1,tgt=0| <--> |ctrl=0,tgt=0><ctrl=1,tgt=1| */
                                {
                                    cplx_t * restrict tA = data + tile_off(tr00, tc10);
                                    cplx_t * restrict tB = data + tile_off(tc11, tr00);
                                    for (gate_idx_t lr = 0; lr < TILE_DIM; ++lr) {
                                        for (gate_idx_t lc = 0; lc < TILE_DIM; ++lc) {
                                            cplx_t tmp = tA[lr * TILE_DIM + lc];
                                            tA[lr * TILE_DIM + lc] = conj(tB[lc * TILE_DIM + lr]);
                                            tB[lc * TILE_DIM + lr] = conj(tmp);
                                        }
                                    }
                                }
                            }
                            else {
                                /* Pair E: |ctrl=0,tgt=1><ctrl=1,tgt=0| <--> |ctrl=0,tgt=1><ctrl=1,tgt=1| */
                                {
                                    cplx_t * restrict tA = data + tile_off(tc11, tr01);
                                    if (tr01 > tc10) {
                                        cplx_t * restrict tB = data + tile_off(tr01, tc10);
                                        for (gate_idx_t lr = 0; lr < TILE_DIM; ++lr) {
                                            for (gate_idx_t lc = 0; lc < TILE_DIM; ++lc) {
                                                cplx_t tmp = tA[lc * TILE_DIM + lr];
                                                tA[lc * TILE_DIM + lr] = conj(tB[lr * TILE_DIM + lc]);
                                                tB[lr * TILE_DIM + lc] = conj(tmp);
                                            }
                                        }
                                    }
                                    else {
                                        cplx_t * restrict tB = data + tile_off(tc10, tr01);
                                        for (gate_idx_t i = 0; i < TILE_SIZE; ++i) {
                                            cplx_t tmp = tA[i]; tA[i] = tB[i]; tB[i] = tmp;
                                        }
                                    }
                                }

                                /* Pair F: |ctrl=0,tgt=0><ctrl=1,tgt=0| <--> |ctrl=0,tgt=0><ctrl=1,tgt=1| */
                                {
                                    cplx_t * restrict tA = data + tile_off(tc10, tr00);
                                    cplx_t * restrict tB = data + tile_off(tc11, tr00);
                                    for (gate_idx_t i = 0; i < TILE_SIZE; ++i) {
                                        cplx_t tmp = tA[i]; tA[i] = tB[i]; tB[i] = tmp;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
