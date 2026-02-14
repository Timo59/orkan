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
    cplx_t * restrict data = state->data;

    const qubit_t lo = (control < target) ? control : target;
    const qubit_t hi = (control > target) ? control : target;
    const gate_idx_t incr_ctrl = (gate_idx_t)1 << control;
    const gate_idx_t incr_tgt  = (gate_idx_t)1 << target;

    if (hi < LOG_TILE_DIM) {
        /*
         * Case 1: Both qubits within a tile.
         *
         * Both bits lo and hi are in local coordinates (0..LOG_TILE_DIM-1).
         * All elements of each swap pair are in the same tile.
         * Tiles store full TILE_DIM x TILE_DIM elements, so no conjugation needed.
         */
        const gate_idx_t n_tiles = (dim + TILE_DIM - 1) >> LOG_TILE_DIM;
        const gate_idx_t tile_dim = (dim < TILE_DIM) ? dim : TILE_DIM;
        const gate_idx_t quarter_td = tile_dim >> 2;

        #pragma omp parallel for schedule(dynamic) if(dim >= OMP_THRESHOLD)
        for (gate_idx_t tr = 0; tr < n_tiles; ++tr) {
            for (gate_idx_t tc = 0; tc <= tr; ++tc) {
                cplx_t *tile = data + tile_off(tr, tc);

                for (gate_idx_t br = 0; br < quarter_td; ++br) {
                    gate_idx_t lr00 = insertBits2_0(br, lo, hi);
                    gate_idx_t lr01 = lr00 | incr_tgt;
                    gate_idx_t lr10 = lr00 | incr_ctrl;
                    gate_idx_t lr11 = lr10 | incr_tgt;

                    for (gate_idx_t bc = 0; bc < quarter_td; ++bc) {
                        gate_idx_t lc00 = insertBits2_0(bc, lo, hi);
                        gate_idx_t lc01 = lc00 | incr_tgt;
                        gate_idx_t lc10 = lc00 | incr_ctrl;
                        gate_idx_t lc11 = lc10 | incr_tgt;

                        cplx_t tmp;

                        /* Pair A: (r10,c00) <-> (r11,c00) */
                        tmp = tile[elem_off(lr10, lc00)];
                        tile[elem_off(lr10, lc00)] = tile[elem_off(lr11, lc00)];
                        tile[elem_off(lr11, lc00)] = tmp;

                        /* Pair E: (r10,c10) <-> (r11,c11) */
                        tmp = tile[elem_off(lr10, lc10)];
                        tile[elem_off(lr10, lc10)] = tile[elem_off(lr11, lc11)];
                        tile[elem_off(lr11, lc11)] = tmp;

                        /* Pair B: (r11,c01) <-> (r10,c01) */
                        tmp = tile[elem_off(lr11, lc01)];
                        tile[elem_off(lr11, lc01)] = tile[elem_off(lr10, lc01)];
                        tile[elem_off(lr10, lc01)] = tmp;

                        /* Pair F: (r11,c10) <-> (r10,c11) */
                        tmp = tile[elem_off(lr11, lc10)];
                        tile[elem_off(lr11, lc10)] = tile[elem_off(lr10, lc11)];
                        tile[elem_off(lr10, lc11)] = tmp;

                        /* Pair C: (r00,c10) <-> (r00,c11) */
                        tmp = tile[elem_off(lr00, lc10)];
                        tile[elem_off(lr00, lc10)] = tile[elem_off(lr00, lc11)];
                        tile[elem_off(lr00, lc11)] = tmp;

                        /* Pair D: (r01,c10) <-> (r01,c11) */
                        tmp = tile[elem_off(lr01, lc10)];
                        tile[elem_off(lr01, lc10)] = tile[elem_off(lr01, lc11)];
                        tile[elem_off(lr01, lc11)] = tmp;
                    }
                }
            }
        }
    }
    else if (lo < LOG_TILE_DIM) {
        /*
         * Case 2: lo within tile, hi across tiles.
         *
         * Bit lo is in local coordinates, bit hi is in tile coordinates
         * (tile_bit = hi - LOG_TILE_DIM). The CX permutation mixes local
         * and tile coordinates depending on which qubit is control vs target.
         *
         * Row/col quartets split across tiles:
         *   r00 = (tr0, lr0), r01 = (tr0, lr1) or (tr1, lr0) depending on
         *         whether lo==target or lo==control.
         *
         * If control > target (hi==control, lo==target):
         *   incr_ctrl = 1<<control affects tile coord, incr_tgt = 1<<target affects local coord
         *   r00 = (tr0, lr0), r01 = (tr0, lr1), r10 = (tr1, lr0), r11 = (tr1, lr1)
         *   c00 = (tc0, lc0), c01 = (tc0, lc1), c10 = (tc1, lc0), c11 = (tc1, lc1)
         *
         * If control < target (hi==target, lo==control):
         *   incr_tgt = 1<<target affects tile coord, incr_ctrl = 1<<control affects local coord
         *   r00 = (tr0, lr0), r01 = (tr1, lr0), r10 = (tr0, lr1), r11 = (tr1, lr1)
         *   c00 = (tc0, lc0), c01 = (tc1, lc0), c10 = (tc0, lc1), c11 = (tc1, lc1)
         *
         * We use _tile and _loc suffixes to disambiguate which increment
         * goes to tile coordinates and which goes to local coordinates.
         */
        const qubit_t tile_bit = hi - LOG_TILE_DIM;
        const gate_idx_t tile_stride = (gate_idx_t)1 << tile_bit;
        const gate_idx_t n_tiles = (dim + TILE_DIM - 1) >> LOG_TILE_DIM;
        const gate_idx_t half_nt = n_tiles >> 1;
        const gate_idx_t tile_dim = (dim < TILE_DIM) ? dim : TILE_DIM;
        const gate_idx_t half_td = tile_dim >> 1;
        const gate_idx_t incr_loc = (gate_idx_t)1 << lo;

        /*
         * Determine whether hi is control or target.
         * hi_is_ctrl: hi==control, lo==target
         *   tile0 = ctrl=0, tile1 = ctrl=1
         *   loc0  = tgt=0,  loc1  = tgt=1
         * !hi_is_ctrl: hi==target, lo==control
         *   tile0 = tgt=0, tile1 = tgt=1
         *   loc0  = ctrl=0, loc1 = ctrl=1
         */
        const int hi_is_ctrl = (control > target);

        #pragma omp parallel for schedule(dynamic) if(dim >= OMP_THRESHOLD)
        for (gate_idx_t i_tr = 0; i_tr < half_nt; ++i_tr) {
            gate_idx_t tr0 = insertBit0(i_tr, tile_bit);
            gate_idx_t tr1 = tr0 | tile_stride;

            for (gate_idx_t i_tc = 0; i_tc <= i_tr; ++i_tc) {
                gate_idx_t tc0 = insertBit0(i_tc, tile_bit);
                gate_idx_t tc1 = tc0 | tile_stride;

                /* Tile pointers: T00, T11, T10 always in lower triangle */
                cplx_t *t00 = data + tile_off(tr0, tc0);
                cplx_t *t11 = data + tile_off(tr1, tc1);
                cplx_t *t10 = data + tile_off(tr1, tc0);

                /* T01 = T(tr0, tc1): only stored if tr0 >= tc1 */
                int lower_01 = (tr0 >= tc1);
                cplx_t *t01 = lower_01 ? data + tile_off(tr0, tc1)
                                       : data + tile_off(tc1, tr0);

                int is_diag = (i_tc == i_tr);

                /*
                 * Map the CX 6 swap pairs to (tile, local) coordinates.
                 *
                 * Notation: T[tile_row, tile_col][local_row, local_col]
                 *
                 * For hi_is_ctrl (tile bit = ctrl bit, local bit = tgt bit):
                 *   r_ab = (tr_a, lr_b), c_ab = (tc_a, lc_b) where a=ctrl, b=tgt
                 *   Pair A: T10[lr0,lc0] <-> T10[lr1,lc0]  -- t10 <-> T(tr1,tc0+tgt)
                 *     But wait: r10=(tr1,lr0), c00=(tc0,lc0) => T(tr1,tc0)[lr0,lc0]
                 *               r11=(tr1,lr1), c00=(tc0,lc0) => T(tr1,tc0)[lr1,lc0]
                 *   So Pair A: T(tr1,tc0)[lr0,lc0] <-> T(tr1,tc0)[lr1,lc0]  (within t10)
                 *   Pair E: T(tr1,tc1)[lr0,lc0] <-> T(tr1,tc1)[lr1,lc1]     (within t11... wait)
                 *     r10=(tr1,lr0), c10=(tc1,lc0) => T(tr1,tc1)[lr0,lc0]
                 *     r11=(tr1,lr1), c11=(tc1,lc1) => T(tr1,tc1)[lr1,lc1]
                 *   Pair B: T(tr1,tc0)[lr1,lc1] <-> T(tr1,tc0)[lr0,lc1]     (within t10)
                 *   Pair F: T(tr1,tc1)[lr1,lc0] <-> T(tr1,tc1)[lr0,lc1]     (within t11... wait)
                 *     r11=(tr1,lr1), c10=(tc1,lc0) => T(tr1,tc1)[lr1,lc0]
                 *     r10=(tr1,lr0), c11=(tc1,lc1) => T(tr1,tc1)[lr0,lc1]
                 *   Pair C: T(tr0,tc1)[lr0,lc0] <-> T(tr0,tc1)[lr0,lc1]     (t01)
                 *     r00=(tr0,lr0), c10=(tc1,lc0) => T(tr0,tc1)[lr0,lc0]
                 *     r00=(tr0,lr0), c11=(tc1,lc1) => T(tr0,tc1)[lr0,lc1]
                 *   Pair D: T(tr0,tc1)[lr1,lc0] <-> T(tr0,tc1)[lr1,lc1]     (t01)
                 *     r01=(tr0,lr1), c10=(tc1,lc0) => T(tr0,tc1)[lr1,lc0]
                 *     r01=(tr0,lr1), c11=(tc1,lc1) => T(tr0,tc1)[lr1,lc1]
                 *
                 * So for hi_is_ctrl: Pairs A,B within t10; Pairs E,F within t11;
                 *   Pairs C,D within t01 (which may be upper-tri).
                 *   t00 is UNTOUCHED.
                 *
                 * For !hi_is_ctrl (tile bit = tgt bit, local bit = ctrl bit):
                 *   r_ab = (tr_b, lr_a), c_ab = (tc_b, lc_a)
                 *   r00=(tr0,lr0), r01=(tr1,lr0), r10=(tr0,lr1), r11=(tr1,lr1)
                 *   c00=(tc0,lc0), c01=(tc1,lc0), c10=(tc0,lc1), c11=(tc1,lc1)
                 *
                 *   Pair A: (r10,c00)<->(r11,c00) = T(tr0,tc0)[lr1,lc0]<->T(tr1,tc0)[lr1,lc0]
                 *     i.e. t00[lr1,lc0] <-> t10[lr1,lc0]
                 *   Pair E: (r10,c10)<->(r11,c11) = T(tr0,tc0)[lr1,lc1]<->T(tr1,tc1)[lr1,lc1]
                 *     i.e. t00[lr1,lc1] <-> t11[lr1,lc1]
                 *   Pair B: (r11,c01)<->(r10,c01) = T(tr1,tc1)[lr1,lc0]<->T(tr0,tc1)[lr1,lc0]
                 *     i.e. t11[lr1,lc0] <-> t01[lr1,lc0]
                 *   Pair F: (r11,c10)<->(r10,c11) = T(tr1,tc0)[lr1,lc1]<->T(tr0,tc1)[lr1,lc1]
                 *     i.e. t10[lr1,lc1] <-> t01[lr1,lc1]
                 *   Pair C: (r00,c10)<->(r00,c11) = T(tr0,tc0)[lr0,lc1]<->T(tr0,tc1)[lr0,lc1]
                 *     i.e. t00[lr0,lc1] <-> t01[lr0,lc1]
                 *   Pair D: (r01,c10)<->(r01,c11) = T(tr1,tc0)[lr0,lc1]<->T(tr1,tc1)[lr0,lc1]
                 *     i.e. t10[lr0,lc1] <-> t11[lr0,lc1]
                 *
                 * So for !hi_is_ctrl: Pairs A,E touch t00; Pairs B,F touch t01;
                 *   all 6 pairs span 2 different tiles each.
                 */

                if (hi_is_ctrl) {
                    /*
                     * hi == control, lo == target.
                     * Pairs A,B are within t10 (always lower-tri off-diag tile).
                     * Pairs E,F are within t11.
                     * Pairs C,D are within t01 (may be upper-tri).
                     *
                     * When is_diag: t00 and t11 are diagonal tiles.
                     * CRITICAL: t01 and t10 point to the SAME physical tile
                     * (T(tc1,tr0) = T(tr1,tr0) = t10 when tr0==tc0).
                     * Pair C's swap t01[lc0,lr0]<->t01[lc1,lr0] is the
                     * transpose of Pair A's swap t10[lr0,lc0]<->t10[lr1,lc0]
                     * from the (b_lc,b_lr) iteration. Iterating all (b_lr,b_lc)
                     * for Pairs A,B already covers the work of Pairs C,D.
                     * So: skip Pairs C,D entirely when is_diag.
                     *
                     * Pairs E,F operate on t11 (diagonal tile): use
                     * dtile_read/dtile_write. Because t11 is diagonal,
                     * dtile helpers canonicalize through the lower triangle.
                     * Iterations (b_lr,b_lc) and (b_lc,b_lr) then access the
                     * SAME physical positions, causing each swap to be done
                     * twice (i.e. undone). Guard with b_lc <= b_lr so each
                     * canonical pair is processed exactly once.
                     */
                    for (gate_idx_t b_lr = 0; b_lr < half_td; ++b_lr) {
                        gate_idx_t lr0 = insertBit0(b_lr, lo);
                        gate_idx_t lr1 = lr0 | incr_loc;

                        for (gate_idx_t b_lc = 0; b_lc < half_td; ++b_lc) {
                            gate_idx_t lc0 = insertBit0(b_lc, lo);
                            gate_idx_t lc1 = lc0 | incr_loc;

                            cplx_t tmp;

                            /* Pair A: t10[lr0,lc0] <-> t10[lr1,lc0] */
                            tmp = t10[elem_off(lr0, lc0)];
                            t10[elem_off(lr0, lc0)] = t10[elem_off(lr1, lc0)];
                            t10[elem_off(lr1, lc0)] = tmp;

                            /* Pair B: t10[lr1,lc1] <-> t10[lr0,lc1] */
                            tmp = t10[elem_off(lr1, lc1)];
                            t10[elem_off(lr1, lc1)] = t10[elem_off(lr0, lc1)];
                            t10[elem_off(lr0, lc1)] = tmp;

                            if (is_diag && b_lc <= b_lr) {
                                /* t11 is a diagonal tile; use dtile helpers.
                                 * Restricted to b_lc <= b_lr: see comment above. */

                                /* Pair E: t11[lr0,lc0] <-> t11[lr1,lc1] */
                                cplx_t vE_a = dtile_read(t11, lr0, lc0);
                                cplx_t vE_b = dtile_read(t11, lr1, lc1);
                                dtile_write(t11, lr0, lc0, vE_b);
                                dtile_write(t11, lr1, lc1, vE_a);

                                /* Pair F: t11[lr1,lc0] <-> t11[lr0,lc1] */
                                cplx_t vF_a = dtile_read(t11, lr1, lc0);
                                cplx_t vF_b = dtile_read(t11, lr0, lc1);
                                dtile_write(t11, lr1, lc0, vF_b);
                                dtile_write(t11, lr0, lc1, vF_a);

                                /* Pairs C,D: SKIPPED when is_diag.
                                 * t01 aliases t10, so Pair C's transposed swap
                                 * on t01 is already covered by Pair A operating
                                 * on t10 across all (b_lr,b_lc) iterations. */
                            } else if (!is_diag) {
                                /* Pair E: t11[lr0,lc0] <-> t11[lr1,lc1] */
                                tmp = t11[elem_off(lr0, lc0)];
                                t11[elem_off(lr0, lc0)] = t11[elem_off(lr1, lc1)];
                                t11[elem_off(lr1, lc1)] = tmp;

                                /* Pair F: t11[lr1,lc0] <-> t11[lr0,lc1] */
                                tmp = t11[elem_off(lr1, lc0)];
                                t11[elem_off(lr1, lc0)] = t11[elem_off(lr0, lc1)];
                                t11[elem_off(lr0, lc1)] = tmp;

                                /* Pairs C,D: within t01 (may be upper-tri) */
                                if (lower_01) {
                                    /* Pair C: t01[lr0,lc0] <-> t01[lr0,lc1] */
                                    tmp = t01[elem_off(lr0, lc0)];
                                    t01[elem_off(lr0, lc0)] = t01[elem_off(lr0, lc1)];
                                    t01[elem_off(lr0, lc1)] = tmp;

                                    /* Pair D: t01[lr1,lc0] <-> t01[lr1,lc1] */
                                    tmp = t01[elem_off(lr1, lc0)];
                                    t01[elem_off(lr1, lc0)] = t01[elem_off(lr1, lc1)];
                                    t01[elem_off(lr1, lc1)] = tmp;
                                } else {
                                    /* t01 = T(tc1, tr0), access via transposed+conj.
                                     * T(tr0,tc1)[lr,lc] = conj(t01[lc,lr])
                                     * Swap conj(t01[lc0,lr0]) <-> conj(t01[lc1,lr0])
                                     * => swap t01[lc0,lr0] <-> t01[lc1,lr0] */
                                    /* Pair C: t01[lc0,lr0] <-> t01[lc1,lr0] */
                                    tmp = t01[elem_off(lc0, lr0)];
                                    t01[elem_off(lc0, lr0)] = t01[elem_off(lc1, lr0)];
                                    t01[elem_off(lc1, lr0)] = tmp;

                                    /* Pair D: t01[lc0,lr1] <-> t01[lc1,lr1] */
                                    tmp = t01[elem_off(lc0, lr1)];
                                    t01[elem_off(lc0, lr1)] = t01[elem_off(lc1, lr1)];
                                    t01[elem_off(lc1, lr1)] = tmp;
                                }
                            }
                        }
                    }

                    /* Mirror upper triangle of diagonal tiles. When is_diag,
                     * t11 is a diagonal tile whose lower triangle was modified
                     * by pairs E,F. t00 is untouched by CX when hi==control.
                     * dtile_write maintained both triangles, but do a full
                     * mirror pass for consistency. */
                    if (is_diag) {
                        for (gate_idx_t lr = 0; lr < TILE_DIM; ++lr) {
                            for (gate_idx_t lc = lr + 1; lc < TILE_DIM; ++lc) {
                                t11[elem_off(lr, lc)] = conj(t11[elem_off(lc, lr)]);
                            }
                        }
                    }
                } else {
                    /*
                     * hi == target, lo == control.
                     * Pair A: t00[lr1,lc0] <-> t10[lr1,lc0]
                     * Pair E: t00[lr1,lc1] <-> t11[lr1,lc1]
                     * Pair B: t11[lr1,lc0] <-> t01[lr1,lc0]
                     * Pair F: t10[lr1,lc1] <-> t01[lr1,lc1]
                     * Pair C: t00[lr0,lc1] <-> t01[lr0,lc1]
                     * Pair D: t10[lr0,lc1] <-> t11[lr0,lc1]
                     *
                     * Here lr0/lr1 refer to ctrl=0/1 in local coord.
                     * Pairs A,E modify t00 which may be diagonal.
                     * Pairs B,C,F involve t01 which may be upper-tri.
                     */
                    if (is_diag) {
                        /* t00 and t11 are diagonal tiles. t01==t10 by storage
                         * (when is_diag, lower_01==0, t01=T(tc1,tr0)=T(tr1,tr0)=t10).
                         *
                         * Blocks (b_lr, b_lc) and (b_lc, b_lr) produce the same
                         * physical swaps due to Hermiticity, so process only b_lc<=b_lr. */
                        for (gate_idx_t b_lr = 0; b_lr < half_td; ++b_lr) {
                            gate_idx_t lr0 = insertBit0(b_lr, lo);
                            gate_idx_t lr1 = lr0 | incr_loc;

                            for (gate_idx_t b_lc = 0; b_lc <= b_lr; ++b_lc) {
                                gate_idx_t lc0 = insertBit0(b_lc, lo);
                                gate_idx_t lc1 = lc0 | incr_loc;

                                /* Read all values needed before writing.
                                 * t01 = T(tr1,tr0) = t10, accessed as conj(t10[lc,lr])
                                 * for logical element T(tr0,tc1)[lr,lc]. */
                                cplx_t pA_a = dtile_read(t00, lr1, lc0);
                                cplx_t pA_b = t10[elem_off(lr1, lc0)];

                                cplx_t pE_a = dtile_read(t00, lr1, lc1);
                                cplx_t pE_b = dtile_read(t11, lr1, lc1);

                                /* t01 logical [lr1,lc0] = conj(t10[lc0,lr1]) */
                                cplx_t pB_a = dtile_read(t11, lr1, lc0);
                                cplx_t pB_b = conj(t10[elem_off(lc0, lr1)]);

                                /* t01 logical [lr1,lc1] = conj(t10[lc1,lr1]) */
                                cplx_t pF_a = t10[elem_off(lr1, lc1)];
                                cplx_t pF_b = conj(t10[elem_off(lc1, lr1)]);

                                /* t01 logical [lr0,lc1] = conj(t10[lc1,lr0]) */
                                cplx_t pC_a = dtile_read(t00, lr0, lc1);
                                cplx_t pC_b = conj(t10[elem_off(lc1, lr0)]);

                                cplx_t pD_a = t10[elem_off(lr0, lc1)];
                                cplx_t pD_b = dtile_read(t11, lr0, lc1);

                                /* Write back */
                                /* Pair A */
                                dtile_write(t00, lr1, lc0, pA_b);
                                t10[elem_off(lr1, lc0)] = pA_a;

                                /* Pair E */
                                dtile_write(t00, lr1, lc1, pE_b);
                                dtile_write(t11, lr1, lc1, pE_a);

                                /* Pair B: t11[lr1,lc0] <-> t01[lr1,lc0] */
                                dtile_write(t11, lr1, lc0, pB_b);
                                t10[elem_off(lc0, lr1)] = conj(pB_a);

                                /* Pair F: t10[lr1,lc1] <-> t01[lr1,lc1] */
                                t10[elem_off(lr1, lc1)] = pF_b;
                                t10[elem_off(lc1, lr1)] = conj(pF_a);

                                /* Pair C: t00[lr0,lc1] <-> t01[lr0,lc1] */
                                dtile_write(t00, lr0, lc1, pC_b);
                                t10[elem_off(lc1, lr0)] = conj(pC_a);

                                /* Pair D */
                                t10[elem_off(lr0, lc1)] = pD_b;
                                dtile_write(t11, lr0, lc1, pD_a);
                            }
                        }

                        /* Mirror upper triangle of diagonal tiles */
                        for (gate_idx_t lr = 0; lr < TILE_DIM; ++lr) {
                            for (gate_idx_t lc = lr + 1; lc < TILE_DIM; ++lc) {
                                t00[elem_off(lr, lc)] = conj(t00[elem_off(lc, lr)]);
                                t11[elem_off(lr, lc)] = conj(t11[elem_off(lc, lr)]);
                            }
                        }
                    } else {
                        for (gate_idx_t b_lr = 0; b_lr < half_td; ++b_lr) {
                            gate_idx_t lr0 = insertBit0(b_lr, lo);
                            gate_idx_t lr1 = lr0 | incr_loc;

                            for (gate_idx_t b_lc = 0; b_lc < half_td; ++b_lc) {
                                gate_idx_t lc0 = insertBit0(b_lc, lo);
                                gate_idx_t lc1 = lc0 | incr_loc;

                                cplx_t tmp;

                                /* Pair A: t00[lr1,lc0] <-> t10[lr1,lc0] */
                                tmp = t00[elem_off(lr1, lc0)];
                                t00[elem_off(lr1, lc0)] = t10[elem_off(lr1, lc0)];
                                t10[elem_off(lr1, lc0)] = tmp;

                                /* Pair E: t00[lr1,lc1] <-> t11[lr1,lc1] */
                                tmp = t00[elem_off(lr1, lc1)];
                                t00[elem_off(lr1, lc1)] = t11[elem_off(lr1, lc1)];
                                t11[elem_off(lr1, lc1)] = tmp;

                                /* Pair D: t10[lr0,lc1] <-> t11[lr0,lc1] */
                                tmp = t10[elem_off(lr0, lc1)];
                                t10[elem_off(lr0, lc1)] = t11[elem_off(lr0, lc1)];
                                t11[elem_off(lr0, lc1)] = tmp;

                                /* Pairs B, C, F involve t01 */
                                if (lower_01) {
                                    /* Pair B: t11[lr1,lc0] <-> t01[lr1,lc0] */
                                    tmp = t11[elem_off(lr1, lc0)];
                                    t11[elem_off(lr1, lc0)] = t01[elem_off(lr1, lc0)];
                                    t01[elem_off(lr1, lc0)] = tmp;

                                    /* Pair F: t10[lr1,lc1] <-> t01[lr1,lc1] */
                                    tmp = t10[elem_off(lr1, lc1)];
                                    t10[elem_off(lr1, lc1)] = t01[elem_off(lr1, lc1)];
                                    t01[elem_off(lr1, lc1)] = tmp;

                                    /* Pair C: t00[lr0,lc1] <-> t01[lr0,lc1] */
                                    tmp = t00[elem_off(lr0, lc1)];
                                    t00[elem_off(lr0, lc1)] = t01[elem_off(lr0, lc1)];
                                    t01[elem_off(lr0, lc1)] = tmp;
                                } else {
                                    /* t01 = T(tc1, tr0): element [lr,lc] of
                                     * logical T(tr0,tc1) stored as conj(t01[lc,lr]) */

                                    /* Pair B: t11[lr1,lc0] <-> conj(t01[lc0,lr1]) */
                                    {
                                        cplx_t val_a = t11[elem_off(lr1, lc0)];
                                        cplx_t val_b = conj(t01[elem_off(lc0, lr1)]);
                                        t11[elem_off(lr1, lc0)] = val_b;
                                        t01[elem_off(lc0, lr1)] = conj(val_a);
                                    }

                                    /* Pair F: t10[lr1,lc1] <-> conj(t01[lc1,lr1]) */
                                    {
                                        cplx_t val_a = t10[elem_off(lr1, lc1)];
                                        cplx_t val_b = conj(t01[elem_off(lc1, lr1)]);
                                        t10[elem_off(lr1, lc1)] = val_b;
                                        t01[elem_off(lc1, lr1)] = conj(val_a);
                                    }

                                    /* Pair C: t00[lr0,lc1] <-> conj(t01[lc1,lr0]) */
                                    {
                                        cplx_t val_a = t00[elem_off(lr0, lc1)];
                                        cplx_t val_b = conj(t01[elem_off(lc1, lr0)]);
                                        t00[elem_off(lr0, lc1)] = val_b;
                                        t01[elem_off(lc1, lr0)] = conj(val_a);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    else {
        /*
         * Case 3: Both qubits across tiles.
         *
         * Both bits lo and hi are in tile coordinates. The CX permutation
         * only affects tile indices, not local coordinates within tiles.
         *
         * tile_ctrl = control - LOG_TILE_DIM, tile_tgt = target - LOG_TILE_DIM.
         * tile_lo = lo - LOG_TILE_DIM, tile_hi = hi - LOG_TILE_DIM.
         *
         * We iterate tile-row quartets and tile-col quartets.
         * Row/col tile quartets: tr00, tr01, tr10, tr11 (bit pattern = ctrl,tgt).
         *
         * The 6 swap pairs at tile level:
         *   Pair A: T(tr10,tc00) <-> T(tr11,tc00)
         *   Pair E: T(tr10,tc10) <-> T(tr11,tc11)
         *   Pair B: T(tr11,tc01) <-> T(tr10,tc01)
         *   Pair F: T(tr11,tc10) <-> T(tr10,tc11)
         *   Pair C: T(tr00,tc10) <-> T(tr00,tc11)
         *   Pair D: T(tr01,tc10) <-> T(tr01,tc11)
         */
        const qubit_t tile_lo = lo - LOG_TILE_DIM;
        const qubit_t tile_hi = hi - LOG_TILE_DIM;
        const qubit_t tile_ctrl = control - LOG_TILE_DIM;
        const qubit_t tile_tgt  = target - LOG_TILE_DIM;
        const gate_idx_t incr_tile_ctrl = (gate_idx_t)1 << tile_ctrl;
        const gate_idx_t incr_tile_tgt  = (gate_idx_t)1 << tile_tgt;
        const gate_idx_t n_tiles = (dim + TILE_DIM - 1) >> LOG_TILE_DIM;
        const gate_idx_t quarter_nt = n_tiles >> 2;

        #pragma omp parallel for schedule(dynamic) if(dim >= OMP_THRESHOLD)
        for (gate_idx_t b_tr = 0; b_tr < quarter_nt; ++b_tr) {
            gate_idx_t tr00 = insertBits2_0(b_tr, tile_lo, tile_hi);
            gate_idx_t tr01 = tr00 | incr_tile_tgt;
            gate_idx_t tr10 = tr00 | incr_tile_ctrl;
            gate_idx_t tr11 = tr10 | incr_tile_tgt;

            for (gate_idx_t b_tc = 0; b_tc <= b_tr; ++b_tc) {
                gate_idx_t tc00 = insertBits2_0(b_tc, tile_lo, tile_hi);
                gate_idx_t tc01 = tc00 | incr_tile_tgt;
                gate_idx_t tc10 = tc00 | incr_tile_ctrl;
                gate_idx_t tc11 = tc10 | incr_tile_tgt;

                int is_diag = (b_tc == b_tr);

                /* ---- Pair A: T(tr10,tc00) <-> T(tr11,tc00) ----
                 * Both always lower-tri (tr10 > tr00 >= tc00, tr11 > tr00 >= tc00). */
                {
                    cplx_t *tA = data + tile_off(tr10, tc00);
                    cplx_t *tB = data + tile_off(tr11, tc00);
                    for (gate_idx_t i = 0; i < TILE_SIZE; ++i) {
                        cplx_t tmp = tA[i]; tA[i] = tB[i]; tB[i] = tmp;
                    }
                }

                /* ---- Pair E: T(tr10,tc10) <-> T(tr11,tc11) ----
                 * Both always lower-tri. */
                {
                    cplx_t *tA = data + tile_off(tr10, tc10);
                    cplx_t *tB = data + tile_off(tr11, tc11);
                    for (gate_idx_t i = 0; i < TILE_SIZE; ++i) {
                        cplx_t tmp = tA[i]; tA[i] = tB[i]; tB[i] = tmp;
                    }
                }

                /* ---- Pair B: T(tr11,tc01) <-> T(tr10,tc01) ----
                 * T(tr11,tc01): always lower-tri (tr11 has both ctrl+tgt bits set).
                 * T(tr10,tc01): lower-tri iff tr10 >= tc01.
                 *   When tile_ctrl < tile_tgt: tr10 = tr00|ctrl_bit (smaller),
                 *   tc01 = tc00|tgt_bit (larger), so tr10 < tc01 is possible.
                 *
                 * When is_diag and !lower_b: T(tr10,tc01) stored as T(tc01,tr10).
                 * Unlike Pair F, these are two distinct tiles (not the same physical
                 * tile), so this is NOT hermitian transpose in place. */
                {
                    cplx_t *tA = data + tile_off(tr11, tc01);
                    int lower_b = (tr10 >= tc01);
                    if (lower_b) {
                        cplx_t *tB = data + tile_off(tr10, tc01);
                        for (gate_idx_t i = 0; i < TILE_SIZE; ++i) {
                            cplx_t tmp = tA[i]; tA[i] = tB[i]; tB[i] = tmp;
                        }
                    } else {
                        /* T(tr10,tc01) upper-tri, stored as T(tc01, tr10) */
                        cplx_t *tB = data + tile_off(tc01, tr10);
                        for (gate_idx_t lr = 0; lr < TILE_DIM; ++lr) {
                            for (gate_idx_t lc = 0; lc < TILE_DIM; ++lc) {
                                cplx_t val_a = tA[elem_off(lr, lc)];
                                cplx_t val_b = conj(tB[elem_off(lc, lr)]);
                                tA[elem_off(lr, lc)] = val_b;
                                tB[elem_off(lc, lr)] = conj(val_a);
                            }
                        }
                    }
                }

                /* ---- Pair F: T(tr11,tc10) <-> T(tr10,tc11) ----
                 * T(tr11,tc10): always lower-tri.
                 * T(tr10,tc11): lower-tri iff tr10 >= tc11.
                 *
                 * When is_diag (tr00==tc00): tr10=tc10, tc11=tr11, so
                 * T(tr10,tc11)=T(tc10,tr11). Since tc10<tr11, T(tr10,tr11)
                 * is upper-tri, stored as T(tr11,tr10). The swap pair
                 * becomes T(tr11,tc10) <-> T(tc10,tr11) where both
                 * physically point to T(tr11,tr10). This is Hermitian
                 * transpose in place. */
                {
                    cplx_t *tA = data + tile_off(tr11, tc10);
                    int lower_f = (tr10 >= tc11);
                    if (lower_f) {
                        cplx_t *tB = data + tile_off(tr10, tc11);
                        for (gate_idx_t i = 0; i < TILE_SIZE; ++i) {
                            cplx_t tmp = tA[i]; tA[i] = tB[i]; tB[i] = tmp;
                        }
                    } else if (is_diag) {
                        /* tA = T(tr11, tc10) = T(tr11, tr10). Hermitian transpose in place.
                         * Logical swap: tA[lr,lc] <-> conj(tA[lc,lr]) */
                        for (gate_idx_t lr = 0; lr < TILE_DIM; ++lr) {
                            tA[elem_off(lr, lr)] = conj(tA[elem_off(lr, lr)]);
                            for (gate_idx_t lc = lr + 1; lc < TILE_DIM; ++lc) {
                                cplx_t val_a = tA[elem_off(lr, lc)];
                                cplx_t val_b = tA[elem_off(lc, lr)];
                                tA[elem_off(lr, lc)] = conj(val_b);
                                tA[elem_off(lc, lr)] = conj(val_a);
                            }
                        }
                    } else {
                        /* T(tr10,tc11) upper-tri, stored as T(tc11, tr10) */
                        cplx_t *tB = data + tile_off(tc11, tr10);
                        for (gate_idx_t lr = 0; lr < TILE_DIM; ++lr) {
                            for (gate_idx_t lc = 0; lc < TILE_DIM; ++lc) {
                                cplx_t val_a = tA[elem_off(lr, lc)];
                                cplx_t val_b = conj(tB[elem_off(lc, lr)]);
                                tA[elem_off(lr, lc)] = val_b;
                                tB[elem_off(lc, lr)] = conj(val_a);
                            }
                        }
                    }
                }

                /* ---- Pairs C+D ----
                 * Pair C: T(tr00,tc10) <-> T(tr00,tc11)
                 * Pair D: T(tr01,tc10) <-> T(tr01,tc11)
                 *
                 * Lower-tri conditions:
                 *   T(tr00,tc10): iff tr00 >= tc10    [lo_10]
                 *   T(tr00,tc11): iff tr00 >= tc11    [lo_11]
                 *   T(tr01,tc10): iff tr01 >= tc10    [lo2_10]
                 *   T(tr01,tc11): iff tr01 >= tc11    [lo2_11]
                 *
                 * When is_diag: tr00==tc00, so tc10>tr00 and tc11>tr00,
                 * meaning lo_10==0 and lo_11==0.
                 *
                 * Pairs C+D, when is_diag, map to the same physical
                 * storage as Pairs A+B (the stored lower-tri reps coincide),
                 * so executing them would undo A+B. Skip.
                 */
                if (!is_diag) {
                    int lo_10 = (tr00 >= tc10);
                    int lo_11 = (tr00 >= tc11);

                    /* Pair C: T(tr00,tc10) <-> T(tr00,tc11) */
                    if (lo_10 && lo_11) {
                        cplx_t *tCa = data + tile_off(tr00, tc10);
                        cplx_t *tCb = data + tile_off(tr00, tc11);
                        for (gate_idx_t i = 0; i < TILE_SIZE; ++i) {
                            cplx_t tmp = tCa[i]; tCa[i] = tCb[i]; tCb[i] = tmp;
                        }
                    }
                    else if (lo_10 && !lo_11) {
                        /* T(tr00,tc10) lower, T(tr00,tc11) upper -> T(tc11,tr00) */
                        cplx_t *tCa = data + tile_off(tr00, tc10);
                        cplx_t *tCb = data + tile_off(tc11, tr00);
                        for (gate_idx_t lr = 0; lr < TILE_DIM; ++lr) {
                            for (gate_idx_t lc = 0; lc < TILE_DIM; ++lc) {
                                cplx_t val_a = tCa[elem_off(lr, lc)];
                                cplx_t val_b = conj(tCb[elem_off(lc, lr)]);
                                tCa[elem_off(lr, lc)] = val_b;
                                tCb[elem_off(lc, lr)] = conj(val_a);
                            }
                        }
                    }
                    else {
                        /* Both upper-tri: swap stored reps directly.
                         * T(tr00,tc10) stored as T(tc10,tr00).
                         * T(tr00,tc11) stored as T(tc11,tr00). */
                        cplx_t *tCa = data + tile_off(tc10, tr00);
                        cplx_t *tCb = data + tile_off(tc11, tr00);
                        for (gate_idx_t i = 0; i < TILE_SIZE; ++i) {
                            cplx_t tmp = tCa[i]; tCa[i] = tCb[i]; tCb[i] = tmp;
                        }
                    }

                    /* Pair D: T(tr01,tc10) <-> T(tr01,tc11) */
                    int lo2_10 = (tr01 >= tc10);
                    int lo2_11 = (tr01 >= tc11);

                    if (lo2_10 && lo2_11) {
                        cplx_t *tDa = data + tile_off(tr01, tc10);
                        cplx_t *tDb = data + tile_off(tr01, tc11);
                        for (gate_idx_t i = 0; i < TILE_SIZE; ++i) {
                            cplx_t tmp = tDa[i]; tDa[i] = tDb[i]; tDb[i] = tmp;
                        }
                    }
                    else if (lo2_10 && !lo2_11) {
                        cplx_t *tDa = data + tile_off(tr01, tc10);
                        cplx_t *tDb = data + tile_off(tc11, tr01);
                        for (gate_idx_t lr = 0; lr < TILE_DIM; ++lr) {
                            for (gate_idx_t lc = 0; lc < TILE_DIM; ++lc) {
                                cplx_t val_a = tDa[elem_off(lr, lc)];
                                cplx_t val_b = conj(tDb[elem_off(lc, lr)]);
                                tDa[elem_off(lr, lc)] = val_b;
                                tDb[elem_off(lc, lr)] = conj(val_a);
                            }
                        }
                    }
                    else {
                        cplx_t *tDa = data + tile_off(tc10, tr01);
                        cplx_t *tDb = data + tile_off(tc11, tr01);
                        for (gate_idx_t i = 0; i < TILE_SIZE; ++i) {
                            cplx_t tmp = tDa[i]; tDa[i] = tDb[i]; tDb[i] = tmp;
                        }
                    }
                }
            }
        }
    }
}
