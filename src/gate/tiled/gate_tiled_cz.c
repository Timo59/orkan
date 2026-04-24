/**
 * @file gate_tiled_cz.c
 * @brief CZ gate for tiled density matrices
 *
 * CZ = diag(1, 1, 1, -1) is a zero-flop diagonal gate (self-adjoint): rho' = CZ * rho * CZ.
 * phi(x) = -1 if both control and target bits of x are set, else +1.
 * rho'[r,c] = phi(r) * conj(phi(c)) * rho[r,c] = phi(r) * phi(c) * rho[r,c].
 *
 * Since phi is real-valued (+/-1), the element is negated when exactly one of
 * row or col has both control+target bits set:
 *
 *   Group A (row=11): negate (r11,c00), (r11,c01), (r11,c10) -- always lower-tri
 *   Group B (col=11): negate (r00,c11), (r01,c11), (r10,c11) -- may cross diagonal
 *   Leave:            (r11,c11) -- phases cancel: (-1)*(-1) = +1
 *   Leave:            all other (rXY,cXY) where neither XY=11
 *
 * In packed lower-triangular storage, Group A and Group B elements occupy
 * distinct memory locations even within a diagonal tile. Both must be negated
 * explicitly. (The Hermitian constraint rho[r,c] = conj(rho[c,r]) means both
 * a position and its mirror must be negated for consistency -- negating one
 * stored element does NOT automatically negate the other.)
 *
 * Four cases based on qubit placement relative to LOG_TILE_DIM:
 *   Case 1: both qubits local (< LOG_TILE_DIM)
 *   Case 2: control local, target tile-level
 *   Case 3: control tile-level, target local
 *   Case 4: both qubits tile-level
 */

#include "gate_tiled.h"

void cz_tiled(state_t *state, const qubit_t control, const qubit_t target) {
    const idx_t dim = (idx_t)1 << state->qubits;
    const idx_t dim_tile = (dim + TILE_DIM - 1) >> LOG_TILE_DIM;
    cplx_t * restrict data = state->data;

    if (control < LOG_TILE_DIM) {

        /*
         * =====================================================================================
         * Case 1: Both control AND target < LOG_TILE_DIM
         * =====================================================================================
         *
         * All affected index pairs are local within each tile. We iterate over
         * every lower-triangular tile (btr, btc) and within each, over all 4x4
         * blocks defined by the two qubit positions.
         *
         * For each block (blr, blc):
         *   r00, r01, r10, r11 are the 4 row indices (ctrl=a, tgt=b)
         *   c00, c01, c10, c11 are the 4 col indices
         *
         * Group A: negate tile[r11, c00], tile[r11, c01], tile[r11, c10]
         *          These are always in the lower triangle (r11 is the largest row index)
         * Group B: negate tile[r00, c11], tile[r01, c11], tile[r10, c11]
         */
        if (target < LOG_TILE_DIM) {
            const idx_t n_base = (dim < TILE_DIM) ? dim >> 2 : TILE_DIM >> 2;
            const qubit_t lo = (control < target) ? control : target;
            const qubit_t hi = (control > target) ? control : target;
            const idx_t incr_ctrl = (idx_t)1 << control;
            const idx_t incr_tgt = (idx_t)1 << target;

            #pragma omp parallel for schedule(static) if(dim >= OMP_THRESHOLD_MIXED)
            for (idx_t btr = 0; btr < dim_tile; ++btr) {
                for (idx_t btc = 0; btc <= btr; ++btc) {
                    cplx_t * restrict tile = data + tile_off(btr, btc);

                    for (idx_t blr = 0; blr < n_base; ++blr) {
                        const idx_t r00 = insertBits2_0(blr, lo, hi);
                        const idx_t r01 = r00 | incr_tgt;
                        const idx_t r10 = r00 | incr_ctrl;
                        const idx_t r11 = r10 | incr_tgt;

                        for (idx_t blc = 0; blc < n_base; ++blc) {
                            const idx_t c00 = insertBits2_0(blc, lo, hi);
                            const idx_t c01 = c00 | incr_tgt;
                            const idx_t c10 = c00 | incr_ctrl;
                            const idx_t c11 = c10 | incr_tgt;

                            /* Group A: row=11, col!=11 -- always lower-tri */
                            tile[elem_off(r11, c00)] = -tile[elem_off(r11, c00)];
                            tile[elem_off(r11, c01)] = -tile[elem_off(r11, c01)];
                            tile[elem_off(r11, c10)] = -tile[elem_off(r11, c10)];

                            /* Group B: col=11, row!=11 */
                            tile[elem_off(r00, c11)] = -tile[elem_off(r00, c11)];
                            tile[elem_off(r01, c11)] = -tile[elem_off(r01, c11)];
                            tile[elem_off(r10, c11)] = -tile[elem_off(r10, c11)];
                        }
                    }
                }
            }
        }

        /*
         * =====================================================================================
         * Case 2: control < LOG_TILE_DIM, target >= LOG_TILE_DIM
         * =====================================================================================
         *
         * The target bit is a tile-level bit; the control bit is tile-local.
         * For each outer pair (btr, btc) we have 4 tiles:
         *   t00 = (tr0, tc0): row target=0, col target=0
         *   t10 = (tr1, tc0): row target=1, col target=0
         *   t11 = (tr1, tc1): row target=1, col target=1
         *   t01 = (tr0, tc1): row target=0, col target=1
         *
         * "row=11" means tile-row has target=1 AND local-row has control=1.
         * "col=11" means tile-col has target=1 AND local-col has control=1.
         *
         * t00: no element can be row=11 or col=11 -> skip entirely
         *
         * t10 (row target=1, col target=0):
         *   row=11 when local row has control bit -> Group A
         *   col=11 never (col target=0) -> no Group B
         *   Action: negate ALL elements where local row has control bit set,
         *           i.e., negate tile[r1, c0] and tile[r1, c1] for each block
         *
         * t11 (row target=1, col target=1):
         *   row=11 when local row has control bit
         *   col=11 when local col has control bit
         *   Negate when exactly one is "11":
         *     tile[r1, c0]: row=11, col!=11 -> negate (Group A in t11)
         *     tile[r0, c1]: row!=11, col=11 -> negate (Group B in t11)
         *     tile[r0, c0]: neither -> leave
         *     tile[r1, c1]: both -> leave (phases cancel)
         *
         * t01 (row target=0, col target=1):
         *   row=11 never (row target=0) -> no Group A
         *   col=11 when local col has control bit -> Group B
         *   Action: negate ALL elements where local col has control bit set,
         *           i.e., negate tile[r0, c1] and tile[r1, c1] for each block
         *
         * Triangle handling for t01 follows the 3-branch pattern from cx_tiled.
         */
        else {
            const idx_t n_bt = dim_tile >> 1;
            const qubit_t tile_bit = target - LOG_TILE_DIM;
            const idx_t incr_tile = (idx_t)1 << tile_bit;
            const idx_t n_bl = TILE_DIM >> 1;
            const idx_t incr_local = (idx_t)1 << control;

            #pragma omp parallel for schedule(static, 1) if(dim >= OMP_THRESHOLD_MIXED)
            for (idx_t btr = 0; btr < n_bt; ++btr) {
                const idx_t tr0 = insertBit0(btr, tile_bit);
                const idx_t tr1 = tr0 | incr_tile;

                for (idx_t btc = 0; btc <= btr; ++btc) {
                    const idx_t tc0 = insertBit0(btc, tile_bit);
                    const idx_t tc1 = tc0 | incr_tile;

                    cplx_t * restrict t10 = data + tile_off(tr1, tc0);
                    cplx_t * restrict t11 = data + tile_off(tr1, tc1);

                    /* Fast path: all 4 tiles are in the lower triangle (tr0 > tc1) */
                    if (tr0 > tc1) {
                        cplx_t * restrict t01 = data + tile_off(tr0, tc1);

                        for (idx_t blr = 0; blr < n_bl; ++blr) {
                            const idx_t r0 = insertBit0(blr, control);
                            const idx_t r1 = r0 | incr_local;

                            for (idx_t blc = 0; blc < n_bl; ++blc) {
                                const idx_t c0 = insertBit0(blc, control);
                                const idx_t c1 = c0 | incr_local;

                                /* t10: negate where local row has control bit (Group A) */
                                t10[elem_off(r1, c0)] = -t10[elem_off(r1, c0)];
                                t10[elem_off(r1, c1)] = -t10[elem_off(r1, c1)];

                                /* t11: negate where exactly one of row/col has control bit */
                                t11[elem_off(r1, c0)] = -t11[elem_off(r1, c0)];
                                t11[elem_off(r0, c1)] = -t11[elem_off(r0, c1)];

                                /* t01: negate where local col has control bit (Group B) */
                                t01[elem_off(r0, c1)] = -t01[elem_off(r0, c1)];
                                t01[elem_off(r1, c1)] = -t01[elem_off(r1, c1)];
                            }
                        }
                    }

                    /* t01 is upper-tri; its stored adjoint is at (tc1, tr0) in lower-tri */
                    else if (btr != btc) {
                        cplx_t * restrict t01_adj = data + tile_off(tc1, tr0);

                        for (idx_t blr = 0; blr < n_bl; ++blr) {
                            const idx_t r0 = insertBit0(blr, control);
                            const idx_t r1 = r0 | incr_local;

                            for (idx_t blc = 0; blc < n_bl; ++blc) {
                                const idx_t c0 = insertBit0(blc, control);
                                const idx_t c1 = c0 | incr_local;

                                /* t10: negate where local row has control bit (Group A) */
                                t10[elem_off(r1, c0)] = -t10[elem_off(r1, c0)];
                                t10[elem_off(r1, c1)] = -t10[elem_off(r1, c1)];

                                /* t11: negate where exactly one of row/col has control bit */
                                t11[elem_off(r1, c0)] = -t11[elem_off(r1, c0)];
                                t11[elem_off(r0, c1)] = -t11[elem_off(r0, c1)];

                                /* t01 stored as adjoint at (tc1, tr0):
                                 * rho[r, c] at t01 = conj(stored[c_local, r_local]) at t01_adj.
                                 * Negating rho[r, c] means negating stored[c_local, r_local].
                                 * Since -conj(z) = conj(-z), just negate the stored element.
                                 *
                                 * Group B in t01: negate where local col has control bit.
                                 * Element t01[r0, c1] -> stored at t01_adj[c1, r0]
                                 * Element t01[r1, c1] -> stored at t01_adj[c1, r1]
                                 */
                                t01_adj[elem_off(c1, r0)] = -t01_adj[elem_off(c1, r0)];
                                t01_adj[elem_off(c1, r1)] = -t01_adj[elem_off(c1, r1)];
                            }
                        }
                    }

                    /* Diagonal block: btr == btc, so t01 = (tr0, tr1) is in the upper
                     * triangle and its stored adjoint is (tr1, tr0) = t10.
                     *
                     * For t10 = (tr1, tr0): every element t10[lr, lc] represents
                     * rho[tr1*TD+lr, tr0*TD+lc]. Since tr0 has target=0, phi_col = +1
                     * always. Since tr1 has target=1, phi_row = -1 iff lr has control bit.
                     * So negate all elements where LOCAL ROW has control bit.
                     *
                     * The Group B elements (from t01's adjoint) require negating t10[lc, lr]
                     * when lc has control bit. But "negate t10[a, b] when a has control bit"
                     * already covers this since both conditions select the same set of
                     * elements (those with first index having control bit set).
                     *
                     * For t11 = (tr1, tr1), a diagonal tile: negate where exactly one of
                     * (lr, lc) has control bit. Both t11[r1, c0] and t11[r0, c1] are
                     * independent positions in raw storage, and their Hermitian mirrors
                     * are handled automatically since we negate raw storage directly.
                     */
                    else {
                        for (idx_t blr = 0; blr < n_bl; ++blr) {
                            const idx_t r0 = insertBit0(blr, control);
                            const idx_t r1 = r0 | incr_local;

                            for (idx_t blc = 0; blc < n_bl; ++blc) {
                                const idx_t c0 = insertBit0(blc, control);
                                const idx_t c1 = c0 | incr_local;

                                /* t10: negate all elements where local row has control bit.
                                 * This handles both Group A (direct) and Group B (adjoint). */
                                t10[elem_off(r1, c0)] = -t10[elem_off(r1, c0)];
                                t10[elem_off(r1, c1)] = -t10[elem_off(r1, c1)];

                                /* t11 (diagonal tile): negate where exactly one of
                                 * row/col has control bit */
                                t11[elem_off(r1, c0)] = -t11[elem_off(r1, c0)];
                                t11[elem_off(r0, c1)] = -t11[elem_off(r0, c1)];
                            }
                        }
                    }
                }
            }
        }
    }

    else {
        /*
         * =====================================================================================
         * Case 3: control >= LOG_TILE_DIM, target < LOG_TILE_DIM
         * =====================================================================================
         *
         * Symmetric to Case 2 with control and target roles swapped at the
         * tile/local level: the control bit is now the tile-level bit, and
         * the target bit is the local bit.
         *
         * "row=11" means tile-row has control=1 AND local-row has target=1.
         * "col=11" means tile-col has control=1 AND local-col has target=1.
         *
         * t00 = (tr0, tc0): neither row nor col can be "11" -> skip
         * t10 = (tr1, tc0): row has control=1, col has control=0
         *   row=11 when local row has target bit -> Group A
         *   col=11 never -> no Group B
         *   Action: negate where local row has target bit: t10[lr1, lc0], t10[lr1, lc1]
         *
         * t11 = (tr1, tc1): both row and col have control=1
         *   row=11 when local row has target bit; col=11 when local col has target bit
         *   Negate when exactly one is "11":
         *     t11[lr1, lc0]: Group A in t11
         *     t11[lr0, lc1]: Group B in t11
         *     t11[lr0, lc0]: neither -> leave
         *     t11[lr1, lc1]: both -> leave (cancel)
         *
         * t01 = (tr0, tc1): row has control=0, col has control=1
         *   col=11 when local col has target bit -> Group B
         *   Action: negate where local col has target bit: t01[lr0, lc1], t01[lr1, lc1]
         */
        if (target < LOG_TILE_DIM) {
            const idx_t n_bt = dim_tile >> 1;
            const qubit_t tile_bit = control - LOG_TILE_DIM;
            const idx_t incr_tile = (idx_t)1 << tile_bit;
            const idx_t n_bl = TILE_DIM >> 1;
            const idx_t incr_local = (idx_t)1 << target;

            #pragma omp parallel for schedule(static, 1) if(dim >= OMP_THRESHOLD_MIXED)
            for (idx_t btr = 0; btr < n_bt; ++btr) {
                const idx_t tr0 = insertBit0(btr, tile_bit);
                const idx_t tr1 = tr0 | incr_tile;

                for (idx_t btc = 0; btc <= btr; ++btc) {
                    const idx_t tc0 = insertBit0(btc, tile_bit);
                    const idx_t tc1 = tc0 | incr_tile;

                    cplx_t * restrict t10 = data + tile_off(tr1, tc0);
                    cplx_t * restrict t11 = data + tile_off(tr1, tc1);

                    /* Fast path: all 4 tiles in lower triangle */
                    if (tr0 > tc1) {
                        cplx_t * restrict t01 = data + tile_off(tr0, tc1);

                        for (idx_t blr = 0; blr < n_bl; ++blr) {
                            const idx_t lr0 = insertBit0(blr, target);
                            const idx_t lr1 = lr0 | incr_local;

                            for (idx_t blc = 0; blc < n_bl; ++blc) {
                                const idx_t lc0 = insertBit0(blc, target);
                                const idx_t lc1 = lc0 | incr_local;

                                /* t10: negate where local row has target bit (Group A) */
                                t10[elem_off(lr1, lc0)] = -t10[elem_off(lr1, lc0)];
                                t10[elem_off(lr1, lc1)] = -t10[elem_off(lr1, lc1)];

                                /* t11: negate where exactly one has target bit */
                                t11[elem_off(lr1, lc0)] = -t11[elem_off(lr1, lc0)];
                                t11[elem_off(lr0, lc1)] = -t11[elem_off(lr0, lc1)];

                                /* t01: negate where local col has target bit (Group B) */
                                t01[elem_off(lr0, lc1)] = -t01[elem_off(lr0, lc1)];
                                t01[elem_off(lr1, lc1)] = -t01[elem_off(lr1, lc1)];
                            }
                        }
                    }

                    /* t01 upper-tri; adjoint stored at (tc1, tr0) */
                    else if (btr != btc) {
                        cplx_t * restrict t01_adj = data + tile_off(tc1, tr0);

                        for (idx_t blr = 0; blr < n_bl; ++blr) {
                            const idx_t lr0 = insertBit0(blr, target);
                            const idx_t lr1 = lr0 | incr_local;

                            for (idx_t blc = 0; blc < n_bl; ++blc) {
                                const idx_t lc0 = insertBit0(blc, target);
                                const idx_t lc1 = lc0 | incr_local;

                                /* t10: Group A */
                                t10[elem_off(lr1, lc0)] = -t10[elem_off(lr1, lc0)];
                                t10[elem_off(lr1, lc1)] = -t10[elem_off(lr1, lc1)];

                                /* t11: exactly one "11" */
                                t11[elem_off(lr1, lc0)] = -t11[elem_off(lr1, lc0)];
                                t11[elem_off(lr0, lc1)] = -t11[elem_off(lr0, lc1)];

                                /* t01 stored as adjoint: negate transposed positions.
                                 * t01[lr0, lc1] -> stored at t01_adj[lc1, lr0]
                                 * t01[lr1, lc1] -> stored at t01_adj[lc1, lr1] */
                                t01_adj[elem_off(lc1, lr0)] = -t01_adj[elem_off(lc1, lr0)];
                                t01_adj[elem_off(lc1, lr1)] = -t01_adj[elem_off(lc1, lr1)];
                            }
                        }
                    }

                    /* Diagonal block: btr == btc.
                     * t01 = (tr0, tr1), adjoint = (tr1, tr0) = t10.
                     *
                     * Same logic as Case 2 diagonal: for t10 = (tr1, tr0), every element
                     * rho[tr1*TD+lr, tr0*TD+lc] has phi_col = +1 (tr0 has control=0)
                     * and phi_row = -1 iff lr has target bit (tr1 has control=1).
                     * So negate all elements where local row has target bit.
                     *
                     * For t11 (diagonal tile): negate where exactly one of (lr, lc) has
                     * target bit. Raw storage negation preserves Hermitian consistency. */
                    else {
                        for (idx_t blr = 0; blr < n_bl; ++blr) {
                            const idx_t lr0 = insertBit0(blr, target);
                            const idx_t lr1 = lr0 | incr_local;

                            for (idx_t blc = 0; blc < n_bl; ++blc) {
                                const idx_t lc0 = insertBit0(blc, target);
                                const idx_t lc1 = lc0 | incr_local;

                                /* t10: negate all elements where local row has target bit */
                                t10[elem_off(lr1, lc0)] = -t10[elem_off(lr1, lc0)];
                                t10[elem_off(lr1, lc1)] = -t10[elem_off(lr1, lc1)];

                                /* t11 (diagonal tile): negate where exactly one has target bit */
                                t11[elem_off(lr1, lc0)] = -t11[elem_off(lr1, lc0)];
                                t11[elem_off(lr0, lc1)] = -t11[elem_off(lr0, lc1)];
                            }
                        }
                    }
                }
            }
        }

        /*
         * =====================================================================================
         * Case 4: Both control AND target >= LOG_TILE_DIM
         * =====================================================================================
         *
         * Both qubits are tile-level bits. "row=11" means tile-row has both
         * control and target bits set. "col=11" means tile-col has both set.
         * Entire tiles are negated (no element-level distinction needed).
         *
         * From the 16 tile combinations (trXY, tcAB):
         *   Group A tiles (row=11, col!=11): (tr11,tc00), (tr11,tc01), (tr11,tc10)
         *   Group B tiles (col=11, row!=11): (tr00,tc11), (tr01,tc11), (tr10,tc11)
         *   Unchanged: (tr11,tc11) and all tiles with neither row=11 nor col=11
         *
         * Group A tiles are always in lower triangle (tr11 > tcXY for XY != 11).
         * Group B tiles may cross the diagonal.
         */
        else {
            const idx_t n_bt = dim_tile >> 2;
            const qubit_t lo = (control < target) ? control - LOG_TILE_DIM : target - LOG_TILE_DIM;
            const qubit_t hi = (control > target) ? control - LOG_TILE_DIM : target - LOG_TILE_DIM;
            const idx_t incr_ctrl = (idx_t)1 << (control - LOG_TILE_DIM);
            const idx_t incr_tgt = (idx_t)1 << (target - LOG_TILE_DIM);

            #pragma omp parallel for schedule(static, 1) if(dim >= OMP_THRESHOLD_MIXED)
            for (idx_t btr = 0; btr < n_bt; ++btr) {
                const idx_t tr00 = insertBits2_0(btr, lo, hi);
                const idx_t tr01 = tr00 | incr_tgt;
                const idx_t tr10 = tr00 | incr_ctrl;
                const idx_t tr11 = tr10 | incr_tgt;

                for (idx_t btc = 0; btc <= btr; ++btc) {
                    const idx_t tc00 = insertBits2_0(btc, lo, hi);
                    const idx_t tc01 = tc00 | incr_tgt;
                    const idx_t tc10 = tc00 | incr_ctrl;
                    const idx_t tc11 = tc10 | incr_tgt;

                    /* Group A: negate entire tiles (tr11, tc00), (tr11, tc01), (tr11, tc10)
                     * These are always in lower triangle since tr11 > tc00, tc01, tc10 */
                    {
                        cplx_t * restrict tA0 = data + tile_off(tr11, tc00);
                        for (idx_t i = 0; i < TILE_SIZE; ++i)
                            tA0[i] = -tA0[i];
                    }
                    {
                        cplx_t * restrict tA1 = data + tile_off(tr11, tc01);
                        for (idx_t i = 0; i < TILE_SIZE; ++i)
                            tA1[i] = -tA1[i];
                    }
                    {
                        cplx_t * restrict tA2 = data + tile_off(tr11, tc10);
                        for (idx_t i = 0; i < TILE_SIZE; ++i)
                            tA2[i] = -tA2[i];
                    }

                    /* Group B: negate tiles (tr00, tc11), (tr01, tc11), (tr10, tc11)
                     * These may be in upper triangle depending on tile indices. */

                    /* Fast path: all Group B tiles in lower triangle */
                    if (tr00 > tc11) {
                        {
                            cplx_t * restrict tB0 = data + tile_off(tr00, tc11);
                            for (idx_t i = 0; i < TILE_SIZE; ++i)
                                tB0[i] = -tB0[i];
                        }
                        {
                            cplx_t * restrict tB1 = data + tile_off(tr01, tc11);
                            for (idx_t i = 0; i < TILE_SIZE; ++i)
                                tB1[i] = -tB1[i];
                        }
                        {
                            cplx_t * restrict tB2 = data + tile_off(tr10, tc11);
                            for (idx_t i = 0; i < TILE_SIZE; ++i)
                                tB2[i] = -tB2[i];
                        }
                    }

                    else if (btr != btc) {
                        /* Some Group B tiles are in upper triangle.
                         * Their stored adjoints are in lower triangle.
                         * Since -conj(z) = conj(-z), we just negate the adjoint tile. */

                        /* (tr00, tc11): always upper-tri here (tr00 <= tc11).
                         * Adjoint at (tc11, tr00) -- always in lower triangle */
                        {
                            cplx_t * restrict tB0 = data + tile_off(tc11, tr00);
                            for (idx_t i = 0; i < TILE_SIZE; ++i)
                                tB0[i] = -tB0[i];
                        }

                        /* (tr01, tc11): check triangle position */
                        if (tr01 > tc11) {
                            cplx_t * restrict tB1 = data + tile_off(tr01, tc11);
                            for (idx_t i = 0; i < TILE_SIZE; ++i)
                                tB1[i] = -tB1[i];
                        } else {
                            cplx_t * restrict tB1 = data + tile_off(tc11, tr01);
                            for (idx_t i = 0; i < TILE_SIZE; ++i)
                                tB1[i] = -tB1[i];
                        }

                        /* (tr10, tc11): check triangle position */
                        if (tr10 > tc11) {
                            cplx_t * restrict tB2 = data + tile_off(tr10, tc11);
                            for (idx_t i = 0; i < TILE_SIZE; ++i)
                                tB2[i] = -tB2[i];
                        } else {
                            cplx_t * restrict tB2 = data + tile_off(tc11, tr10);
                            for (idx_t i = 0; i < TILE_SIZE; ++i)
                                tB2[i] = -tB2[i];
                        }
                    }

                    /* Diagonal block: btr == btc.
                     * Group B tile adjoints coincide with Group A tiles:
                     *   (tr00, tc11) adjoint = (tc11, tr00) = (tr11, tc00) -- Group A tile!
                     *   (tr01, tc11) adjoint = (tc11, tr01) = (tr11, tc01) -- Group A tile!
                     *   (tr10, tc11) adjoint = (tc11, tr10) = (tr11, tc10) -- Group A tile!
                     *
                     * Group A already negated these tiles. Negating Group B's adjoint
                     * would double-negate them. So skip Group B entirely on diagonal. */
                    /* else: btr == btc -- nothing to do for Group B */
                }
            }
        }
    }
}
