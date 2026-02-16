/**
 * @file gate_tiled_swap.c
 * @brief SWAP gate for tiled density matrices
 *
 * SWAP on qubits (lo, hi) permutes basis states by swapping the bit values at
 * positions lo and hi. For density matrix element rho[r,c]:
 *   rho'[perm(r), perm(c)] = rho[r,c]
 * where perm(x) exchanges bits lo and hi in x. SWAP is Hermitian and unitary.
 *
 * Since SWAP is a self-inverse permutation, every element is either a fixed point
 * (bits lo and hi have the same value in the index) or part of a 2-cycle. We only
 * need to swap each 2-cycle pair once.
 *
 * The SWAP permutation acts on 4x4 blocks indexed by the 2-bit pattern at positions
 * (lo, hi). Of the 16 element positions, 4 are fixed (00->00, 11->11) and 12 form
 * 6 swap pairs (see swap_packed for the derivation):
 *   Pair 1: (r01, c00) <-> (r10, c00)
 *   Pair 2: (r01, c01) <-> (r10, c10)
 *   Pair 3: (r11, c01) <-> (r11, c10)
 *   Pair 4: (r10, c01) <-> (r01, c10)  -- may cross the diagonal
 *   Pair 5: (r00, c01) <-> (r00, c10)  -- may cross the diagonal
 *   Pair 6: (r10, c11) <-> (r01, c11)  -- may cross the diagonal
 *
 * Three cases based on where lo and hi fall relative to LOG_TILE_DIM:
 *   Case 1: hi < LOG_TILE_DIM   -- both qubits within a tile (local coordinates)
 *   Case 2: lo < LOG_TILE_DIM, hi >= LOG_TILE_DIM -- lo within tile, hi across tiles
 *   Case 3: lo >= LOG_TILE_DIM  -- both qubits across tiles
 */

#include "gate_tiled.h"

void swap_tiled(state_t *state, const qubit_t q1, const qubit_t q2) {
    const gate_idx_t dim = (gate_idx_t)1 << state->qubits;
    const gate_idx_t n_tiles = (dim + TILE_DIM - 1) >> LOG_TILE_DIM;
    cplx_t * restrict data = state->data;

    const qubit_t lo = (q1 < q2) ? q1 : q2;
    const qubit_t hi = (q1 > q2) ? q1 : q2;

    if (hi < LOG_TILE_DIM) {
        /*
         * Case 1: Both qubits within a tile.
         *
         * Both bits lo and hi are in local coordinates (0..LOG_TILE_DIM-1).
         * The permutation only affects local row/col within each tile.
         * All 4 elements of each swap pair are in the same tile.
         * Tiles store full TILE_DIM x TILE_DIM elements, so no conjugation needed.
         *
         * For diagonal tiles: the swaps modify both the lower and upper triangle
         * positions. Since SWAP preserves Hermiticity and both conjugate positions
         * are swapped consistently (both are within the same tile), the result is
         * automatically Hermitian-consistent.
         */
        const gate_idx_t dim_tile = (dim < TILE_DIM) ? dim : TILE_DIM;
        const gate_idx_t n_base = dim_tile >> 2;
        const gate_idx_t incr_lo = (gate_idx_t)1 << lo;
        const gate_idx_t incr_hi = (gate_idx_t)1 << hi;

        #pragma omp parallel for schedule(static, 1) if(dim >= OMP_THRESHOLD)
        for (gate_idx_t tr = 0; tr < n_tiles; ++tr) {
            for (gate_idx_t tc = 0; tc <= tr; ++tc) {
                cplx_t *tile = data + tile_off(tr, tc);

                for (gate_idx_t br = 0; br < n_base; ++br) {
                    const gate_idx_t r00 = insertBits2_0(br, lo, hi);
                    const gate_idx_t r01 = r00 | incr_lo;
                    const gate_idx_t r10 = r00 | incr_hi;
                    const gate_idx_t r11 = r10 | incr_lo;

                    for (gate_idx_t bc = 0; bc < n_base; ++bc) {
                        const gate_idx_t c00 = insertBits2_0(bc, lo, hi);
                        const gate_idx_t c01 = c00 | incr_lo;
                        const gate_idx_t c10 = c00 | incr_hi;
                        const gate_idx_t c11 = c10 | incr_lo;

                        /* Pair 1: (r01,c00) <-> (r10,c00) */
                        cplx_t tmp = tile[r01 * TILE_DIM + c00];
                        tile[r01 * TILE_DIM + c00] = tile[r10 * TILE_DIM + c00];
                        tile[r10 * TILE_DIM + c00] = tmp;

                        /* Pair 2: (r01,c01) <-> (r10,c10) */
                        tmp = tile[r01 * TILE_DIM + c01];
                        tile[r01 * TILE_DIM + c01] = tile[r10 * TILE_DIM + c10];
                        tile[r10 * TILE_DIM + c10] = tmp;

                        /* Pair 3: (r11,c01) <-> (r11,c10) */
                        tmp = tile[r11 * TILE_DIM + c01];
                        tile[r11 * TILE_DIM + c01] = tile[r11 * TILE_DIM + c10];
                        tile[r11 * TILE_DIM + c10] = tmp;

                        /* Pair 4: (r10,c01) <-> (r01,c10) */
                        tmp = tile[r10 * TILE_DIM + c01];
                        tile[r10 * TILE_DIM + c01] = tile[r01 * TILE_DIM + c10];
                        tile[r01 * TILE_DIM + c10] = tmp;

                        /* Pair 5: (r00,c01) <-> (r00,c10) */
                        tmp = tile[r00 * TILE_DIM + c01];
                        tile[r00 * TILE_DIM + c01] = tile[r00 * TILE_DIM + c10];
                        tile[r00 * TILE_DIM + c10] = tmp;

                        /* Pair 6: (r01,c11) <-> (r10,c11) */
                        tmp = tile[r01 * TILE_DIM + c11];
                        tile[r01 * TILE_DIM + c11] = tile[r10 * TILE_DIM + c11];
                        tile[r10 * TILE_DIM + c11] = tmp;
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
         * (tile_bit = hi - LOG_TILE_DIM). The permutation swaps a local bit
         * with a tile-coordinate bit, so elements move between tiles.
         *
         * For row index r with tile-row tr and local-row lr:
         *   perm(r) swaps bit lo in lr with bit tile_bit in tr
         *
         * Non-trivial action: elements in (tr0, lr|lo=1) swap with
         * elements in (tr1, lr&~lo), and similarly for columns.
         *
         * We iterate tile-row pairs (tr0, tr1) and tile-col pairs (tc0, tc1)
         * based on tile_bit, and within each 4-tile group, iterate local
         * coordinate half-pairs based on lo.
         *
         * The 6 swap pairs, with row/col split into (tile, local):
         *   Pair 1: T(tr0,tc0)[lr1,lc0] <-> T(tr1,tc0)[lr0,lc0]
         *   Pair 2: T(tr0,tc0)[lr1,lc1] <-> T(tr1,tc1)[lr0,lc0]
         *   Pair 3: T(tr1,tc0)[lr1,lc1] <-> T(tr1,tc1)[lr1,lc0]
         *   Pair 4: T(tr1,tc0)[lr0,lc1] <-> T(tr0,tc1)[lr1,lc0]
         *   Pair 5: T(tr0,tc0)[lr0,lc1] <-> T(tr0,tc1)[lr0,lc0]
         *   Pair 6: T(tr1,tc1)[lr0,lc1] <-> T(tr0,tc1)[lr1,lc1]
         *
         * Pairs 1-3: both tiles always in lower triangle.
         * Pairs 4-6: T(tr0,tc1) may be upper-tri, accessed via T(tc1,tr0)+conj.
         *
         * Diagonal tile-block note (is_diag, i.e., tr0==tc0, tr1==tc1):
         *   T(tr0,tr0) and T(tr1,tr1) are diagonal tiles where elem(lr,lc)
         *   and elem(lc,lr) are conjugates. Swaps modify some positions in
         *   these tiles. To avoid reading stale upper-triangle values, we
         *   use dtile_read/dtile_write which canonicalize through the
         *   lower-triangle position. After all swaps, we mirror the upper
         *   triangle of diagonal tiles to restore consistency.
         */
        const gate_idx_t half_nt = n_tiles >> 1;
        const qubit_t tile_bit = hi - LOG_TILE_DIM;
        const gate_idx_t incr_tile = (gate_idx_t)1 << tile_bit;
        const gate_idx_t tile_dim = (dim < TILE_DIM) ? dim : TILE_DIM;
        const gate_idx_t half_td = tile_dim >> 1;
        const gate_idx_t incr_lo = (gate_idx_t)1 << lo;

        #pragma omp parallel for schedule(static, 1) if(dim >= OMP_THRESHOLD)
        for (gate_idx_t btr = 0; btr < half_nt; ++btr) {
            const gate_idx_t tr0 = insertBit0(btr, tile_bit);
            const gate_idx_t tr1 = tr0 | incr_tile;

            for (gate_idx_t btc = 0; btc <= btr; ++btc) {
                const gate_idx_t tc0 = insertBit0(btc, tile_bit);
                const gate_idx_t tc1 = tc0 | incr_tile;

                /* Tile pointers: T00, T11, T10 always in lower triangle */
                cplx_t * restrict t00 = data + tile_off(tr0, tc0);
                cplx_t * restrict t10 = data + tile_off(tr1, tc0);
                cplx_t * restrict t11 = data + tile_off(tr1, tc1);
                
                if (tr0 > tc1) {
                    cplx_t * restrict t01 = data + tile_off(tr0, tc1);

                    for (gate_idx_t blr = 0; blr < half_td; ++blr) {
                        const gate_idx_t lr0 = insertBit0(blr, lo);
                        const gate_idx_t lr1 = lr0 | incr_lo;

                        for (gate_idx_t blc = 0; blc < half_td; ++blc) {
                            const gate_idx_t lc0 = insertBit0(blc, lo);
                            const gate_idx_t lc1 = lc0 | incr_lo;

                            /* |hi=0,lo=1><hi=0,lo=0| <--> |hi=1,lo=0><hi=0,lo=0| */
                            cplx_t tmp = t00[lr1 * TILE_DIM + lc0];
                            t00[lr1 * TILE_DIM + lc0] = t10[lr0 * TILE_DIM + lc0];
                            t10[lr0 * TILE_DIM + lc0] = tmp;

                            /* |hi=0,lo=1><hi=0,lo=1| <--> |hi=1,lo=0><hi=1,lo=0| */
                            tmp = t00[lr1 * TILE_DIM + lc1];
                            t00[lr1 * TILE_DIM + lc1] = t11[lr0 * TILE_DIM + lc0];
                            t11[lr0 * TILE_DIM + lc0] = tmp;
 
                            /* |hi=1,lo=1><hi=0,lo=1| <--> |hi=1,lo=1><hi=1,lo=0| */
                            tmp = t10[lr1 * TILE_DIM + lc1];
                            t10[lr1 * TILE_DIM + lc1] = t11[lr1 * TILE_DIM + lc0];
                            t11[lr1 * TILE_DIM + lc0] = tmp;                           

                            /* |hi=1,lo=0><hi=0,lo=1| <--> |hi=0,lo=1><hi=1,lo=0| */
                            tmp = t10[lr0 * TILE_DIM + lc1];
                            t10[lr0 * TILE_DIM + lc1] = t01[lr1 * TILE_DIM + lc0];
                            t01[lr1 * TILE_DIM + lc0] = tmp; 

                            /* |hi=1,lo=0><hi=1,lo=1| <--> |hi=0,lo=1><hi=1,lo=1| */
                            tmp = t11[lr0 * TILE_DIM + lc1];
                            t11[lr0 * TILE_DIM + lc1] = t01[lr1 * TILE_DIM + lc1];
                            t01[lr1 * TILE_DIM + lc1] = tmp; 

                            /* |hi=0,lo=0><hi=0,lo=1| <--> |hi=0,lo=0><hi=1,lo=0| */
                            tmp = t00[lr0 * TILE_DIM + lc1];
                            t00[lr0 * TILE_DIM + lc1] = t01[lr0 * TILE_DIM + lc0];
                            t01[lr0 * TILE_DIM + lc0] = tmp; 
                        }
                    }
                }
                else if (btr != btc) {
                    cplx_t * restrict t01 = data + tile_off(tc1, tr0);

                    for (gate_idx_t blr = 0; blr < half_td; ++blr) {
                        const gate_idx_t lr0 = insertBit0(blr, lo);
                        const gate_idx_t lr1 = lr0 | incr_lo;

                        for (gate_idx_t blc = 0; blc < half_td; ++blc) {
                            const gate_idx_t lc0 = insertBit0(blc, lo);
                            const gate_idx_t lc1 = lc0 | incr_lo;

                            /* |hi=0,lo=1><hi=0,lo=0| <--> |hi=1,lo=0><hi=0,lo=0| */
                            cplx_t tmp = t00[lr1 * TILE_DIM + lc0];
                            t00[lr1 * TILE_DIM + lc0] = t10[lr0 * TILE_DIM + lc0];
                            t10[lr0 * TILE_DIM + lc0] = tmp;

                            /* |hi=0,lo=1><hi=0,lo=1| <--> |hi=1,lo=0><hi=1,lo=0| */
                            tmp = t00[lr1 * TILE_DIM + lc1];
                            t00[lr1 * TILE_DIM + lc1] = t11[lr0 * TILE_DIM + lc0];
                            t11[lr0 * TILE_DIM + lc0] = tmp;
 
                            /* |hi=1,lo=1><hi=0,lo=1| <--> |hi=1,lo=1><hi=1,lo=0| */
                            tmp = t10[lr1 * TILE_DIM + lc1];
                            t10[lr1 * TILE_DIM + lc1] = t11[lr1 * TILE_DIM + lc0];
                            t11[lr1 * TILE_DIM + lc0] = tmp;                           

                            /* |hi=1,lo=0><hi=0,lo=1| <--> |hi=0,lo=1><hi=1,lo=0| */
                            tmp = t10[lr0 * TILE_DIM + lc1];
                            t10[lr0 * TILE_DIM + lc1] = conj(t01[lc0 * TILE_DIM + lr1]);
                            t01[lc0 * TILE_DIM + lr1] = conj(tmp);

                            /* |hi=1,lo=0><hi=1,lo=1| <--> |hi=0,lo=1><hi=1,lo=1| */
                            tmp = t11[lr0 * TILE_DIM + lc1];
                            t11[lr0 * TILE_DIM + lc1] = conj(t01[lc1 * TILE_DIM + lr1]);
                            t01[lc1 * TILE_DIM + lr1] = conj(tmp);

                            /* |hi=0,lo=0><hi=0,lo=1| <--> |hi=0,lo=0><hi=1,lo=0| */
                            tmp = t00[lr0 * TILE_DIM + lc1];
                            t00[lr0 * TILE_DIM + lc1] = conj(t01[lc0 * TILE_DIM + lr0]);
                            t01[lc0 * TILE_DIM + lr0] = conj(tmp);
                        }
                    }
                }
                else {
                    for (gate_idx_t blr = 0; blr < half_td; ++blr) {
                        const gate_idx_t lr0 = insertBit0(blr, lo);
                        const gate_idx_t lr1 = lr0 | incr_lo;

                        for (gate_idx_t blc = 0; blc < half_td; ++blc) {
                            const gate_idx_t lc0 = insertBit0(blc, lo);
                            const gate_idx_t lc1 = lc0 | incr_lo;

                            /* |hi=0,lo=0><hi=0,lo=1| <--> |hi=0,lo=0><hi=1,lo=0|
                             * T(1,0) = conj(T(0,1))*/
                            t00[lr0 * TILE_DIM + lc1] = conj(t10[lc0 * TILE_DIM + lr0]);

                            /* |hi=1,lo=0><hi=1,lo=1| <--> |hi=0,lo=1><hi=1,lo=1| */
                            t11[lr0 * TILE_DIM + lc1] = conj(t10[lc1 * TILE_DIM + lr1]);
                        }
                    }

                    for (gate_idx_t blr = 0; blr < half_td; ++blr) {
                        const gate_idx_t lr0 = insertBit0(blr, lo);
                        const gate_idx_t lr1 = lr0 | incr_lo;

                        for (gate_idx_t blc = 0; blc < half_td; ++blc) {
                            const gate_idx_t lc0 = insertBit0(blc, lo);
                            const gate_idx_t lc1 = lc0 | incr_lo;

                            /* |hi=0,lo=1><hi=0,lo=0| <--> |hi=1,lo=0><hi=0,lo=0| */
                            cplx_t tmp = t00[lr1 * TILE_DIM + lc0];
                            t00[lr1 * TILE_DIM + lc0] = t10[lr0 * TILE_DIM + lc0];
                            t10[lr0 * TILE_DIM + lc0] = tmp;

                            /* |hi=0,lo=1><hi=0,lo=1| <--> |hi=1,lo=0><hi=1,lo=0| */
                            tmp = t00[lr1 * TILE_DIM + lc1];
                            t00[lr1 * TILE_DIM + lc1] = t11[lr0 * TILE_DIM + lc0];
                            t11[lr0 * TILE_DIM + lc0] = tmp;
 
                            /* |hi=1,lo=1><hi=0,lo=1| <--> |hi=1,lo=1><hi=1,lo=0| */
                            tmp = t10[lr1 * TILE_DIM + lc1];
                            t10[lr1 * TILE_DIM + lc1] = t11[lr1 * TILE_DIM + lc0];
                            t11[lr1 * TILE_DIM + lc0] = tmp;                           

                            /* |hi=1,lo=0><hi=0,lo=1| <--> |hi=0,lo=1><hi=1,lo=0| */
                            if (blc <= blr) {
                                tmp = t10[lr0 * TILE_DIM + lc1];
                                t10[lr0 * TILE_DIM + lc1] = conj(t10[lc0 * TILE_DIM + lr1]);
                                t10[lc0 * TILE_DIM + lr1] = conj(tmp);
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
         * Both bits lo and hi are in tile coordinates. The permutation only
         * affects tile indices, not local coordinates within tiles. For each
         * local position (lr, lc), the swap operates on tile-level quartets.
         *
         * tile_lo = lo - LOG_TILE_DIM, tile_hi = hi - LOG_TILE_DIM.
         *
         * We iterate tile-row quartets and tile-col quartets (lower triangle
         * at the quartet base level). Since local coordinates are unaffected,
         * the swaps operate on entire tiles (TILE_SIZE elements at a time).
         *
         * The 6 swap pairs (at the tile level):
         *   Pair 1: T(tr01,tc00) <-> T(tr10,tc00)
         *   Pair 2: T(tr01,tc01) <-> T(tr10,tc10)
         *   Pair 3: T(tr11,tc01) <-> T(tr11,tc10)
         *   Pair 4: T(tr10,tc01) <-> T(tr01,tc10)
         *   Pair 5: T(tr00,tc01) <-> T(tr00,tc10)
         *   Pair 6: T(tr10,tc11) <-> T(tr01,tc11)
         *
         * Tile lower-triangle status:
         *   Pairs 1-3: both tiles always in lower triangle.
         *   Pair 4: T(tr01,tc10) may be upper-tri if tr01 < tc10.
         *   Pairs 5,6: may involve upper-tri tiles depending on whether
         *     tr00 >= tc01 and tr00 >= tc10.
         */
        const qubit_t tile_lo = lo - LOG_TILE_DIM;
        const qubit_t tile_hi = hi - LOG_TILE_DIM;
        const gate_idx_t incr_lo = (gate_idx_t)1 << tile_lo;
        const gate_idx_t incr_hi = (gate_idx_t)1 << tile_hi;
        const gate_idx_t quarter_nt = n_tiles >> 2;

        #pragma omp parallel for schedule(static, 1) if(dim >= OMP_THRESHOLD)
        for (gate_idx_t btr = 0; btr < quarter_nt; ++btr) {
            const gate_idx_t tr00 = insertBits2_0(btr, tile_lo, tile_hi);
            const gate_idx_t tr01 = tr00 | incr_lo, tr10 = tr00 | incr_hi, tr11 = tr10 | incr_lo;

            for (gate_idx_t btc = 0; btc <= btr; ++btc) {
                const gate_idx_t tc00 = insertBits2_0(btc, tile_lo, tile_hi);
                const gate_idx_t tc01 = tc00 | incr_lo, tc10 = tc00 | incr_hi, tc11 = tc10 | incr_lo;
                cplx_t tmp;

                /* ---- Pair 1: T(tr01,tc00) <-> T(tr10,tc00) ----
                 * Both always lower-tri (tr01 >= tr00 >= tc00, tr10 > tr00 >= tc00). */
                {
                    cplx_t * restrict tA = data + tile_off(tr01, tc00);
                    cplx_t * restrict tB = data + tile_off(tr10, tc00);
                    for (gate_idx_t i = 0; i < TILE_SIZE; ++i) {
                        tmp = tA[i]; tA[i] = tB[i]; tB[i] = tmp;
                    }
                }

                /* ---- Pair 2: T(tr01,tc01) <-> T(tr10,tc10) ----
                 * Both always lower-tri (differences equal tr00-tc00 >= 0). */
                {
                    cplx_t * restrict tA = data + tile_off(tr01, tc01);
                    cplx_t * restrict tB = data + tile_off(tr10, tc10);
                    for (gate_idx_t i = 0; i < TILE_SIZE; ++i) {
                        tmp = tA[i]; tA[i] = tB[i]; tB[i] = tmp;
                    }
                }

                /* ---- Pair 3: T(tr11,tc01) <-> T(tr11,tc10) ----
                 * Both always lower-tri. */
                {
                    cplx_t * restrict tA = data + tile_off(tr11, tc01);
                    cplx_t * restrict tB = data + tile_off(tr11, tc10);
                    for (gate_idx_t i = 0; i < TILE_SIZE; ++i) {
                        tmp = tA[i]; tA[i] = tB[i]; tB[i] = tmp;
                    }
                }

                if (tr00 > tc10) {
                    /* ---- Pair 4: T(tr10,tc01) <-> T(tr01,tc10) ----
                     * T(tr10,tc01): always lower-tri.
                     * T(tr01,tc10): lower-tri iff tr01 >= tc10. */
                    {
                        cplx_t * restrict tA = data + tile_off(tr10, tc01);
                        cplx_t * restrict tB = data + tile_off(tr01, tc10);
                        for (gate_idx_t i = 0; i < TILE_SIZE; ++i) {
                            tmp = tA[i]; tA[i] = tB[i]; tB[i] = tmp;
                        }
                    }

                    {
                        cplx_t * restrict tA = data + tile_off(tr10, tc11);
                        cplx_t * restrict tB = data + tile_off(tr01, tc11);
                        for (gate_idx_t i = 0; i < TILE_SIZE; ++i) {
                            tmp = tA[i]; tA[i] = tB[i]; tB[i] = tmp;
                        }
                    }

                    {
                        cplx_t * restrict tA = data + tile_off(tr00, tc01);
                        cplx_t * restrict tB = data + tile_off(tr00, tc10);
                        for (gate_idx_t i = 0; i < TILE_SIZE; ++i) {
                            tmp = tA[i]; tA[i] = tB[i]; tB[i] = tmp;
                        }
                    }
                }
                else {
                    /* ---- Pair 4: T(tr10,tc01) <-> T(tr01,tc10) ----
                     * T(tr10,tc01): always lower-tri.
                     * T(tr01,tc10): lower-tri iff tr01 >= tc10. */
                    {
                        cplx_t * restrict tA = data + tile_off(tr10, tc01);
                        if (tr01 > tc10) {
                            cplx_t * restrict tB = data + tile_off(tr01, tc10);
                            for (gate_idx_t i = 0; i < TILE_SIZE; ++i) {
                                tmp = tA[i]; tA[i] = tB[i]; tB[i] = tmp;
                            }
                        }
                        else {
                            cplx_t * restrict tB = data + tile_off(tc10, tr01);
                            for (gate_idx_t lr = 0; lr < TILE_DIM; ++lr) {
                                for (gate_idx_t lc = 0; lc < TILE_DIM; ++lc) {
                                    tmp = tA[lr * TILE_DIM + lc];
                                    tA[lr * TILE_DIM + lc] = conj(tB[lc * TILE_DIM + lr]);
                                    tB[lc * TILE_DIM + lr] = conj(tmp);
                                }
                            }
                        }
                    }

                    if (tr00 > tc01) {
                        {
                            cplx_t * restrict tA = data + tile_off(tr00, tc01);
                            cplx_t * restrict tB = data + tile_off(tc10, tr00);
                            for (gate_idx_t lr = 0; lr < TILE_DIM; ++lr) {
                                for (gate_idx_t lc = 0; lc < TILE_DIM; ++lc) {
                                    tmp = tA[lr * TILE_DIM + lc];
                                    tA[lr * TILE_DIM + lc] = conj(tB[lc * TILE_DIM + lr]);
                                    tB[lc * TILE_DIM + lr] = conj(tmp);
                                }
                            }
                        }

                        {
                            cplx_t * restrict tA = data + tile_off(tr10, tc11);
                            cplx_t * restrict tB = data + tile_off(tc11, tr01);
                            for (gate_idx_t lr = 0; lr < TILE_DIM; ++lr) {
                                for (gate_idx_t lc = 0; lc < TILE_DIM; ++lc) {
                                    tmp = tA[lr * TILE_DIM + lc];
                                    tA[lr * TILE_DIM + lc] = conj(tB[lc * TILE_DIM + lr]);
                                    tB[lc * TILE_DIM + lr] = conj(tmp);
                                }
                            }
                        }
                    }
                    else if (btc != btr) {
                        {
                            cplx_t * restrict tA = data + tile_off(tc01, tr00);
                            cplx_t * restrict tB = data + tile_off(tc10, tr00);
                            for (gate_idx_t i = 0; i < TILE_SIZE; ++i) {
                                tmp = tA[i]; tA[i] = tB[i]; tB[i] = tmp;
                            }
                        }

                        {
                            cplx_t * restrict tA = data + tile_off(tc11, tr10);
                            cplx_t * restrict tB = data + tile_off(tc11, tr01);
                            for (gate_idx_t i = 0; i < TILE_SIZE; ++i) {
                                tmp = tA[i]; tA[i] = tB[i]; tB[i] = tmp;
                            }
                        }
                    }
                }
            }
        }
    }
}
