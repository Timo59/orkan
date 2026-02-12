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
        const gate_idx_t n_tiles = (dim + TILE_DIM - 1) >> LOG_TILE_DIM;
        const gate_idx_t tile_dim = (dim < TILE_DIM) ? dim : TILE_DIM;
        const gate_idx_t quarter_td = tile_dim >> 2;
        const gate_idx_t incr_lo = (gate_idx_t)1 << lo;
        const gate_idx_t incr_hi = (gate_idx_t)1 << hi;

        #pragma omp parallel for schedule(dynamic) if(dim >= OMP_THRESHOLD)
        for (gate_idx_t tr = 0; tr < n_tiles; ++tr) {
            for (gate_idx_t tc = 0; tc <= tr; ++tc) {
                cplx_t *tile = data + tile_off(tr, tc);

                for (gate_idx_t br = 0; br < quarter_td; ++br) {
                    gate_idx_t lr00 = insertBits2_0(br, lo, hi);
                    gate_idx_t lr01 = lr00 | incr_lo;
                    gate_idx_t lr10 = lr00 | incr_hi;
                    gate_idx_t lr11 = lr10 | incr_lo;

                    for (gate_idx_t bc = 0; bc < quarter_td; ++bc) {
                        gate_idx_t lc00 = insertBits2_0(bc, lo, hi);
                        gate_idx_t lc01 = lc00 | incr_lo;
                        gate_idx_t lc10 = lc00 | incr_hi;
                        gate_idx_t lc11 = lc10 | incr_lo;

                        cplx_t tmp;

                        /* Pair 1: (r01,c00) <-> (r10,c00) */
                        tmp = tile[elem_off(lr01, lc00)];
                        tile[elem_off(lr01, lc00)] = tile[elem_off(lr10, lc00)];
                        tile[elem_off(lr10, lc00)] = tmp;

                        /* Pair 2: (r01,c01) <-> (r10,c10) */
                        tmp = tile[elem_off(lr01, lc01)];
                        tile[elem_off(lr01, lc01)] = tile[elem_off(lr10, lc10)];
                        tile[elem_off(lr10, lc10)] = tmp;

                        /* Pair 3: (r11,c01) <-> (r11,c10) */
                        tmp = tile[elem_off(lr11, lc01)];
                        tile[elem_off(lr11, lc01)] = tile[elem_off(lr11, lc10)];
                        tile[elem_off(lr11, lc10)] = tmp;

                        /* Pair 4: (r10,c01) <-> (r01,c10) */
                        tmp = tile[elem_off(lr10, lc01)];
                        tile[elem_off(lr10, lc01)] = tile[elem_off(lr01, lc10)];
                        tile[elem_off(lr01, lc10)] = tmp;

                        /* Pair 5: (r00,c01) <-> (r00,c10) */
                        tmp = tile[elem_off(lr00, lc01)];
                        tile[elem_off(lr00, lc01)] = tile[elem_off(lr00, lc10)];
                        tile[elem_off(lr00, lc10)] = tmp;

                        /* Pair 6: (r01,c11) <-> (r10,c11) */
                        tmp = tile[elem_off(lr01, lc11)];
                        tile[elem_off(lr01, lc11)] = tile[elem_off(lr10, lc11)];
                        tile[elem_off(lr10, lc11)] = tmp;
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
        const qubit_t tile_bit = hi - LOG_TILE_DIM;
        const gate_idx_t tile_stride = (gate_idx_t)1 << tile_bit;
        const gate_idx_t n_tiles = (dim + TILE_DIM - 1) >> LOG_TILE_DIM;
        const gate_idx_t half_nt = n_tiles >> 1;
        const gate_idx_t tile_dim = (dim < TILE_DIM) ? dim : TILE_DIM;
        const gate_idx_t half_td = tile_dim >> 1;
        const gate_idx_t incr_lo = (gate_idx_t)1 << lo;

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

                for (gate_idx_t b_lr = 0; b_lr < half_td; ++b_lr) {
                    gate_idx_t lr0 = insertBit0(b_lr, lo);
                    gate_idx_t lr1 = lr0 | incr_lo;

                    for (gate_idx_t b_lc = 0; b_lc < half_td; ++b_lc) {
                        gate_idx_t lc0 = insertBit0(b_lc, lo);
                        gate_idx_t lc1 = lc0 | incr_lo;

                        if (is_diag) {
                            /* t01 == t10 (same memory), and t00/t11 are diagonal
                             * tiles where dtile canonicalization creates further
                             * aliasing. Blocks (b_lr, b_lc) and (b_lc, b_lr)
                             * write identical values to the same physical
                             * locations, so process only b_lc <= b_lr. */
                            if (b_lc > b_lr) break;
                            cplx_t p1a = dtile_read(t00, lr1, lc0);
                            cplx_t p2a = dtile_read(t00, lr1, lc1);
                            cplx_t p5a = dtile_read(t00, lr0, lc1);
                            cplx_t p2b = dtile_read(t11, lr0, lc0);
                            cplx_t p3b = dtile_read(t11, lr1, lc0);
                            cplx_t p6a = dtile_read(t11, lr0, lc1);
                            cplx_t p1b = t10[elem_off(lr0, lc0)];
                            cplx_t p3a = t10[elem_off(lr1, lc1)];
                            cplx_t p4a = t10[elem_off(lr0, lc1)];
                            cplx_t p4b = conj(t10[elem_off(lc0, lr1)]);
                            cplx_t p5b = conj(t10[elem_off(lc0, lr0)]);
                            cplx_t p6b = conj(t10[elem_off(lc1, lr1)]);

                            /* Pair 1 */ dtile_write(t00, lr1, lc0, p1b);
                                         t10[elem_off(lr0, lc0)] = p1a;
                            /* Pair 2 */ dtile_write(t00, lr1, lc1, p2b);
                                         dtile_write(t11, lr0, lc0, p2a);
                            /* Pair 3 */ t10[elem_off(lr1, lc1)] = p3b;
                                         dtile_write(t11, lr1, lc0, p3a);
                            /* Pair 4 */ t10[elem_off(lr0, lc1)] = p4b;
                                         t10[elem_off(lc0, lr1)] = conj(p4a);
                            /* Pair 5 */ dtile_write(t00, lr0, lc1, p5b);
                                         t10[elem_off(lc0, lr0)] = conj(p5a);
                            /* Pair 6 */ dtile_write(t11, lr0, lc1, p6b);
                                         t10[elem_off(lc1, lr1)] = conj(p6a);
                        } else {
                            cplx_t tmp;

                            /* Pair 1: T(tr0,tc0)[lr1,lc0] <-> T(tr1,tc0)[lr0,lc0] */
                            tmp = t00[elem_off(lr1, lc0)];
                            t00[elem_off(lr1, lc0)] = t10[elem_off(lr0, lc0)];
                            t10[elem_off(lr0, lc0)] = tmp;

                            /* Pair 2: T(tr0,tc0)[lr1,lc1] <-> T(tr1,tc1)[lr0,lc0] */
                            tmp = t00[elem_off(lr1, lc1)];
                            t00[elem_off(lr1, lc1)] = t11[elem_off(lr0, lc0)];
                            t11[elem_off(lr0, lc0)] = tmp;

                            /* Pair 3: T(tr1,tc0)[lr1,lc1] <-> T(tr1,tc1)[lr1,lc0] */
                            tmp = t10[elem_off(lr1, lc1)];
                            t10[elem_off(lr1, lc1)] = t11[elem_off(lr1, lc0)];
                            t11[elem_off(lr1, lc0)] = tmp;

                            /* Pair 4: T(tr1,tc0)[lr0,lc1] <-> T(tr0,tc1)[lr1,lc0] */
                            {
                                cplx_t val_a = t10[elem_off(lr0, lc1)];
                                cplx_t val_b;
                                if (lower_01) {
                                    val_b = t01[elem_off(lr1, lc0)];
                                    t10[elem_off(lr0, lc1)] = val_b;
                                    t01[elem_off(lr1, lc0)] = val_a;
                                } else {
                                    val_b = conj(t01[elem_off(lc0, lr1)]);
                                    t10[elem_off(lr0, lc1)] = val_b;
                                    t01[elem_off(lc0, lr1)] = conj(val_a);
                                }
                            }

                            /* Pair 5: T(tr0,tc0)[lr0,lc1] <-> T(tr0,tc1)[lr0,lc0] */
                            {
                                cplx_t val_a = t00[elem_off(lr0, lc1)];
                                cplx_t val_b;
                                if (lower_01) {
                                    val_b = t01[elem_off(lr0, lc0)];
                                    t01[elem_off(lr0, lc0)] = val_a;
                                } else {
                                    val_b = conj(t01[elem_off(lc0, lr0)]);
                                    t01[elem_off(lc0, lr0)] = conj(val_a);
                                }
                                t00[elem_off(lr0, lc1)] = val_b;
                            }

                            /* Pair 6: T(tr1,tc1)[lr0,lc1] <-> T(tr0,tc1)[lr1,lc1] */
                            {
                                cplx_t val_a = t11[elem_off(lr0, lc1)];
                                cplx_t val_b;
                                if (lower_01) {
                                    val_b = t01[elem_off(lr1, lc1)];
                                    t01[elem_off(lr1, lc1)] = val_a;
                                } else {
                                    val_b = conj(t01[elem_off(lc1, lr1)]);
                                    t01[elem_off(lc1, lr1)] = conj(val_a);
                                }
                                t11[elem_off(lr0, lc1)] = val_b;
                            }
                        }
                    }
                }

                /* Mirror upper triangle of diagonal tiles from lower triangle.
                 * When is_diag, T(tr0,tr0) and T(tr1,tr1) are diagonal tiles
                 * whose lower-triangle elements are now correct. Propagate to
                 * upper triangle for storage consistency. */
                if (is_diag) {
                    for (gate_idx_t lr = 0; lr < TILE_DIM; ++lr) {
                        for (gate_idx_t lc = lr + 1; lc < TILE_DIM; ++lc) {
                            t00[elem_off(lr, lc)] = conj(t00[elem_off(lc, lr)]);
                            t11[elem_off(lr, lc)] = conj(t11[elem_off(lc, lr)]);
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
        const gate_idx_t n_tiles = (dim + TILE_DIM - 1) >> LOG_TILE_DIM;
        const gate_idx_t quarter_nt = n_tiles >> 2;

        #pragma omp parallel for schedule(dynamic) if(dim >= OMP_THRESHOLD)
        for (gate_idx_t b_tr = 0; b_tr < quarter_nt; ++b_tr) {
            gate_idx_t tr00 = insertBits2_0(b_tr, tile_lo, tile_hi);
            gate_idx_t tr01 = tr00 | incr_lo;
            gate_idx_t tr10 = tr00 | incr_hi;
            gate_idx_t tr11 = tr10 | incr_lo;

            for (gate_idx_t b_tc = 0; b_tc <= b_tr; ++b_tc) {
                gate_idx_t tc00 = insertBits2_0(b_tc, tile_lo, tile_hi);
                gate_idx_t tc01 = tc00 | incr_lo;
                gate_idx_t tc10 = tc00 | incr_hi;
                gate_idx_t tc11 = tc10 | incr_lo;

                int is_diag = (b_tc == b_tr);

                /* ---- Pair 1: T(tr01,tc00) <-> T(tr10,tc00) ----
                 * Both always lower-tri (tr01 >= tr00 >= tc00, tr10 > tr00 >= tc00). */
                {
                    cplx_t *tA = data + tile_off(tr01, tc00);
                    cplx_t *tB = data + tile_off(tr10, tc00);
                    for (gate_idx_t i = 0; i < TILE_SIZE; ++i) {
                        cplx_t tmp = tA[i]; tA[i] = tB[i]; tB[i] = tmp;
                    }
                }

                /* ---- Pair 2: T(tr01,tc01) <-> T(tr10,tc10) ----
                 * Both always lower-tri (differences equal tr00-tc00 >= 0). */
                {
                    cplx_t *tA = data + tile_off(tr01, tc01);
                    cplx_t *tB = data + tile_off(tr10, tc10);
                    for (gate_idx_t i = 0; i < TILE_SIZE; ++i) {
                        cplx_t tmp = tA[i]; tA[i] = tB[i]; tB[i] = tmp;
                    }
                }

                /* ---- Pair 3: T(tr11,tc01) <-> T(tr11,tc10) ----
                 * Both always lower-tri. */
                {
                    cplx_t *tA = data + tile_off(tr11, tc01);
                    cplx_t *tB = data + tile_off(tr11, tc10);
                    for (gate_idx_t i = 0; i < TILE_SIZE; ++i) {
                        cplx_t tmp = tA[i]; tA[i] = tB[i]; tB[i] = tmp;
                    }
                }

                /* ---- Pair 4: T(tr10,tc01) <-> T(tr01,tc10) ----
                 * T(tr10,tc01): always lower-tri.
                 * T(tr01,tc10): lower-tri iff tr01 >= tc10.
                 *
                 * When is_diag (tr00==tc00): tr01=tc01, tc10=tr10, so this becomes
                 * T(tr10,tr01) <-> T(tr01,tr10). Since tr10 > tr01, T(tr01,tr10) is
                 * upper-tri, stored as T(tr10,tr01) = tA. So tA == tB and the "swap"
                 * is the in-place Hermitian transpose of a single tile. */
                {
                    cplx_t *tA = data + tile_off(tr10, tc01);
                    int lower_b = (tr01 >= tc10);
                    if (lower_b) {
                        cplx_t *tB = data + tile_off(tr01, tc10);
                        for (gate_idx_t i = 0; i < TILE_SIZE; ++i) {
                            cplx_t tmp = tA[i]; tA[i] = tB[i]; tB[i] = tmp;
                        }
                    } else if (is_diag) {
                        /* tA points to T(tr10, tr01). The logical swap maps
                         * rho[tr10+lr, tr01+lc] <-> rho[tr01+lr, tr10+lc].
                         * In storage: tA[lr,lc] <-> conj(tA[lc,lr]).
                         * Hermitian-transpose the tile in place. */
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
                        /* tA = T(tr10,tc01), tB = T(tc10,tr01). Distinct tiles.
                         * Logical: tA[lr,lc] <-> conj(tB[lc,lr]). */
                        cplx_t *tB = data + tile_off(tc10, tr01);
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

                /* ---- Pairs 5+6 ----
                 * Pair 5: T(tr00,tc01) <-> T(tr00,tc10)
                 * Pair 6: T(tr10,tc11) <-> T(tr01,tc11)
                 *
                 * Lower-tri conditions (using incr_hi > incr_lo):
                 *   T(tr00,tc01): iff tr00 >= tc01     [lo_01]
                 *   T(tr00,tc10): iff tr00 >= tc10     [lo_10]
                 *   T(tr10,tc11): iff tr10 >= tc11, equiv tr00 >= tc01  [lo_01]
                 *   T(tr01,tc11): iff tr01 >= tc11, equiv tr00 >= tc10  [lo_10]
                 *
                 * Three sub-cases for !is_diag:
                 * (a) lo_01 && lo_10: all 4 tiles lower-tri -> plain swap.
                 * (b) lo_01 && !lo_10: mixed -> conjugate access for 2 tiles.
                 * (c) !lo_01: all 4 upper-tri -> swap stored reps directly.
                 * For is_diag: tr00==tc00 implies tc01>tr00, so !lo_01 always.
                 * Handled separately after the if/else chain.
                 */
                {
                    int lo_01 = (tr00 >= tc01);
                    int lo_10 = (tr00 >= tc10);

                    if (lo_01 && lo_10) {
                        /* (a) All tiles in lower triangle: plain tile swaps */
                        cplx_t *t5a = data + tile_off(tr00, tc01);
                        cplx_t *t5b = data + tile_off(tr00, tc10);
                        for (gate_idx_t i = 0; i < TILE_SIZE; ++i) {
                            cplx_t tmp = t5a[i]; t5a[i] = t5b[i]; t5b[i] = tmp;
                        }
                        cplx_t *t6a = data + tile_off(tr10, tc11);
                        cplx_t *t6b = data + tile_off(tr01, tc11);
                        for (gate_idx_t i = 0; i < TILE_SIZE; ++i) {
                            cplx_t tmp = t6a[i]; t6a[i] = t6b[i]; t6b[i] = tmp;
                        }
                    }
                    else if (lo_01) {
                        /* (b) T(tr00,tc01) lower, T(tr00,tc10) upper -> T(tc10,tr00)
                         *     T(tr10,tc11) lower, T(tr01,tc11) upper -> T(tc11,tr01) */
                        cplx_t *t5a = data + tile_off(tr00, tc01);
                        cplx_t *t5b = data + tile_off(tc10, tr00);
                        for (gate_idx_t lr = 0; lr < TILE_DIM; ++lr) {
                            for (gate_idx_t lc = 0; lc < TILE_DIM; ++lc) {
                                cplx_t val_a = t5a[elem_off(lr, lc)];
                                cplx_t val_b = conj(t5b[elem_off(lc, lr)]);
                                t5a[elem_off(lr, lc)] = val_b;
                                t5b[elem_off(lc, lr)] = conj(val_a);
                            }
                        }
                        cplx_t *t6a = data + tile_off(tr10, tc11);
                        cplx_t *t6b = data + tile_off(tc11, tr01);
                        for (gate_idx_t lr = 0; lr < TILE_DIM; ++lr) {
                            for (gate_idx_t lc = 0; lc < TILE_DIM; ++lc) {
                                cplx_t val_a = t6a[elem_off(lr, lc)];
                                cplx_t val_b = conj(t6b[elem_off(lc, lr)]);
                                t6a[elem_off(lr, lc)] = val_b;
                                t6b[elem_off(lc, lr)] = conj(val_a);
                            }
                        }
                    }
                    else if (!is_diag) {
                        /* (c) All 4 tiles upper-tri. Swap their stored lower-tri
                         * representatives directly (conjugation on both sides
                         * cancels for a plain swap between transposed tiles). */
                        cplx_t *t5a = data + tile_off(tc01, tr00);
                        cplx_t *t5b = data + tile_off(tc10, tr00);
                        for (gate_idx_t i = 0; i < TILE_SIZE; ++i) {
                            cplx_t tmp = t5a[i]; t5a[i] = t5b[i]; t5b[i] = tmp;
                        }
                        cplx_t *t6a = data + tile_off(tc11, tr10);
                        cplx_t *t6b = data + tile_off(tc11, tr01);
                        for (gate_idx_t i = 0; i < TILE_SIZE; ++i) {
                            cplx_t tmp = t6a[i]; t6a[i] = t6b[i]; t6b[i] = tmp;
                        }
                    }
                    /* is_diag && !lo_01: When is_diag, tr==tc so Pairs 5+6
                     * operate on the same physical tiles as Pairs 1+3
                     * (the stored lower-tri representatives coincide).
                     * Executing them would undo the Pair 1+3 swaps. Skip. */
                }
            }
        }
    }
}
