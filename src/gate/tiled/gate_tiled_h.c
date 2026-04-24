/**
 * @file gate_tiled_h.c
 * @brief Hadamard gate for tiled density matrices
 *
 * H = (1/sqrt(2)) [[1, 1], [1, -1]].  Self-adjoint (H = H†), so ρ' = HρH.
 *
 * Conjugation on each 2x2 block of the density matrix:
 *   ρ'_00 = 0.5 * (ρ00 + ρ01 + ρ10 + ρ11)
 *   ρ'_01 = 0.5 * (ρ00 - ρ01 + ρ10 - ρ11)
 *   ρ'_10 = 0.5 * (ρ00 + ρ01 - ρ10 - ρ11)
 *   ρ'_11 = 0.5 * (ρ00 - ρ01 - ρ10 + ρ11)
 *
 * All 4 elements mix -- this is a non-permutation, non-diagonal gate.
 *
 * Two traversal modes:
 *   Within-tile (target < LOG_TILE_DIM):
 *     All 4 butterfly elements lie within the same tile.
 *   Cross-tile (target >= LOG_TILE_DIM):
 *     The 4 elements span 4 tiles: T(tr0,tc0), T(tr0,tc1), T(tr1,tc0), T(tr1,tc1).
 *     Three sub-cases based on lower-triangle storage:
 *       tr0 > tc1:  all 4 tiles stored directly
 *       btr != btc (off-diag): T(tr0,tc1) accessed via T(tc1,tr0) with conjugation
 *       btr == btc (diagonal):  t00/t11 are diagonal tiles, t01 shares storage with t10
 */

#include "gate_tiled.h"

void h_tiled(state_t *state, const qubit_t target) {
    const idx_t dim = (idx_t)1 << state->qubits;
    const idx_t n_tiles = (dim + TILE_DIM - (idx_t)1) >> LOG_TILE_DIM;
    const idx_t dim_tile = dim < TILE_DIM ? dim : TILE_DIM;
    const idx_t incr = (idx_t)1 << target;
    cplx_t * restrict data = state->data;

    if (target < LOG_TILE_DIM) {
        /*
         * Within-tile path: all 4 butterfly elements are in the same tile.
         * Tiles store full TILE_DIM x TILE_DIM data (including diagonal tiles),
         * so all elements are directly accessible without conjugation logic.
         */
        const idx_t n_base = dim_tile >> 1;

        #pragma omp parallel for schedule(static, 1) if(dim >= OMP_THRESHOLD_MIXED)
        for (idx_t tr = 0; tr < n_tiles; ++tr) {
            for (idx_t tc = 0; tc <= tr; ++tc) {
                idx_t offset_tile = (tr * (tr + 1) / 2 + tc) * TILE_SIZE;

                for (idx_t br = 0; br < n_base; ++br) {
                    const idx_t r0 = insertBit0(br, target);
                    const idx_t r1 = r0 | incr;

                    const idx_t offset_r0 = r0 * TILE_DIM + offset_tile;
                    const idx_t offset_r1 = r1 * TILE_DIM + offset_tile;

                    for (idx_t bc = 0; bc < n_base; ++bc) {
                        const idx_t c0 = insertBit0(bc, target);
                        const idx_t c1 = c0 | incr;

                        cplx_t v00 = data[c0 + offset_r0];
                        cplx_t v01 = data[c1 + offset_r0];
                        cplx_t v10 = data[c0 + offset_r1];
                        cplx_t v11 = data[c1 + offset_r1];

                        data[c0 + offset_r0] = 0.5 * (v00 + v01 + v10 + v11);
                        data[c1 + offset_r0] = 0.5 * (v00 - v01 + v10 - v11);
                        data[c0 + offset_r1] = 0.5 * (v00 + v01 - v10 - v11);
                        data[c1 + offset_r1] = 0.5 * (v00 - v01 - v10 + v11);
                    }
                }
            }
        }
    }

    else {
        /*
         * Cross-tile path: the 4 elements span 4 tiles.
         * Enumerate tile-row pairs (tr0, tr1) and tile-col pairs (tc0, tc1)
         * where bit target_tile differs.
         */
        const idx_t n_base = n_tiles >> 1;
        const qubit_t target_tile = target - LOG_TILE_DIM;
        const idx_t incr_tile = (idx_t)1 << target_tile;

        #pragma omp parallel for schedule(static, 1) if(dim >= OMP_THRESHOLD_MIXED)
        for (idx_t btr = 0; btr < n_base; ++btr) {
            const idx_t tr0 = insertBit0(btr, target_tile);
            const idx_t tr1 = tr0 | incr_tile;

            for (idx_t btc = 0; btc <= btr; ++btc) {
                const idx_t tc0 = insertBit0(btc, target_tile);
                const idx_t tc1 = tc0 | incr_tile;

                const idx_t offset_t00 = (tr0 * (tr0 + 1) / 2 + tc0) * TILE_SIZE;
                const idx_t offset_t10 = (tr1 * (tr1 + 1) / 2 + tc0) * TILE_SIZE;
                const idx_t offset_t11 = (tr1 * (tr1 + 1) / 2 + tc1) * TILE_SIZE;

                if (tr0 > tc1) {
                    /*
                     * Sub-case a: all 4 tiles stored directly in lower triangle.
                     * T(tr0,tc1) is stored as-is since tr0 > tc1.
                     */
                    const idx_t offset_t01 = (tr0 * (tr0 + 1) / 2 + tc1) * TILE_SIZE;

                    for (idx_t lr = 0; lr < TILE_DIM; ++lr) {
                        cplx_t * restrict row00 = data + lr * TILE_DIM + offset_t00;
                        cplx_t * restrict row01 = data + lr * TILE_DIM + offset_t01;
                        cplx_t * restrict row10 = data + lr * TILE_DIM + offset_t10;
                        cplx_t * restrict row11 = data + lr * TILE_DIM + offset_t11;

                        #pragma clang loop vectorize(enable)
                        for (idx_t lc = 0; lc < TILE_DIM; ++lc) {
                            cplx_t v00 = row00[lc];
                            cplx_t v01 = row01[lc];
                            cplx_t v10 = row10[lc];
                            cplx_t v11 = row11[lc];

                            row00[lc] = 0.5 * (v00 + v01 + v10 + v11);
                            row01[lc] = 0.5 * (v00 - v01 + v10 - v11);
                            row10[lc] = 0.5 * (v00 + v01 - v10 - v11);
                            row11[lc] = 0.5 * (v00 - v01 - v10 + v11);
                        }
                    }
                }
                else if (btr != btc) {
                    /*
                     * Sub-case b: off-diagonal (btr != btc), tr0 <= tc1.
                     * T(tr0,tc1) not in lower triangle; access via T(tc1,tr0).
                     * stored[lc,lr] = conj(rho01[lr,lc]).
                     *   Read:  rho01 = conj(stored[lc,lr])
                     *   Write: stored[lc,lr] = conj(new_rho01)
                     */
                    const idx_t offset_t01 = (tc1 * (tc1 + 1) / 2 + tr0) * TILE_SIZE;

                    for (idx_t lr = 0; lr < TILE_DIM; ++lr) {
                        const idx_t idx00 = lr * TILE_DIM + offset_t00;
                        const idx_t idx10 = lr * TILE_DIM + offset_t10;
                        const idx_t idx11 = lr * TILE_DIM + offset_t11;

                        for (idx_t lc = 0; lc < TILE_DIM; ++lc) {
                            cplx_t v00 = data[idx00 + lc];
                            cplx_t v01 = conj(data[lr + lc * TILE_DIM + offset_t01]);
                            cplx_t v10 = data[idx10 + lc];
                            cplx_t v11 = data[idx11 + lc];

                            data[idx00 + lc] = 0.5 * (v00 + v01 + v10 + v11);
                            data[idx11 + lc] = 0.5 * (v00 - v01 - v10 + v11);
                            data[idx10 + lc] = 0.5 * (v00 + v01 - v10 - v11);

                            cplx_t new01 = 0.5 * (v00 - v01 + v10 - v11);
                            data[lr + lc * TILE_DIM + offset_t01] = conj(new01);
                        }
                    }
                }
                else {
                    /*
                     * Sub-case c: diagonal (btr == btc), so tr0==tc0, tr1==tc1.
                     * t00 and t11 are diagonal tiles (only lower triangle valid).
                     * t01 = T(tc1,tr0) = T(tr1,tc0) = t10 storage (same memory).
                     * lower_01 = 0 (tr0 < tc1 on diagonal).
                     *
                     * Iterate only lower triangle within-tile: lc from 0 to lr.
                     * After computing new t00/t11, mirror upper triangle.
                     *
                     * For t01 access: a01 = t10 + elem_off(lc,lr).
                     * On diagonal (lr==lc): a01 == a10, so the last write wins.
                     * The macro writes *a10 = new10 then *a01 = new10 (lower_01=0),
                     * so final value is new10 (no conjugation on diag_block).
                     */

                    for (idx_t lr = 0; lr < TILE_DIM; ++lr) {
                        const idx_t idx00 = lr * TILE_DIM + offset_t00;
                        const idx_t idx10 = lr * TILE_DIM + offset_t10;
                        const idx_t idx11 = lr * TILE_DIM + offset_t11;

                        for (idx_t lc = 0; lc <= lr; ++lc) {
                            cplx_t v00 = data[idx00 + lc];
                            cplx_t v10 = data[idx10 + lc];
                            cplx_t v11 = data[idx11 + lc];

                            cplx_t v01;
                            if (lr == lc) {
                                /* Diagonal within-tile element: rho01 = conj(rho10)
                                 * by Hermiticity of the density matrix. */
                                v01 = conj(v10);
                            } else {
                                /* Off-diagonal: t01 stored as T(tc1,tr0) = T(tr1,tc0)
                                 * = t10 storage. Stored at elem_off(lc,lr) and
                                 * stored value = conj(rho01). */
                                v01 = conj(data[lc * TILE_DIM + lr + offset_t10]);
                            }

                            cplx_t n00 = 0.5 * (v00 + v01 + v10 + v11);
                            cplx_t n01 = 0.5 * (v00 - v01 + v10 - v11);
                            cplx_t n10 = 0.5 * (v00 + v01 - v10 - v11);
                            cplx_t n11 = 0.5 * (v00 - v01 - v10 + v11);

                            /* Write t00 and t11 (diagonal tiles): lower triangle + mirror */
                            data[idx00 + lc] = n00;
                            data[idx11 + lc] = n11;
                            if (lr != lc) {
                                data[lc * TILE_DIM + lr + offset_t00] = conj(n00);
                                data[lc * TILE_DIM + lr + offset_t11] = conj(n11);
                            }

                            /* Write t10 and t01.
                             * t10[lr,lc] gets n10. t01 is stored at t10[lc,lr].
                             * When lr==lc: a10 and a01 alias the same address.
                             *   Macro behavior: write n10 then overwrite with n10
                             *   (lower_01=0, diag_block=1 => *a01 = new10).
                             *   Final value = n10.
                             * When lr!=lc: t10[lc,lr] = conj(n01)
                             *   (lower_01=0, diag_block=0 => stored = conj(new01)). */
                            data[idx10 + lc] = n10;
                            if (lr != lc) {
                                data[lc * TILE_DIM + lr + offset_t10] = conj(n01);
                            }
                            /* When lr==lc: no separate t01 write needed;
                             * a01 aliases a10, final value n10 is already correct. */
                        }
                    }
                }
            }
        }
    }
}
