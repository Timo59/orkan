/**
 * @file gate_tiled_y.c
 * @brief Pauli-Y gate for tiled density matrices
 *
 * Y is a self-adjoint permutation-with-sign gate: ρ' = YρY.
 * Y = [[0,-i],[i,0]], so:
 *   ρ'(r0,c0) = ρ(r1,c1)      (diagonal swap, same as X)
 *   ρ'(r1,c1) = ρ(r0,c0)
 *   ρ'(r0,c1) = -ρ(r1,c0)     (off-diagonal negate-swap)
 *   ρ'(r1,c0) = -ρ(r0,c1)
 *
 * Two traversal modes:
 *   Within-tile (target < LOG_TILE_DIM):
 *     All 4 butterfly elements lie within the same tile.
 *   Cross-tile (target >= LOG_TILE_DIM):
 *     The 4 elements span 4 tiles: T(tr0,tc0), T(tr0,tc1), T(tr1,tc0), T(tr1,tc1).
 *     Three sub-cases based on lower-triangle storage:
 *       tr0 > tc1:  all 4 tiles stored directly
 *       tr0 == tc0 (diagonal): t00/t11 are diagonal tiles, t01 accessed via T(tc1,tr0)
 *       tr0 < tc1 (off-diag):  t01 accessed via T(tc1,tr0) with conjugation
 */

#include "gate_tiled.h"

void y_tiled(state_t *state, const qubit_t target) {
    const idx_t dim = (idx_t)1 << state->qubits;
    const idx_t n_tiles = (dim + TILE_DIM - (idx_t)1) >> LOG_TILE_DIM;
    const idx_t incr = (idx_t)1 << target;
    cplx_t * restrict data = state->data;

    if (target < LOG_TILE_DIM) {
        const idx_t dim_tile = dim < TILE_DIM ? dim : TILE_DIM;
        const idx_t n_base = dim_tile >> 1;

        #pragma omp parallel for schedule(static, 1) if(dim >= OMP_THRESHOLD)
        for (idx_t tr = 0; tr < n_tiles; ++tr) {
            for (idx_t tc = 0; tc <= tr; ++tc) {
                idx_t offset_tile = tile_off(tr, tc);

                for (idx_t br = 0; br < n_base; ++br) {
                    const idx_t r0 = insertBit0(br, target);
                    const idx_t r1 = r0 | incr;

                    const idx_t offset_r0 = r0 * TILE_DIM + offset_tile;
                    const idx_t offset_r1 = r1 * TILE_DIM + offset_tile;

                    for (idx_t bc = 0; bc < n_base; ++bc) {
                        const idx_t c0 = insertBit0(bc, target);
                        const idx_t c1 = c0 | incr;

                        /* swap ρ(r0,c0) <-> ρ(r1,c1) */
                        cplx_t tmp = data[c0 + offset_r0];
                        data[c0 + offset_r0] = data[c1 + offset_r1];
                        data[c1 + offset_r1] = tmp;

                        /* negate-swap ρ(r0,c1) <-> ρ(r1,c0) */
                        tmp = data[c0 + offset_r1];
                        data[c0 + offset_r1] = -data[c1 + offset_r0];
                        data[c1 + offset_r0] = -tmp;
                    }
                }
            }
        }
    }

    else {
        const idx_t n_base = n_tiles >> 1;
        const qubit_t target_tile = target - LOG_TILE_DIM;
        const idx_t incr_tile = (idx_t)1 << target_tile;

        #pragma omp parallel for schedule(static, 1) if(dim >= OMP_THRESHOLD)
        for (idx_t btr = 0; btr < n_base; ++btr) {
            const idx_t tr0 = insertBit0(btr, target_tile);
            const idx_t tr1 = tr0 | incr_tile;

            for (idx_t btc = 0; btc <= btr; ++btc) {
                const idx_t tc0 = insertBit0(btc, target_tile);
                const idx_t tc1 = tc0 | incr_tile;

                const idx_t offset_t00 = tile_off(tr0, tc0);
                const idx_t offset_t10 = tile_off(tr1, tc0);
                const idx_t offset_t11 = tile_off(tr1, tc1);

                if (tr0 > tc1) {
                    /* All 4 tiles stored directly */
                    const idx_t offset_t01 = tile_off(tr0, tc1);

                    for (idx_t lr = 0; lr < TILE_DIM; ++lr) {
                        cplx_t * restrict row00 = data + lr * TILE_DIM + offset_t00;
                        cplx_t * restrict row01 = data + lr * TILE_DIM + offset_t01;
                        cplx_t * restrict row10 = data + lr * TILE_DIM + offset_t10;
                        cplx_t * restrict row11 = data + lr * TILE_DIM + offset_t11;

                        #pragma clang loop vectorize(enable)
                        for (idx_t lc = 0; lc < TILE_DIM; ++lc) {
                            /* swap ρ(r0,c0) <-> ρ(r1,c1) */
                            cplx_t tmp = row00[lc];
                            row00[lc] = row11[lc];
                            row11[lc] = tmp;

                            /* negate-swap ρ(r0,c1) <-> ρ(r1,c0) */
                            tmp = row01[lc];
                            row01[lc] = -row10[lc];
                            row10[lc] = -tmp;
                        }
                    }
                }
                else if (btr != btc) {
                    /* Off-diagonal: t01 accessed via T(tc1,tr0) with conjugation.
                     * stored_t01[lc,lr] = conj(ρ[r0,c1]).
                     * ρ'[r1,c0] = -ρ[r0,c1] = -conj(stored)
                     * new stored = conj(ρ'[r0,c1]) = conj(-ρ[r1,c0]) = -conj(ρ[r1,c0]) */
                    const idx_t offset_t01 = tile_off(tc1, tr0);

                    for (idx_t lr = 0; lr < TILE_DIM; ++lr) {
                        const idx_t idx00 = lr * TILE_DIM + offset_t00;
                        const idx_t idx10 = lr * TILE_DIM + offset_t10;
                        const idx_t idx11 = lr * TILE_DIM + offset_t11;

                        for (idx_t lc = 0; lc < TILE_DIM; ++lc) {
                            cplx_t tmp = data[idx00 + lc];
                            data[idx00 + lc] = data[idx11 + lc];
                            data[idx11 + lc] = tmp;

                            tmp = data[lr + lc * TILE_DIM + offset_t01];
                            data[lr + lc * TILE_DIM + offset_t01] = -conj(data[idx10 + lc]);
                            data[idx10 + lc] = -conj(tmp);
                        }
                    }
                }
                else {
                    /* Diagonal: btr == btc, so tr0==tc0, tr1==tc1.
                     * t00/t11 are diagonal tiles; t01 = T(tc1,tr0) = T(tr1,tc0) = t10.
                     * On the diagonal element (lr==lc): ρ[r0,c1] = conj(ρ[r1,c0]),
                     * so ρ'[r1,c0] = -conj(ρ[r1,c0]) → t10[lr,lr] = -conj(t10[lr,lr]). */

                    for (idx_t lr = 0; lr < TILE_DIM; ++lr) {
                        cplx_t * restrict row00 = data + lr * TILE_DIM + offset_t00;
                        cplx_t * restrict row11 = data + lr * TILE_DIM + offset_t11;
                        const idx_t idx10 = lr * TILE_DIM + offset_t10;

                        #pragma clang loop vectorize(enable)
                        for (idx_t lc = 0; lc < TILE_DIM; ++lc) {
                            cplx_t tmp = row00[lc];
                            row00[lc] = row11[lc];
                            row11[lc] = tmp;
                        }

                        for (idx_t lc = 0; lc <= lr; ++lc) {
                            cplx_t tmp = data[idx10 + lc];
                            data[idx10 + lc] = -conj(data[lr + lc * TILE_DIM + offset_t10]);
                            data[lr + lc * TILE_DIM + offset_t10] = -conj(tmp);
                        }
                    }
                }
            }
        }
    }
}
