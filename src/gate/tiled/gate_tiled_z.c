/**
 * @file gate_tiled_z.c
 * @brief Pauli-Z gate for tiled density matrices
 *
 * Z is a diagonal gate (self-adjoint): ρ' = ZρZ.
 * Z = diag(1, -1), so ρ'[r,c] = (-1)^(r_t XOR c_t) * ρ[r,c]
 * where r_t, c_t are the target bits of r, c.
 *
 * Only off-diagonal elements (where target bit differs) are negated.
 * Diagonal elements (r_t == c_t) are unchanged.
 *
 * Two traversal modes:
 *   Within-tile (target < LOG_TILE_DIM):
 *     Both affected elements lie within the same tile.
 *   Cross-tile (target >= LOG_TILE_DIM):
 *     The off-diagonal elements are in tiles T(tr1,tc0) and T(tr0,tc1).
 *     T(tr0,tc0) and T(tr1,tc1) are unchanged.
 */

#include "gate_tiled.h"

void z_tiled(state_t *state, const qubit_t target) {
    const idx_t dim = (idx_t)1 << state->qubits;
    const idx_t n_tiles = (dim + TILE_DIM - (idx_t)1) >> LOG_TILE_DIM;
    const idx_t incr = (idx_t)1 << target;
    cplx_t * restrict data = state->data;

    if (target < LOG_TILE_DIM) {
        const idx_t dim_tile = dim < TILE_DIM ? dim : TILE_DIM;
        const idx_t n_base = dim_tile >> 1;

        #pragma omp parallel for schedule(static, 1) if(dim >= OMP_THRESHOLD_MIXED)
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

                        /* negate ρ(r0,c1) and ρ(r1,c0) — target bit differs */
                        data[c1 + offset_r0] = -data[c1 + offset_r0];
                        data[c0 + offset_r1] = -data[c0 + offset_r1];
                    }
                }
            }
        }
    }

    else {
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

                /* T(tr1,tc0): always in lower triangle — negate all elements */
                const idx_t offset_t10 = tile_off(tr1, tc0);

                for (idx_t lr = 0; lr < TILE_DIM; ++lr) {
                    cplx_t * restrict row10 = data + lr * TILE_DIM + offset_t10;
                    #pragma clang loop vectorize(enable)
                    for (idx_t lc = 0; lc < TILE_DIM; ++lc) {
                        row10[lc] = -row10[lc];
                    }
                }

                /* T(tr0,tc1): negate, but skip on diagonal blocks (a01 = conj(a10)) */
                if (btr != btc) {
                    if (tr0 > tc1) {
                        /* stored directly as T(tr0,tc1) */
                        const idx_t offset_t01 = tile_off(tr0, tc1);

                        for (idx_t lr = 0; lr < TILE_DIM; ++lr) {
                            cplx_t * restrict row01 = data + lr * TILE_DIM + offset_t01;
                            #pragma clang loop vectorize(enable)
                            for (idx_t lc = 0; lc < TILE_DIM; ++lc) {
                                row01[lc] = -row01[lc];
                            }
                        }
                    } else {
                        /* stored as T(tc1,tr0) — negate stored conjugate
                         * -conj(z) == conj(-z), so just negate in place */
                        const idx_t offset_t01 = tile_off(tc1, tr0);

                        for (idx_t lr = 0; lr < TILE_DIM; ++lr) {
                            cplx_t * restrict row01 = data + lr * TILE_DIM + offset_t01;
                            #pragma clang loop vectorize(enable)
                            for (idx_t lc = 0; lc < TILE_DIM; ++lc) {
                                row01[lc] = -row01[lc];
                            }
                        }
                    }
                }
                /* When btr == btc (diagonal block): t01's stored value is
                 * conj(t10), so negating t10 implicitly negates t01.
                 * No separate write needed. */
            }
        }
    }
}
