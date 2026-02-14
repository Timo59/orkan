/**
 * @file gate_tiled_sdg.c
 * @brief S-dagger gate for tiled density matrices
 *
 * Sdg = diag(1, -i).  Diagonal gate: rho' = Sdg rho Sdg^dagger.
 *
 * Only off-diagonal density matrix elements change:
 *   rho'_00 = rho_00   (unchanged)
 *   rho'_11 = rho_11   (unchanged)
 *   rho'_10 = -i * rho_10  (multiply by -i)
 *   rho'_01 = i * rho_01   (multiply by i)
 *
 * Multiply by i:  CMPLX(-cimag(z), creal(z))
 * Multiply by -i: CMPLX(cimag(z), -creal(z))
 *
 * Two traversal modes:
 *   Within-tile (target < LOG_TILE_DIM):
 *     Both affected elements lie within the same tile.
 *   Cross-tile (target >= LOG_TILE_DIM):
 *     The off-diagonal elements are in tiles T(tr1,tc0) and T(tr0,tc1).
 *     T(tr0,tc0) and T(tr1,tc1) are unchanged.
 *     Three sub-cases for T(tr0,tc1) based on lower-triangle storage:
 *       tr0 > tc1:  stored directly, apply *= i
 *       btr != btc (off-diag): stored as conj(rho01) in T(tc1,tr0),
 *                   new stored = conj(rho01 * i) = conj(rho01) * (-i) = stored * (-i)
 *       btr == btc (diagonal): t01 is conj of t10, skip (handled by Hermiticity)
 */

#include "gate_tiled.h"

void sdg_tiled(state_t *state, const qubit_t target) {
    const gate_idx_t dim = (gate_idx_t)1 << state->qubits;
    const gate_idx_t n_tiles = (dim + TILE_DIM - (gate_idx_t)1) >> LOG_TILE_DIM;
    const gate_idx_t incr = (gate_idx_t)1 << target;
    cplx_t * restrict data = state->data;

    if (target < LOG_TILE_DIM) {
        /*
         * Within-tile path: both rho01 and rho10 elements are in the same tile.
         * Tiles store full TILE_DIM x TILE_DIM data, so all elements are
         * directly accessible with no conjugation logic needed.
         * Only a01 and a10 are modified; a00 and a11 are unchanged.
         */
        const gate_idx_t dim_tile = dim < TILE_DIM ? dim : TILE_DIM;
        const gate_idx_t n_base = dim_tile >> 1;

        #pragma omp parallel for schedule(static, 1) if(dim >= OMP_THRESHOLD)
        for (gate_idx_t tr = 0; tr < n_tiles; ++tr) {
            for (gate_idx_t tc = 0; tc <= tr; ++tc) {
                gate_idx_t offset_tile = tile_off(tr, tc);

                for (gate_idx_t br = 0; br < n_base; ++br) {
                    const gate_idx_t r0 = insertBit0(br, target);
                    const gate_idx_t r1 = r0 | incr;

                    const gate_idx_t offset_r0 = r0 * TILE_DIM + offset_tile;
                    const gate_idx_t offset_r1 = r1 * TILE_DIM + offset_tile;

                    for (gate_idx_t bc = 0; bc < n_base; ++bc) {
                        const gate_idx_t c0 = insertBit0(bc, target);
                        const gate_idx_t c1 = c0 | incr;

                        /* rho01 *= i */
                        cplx_t *p01 = data + c1 + offset_r0;
                        *p01 = CMPLX(-cimag(*p01), creal(*p01));

                        /* rho10 *= -i */
                        cplx_t *p10 = data + c0 + offset_r1;
                        *p10 = CMPLX(cimag(*p10), -creal(*p10));
                    }
                }
            }
        }
    }

    else {
        /*
         * Cross-tile path: rho10 is in T(tr1,tc0) and rho01 is in T(tr0,tc1).
         * T(tr0,tc0) and T(tr1,tc1) are unchanged (diagonal elements).
         * Three sub-cases handle T(tr0,tc1) storage.
         */
        const gate_idx_t n_base = n_tiles >> 1;
        const qubit_t target_tile = target - LOG_TILE_DIM;
        const gate_idx_t incr_tile = (gate_idx_t)1 << target_tile;

        #pragma omp parallel for schedule(static, 1) if(dim >= OMP_THRESHOLD)
        for (gate_idx_t btr = 0; btr < n_base; ++btr) {
            const gate_idx_t tr0 = insertBit0(btr, target_tile);
            const gate_idx_t tr1 = tr0 | incr_tile;

            for (gate_idx_t btc = 0; btc <= btr; ++btc) {
                const gate_idx_t tc0 = insertBit0(btc, target_tile);
                const gate_idx_t tc1 = tc0 | incr_tile;

                /* T(tr1,tc0): always in lower triangle -- rho10 *= -i */
                const gate_idx_t offset_t10 = tile_off(tr1, tc0);

                for (gate_idx_t lr = 0; lr < TILE_DIM; ++lr) {
                    cplx_t * restrict row10 = data + lr * TILE_DIM + offset_t10;
                    #pragma clang loop vectorize(enable)
                    for (gate_idx_t lc = 0; lc < TILE_DIM; ++lc) {
                        row10[lc] = CMPLX(cimag(row10[lc]), -creal(row10[lc]));
                    }
                }

                /* T(tr0,tc1): rho01 *= i, but skip on diagonal blocks */
                if (btr != btc) {
                    if (tr0 > tc1) {
                        /* stored directly as T(tr0,tc1): apply *= i */
                        const gate_idx_t offset_t01 = tile_off(tr0, tc1);

                        for (gate_idx_t lr = 0; lr < TILE_DIM; ++lr) {
                            cplx_t * restrict row01 = data + lr * TILE_DIM + offset_t01;
                            #pragma clang loop vectorize(enable)
                            for (gate_idx_t lc = 0; lc < TILE_DIM; ++lc) {
                                row01[lc] = CMPLX(-cimag(row01[lc]), creal(row01[lc]));
                            }
                        }
                    } else {
                        /* stored as T(tc1,tr0): stored = conj(rho01)
                         * new stored = conj(rho01 * i) = conj(rho01) * (-i) = stored * (-i) */
                        const gate_idx_t offset_t01 = tile_off(tc1, tr0);

                        for (gate_idx_t lr = 0; lr < TILE_DIM; ++lr) {
                            cplx_t * restrict row01 = data + lr * TILE_DIM + offset_t01;
                            #pragma clang loop vectorize(enable)
                            for (gate_idx_t lc = 0; lc < TILE_DIM; ++lc) {
                                row01[lc] = CMPLX(cimag(row01[lc]), -creal(row01[lc]));
                            }
                        }
                    }
                }
                /* When btr == btc (diagonal block): t01's stored value is
                 * conj(t10), so applying *(-i) to t10 implicitly applies *(i) to
                 * the logical t01.  No separate write needed. */
            }
        }
    }
}
