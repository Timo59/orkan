/**
 * @file gate_tiled_rz.c
 * @brief Rz rotation gate for tiled density matrices
 *
 * Rz(theta) = diag(e^(-i*theta/2), e^(i*theta/2)).  Diagonal gate: rho' = Rz rho Rz+.
 *
 * Only off-diagonal density matrix elements change:
 *   rho'_00 = rho_00                  (unchanged)
 *   rho'_11 = rho_11                  (unchanged)
 *   rho'_10 *= e^(i*theta):
 *     given r = creal, i = cimag -> CMPLX(r*c - i*sn, i*c + r*sn)
 *   rho'_01 *= e^(-i*theta):
 *     given r = creal, i = cimag -> CMPLX(r*c + i*sn, i*c - r*sn)
 *
 * where c = cos(theta), sn = sin(theta).
 *
 * Two traversal modes:
 *   Within-tile (target < LOG_TILE_DIM):
 *     Both affected elements (a01, a10) lie within the same tile.
 *     a00 and a11 are skipped (diagonal gate).
 *   Cross-tile (target >= LOG_TILE_DIM):
 *     The off-diagonal elements are in tiles T(tr1,tc0) and T(tr0,tc1).
 *     T(tr0,tc0) and T(tr1,tc1) are unchanged (diagonal elements).
 *     Three sub-cases based on lower-triangle storage:
 *       tr0 > tc1:  t01 stored directly as T(tr0,tc1)
 *       btr != btc (off-diag): t01 accessed via T(tc1,tr0) with conjugation
 *       btr == btc (diagonal):  t01 shares storage with t10, skip a01
 */

#include "gate_tiled.h"

void rz_tiled(state_t *state, const qubit_t target, const double theta) {
    const double c = cos(theta), sn = sin(theta);

    const idx_t dim = (idx_t)1 << state->qubits;
    const idx_t n_tiles = (dim + TILE_DIM - (idx_t)1) >> LOG_TILE_DIM;
    const idx_t incr = (idx_t)1 << target;
    cplx_t * restrict data = state->data;

    if (target < LOG_TILE_DIM) {
        /*
         * Within-tile path: all 4 butterfly elements are in the same tile.
         * Diagonal gate: only a01 and a10 are modified; a00/a11 unchanged.
         */
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

                        /* rho01 *= e^(-i*theta) */
                        double r01 = creal(data[c1 + offset_r0]);
                        double i01 = cimag(data[c1 + offset_r0]);
                        data[c1 + offset_r0] = CMPLX(r01 * c + i01 * sn,
                                                     i01 * c - r01 * sn);

                        /* rho10 *= e^(i*theta) */
                        double r10 = creal(data[c0 + offset_r1]);
                        double i10 = cimag(data[c0 + offset_r1]);
                        data[c0 + offset_r1] = CMPLX(r10 * c - i10 * sn,
                                                     i10 * c + r10 * sn);
                    }
                }
            }
        }
    }

    else {
        /*
         * Cross-tile path: the off-diagonal elements span tiles T(tr1,tc0)
         * and T(tr0,tc1). Diagonal tiles T(tr0,tc0) and T(tr1,tc1) are unchanged.
         */
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

                /* T(tr1,tc0): always in lower triangle.
                 * Apply rho10 *= e^(i*theta) to all elements. */
                const idx_t offset_t10 = tile_off(tr1, tc0);

                for (idx_t lr = 0; lr < TILE_DIM; ++lr) {
                    cplx_t * restrict row10 = data + lr * TILE_DIM + offset_t10;
                    #pragma clang loop vectorize(enable)
                    for (idx_t lc = 0; lc < TILE_DIM; ++lc) {
                        double r = creal(row10[lc]);
                        double im = cimag(row10[lc]);
                        row10[lc] = CMPLX(r * c - im * sn,
                                          im * c + r * sn);
                    }
                }

                /* T(tr0,tc1): apply rho01 *= e^(-i*theta), but skip on diagonal
                 * blocks where t01 shares storage with t10 (conj pair). */
                if (btr != btc) {
                    if (tr0 > tc1) {
                        /* Sub-case a: t01 stored directly as T(tr0,tc1).
                         * Apply e^(-i*theta) directly. */
                        const idx_t offset_t01 = tile_off(tr0, tc1);

                        for (idx_t lr = 0; lr < TILE_DIM; ++lr) {
                            cplx_t * restrict row01 = data + lr * TILE_DIM + offset_t01;
                            #pragma clang loop vectorize(enable)
                            for (idx_t lc = 0; lc < TILE_DIM; ++lc) {
                                double r = creal(row01[lc]);
                                double im = cimag(row01[lc]);
                                row01[lc] = CMPLX(r * c + im * sn,
                                                  im * c - r * sn);
                            }
                        }
                    } else {
                        /* Sub-case b: off-diagonal, t01 stored as T(tc1,tr0).
                         * stored = conj(rho01).
                         * new rho01 = rho01 * e^(-i*theta)
                         * new stored = conj(rho01 * e^(-i*theta))
                         *            = conj(rho01) * e^(i*theta)
                         *            = stored * e^(i*theta)
                         * So: CMPLX(r*c - im*sn, im*c + r*sn) */
                        const idx_t offset_t01 = tile_off(tc1, tr0);

                        for (idx_t lr = 0; lr < TILE_DIM; ++lr) {
                            cplx_t * restrict row01 = data + lr * TILE_DIM + offset_t01;
                            #pragma clang loop vectorize(enable)
                            for (idx_t lc = 0; lc < TILE_DIM; ++lc) {
                                double r = creal(row01[lc]);
                                double im = cimag(row01[lc]);
                                row01[lc] = CMPLX(r * c - im * sn,
                                                  im * c + r * sn);
                            }
                        }
                    }
                }
                /* When btr == btc (diagonal block): t01's stored value is
                 * conj(t10). Modifying t10 implicitly updates t01 through
                 * the conjugation invariant. No separate write needed. */
            }
        }
    }
}
