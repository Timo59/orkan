/**
 * @file gate_tiled_ry.c
 * @brief Ry rotation gate for tiled density matrices
 *
 * Ry(theta) = [[cos(theta/2), -sin(theta/2)], [sin(theta/2), cos(theta/2)]].
 * Mixing gate -- all 4 quadrants of each 2x2 block participate.
 *
 * With c=cos(theta/2), s=sin(theta/2), c2=c*c, s2=s*s, cs=c*s:
 *   rho'_00 = c2*rho00 + s2*rho11 - cs*(rho01 + rho10)
 *   rho'_11 = s2*rho00 + c2*rho11 + cs*(rho01 + rho10)
 *   rho'_01 = cs*(rho00 - rho11) + c2*rho01 - s2*rho10
 *   rho'_10 = cs*(rho00 - rho11) + c2*rho10 - s2*rho01
 *
 * Unlike Rx, all cs coefficients are real -- no i*cs terms appear.
 * This means the kernel uses cs*sum and cs*dd directly (real scalar
 * times complex), rather than the CMPLX rotation trick needed for Rx.
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
#include <math.h>

void ry_tiled(state_t *state, const qubit_t target, const double theta) {
    const double co = cos(theta / 2.0);
    const double si = sin(theta / 2.0);
    const double c2 = co * co, s2 = si * si, cs = co * si;

    const gate_idx_t dim = (gate_idx_t)1 << state->qubits;
    const gate_idx_t n_tiles = (dim + TILE_DIM - (gate_idx_t)1) >> LOG_TILE_DIM;
    const gate_idx_t dim_tile = dim < TILE_DIM ? dim : TILE_DIM;
    const gate_idx_t incr = (gate_idx_t)1 << target;
    cplx_t * restrict data = state->data;

    if (target < LOG_TILE_DIM) {
        /*
         * Within-tile path: all 4 butterfly elements are in the same tile.
         * Tiles store full TILE_DIM x TILE_DIM data (including diagonal tiles),
         * so all elements are directly accessible without conjugation logic.
         */
        const gate_idx_t n_base = dim_tile >> 1;

        #pragma omp parallel for schedule(static, 1) if(dim >= OMP_THRESHOLD)
        for (gate_idx_t tr = 0; tr < n_tiles; ++tr) {
            for (gate_idx_t tc = 0; tc <= tr; ++tc) {
                gate_idx_t offset_tile = (tr * (tr + 1) / 2 + tc) * TILE_SIZE;

                for (gate_idx_t br = 0; br < n_base; ++br) {
                    const gate_idx_t r0 = insertBit0(br, target);
                    const gate_idx_t r1 = r0 | incr;

                    const gate_idx_t offset_r0 = r0 * TILE_DIM + offset_tile;
                    const gate_idx_t offset_r1 = r1 * TILE_DIM + offset_tile;

                    for (gate_idx_t bc = 0; bc < n_base; ++bc) {
                        const gate_idx_t c0 = insertBit0(bc, target);
                        const gate_idx_t c1 = c0 | incr;

                        cplx_t v00 = data[c0 + offset_r0];
                        cplx_t v01 = data[c1 + offset_r0];
                        cplx_t v10 = data[c0 + offset_r1];
                        cplx_t v11 = data[c1 + offset_r1];

                        cplx_t sum = v01 + v10;
                        cplx_t dd  = v00 - v11;

                        data[c0 + offset_r0] = c2 * v00 + s2 * v11 - cs * sum;
                        data[c1 + offset_r0] = cs * dd + c2 * v01 - s2 * v10;
                        data[c0 + offset_r1] = cs * dd + c2 * v10 - s2 * v01;
                        data[c1 + offset_r1] = s2 * v00 + c2 * v11 + cs * sum;
                    }
                }
            }
        }
    }

    else {
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

                const gate_idx_t offset_t00 = (tr0 * (tr0 + 1) / 2 + tc0) * TILE_SIZE;
                const gate_idx_t offset_t10 = (tr1 * (tr1 + 1) / 2 + tc0) * TILE_SIZE;
                const gate_idx_t offset_t11 = (tr1 * (tr1 + 1) / 2 + tc1) * TILE_SIZE;

                if (tr0 > tc1) {
                    const gate_idx_t offset_t01 = (tr0 * (tr0 + 1) / 2 + tc1) * TILE_SIZE;

                    for (gate_idx_t lr = 0; lr < TILE_DIM; ++lr) {
                        cplx_t * restrict row00 = data + lr * TILE_DIM + offset_t00;
                        cplx_t * restrict row01 = data + lr * TILE_DIM + offset_t01;
                        cplx_t * restrict row10 = data + lr * TILE_DIM + offset_t10;
                        cplx_t * restrict row11 = data + lr * TILE_DIM + offset_t11;

                        #pragma clang loop vectorize(enable)
                        for (gate_idx_t lc = 0; lc < TILE_DIM; ++lc) {
                            cplx_t v00 = row00[lc], v01 = row01[lc];
                            cplx_t v10 = row10[lc], v11 = row11[lc];

                            cplx_t sum = v01 + v10;
                            cplx_t dd  = v00 - v11;

                            row00[lc] = c2 * v00 + s2 * v11 - cs * sum;
                            row01[lc] = cs * dd + c2 * v01 - s2 * v10;
                            row10[lc] = cs * dd + c2 * v10 - s2 * v01;
                            row11[lc] = s2 * v00 + c2 * v11 + cs * sum;
                        }
                    }
                }
                else if (btr != btc) {
                    const gate_idx_t offset_t01 = (tc1 * (tc1 + 1) / 2 + tr0) * TILE_SIZE;

                    for (gate_idx_t lr = 0; lr < TILE_DIM; ++lr) {
                        const gate_idx_t idx00 = lr * TILE_DIM + offset_t00;
                        const gate_idx_t idx10 = lr * TILE_DIM + offset_t10;
                        const gate_idx_t idx11 = lr * TILE_DIM + offset_t11;

                        for (gate_idx_t lc = 0; lc < TILE_DIM; ++lc) {
                            cplx_t v00 = data[idx00 + lc];
                            cplx_t v01 = conj(data[lr + lc * TILE_DIM + offset_t01]);
                            cplx_t v10 = data[idx10 + lc];
                            cplx_t v11 = data[idx11 + lc];

                            cplx_t sum = v01 + v10;
                            cplx_t dd  = v00 - v11;

                            data[idx00 + lc] = c2 * v00 + s2 * v11 - cs * sum;
                            data[idx11 + lc] = s2 * v00 + c2 * v11 + cs * sum;
                            data[idx10 + lc] = cs * dd + c2 * v10 - s2 * v01;

                            cplx_t new01 = cs * dd + c2 * v01 - s2 * v10;
                            data[lr + lc * TILE_DIM + offset_t01] = conj(new01);
                        }
                    }
                }
                else {
                    for (gate_idx_t lr = 0; lr < TILE_DIM; ++lr) {
                        const gate_idx_t idx00 = lr * TILE_DIM + offset_t00;
                        const gate_idx_t idx10 = lr * TILE_DIM + offset_t10;
                        const gate_idx_t idx11 = lr * TILE_DIM + offset_t11;

                        for (gate_idx_t lc = 0; lc <= lr; ++lc) {
                            cplx_t v00 = data[idx00 + lc];
                            cplx_t v10 = data[idx10 + lc];
                            cplx_t v11 = data[idx11 + lc];

                            cplx_t v01;
                            if (lr == lc) {
                                v01 = conj(v10);
                            } else {
                                v01 = conj(data[lc * TILE_DIM + lr + offset_t10]);
                            }

                            cplx_t sum = v01 + v10;
                            cplx_t dd  = v00 - v11;

                            cplx_t n00 = c2 * v00 + s2 * v11 - cs * sum;
                            cplx_t n01 = cs * dd + c2 * v01 - s2 * v10;
                            cplx_t n10 = cs * dd + c2 * v10 - s2 * v01;
                            cplx_t n11 = s2 * v00 + c2 * v11 + cs * sum;

                            data[idx00 + lc] = n00;
                            data[idx11 + lc] = n11;
                            if (lr != lc) {
                                data[lc * TILE_DIM + lr + offset_t00] = conj(n00);
                                data[lc * TILE_DIM + lr + offset_t11] = conj(n11);
                            }

                            data[idx10 + lc] = n10;
                            if (lr != lc) {
                                data[lc * TILE_DIM + lr + offset_t10] = conj(n01);
                            }
                        }
                    }
                }
            }
        }
    }
}
