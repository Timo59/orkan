/**
 * @file gate_tiled_cy.c
 * @brief CY (controlled-Y) gate for tiled density matrices
 *
 * CY has the same permutation as CX (control flips target), plus per-element
 * phase factors: phi(00)=1, phi(01)=1, phi(10)=-i, phi(11)=i.
 * rho'[a,b] = phi(a)*conj(phi(b))*rho[perm(a), perm(b)].
 *
 * The 6 swap pairs (same positions as CX) with their CY phase products:
 *   Pair A: (r10,c00) <-> (r11,c00)   -- phases (-i, +i), always lower-tri
 *   Pair E: (r10,c10) <-> (r11,c11)   -- plain swap (phases cancel), always lower-tri
 *   Pair B: (r11,c01) <-> (r10,c01)   -- phases (+i, -i), may cross diagonal
 *   Pair F: (r11,c10) <-> (r10,c11)   -- negate-swap (phase -1), may cross diagonal
 *   Pair C: (r00,c10) <-> (r00,c11)   -- phases (+i, -i), may cross diagonal
 *   Pair D: (r01,c10) <-> (r01,c11)   -- phases (+i, -i), may cross diagonal
 *
 * Phase arithmetic on complex value v:
 *   +i * v = CMPLX(-cimag(v), creal(v))       (lower-tri, direct)
 *   -i * v = CMPLX(cimag(v), -creal(v))       (lower-tri, direct)
 *
 * Upper-tri pitfall: stored s = conj(rho). Applying phase phi to rho gives
 * new stored = conj(phi) * s. So:
 *   "apply +i" on upper-tri stored s: conj(i)*s = -i*s = CMPLX(cimag(s), -creal(s))
 *   "apply -i" on upper-tri stored s: conj(-i)*s = +i*s = CMPLX(-cimag(s), creal(s))
 *   "apply -1" on upper-tri stored s: -s (same as lower-tri)
 *
 * Structure follows gate_tiled_cx.c with identical case/branch topology.
 * CX code labels map to CY spec pairs as follows (varies by case):
 *   Case 2: CX-A=A(-i,+i) CX-B=E(swap) CX-C=B(+i,-i) CX-D=F(neg) CX-E=D(+i,-i) CX-F=C(+i,-i)
 *   Case 3: CX-A=A(-i,+i) CX-B=E(swap) CX-C=B(+i,-i) CX-D=F(neg) CX-E=D(+i,-i) CX-F=C(+i,-i)
 *   Case 4: CX-A=A(-i,+i) CX-B=E(swap) CX-C=B(+i,-i) CX-D=F(neg) CX-E=D(+i,-i) CX-F=C(+i,-i)
 */

#include "gate_tiled.h"

void cy_tiled(state_t *state, const qubit_t control, const qubit_t target) {
    const gate_idx_t dim = (gate_idx_t)1 << state->qubits;
    const gate_idx_t dim_tile = (dim + TILE_DIM - 1) >> LOG_TILE_DIM;
    cplx_t * restrict data = state->data;

    if (control < LOG_TILE_DIM) {
        /* Case 1: Both control and target smaller than LOG_TILE_DIM
         *
         * All swap pairs are contained within the same tile.
         * All elements are lower-tri -> phases apply directly.
         */
        if (target < LOG_TILE_DIM) {
            const gate_idx_t n_tiles = dim_tile * (dim_tile + 1) / 2;
            const gate_idx_t n_base = (dim < TILE_DIM) ? dim >> 2 : TILE_DIM >> 2;
            const qubit_t lo = (control < target) ? control : target;
            const qubit_t hi = (control > target) ? control : target;
            const gate_idx_t incr_ctrl = (gate_idx_t)1 << control;
            const gate_idx_t incr_tgt = (gate_idx_t)1 << target;

            #pragma omp parallel for schedule(static) if(dim >= OMP_THRESHOLD)
            for (gate_idx_t i = 0; i < n_tiles; ++i) {
                cplx_t * restrict tile = data + i * TILE_SIZE;

                for (gate_idx_t blr = 0; blr < n_base; ++blr) {
                    const gate_idx_t r00 = insertBits2_0(blr, lo, hi);
                    const gate_idx_t r01 = r00 | incr_tgt, r10 = r00 | incr_ctrl, r11 = r10 | incr_tgt;

                    for (gate_idx_t blc = 0; blc < n_base; ++blc) {
                        const gate_idx_t c00 = insertBits2_0(blc, lo, hi);
                        const gate_idx_t c01 = c00 | incr_tgt, c10 = c00 | incr_ctrl, c11 = c10 | incr_tgt;

                        /* Pair A: (r10,c00) <-> (r11,c00); phases (-i, +i) */
                        cplx_t tmp = tile[elem_off(r10, c00)];
                        tile[elem_off(r10, c00)] = CMPLX( cimag(tile[elem_off(r11, c00)]),
                                                         -creal(tile[elem_off(r11, c00)]));
                        tile[elem_off(r11, c00)] = CMPLX(-cimag(tmp), creal(tmp));

                        /* Pair E: (r10,c10) <-> (r11,c11); plain swap */
                        tmp = tile[elem_off(r10, c10)];
                        tile[elem_off(r10, c10)] = tile[elem_off(r11, c11)];
                        tile[elem_off(r11, c11)] = tmp;

                        /* Pair B: (r11,c01) <-> (r10,c01); phases (+i, -i) */
                        tmp = tile[elem_off(r11, c01)];
                        tile[elem_off(r11, c01)] = CMPLX(-cimag(tile[elem_off(r10, c01)]),
                                                          creal(tile[elem_off(r10, c01)]));
                        tile[elem_off(r10, c01)] = CMPLX( cimag(tmp), -creal(tmp));

                        /* Pair F: (r11,c10) <-> (r10,c11); negate-swap */
                        tmp = tile[elem_off(r11, c10)];
                        tile[elem_off(r11, c10)] = -tile[elem_off(r10, c11)];
                        tile[elem_off(r10, c11)] = -tmp;

                        /* Pair D: (r01,c10) <-> (r01,c11); phases (+i, -i) */
                        tmp = tile[elem_off(r01, c10)];
                        tile[elem_off(r01, c10)] = CMPLX(-cimag(tile[elem_off(r01, c11)]),
                                                          creal(tile[elem_off(r01, c11)]));
                        tile[elem_off(r01, c11)] = CMPLX( cimag(tmp), -creal(tmp));

                        /* Pair C: (r00,c10) <-> (r00,c11); phases (+i, -i) */
                        tmp = tile[elem_off(r00, c10)];
                        tile[elem_off(r00, c10)] = CMPLX(-cimag(tile[elem_off(r00, c11)]),
                                                          creal(tile[elem_off(r00, c11)]));
                        tile[elem_off(r00, c11)] = CMPLX( cimag(tmp), -creal(tmp));
                    }
                }
            }
        }

        /* Case 2: Target >= LOG_TILE_DIM, control < LOG_TILE_DIM
         *
         * Swap pairs span four tiles: t00=(tr0,tc0), t10=(tr1,tc0),
         * t11=(tr1,tc1), t01=(tr0,tc1). The target bit distinguishes
         * tile rows/cols; the control bit distinguishes local rows/cols.
         */
        else {
            const gate_idx_t n_bt = dim_tile >> 1;
            const qubit_t tile_bit = target - LOG_TILE_DIM;
            const gate_idx_t incr_tile = (gate_idx_t)1 << tile_bit;
            const gate_idx_t n_bl = TILE_DIM >> 1;
            const gate_idx_t incr_local = (gate_idx_t)1 << control;

            #pragma omp parallel for schedule(static, 1) if(dim >= OMP_THRESHOLD)
            for (gate_idx_t btr = 0; btr < n_bt; ++btr) {
                const gate_idx_t tr0 = insertBit0(btr, tile_bit);
                const gate_idx_t tr1 = tr0 | incr_tile;

                for (gate_idx_t btc = 0; btc <= btr; ++btc) {
                    const gate_idx_t tc0 = insertBit0(btc, tile_bit);
                    const gate_idx_t tc1 = tc0 | incr_tile;

                    cplx_t * restrict t00 = data + tile_off(tr0, tc0);
                    cplx_t * restrict t10 = data + tile_off(tr1, tc0);
                    cplx_t * restrict t11 = data + tile_off(tr1, tc1);

                    /* Fast path: All tiles are in the lower triangle */
                    if (tr0 > tc1) {
                        cplx_t * restrict t01 = data + tile_off(tr0, tc1);

                        for (gate_idx_t blr = 0; blr < n_bl; ++blr) {
                            const gate_idx_t r0 = insertBit0(blr, control);
                            const gate_idx_t r1 = r0 | incr_local;

                            for (gate_idx_t blc = 0; blc < n_bl; ++blc) {
                                const gate_idx_t c0 = insertBit0(blc, control);
                                const gate_idx_t c1 = c0 | incr_local;

                                /* Pair A: phases (-i, +i) */
                                cplx_t tmp = t00[elem_off(r1, c0)];
                                t00[elem_off(r1, c0)] = CMPLX( cimag(t10[elem_off(r1, c0)]),
                                                               -creal(t10[elem_off(r1, c0)]));
                                t10[elem_off(r1, c0)] = CMPLX(-cimag(tmp), creal(tmp));

                                /* Pair E: plain swap */
                                tmp = t00[elem_off(r1, c1)];
                                t00[elem_off(r1, c1)] = t11[elem_off(r1, c1)];
                                t11[elem_off(r1, c1)] = tmp;

                                /* Pair D: phases (+i, -i) */
                                tmp = t10[elem_off(r0, c1)];
                                t10[elem_off(r0, c1)] = CMPLX(-cimag(t11[elem_off(r0, c1)]),
                                                                creal(t11[elem_off(r0, c1)]));
                                t11[elem_off(r0, c1)] = CMPLX( cimag(tmp), -creal(tmp));

                                /* Pair B: phases (+i, -i) */
                                tmp = t11[elem_off(r1, c0)];
                                t11[elem_off(r1, c0)] = CMPLX(-cimag(t01[elem_off(r1, c0)]),
                                                                creal(t01[elem_off(r1, c0)]));
                                t01[elem_off(r1, c0)] = CMPLX( cimag(tmp), -creal(tmp));

                                /* Pair F: negate-swap */
                                tmp = t10[elem_off(r1, c1)];
                                t10[elem_off(r1, c1)] = -t01[elem_off(r1, c1)];
                                t01[elem_off(r1, c1)] = -tmp;

                                /* Pair C: phases (+i, -i) */
                                tmp = t00[elem_off(r0, c1)];
                                t00[elem_off(r0, c1)] = CMPLX(-cimag(t01[elem_off(r0, c1)]),
                                                                creal(t01[elem_off(r0, c1)]));
                                t01[elem_off(r0, c1)] = CMPLX( cimag(tmp), -creal(tmp));
                            }
                        }
                    }

                    /* t01 is in upper triangle but its adjoint is in lower triangle */
                    else if (btr != btc) {
                        cplx_t * restrict t01 = data + tile_off(tc1, tr0);

                        for (gate_idx_t blr = 0; blr < n_bl; ++blr) {
                            const gate_idx_t r0 = insertBit0(blr, control);
                            const gate_idx_t r1 = r0 | incr_local;

                            for (gate_idx_t blc = 0; blc < n_bl; ++blc) {
                                const gate_idx_t c0 = insertBit0(blc, control);
                                const gate_idx_t c1 = c0 | incr_local;

                                /* Pair A: phases (-i, +i) — both lower-tri */
                                cplx_t tmp = t00[elem_off(r1, c0)];
                                t00[elem_off(r1, c0)] = CMPLX( cimag(t10[elem_off(r1, c0)]),
                                                               -creal(t10[elem_off(r1, c0)]));
                                t10[elem_off(r1, c0)] = CMPLX(-cimag(tmp), creal(tmp));

                                /* Pair E: plain swap — both lower-tri */
                                tmp = t00[elem_off(r1, c1)];
                                t00[elem_off(r1, c1)] = t11[elem_off(r1, c1)];
                                t11[elem_off(r1, c1)] = tmp;

                                /* Pair D: phases (+i, -i) — both lower-tri */
                                tmp = t10[elem_off(r0, c1)];
                                t10[elem_off(r0, c1)] = CMPLX(-cimag(t11[elem_off(r0, c1)]),
                                                                creal(t11[elem_off(r0, c1)]));
                                t11[elem_off(r0, c1)] = CMPLX( cimag(tmp), -creal(tmp));

                                /* Pair B: phases (+i, -i) — t11 lower-tri, t01 upper-tri stored
                                 * new rho[r11,c01] = +i * old rho[r10,c01]
                                 *   old rho[r10,c01] = conj(t01[c0,r1])
                                 *   +i * conj(s) = CMPLX(cimag(s), creal(s))
                                 * new stored at t01[c0,r1] = conj(-i * old rho[r11,c01])
                                 *   = conj(-i * tmp) -> CMPLX(cimag(tmp), creal(tmp))
                                 */
                                tmp = t11[elem_off(r1, c0)];
                                cplx_t s = t01[elem_off(c0, r1)];
                                t11[elem_off(r1, c0)] = CMPLX(cimag(s), creal(s));
                                t01[elem_off(c0, r1)] = CMPLX(cimag(tmp), creal(tmp));

                                /* Pair F: negate-swap — t10 lower-tri, t01 upper-tri stored
                                 * new rho[r11,c10] = -old rho[r10,c11] = -conj(t01[c1,r1])
                                 * new stored at t01[c1,r1] = conj(-old rho[r11,c10]) = -conj(tmp)
                                 */
                                tmp = t10[elem_off(r1, c1)];
                                t10[elem_off(r1, c1)] = -conj(t01[elem_off(c1, r1)]);
                                t01[elem_off(c1, r1)] = -conj(tmp);

                                /* Pair C: phases (+i, -i) — t00 lower-tri, t01 upper-tri stored
                                 * new rho[r00,c10] = +i * old rho[r00,c11]
                                 *   old rho[r00,c11] = conj(t01[c1,r0])
                                 *   +i * conj(s) = CMPLX(cimag(s), creal(s))
                                 * new stored at t01[c1,r0] = conj(-i * old rho[r00,c10])
                                 *   = conj(-i * tmp) -> CMPLX(cimag(tmp), creal(tmp))
                                 */
                                tmp = t00[elem_off(r0, c1)];
                                s = t01[elem_off(c1, r0)];
                                t00[elem_off(r0, c1)] = CMPLX(cimag(s), creal(s));
                                t01[elem_off(c1, r0)] = CMPLX(cimag(tmp), creal(tmp));
                            }
                        }
                    }

                    /* t01 is in upper triangle and its adjoint is t10 (btr == btc) */
                    else {
                        /* One-way writes for Pair B and Pair C first.
                         * These read from t10 (= t01's stored adjoint) before
                         * the second loop modifies t10.
                         *
                         * Pair B: new rho[r11,c01] = +i * conj(t10[c0,r1])
                         *   = CMPLX(cimag(t10[c0,r1]), creal(t10[c0,r1]))
                         * Pair C: new rho[r00,c10] = +i * conj(t10[c1,r0])
                         *   = CMPLX(cimag(t10[c1,r0]), creal(t10[c1,r0]))
                         */
                        for (gate_idx_t blr = 0; blr < n_bl; ++blr) {
                            const gate_idx_t r0 = insertBit0(blr, control);
                            const gate_idx_t r1 = r0 | incr_local;

                            for (gate_idx_t blc = 0; blc < n_bl; ++blc) {
                                const gate_idx_t c0 = insertBit0(blc, control);
                                const gate_idx_t c1 = c0 | incr_local;

                                /* Pair B: one-way write to t11 */
                                cplx_t s = t10[elem_off(c0, r1)];
                                t11[elem_off(r1, c0)] = CMPLX(cimag(s), creal(s));

                                /* Pair C: one-way write to t00 */
                                s = t10[elem_off(c1, r0)];
                                t00[elem_off(r0, c1)] = CMPLX(cimag(s), creal(s));
                            }
                        }

                        for (gate_idx_t blr = 0; blr < n_bl; ++blr) {
                            const gate_idx_t r0 = insertBit0(blr, control);
                            const gate_idx_t r1 = r0 | incr_local;

                            for (gate_idx_t blc = 0; blc < n_bl; ++blc) {
                                const gate_idx_t c0 = insertBit0(blc, control);
                                const gate_idx_t c1 = c0 | incr_local;

                                /* Pair A: phases (-i, +i) — both in t00/t10 lower-tri */
                                cplx_t tmp = t00[elem_off(r1, c0)];
                                t00[elem_off(r1, c0)] = CMPLX( cimag(t10[elem_off(r1, c0)]),
                                                               -creal(t10[elem_off(r1, c0)]));
                                t10[elem_off(r1, c0)] = CMPLX(-cimag(tmp), creal(tmp));

                                /* Pair E: plain swap — t00/t11 */
                                tmp = t00[elem_off(r1, c1)];
                                t00[elem_off(r1, c1)] = t11[elem_off(r1, c1)];
                                t11[elem_off(r1, c1)] = tmp;

                                /* Pair D: phases (+i, -i) — t10/t11 */
                                tmp = t10[elem_off(r0, c1)];
                                t10[elem_off(r0, c1)] = CMPLX(-cimag(t11[elem_off(r0, c1)]),
                                                                creal(t11[elem_off(r0, c1)]));
                                t11[elem_off(r0, c1)] = CMPLX( cimag(tmp), -creal(tmp));

                                /* Pair F: negate-swap — within t10 (t01's adjoint = t10)
                                 * new rho[r11,c10] = -rho[r10,c11] = -conj(t10[c1,r1])
                                 * new stored t10[c1,r1] = conj(-rho[r11,c10]) = -conj(tmp)
                                 * Only process lower triangle to avoid double-write.
                                 */
                                if (blc <= blr) {
                                    tmp = t10[elem_off(r1, c1)];
                                    t10[elem_off(r1, c1)] = -conj(t10[elem_off(c1, r1)]);
                                    t10[elem_off(c1, r1)] = -conj(tmp);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    else {
        /* Case 3: Control >= LOG_TILE_DIM, target < LOG_TILE_DIM
         *
         * The control bit distinguishes tile rows/cols; the target bit
         * distinguishes local rows/cols within a tile.
         * Pairs A, B (spec) are within t10; Pair E (spec) is within t11;
         * Pair F (spec) is within t11; Pairs C, D (spec) are in t01.
         */
        if (target < LOG_TILE_DIM) {
            const gate_idx_t n_bt = dim_tile >> 1;
            const gate_idx_t tile_bit = control - LOG_TILE_DIM;
            const gate_idx_t incr_tile = (gate_idx_t)1 << tile_bit;
            const gate_idx_t n_bl = TILE_DIM >> 1;
            const gate_idx_t incr_local = (gate_idx_t)1 << target;

            #pragma omp parallel for schedule(static, 1) if(dim >= OMP_THRESHOLD)
            for (gate_idx_t btr = 0; btr < n_bt; ++btr) {
                const gate_idx_t tr0 = insertBit0(btr, tile_bit);
                const gate_idx_t tr1 = tr0 | incr_tile;

                for (gate_idx_t btc = 0; btc <= btr; ++btc) {
                    const gate_idx_t tc0 = insertBit0(btc, tile_bit);
                    const gate_idx_t tc1 = tc0 | incr_tile;

                    cplx_t * restrict t10 = data + tile_off(tr1, tc0);
                    cplx_t * restrict t11 = data + tile_off(tr1, tc1);

                    if (tr0 > tc1) {
                        cplx_t * restrict t01 = data + tile_off(tr0, tc1);

                        for (gate_idx_t blr = 0; blr < n_bl; ++blr) {
                            const gate_idx_t lr0 = insertBit0(blr, target);
                            const gate_idx_t lr1 = lr0 | incr_local;

                            for (gate_idx_t blc = 0; blc < n_bl; ++blc) {
                                const gate_idx_t lc0 = insertBit0(blc, target);
                                const gate_idx_t lc1 = lc0 | incr_local;

                                /* Pair A: phases (-i, +i) — within t10 */
                                cplx_t tmp = t10[elem_off(lr0, lc0)];
                                t10[elem_off(lr0, lc0)] = CMPLX( cimag(t10[elem_off(lr1, lc0)]),
                                                                 -creal(t10[elem_off(lr1, lc0)]));
                                t10[elem_off(lr1, lc0)] = CMPLX(-cimag(tmp), creal(tmp));

                                /* Pair B: phases (+i, -i) — within t10 */
                                tmp = t10[elem_off(lr0, lc1)];
                                t10[elem_off(lr0, lc1)] = CMPLX( cimag(t10[elem_off(lr1, lc1)]),
                                                                 -creal(t10[elem_off(lr1, lc1)]));
                                t10[elem_off(lr1, lc1)] = CMPLX(-cimag(tmp), creal(tmp));

                                /* Pair E: plain swap — within t11 */
                                tmp = t11[elem_off(lr0, lc0)];
                                t11[elem_off(lr0, lc0)] = t11[elem_off(lr1, lc1)];
                                t11[elem_off(lr1, lc1)] = tmp;

                                /* Pair F: negate-swap — within t11 */
                                tmp = t11[elem_off(lr0, lc1)];
                                t11[elem_off(lr0, lc1)] = -t11[elem_off(lr1, lc0)];
                                t11[elem_off(lr1, lc0)] = -tmp;

                                /* Pair D: phases (+i, -i) — within t01 */
                                tmp = t01[elem_off(lr1, lc0)];
                                t01[elem_off(lr1, lc0)] = CMPLX(-cimag(t01[elem_off(lr1, lc1)]),
                                                                  creal(t01[elem_off(lr1, lc1)]));
                                t01[elem_off(lr1, lc1)] = CMPLX( cimag(tmp), -creal(tmp));

                                /* Pair C: phases (+i, -i) — within t01 */
                                tmp = t01[elem_off(lr0, lc0)];
                                t01[elem_off(lr0, lc0)] = CMPLX(-cimag(t01[elem_off(lr0, lc1)]),
                                                                  creal(t01[elem_off(lr0, lc1)]));
                                t01[elem_off(lr0, lc1)] = CMPLX( cimag(tmp), -creal(tmp));
                            }
                        }
                    }

                    else if (btr != btc) {
                        cplx_t * restrict t01 = data + tile_off(tc1, tr0);

                        for (gate_idx_t blr = 0; blr < n_bl; ++blr) {
                            const gate_idx_t lr0 = insertBit0(blr, target);
                            const gate_idx_t lr1 = lr0 | incr_local;

                            for (gate_idx_t blc = 0; blc < n_bl; ++blc) {
                                const gate_idx_t lc0 = insertBit0(blc, target);
                                const gate_idx_t lc1 = lc0 | incr_local;

                                /* Pair A: phases (-i, +i) — within t10 */
                                cplx_t tmp = t10[elem_off(lr0, lc0)];
                                t10[elem_off(lr0, lc0)] = CMPLX( cimag(t10[elem_off(lr1, lc0)]),
                                                                 -creal(t10[elem_off(lr1, lc0)]));
                                t10[elem_off(lr1, lc0)] = CMPLX(-cimag(tmp), creal(tmp));

                                /* Pair B: phases (+i, -i) — within t10 */
                                tmp = t10[elem_off(lr0, lc1)];
                                t10[elem_off(lr0, lc1)] = CMPLX( cimag(t10[elem_off(lr1, lc1)]),
                                                                 -creal(t10[elem_off(lr1, lc1)]));
                                t10[elem_off(lr1, lc1)] = CMPLX(-cimag(tmp), creal(tmp));

                                /* Pair E: plain swap — within t11 */
                                tmp = t11[elem_off(lr0, lc0)];
                                t11[elem_off(lr0, lc0)] = t11[elem_off(lr1, lc1)];
                                t11[elem_off(lr1, lc1)] = tmp;

                                /* Pair F: negate-swap — within t11 */
                                tmp = t11[elem_off(lr0, lc1)];
                                t11[elem_off(lr0, lc1)] = -t11[elem_off(lr1, lc0)];
                                t11[elem_off(lr1, lc0)] = -tmp;

                                /* Pair D: phases (+i, -i) — t01 is upper-tri stored
                                 * rho[r01,c10] stored as conj at t01[lc0,lr1]
                                 * rho[r01,c11] stored as conj at t01[lc1,lr1]
                                 * Apply conj(+i)=-i to first stored, conj(-i)=+i to second.
                                 */
                                tmp = t01[elem_off(lc0, lr1)];
                                cplx_t s = t01[elem_off(lc1, lr1)];
                                t01[elem_off(lc0, lr1)] = CMPLX( cimag(s), -creal(s));
                                t01[elem_off(lc1, lr1)] = CMPLX(-cimag(tmp), creal(tmp));

                                /* Pair C: phases (+i, -i) — t01 is upper-tri stored
                                 * rho[r00,c10] stored as conj at t01[lc0,lr0]
                                 * rho[r00,c11] stored as conj at t01[lc1,lr0]
                                 */
                                tmp = t01[elem_off(lc0, lr0)];
                                s = t01[elem_off(lc1, lr0)];
                                t01[elem_off(lc0, lr0)] = CMPLX( cimag(s), -creal(s));
                                t01[elem_off(lc1, lr0)] = CMPLX(-cimag(tmp), creal(tmp));
                            }
                        }
                    }

                    /* Diagonal: btr == btc, t01 adjoint = t10.
                     * Pairs D and C (spec) in the upper-tri t01 are automatically
                     * handled by Pairs A and B (spec) at transposed (blr,blc)
                     * iterations within t10, since the phases compose correctly
                     * through the Hermitian conjugate relationship.
                     */
                    else {
                        for (gate_idx_t blr = 0; blr < n_bl; ++blr) {
                            const gate_idx_t lr0 = insertBit0(blr, target);
                            const gate_idx_t lr1 = lr0 | incr_local;

                            for (gate_idx_t blc = 0; blc < n_bl; ++blc) {
                                const gate_idx_t lc0 = insertBit0(blc, target);
                                const gate_idx_t lc1 = lc0 | incr_local;

                                /* Pair A: phases (-i, +i) — within t10 */
                                cplx_t tmp = t10[elem_off(lr0, lc0)];
                                t10[elem_off(lr0, lc0)] = CMPLX( cimag(t10[elem_off(lr1, lc0)]),
                                                                 -creal(t10[elem_off(lr1, lc0)]));
                                t10[elem_off(lr1, lc0)] = CMPLX(-cimag(tmp), creal(tmp));

                                /* Pair B: phases (+i, -i) — within t10 */
                                tmp = t10[elem_off(lr0, lc1)];
                                t10[elem_off(lr0, lc1)] = CMPLX( cimag(t10[elem_off(lr1, lc1)]),
                                                                 -creal(t10[elem_off(lr1, lc1)]));
                                t10[elem_off(lr1, lc1)] = CMPLX(-cimag(tmp), creal(tmp));

                                /* Pair E: plain swap — within t11 */
                                tmp = t11[elem_off(lr0, lc0)];
                                t11[elem_off(lr0, lc0)] = t11[elem_off(lr1, lc1)];
                                t11[elem_off(lr1, lc1)] = tmp;

                                /* Pair F: negate-swap — within t11 */
                                tmp = t11[elem_off(lr0, lc1)];
                                t11[elem_off(lr0, lc1)] = -t11[elem_off(lr1, lc0)];
                                t11[elem_off(lr1, lc0)] = -tmp;
                            }
                        }
                    }
                }
            }
        }

        /* Case 4: Both control and target >= LOG_TILE_DIM
         *
         * Entire tiles are swapped/phase-swapped. The tile-level structure
         * is identical to CX; only the per-element operations change.
         */
        else {
            const gate_idx_t n_bt = dim_tile >> 2;
            const qubit_t lo = (control < target) ? control - LOG_TILE_DIM : target - LOG_TILE_DIM;
            const qubit_t hi = (control > target) ? control - LOG_TILE_DIM : target - LOG_TILE_DIM;
            const gate_idx_t incr_ctrl = (gate_idx_t)1 << (control - LOG_TILE_DIM);
            const gate_idx_t incr_tgt = (gate_idx_t)1 << (target - LOG_TILE_DIM);

            #pragma omp parallel for schedule(static, 1) if(dim >= OMP_THRESHOLD)
            for (gate_idx_t btr = 0; btr < n_bt; ++btr) {
                const gate_idx_t tr00 = insertBits2_0(btr, lo, hi);
                const gate_idx_t tr01 = tr00 | incr_tgt, tr10 = tr00 | incr_ctrl, tr11 = tr10 | incr_tgt;

                for (gate_idx_t btc = 0; btc <= btr; ++btc) {
                    const gate_idx_t tc00 = insertBits2_0(btc, lo, hi);
                    const gate_idx_t tc01 = tc00 | incr_tgt, tc10 = tc00 | incr_ctrl, tc11 = tc10 | incr_tgt;

                    /* Pair A: always lower-tri; phases (-i, +i) */
                    {
                        cplx_t * restrict tA = data + tile_off(tr10, tc00);
                        cplx_t * restrict tB = data + tile_off(tr11, tc00);
                        for (gate_idx_t i = 0; i < TILE_SIZE; ++i) {
                            cplx_t tmp = tA[i];
                            tA[i] = CMPLX( cimag(tB[i]), -creal(tB[i]));
                            tB[i] = CMPLX(-cimag(tmp),    creal(tmp));
                        }
                    }

                    /* Pair E: always lower-tri; plain swap */
                    {
                        cplx_t * restrict tA = data + tile_off(tr10, tc10);
                        cplx_t * restrict tB = data + tile_off(tr11, tc11);
                        for (gate_idx_t i = 0; i < TILE_SIZE; ++i) {
                            cplx_t tmp = tA[i]; tA[i] = tB[i]; tB[i] = tmp;
                        }
                    }

                    /* Fast path: tr00 > tc11, all tiles lower-tri */
                    if (tr00 > tc11) {
                        /* Pair B: phases (+i, -i) */
                        {
                            cplx_t * restrict tA = data + tile_off(tr11, tc01);
                            cplx_t * restrict tB = data + tile_off(tr10, tc01);
                            for (gate_idx_t i = 0; i < TILE_SIZE; ++i) {
                                cplx_t tmp = tA[i];
                                tA[i] = CMPLX(-cimag(tB[i]),  creal(tB[i]));
                                tB[i] = CMPLX( cimag(tmp),   -creal(tmp));
                            }
                        }
                        /* Pair F: negate-swap */
                        {
                            cplx_t * restrict tA = data + tile_off(tr11, tc10);
                            cplx_t * restrict tB = data + tile_off(tr10, tc11);
                            for (gate_idx_t i = 0; i < TILE_SIZE; ++i) {
                                cplx_t tmp = tA[i]; tA[i] = -tB[i]; tB[i] = -tmp;
                            }
                        }
                        /* Pair D: phases (+i, -i) */
                        {
                            cplx_t * restrict tA = data + tile_off(tr01, tc10);
                            cplx_t * restrict tB = data + tile_off(tr01, tc11);
                            for (gate_idx_t i = 0; i < TILE_SIZE; ++i) {
                                cplx_t tmp = tA[i];
                                tA[i] = CMPLX(-cimag(tB[i]),  creal(tB[i]));
                                tB[i] = CMPLX( cimag(tmp),   -creal(tmp));
                            }
                        }
                        /* Pair C: phases (+i, -i) */
                        {
                            cplx_t * restrict tA = data + tile_off(tr00, tc10);
                            cplx_t * restrict tB = data + tile_off(tr00, tc11);
                            for (gate_idx_t i = 0; i < TILE_SIZE; ++i) {
                                cplx_t tmp = tA[i];
                                tA[i] = CMPLX(-cimag(tB[i]),  creal(tB[i]));
                                tB[i] = CMPLX( cimag(tmp),   -creal(tmp));
                            }
                        }
                    }
                    else {
                        if (tr10 > tc11) {
                            {
                                /* Pair B: phases (+i, -i) — both lower-tri */
                                cplx_t * restrict tA = data + tile_off(tr11, tc01);
                                cplx_t * restrict tB = data + tile_off(tr10, tc01);
                                for (gate_idx_t i = 0; i < TILE_SIZE; ++i) {
                                    cplx_t tmp = tA[i];
                                    tA[i] = CMPLX(-cimag(tB[i]),  creal(tB[i]));
                                    tB[i] = CMPLX( cimag(tmp),   -creal(tmp));
                                }
                            }
                            {
                                /* Pair F: negate-swap — both lower-tri */
                                cplx_t * restrict tA = data + tile_off(tr11, tc10);
                                cplx_t * restrict tB = data + tile_off(tr10, tc11);
                                for (gate_idx_t i = 0; i < TILE_SIZE; ++i) {
                                    cplx_t tmp = tA[i]; tA[i] = -tB[i]; tB[i] = -tmp;
                                }
                            }
                        }
                        else {
                            {
                                /* Pair B: (tr11,tc01) lower; (tr10,tc01) may cross diagonal */
                                cplx_t * restrict tA = data + tile_off(tr11, tc01);
                                if (tr10 > tc01) {
                                    /* Both lower-tri */
                                    cplx_t * restrict tB = data + tile_off(tr10, tc01);
                                    for (gate_idx_t lr = 0; lr < TILE_DIM; ++lr) {
                                        for (gate_idx_t lc = 0; lc < TILE_DIM; ++lc) {
                                            cplx_t tmp = tA[lr * TILE_DIM + lc];
                                            tA[lr * TILE_DIM + lc] = CMPLX(-cimag(tB[lr * TILE_DIM + lc]),
                                                                            creal(tB[lr * TILE_DIM + lc]));
                                            tB[lr * TILE_DIM + lc] = CMPLX( cimag(tmp), -creal(tmp));
                                        }
                                    }
                                }
                                else {
                                    /* tB stores adjoint of (tr10,tc01) at (tc01,tr10) */
                                    cplx_t * restrict tB = data + tile_off(tc01, tr10);
                                    for (gate_idx_t lr = 0; lr < TILE_DIM; ++lr) {
                                        for (gate_idx_t lc = 0; lc < TILE_DIM; ++lc) {
                                            /* tA[lr,lc] is lower-tri = rho[r11,c01]
                                             * tB[lc,lr] is stored conj(rho[r10,c01])
                                             * Pair B: new rho[r11,c01] = +i*conj(tB[lc,lr])
                                             *   = CMPLX(cimag(tB[lc,lr]), creal(tB[lc,lr]))
                                             * new stored at tB[lc,lr] = conj(-i*tA[lr,lc])
                                             *   = CMPLX(cimag(tA[lr,lc]), creal(tA[lr,lc]))
                                             */
                                            cplx_t tmp = tA[lr * TILE_DIM + lc];
                                            cplx_t s = tB[lc * TILE_DIM + lr];
                                            tA[lr * TILE_DIM + lc] = CMPLX(cimag(s), creal(s));
                                            tB[lc * TILE_DIM + lr] = CMPLX(cimag(tmp), creal(tmp));
                                        }
                                    }
                                }
                            }

                            if (btr != btc) {
                                /* Pair F: negate-swap; (tr11,tc10) lower, (tr10,tc11) upper */
                                cplx_t * restrict tA = data + tile_off(tr11, tc10);
                                cplx_t * restrict tB = data + tile_off(tc11, tr10);
                                for (gate_idx_t lr = 0; lr < TILE_DIM; ++lr) {
                                    for (gate_idx_t lc = 0; lc < TILE_DIM; ++lc) {
                                        /* tA[lr,lc] = rho[r11,c10] (lower-tri)
                                         * tB[lc,lr] = conj(rho[r10,c11]) (upper-tri stored)
                                         * new rho[r11,c10] = -rho[r10,c11] = -conj(tB[lc,lr])
                                         * new stored tB[lc,lr] = conj(-rho[r11,c10]) = -conj(tA[lr,lc])
                                         */
                                        cplx_t tmp = tA[lr * TILE_DIM + lc];
                                        tA[lr * TILE_DIM + lc] = -conj(tB[lc * TILE_DIM + lr]);
                                        tB[lc * TILE_DIM + lr] = -conj(tmp);
                                    }
                                }
                            }
                            else {
                                /* Pair F: negate-swap; diagonal case (btr == btc) */
                                cplx_t * restrict tA = data + tile_off(tr11, tc10);
                                cplx_t * restrict tB = data + tile_off(tc11, tr10);
                                for (gate_idx_t lr = 0; lr < TILE_DIM; ++lr) {
                                    for (gate_idx_t lc = 0; lc <= lr; ++lc) {
                                        cplx_t tmp = tA[lr * TILE_DIM + lc];
                                        tA[lr * TILE_DIM + lc] = -conj(tB[lc * TILE_DIM + lr]);
                                        tB[lc * TILE_DIM + lr] = -conj(tmp);
                                    }
                                }
                            }
                        }

                        if (btr != btc) {
                            if (tr00 > tc10) {
                                /* Pair D: phases (+i, -i) — both lower-tri */
                                {
                                    cplx_t * restrict tA = data + tile_off(tr01, tc10);
                                    cplx_t * restrict tB = data + tile_off(tr01, tc11);
                                    for (gate_idx_t i = 0; i < TILE_SIZE; ++i) {
                                        cplx_t tmp = tA[i];
                                        tA[i] = CMPLX(-cimag(tB[i]),  creal(tB[i]));
                                        tB[i] = CMPLX( cimag(tmp),   -creal(tmp));
                                    }
                                }
                                /* Pair C: (tr00,tc10) lower, (tr00,tc11) upper (via tc11,tr00) */
                                {
                                    cplx_t * restrict tA = data + tile_off(tr00, tc10);
                                    cplx_t * restrict tB = data + tile_off(tc11, tr00);
                                    for (gate_idx_t lr = 0; lr < TILE_DIM; ++lr) {
                                        for (gate_idx_t lc = 0; lc < TILE_DIM; ++lc) {
                                            /* tA[lr,lc] = rho[r00,c10] (lower-tri)
                                             * tB[lc,lr] = conj(rho[r00,c11]) (upper-tri stored)
                                             * Pair C: new rho[r00,c10] = +i*rho[r00,c11] = +i*conj(tB[lc,lr])
                                             *   = CMPLX(cimag(tB[lc,lr]), creal(tB[lc,lr]))
                                             * new stored tB[lc,lr] = conj(-i*rho[r00,c10]) = conj(-i*tA[lr,lc])
                                             *   = CMPLX(cimag(tA[lr,lc]), creal(tA[lr,lc]))
                                             */
                                            cplx_t tmp = tA[lr * TILE_DIM + lc];
                                            cplx_t s = tB[lc * TILE_DIM + lr];
                                            tA[lr * TILE_DIM + lc] = CMPLX(cimag(s), creal(s));
                                            tB[lc * TILE_DIM + lr] = CMPLX(cimag(tmp), creal(tmp));
                                        }
                                    }
                                }
                            }
                            else {
                                /* Pair D: (tr01,tc11) upper at (tc11,tr01);
                                 *         (tr01,tc10) may be lower or upper */
                                {
                                    cplx_t * restrict tA = data + tile_off(tc11, tr01);
                                    if (tr01 > tc10) {
                                        cplx_t * restrict tB = data + tile_off(tr01, tc10);
                                        for (gate_idx_t lr = 0; lr < TILE_DIM; ++lr) {
                                            for (gate_idx_t lc = 0; lc < TILE_DIM; ++lc) {
                                                /* tA[lc,lr] = conj(rho[r01,c11]) (upper stored)
                                                 * tB[lr,lc] = rho[r01,c10] (lower-tri)
                                                 * Pair D: new rho[r01,c10] = +i*rho[r01,c11] = +i*conj(tA[lc,lr])
                                                 *   = CMPLX(cimag(tA[lc,lr]), creal(tA[lc,lr]))
                                                 * new stored tA[lc,lr] = conj(-i*rho[r01,c10]) = conj(-i*tB[lr,lc])
                                                 *   = CMPLX(cimag(tB[lr,lc]), creal(tB[lr,lc]))
                                                 */
                                                cplx_t tmp = tB[lr * TILE_DIM + lc];
                                                cplx_t s = tA[lc * TILE_DIM + lr];
                                                tB[lr * TILE_DIM + lc] = CMPLX(cimag(s), creal(s));
                                                tA[lc * TILE_DIM + lr] = CMPLX(cimag(tmp), creal(tmp));
                                            }
                                        }
                                    }
                                    else {
                                        /* Both upper-tri: (tr01,tc10) stored at (tc10,tr01),
                                         *                 (tr01,tc11) stored at (tc11,tr01) = tA */
                                        cplx_t * restrict tB = data + tile_off(tc10, tr01);
                                        /* Both stored values are conj of actual rho.
                                         * tA[i] = conj(rho[r01,c11]), tB[i] = conj(rho[r01,c10])
                                         * stored at same (lr,lc) positions (not transposed between each other).
                                         * Pair D: new rho[r01,c10] = +i*rho[r01,c11]
                                         *   new stored tB[i] = conj(+i*rho[r01,c11]) = -i*conj(rho[r01,c11]) = -i*tA[i]
                                         * Pair D: new rho[r01,c11] = -i*rho[r01,c10]
                                         *   new stored tA[i] = conj(-i*rho[r01,c10]) = +i*conj(rho[r01,c10]) = +i*tB[i]
                                         */
                                        for (gate_idx_t i = 0; i < TILE_SIZE; ++i) {
                                            cplx_t tmp = tA[i];
                                            cplx_t s = tB[i];
                                            tA[i] = CMPLX(-cimag(s), creal(s));
                                            tB[i] = CMPLX( cimag(tmp), -creal(tmp));
                                        }
                                    }
                                }
                                /* Pair C: (tr00,tc10) upper at (tc10,tr00);
                                 *         (tr00,tc11) upper at (tc11,tr00) */
                                {
                                    cplx_t * restrict tA = data + tile_off(tc10, tr00);
                                    cplx_t * restrict tB = data + tile_off(tc11, tr00);
                                    /* Both upper-tri stored at same-position tile indices.
                                     * tA[i] = conj(rho[r00,c10]) stored at (tc10,tr00)
                                     * tB[i] = conj(rho[r00,c11]) stored at (tc11,tr00)
                                     * Pair C: new rho[r00,c10] = +i*rho[r00,c11]
                                     *   new stored tA[i] = conj(+i)*tB[i] = -i*tB[i]
                                     * new rho[r00,c11] = -i*rho[r00,c10]
                                     *   new stored tB[i] = conj(-i)*tA[i] = +i*tA[i]
                                     */
                                    for (gate_idx_t i = 0; i < TILE_SIZE; ++i) {
                                        cplx_t tmp = tA[i];
                                        cplx_t s = tB[i];
                                        tA[i] = CMPLX( cimag(s), -creal(s));
                                        tB[i] = CMPLX(-cimag(tmp), creal(tmp));
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
