#include "gate_tiled.h"

void two_from_mat(state_t *state, const qubit_t q1, const qubit_t q2, const cplx_t *mat) {
    const gate_idx_t dim = (gate_idx_t)1 << state->qubits;
    const gate_idx_t dim_tile = (dim + TILE_DIM - 1) >> LOG_TILE_DIM;
    qubit_t lo = (q1 < q2) ? q1 : q2;
    qubit_t hi = (q1 > q2) ? q1 : q2;
    cplx_t * restrict data = state->data;

    if (hi < LOG_TILE_DIM) {
        const gate_idx_t n_tiles = dim_tile * (dim_tile + 1) / 2;
        const gate_idx_t nbl = (dim < TILE_DIM) ? dim >> 2 : TILE_DIM >> 2;
        const gate_idx_t incr_lo = (gate_idx_t)1 << lo;
        const gate_idx_t incr_hi = (gate_idx_t)1 << hi;

        #pragma omp parallel for schedule(static) if(dim >= OMP_THRESHOLD)
        for (gate_idx_t i = 0; i < n_tiles; ++i) {
            cplx_t * restrict tile = data + i * TILE_SIZE;

            for (gate_idx_t blr = 0; blr < nbl; ++blr) {
                const gate_idx_t lr00 = insertBits2_0(blr, lo, hi);
                const gate_idx_t lr01 = lr00 | incr_lo;
                const gate_idx_t lr10 = lr00 | incr_hi;
                const gate_idx_t lr11 = lr10 | incr_lo;

                for (gate_idx_t blc = 0; blc < nbl; ++blc) {
                    const gate_idx_t lc00 = insertBits2_0(blc, lo, hi);
                    const gate_idx_t lc01 = lc00 | incr_lo;
                    const gate_idx_t lc10 = lc00 | incr_hi;
                    const gate_idx_t lc11 = lc10 | incr_lo;
                    cplx_t old[16], new[16] = {0};

                    old[0] = tile[elem_off(lr00, lc00)];
                    old[1] = tile[elem_off(lr01, lc00)];
                    old[2] = tile[elem_off(lr10, lc00)];
                    old[3] = tile[elem_off(lr11, lc00)];
                    old[4] = tile[elem_off(lr00, lc01)];
                    old[5] = tile[elem_off(lr01, lc01)];
                    old[6] = tile[elem_off(lr10, lc01)];
                    old[7] = tile[elem_off(lr11, lc01)];
                    old[8] = tile[elem_off(lr00, lc10)];
                    old[9] = tile[elem_off(lr01, lc10)];
                    old[10] = tile[elem_off(lr10, lc10)];
                    old[11] = tile[elem_off(lr11, lc10)];
                    old[12] = tile[elem_off(lr00, lc11)];
                    old[13] = tile[elem_off(lr01, lc11)];
                    old[14] = tile[elem_off(lr10, lc11)];
                    old[15] = tile[elem_off(lr11, lc11)];

                    for (unsigned char col = 0; col < 4; ++col) {
                        for (unsigned char row = 0; row < 4; ++row) {
                            for (unsigned char k = 0; k < 4; ++k) {
                                for (unsigned char l = 0; l < 4; ++l) {
                                    new[row + col * 4] += mat[row + k * 4] * old[k + l * 4] * conj(mat[col + l * 4]);
                                }
                            }
                        }
                    }

                    tile[elem_off(lr00, lc00)] = new[0];
                    tile[elem_off(lr01, lc00)] = new[1];
                    tile[elem_off(lr10, lc00)] = new[2];
                    tile[elem_off(lr11, lc00)] = new[3];
                    tile[elem_off(lr00, lc01)] = new[4];
                    tile[elem_off(lr01, lc01)] = new[5];
                    tile[elem_off(lr10, lc01)] = new[6];
                    tile[elem_off(lr11, lc01)] = new[7];
                    tile[elem_off(lr00, lc10)] = new[8];
                    tile[elem_off(lr01, lc10)] = new[9];
                    tile[elem_off(lr10, lc10)] = new[10];
                    tile[elem_off(lr11, lc10)] = new[11];
                    tile[elem_off(lr00, lc11)] = new[12];
                    tile[elem_off(lr01, lc11)] = new[13];
                    tile[elem_off(lr10, lc11)] = new[14];
                    tile[elem_off(lr11, lc11)] = new[15];
                }
            }
        }
    }

    else if (lo < LOG_TILE_DIM) {
        const qubit_t idx_tile = hi - LOG_TILE_DIM;
        const gate_idx_t incr_tile = (gate_idx_t)1 << idx_tile;
        const gate_idx_t nbt = dim_tile >> 1;
        const gate_idx_t incr_lo = (gate_idx_t)1 << lo;
        const gate_idx_t nbl = TILE_DIM >> 1;

        #pragma omp parallel for schedule(static, 1) if(dim >= OMP_THRESHOLD)
        for (gate_idx_t btr = 0; btr < nbt; ++btr) {
            const gate_idx_t tr0 = insertBit0(btr, idx_tile);
            const gate_idx_t tr1 = tr0 | incr_tile;

            for (gate_idx_t btc = 0; btc <= btr; ++btc) {
                const gate_idx_t tc0 = insertBit0(btc, idx_tile);
                const gate_idx_t tc1 = tc0 | incr_tile;

                cplx_t * restrict t00 = data + tile_off(tr0, tc0);
                cplx_t * restrict t10 = data + tile_off(tr1, tc0);
                cplx_t * restrict t11 = data + tile_off(tr1, tc1);

                if (tr0 > tc1) {
                    cplx_t * restrict t01 = data + tile_off(tr0, tc1);

                    for (gate_idx_t blr = 0; blr < nbl; ++blr) {
                        const gate_idx_t lr0 = insertBit0(blr, lo);
                        const gate_idx_t lr1 = lr0 | incr_lo;

                        for (gate_idx_t blc = 0; blc < nbl; ++blc) {
                            const gate_idx_t lc0 = insertBit0(blc, lo);
                            const gate_idx_t lc1 = lc0 | incr_lo;
                            cplx_t old[16], new[16] = {0};

                            old[0] = t00[elem_off(lr0, lc0)];
                            old[1] = t00[elem_off(lr1, lc0)];
                            old[2] = t10[elem_off(lr0, lc0)];
                            old[3] = t10[elem_off(lr1, lc0)];
                            old[4] = t00[elem_off(lr0, lc1)];
                            old[5] = t00[elem_off(lr1, lc1)];
                            old[6] = t10[elem_off(lr0, lc1)];
                            old[7] = t10[elem_off(lr1, lc1)];
                            old[8] = t01[elem_off(lr0, lc0)];
                            old[9] = t01[elem_off(lr1, lc0)];
                            old[10] = t11[elem_off(lr0, lc0)];
                            old[11] = t11[elem_off(lr1, lc0)];
                            old[12] = t01[elem_off(lr0, lc1)];
                            old[13] = t01[elem_off(lr1, lc1)];
                            old[14] = t11[elem_off(lr0, lc1)];
                            old[15] = t11[elem_off(lr1, lc1)];

                            for (unsigned char col = 0; col < 4; ++col) {
                                for (unsigned char row = 0; row < 4; ++row) {
                                    for (unsigned char k = 0; k < 4; ++k) {
                                        for (unsigned char l = 0; l < 4; ++l) {
                                            new[row + col * 4] += mat[row + k * 4] * old[k + l * 4] * conj(mat[col + l * 4]);
                                        }
                                    }
                                }
                            }

                            t00[elem_off(lr0, lc0)] = new[0];
                            t00[elem_off(lr1, lc0)] = new[1];
                            t10[elem_off(lr0, lc0)] = new[2];
                            t10[elem_off(lr1, lc0)] = new[3];
                            t00[elem_off(lr0, lc1)] = new[4];
                            t00[elem_off(lr1, lc1)] = new[5];
                            t10[elem_off(lr0, lc1)] = new[6];
                            t10[elem_off(lr1, lc1)] = new[7];
                            t01[elem_off(lr0, lc0)] = new[8];
                            t01[elem_off(lr1, lc0)] = new[9];
                            t11[elem_off(lr0, lc0)] = new[10];
                            t11[elem_off(lr1, lc0)] = new[11];
                            t01[elem_off(lr0, lc1)] = new[12];
                            t01[elem_off(lr1, lc1)] = new[13];
                            t11[elem_off(lr0, lc1)] = new[14];
                            t11[elem_off(lr1, lc1)] = new[15];
                        }
                    }
                }

                else if (btc != btr) {
                    cplx_t * restrict t01 = data + tile_off(tc1, tr0);

                    for (gate_idx_t blr = 0; blr < nbl; ++blr) {
                        const gate_idx_t lr0 = insertBit0(blr, lo);
                        const gate_idx_t lr1 = lr0 | incr_lo;

                        for (gate_idx_t blc = 0; blc < nbl; ++blc) {
                            const gate_idx_t lc0 = insertBit0(blc, lo);
                            const gate_idx_t lc1 = lc0 | incr_lo;
                            cplx_t old[16], new[16] = {0};

                            old[0] = t00[elem_off(lr0, lc0)];
                            old[1] = t00[elem_off(lr1, lc0)];
                            old[2] = t10[elem_off(lr0, lc0)];
                            old[3] = t10[elem_off(lr1, lc0)];
                            old[4] = t00[elem_off(lr0, lc1)];
                            old[5] = t00[elem_off(lr1, lc1)];
                            old[6] = t10[elem_off(lr0, lc1)];
                            old[7] = t10[elem_off(lr1, lc1)];
                            old[8] = conj(t01[elem_off(lc0, lr0)]);
                            old[9] = conj(t01[elem_off(lc0, lr1)]);
                            old[10] = t11[elem_off(lr0, lc0)];
                            old[11] = t11[elem_off(lr1, lc0)];
                            old[12] = conj(t01[elem_off(lc1, lr0)]);
                            old[13] = conj(t01[elem_off(lc1, lr1)]);
                            old[14] = t11[elem_off(lr0, lc1)];
                            old[15] = t11[elem_off(lr1, lc1)];

                            for (unsigned char col = 0; col < 4; ++col) {
                                for (unsigned char row = 0; row < 4; ++row) {
                                    for (unsigned char k = 0; k < 4; ++k) {
                                        for (unsigned char l = 0; l < 4; ++l) {
                                            new[row + col * 4] += mat[row + k * 4] * old[k + l * 4] * conj(mat[col + l * 4]);
                                        }
                                    }
                                }
                            }

                            t00[elem_off(lr0, lc0)] = new[0];
                            t00[elem_off(lr1, lc0)] = new[1];
                            t10[elem_off(lr0, lc0)] = new[2];
                            t10[elem_off(lr1, lc0)] = new[3];
                            t00[elem_off(lr0, lc1)] = new[4];
                            t00[elem_off(lr1, lc1)] = new[5];
                            t10[elem_off(lr0, lc1)] = new[6];
                            t10[elem_off(lr1, lc1)] = new[7];
                            t01[elem_off(lc0, lr0)] = conj(new[8]);
                            t01[elem_off(lc0, lr1)] = conj(new[9]);
                            t11[elem_off(lr0, lc0)] = new[10];
                            t11[elem_off(lr1, lc0)] = new[11];
                            t01[elem_off(lc1, lr0)] = conj(new[12]);
                            t01[elem_off(lc1, lr1)] = conj(new[13]);
                            t11[elem_off(lr0, lc1)] = new[14];
                            t11[elem_off(lr1, lc1)] = new[15];
                        }
                    }
                }

                else {
                    for (gate_idx_t blr = 0; blr < nbl; ++blr) {
                        const gate_idx_t lr0 = insertBit0(blr, lo);
                        const gate_idx_t lr1 = lr0 | incr_lo;

                        for (gate_idx_t blc = 0; blc <= blr; ++blc) {
                            const gate_idx_t lc0 = insertBit0(blc, lo);
                            const gate_idx_t lc1 = lc0 | incr_lo;
                            cplx_t old[16], new[16] = {0};

                            old[0] = t00[elem_off(lr0, lc0)];
                            old[1] = t00[elem_off(lr1, lc0)];
                            old[2] = t10[elem_off(lr0, lc0)];
                            old[3] = t10[elem_off(lr1, lc0)];
                            old[4] = t00[elem_off(lr0, lc1)];
                            old[5] = t00[elem_off(lr1, lc1)];
                            old[6] = t10[elem_off(lr0, lc1)];
                            old[7] = t10[elem_off(lr1, lc1)];
                            old[8] = conj(t10[elem_off(lc0, lr0)]);
                            old[9] = conj(t10[elem_off(lc0, lr1)]);
                            old[10] = t11[elem_off(lr0, lc0)];
                            old[11] = t11[elem_off(lr1, lc0)];
                            old[12] = conj(t10[elem_off(lc1, lr0)]);
                            old[13] = conj(t10[elem_off(lc1, lr1)]);
                            old[14] = t11[elem_off(lr0, lc1)];
                            old[15] = t11[elem_off(lr1, lc1)];

                            for (unsigned char col = 0; col < 4; ++col) {
                                for (unsigned char row = 0; row < 4; ++row) {
                                    for (unsigned char k = 0; k < 4; ++k) {
                                        for (unsigned char l = 0; l < 4; ++l) {
                                            new[row + col * 4] += mat[row + k * 4] * old[k + l * 4] * conj(mat[col + l * 4]);
                                        }
                                    }
                                }
                            }

                            t00[elem_off(lr0, lc0)] = new[0];
                            t00[elem_off(lc0, lr0)] = conj(new[0]);
                            t00[elem_off(lr1, lc0)] = new[1];
                            t00[elem_off(lc0, lr1)] = conj(new[1]);
                            t10[elem_off(lr0, lc0)] = new[2];
                            t10[elem_off(lr1, lc0)] = new[3];
                            t00[elem_off(lr0, lc1)] = new[4];
                            t00[elem_off(lc1, lr0)] = conj(new[4]);
                            t00[elem_off(lr1, lc1)] = new[5];
                            t00[elem_off(lc1, lr1)] = conj(new[5]);
                            t10[elem_off(lr0, lc1)] = new[6];
                            t10[elem_off(lr1, lc1)] = new[7];
                            t10[elem_off(lc0, lr0)] = conj(new[8]);
                            t10[elem_off(lc0, lr1)] = conj(new[9]);
                            t11[elem_off(lr0, lc0)] = new[10];
                            t11[elem_off(lc0, lr0)] = conj(new[10]);
                            t11[elem_off(lr1, lc0)] = new[11];
                            t11[elem_off(lc0, lr1)] = conj(new[11]);
                            t10[elem_off(lc1, lr0)] = conj(new[12]);
                            t10[elem_off(lc1, lr1)] = conj(new[13]);
                            t11[elem_off(lr0, lc1)] = new[14];
                            t11[elem_off(lc1, lr0)] = conj(new[14]);
                            t11[elem_off(lr1, lc1)] = new[15];
                            t11[elem_off(lc1, lr1)] = conj(new[15]);
                        }
                    }
                }
            }
        }
    }

    else {
        const gate_idx_t nbt = dim_tile >> 2;
        lo -= LOG_TILE_DIM;
        hi -= LOG_TILE_DIM;
        const gate_idx_t incr_lo = (gate_idx_t)1 << lo;
        const gate_idx_t incr_hi = (gate_idx_t)1 << hi;

        #pragma omp parallel for schedule(static, 1) if(dim >= OMP_THRESHOLD)
        for (gate_idx_t btr = 0; btr < nbt; ++btr) {
            const gate_idx_t tr00 = insertBits2_0(btr, lo, hi);
            const gate_idx_t tr01 = tr00 | incr_lo, tr10 = tr00 | incr_hi, tr11 = tr10 | incr_lo;

            for (gate_idx_t btc = 0; btc <= btr; ++btc) {
                const gate_idx_t tc00 = insertBits2_0(btc, lo, hi);
                const gate_idx_t tc01 = tc00 | incr_lo, tc10 = tc00 | incr_hi, tc11 = tc10 | incr_lo;

                if (tr00 > tc11) {

                    for (gate_idx_t i = 0; i < TILE_SIZE; ++i) {
                        cplx_t old[16], new[16] = {0};

                        old[0] = data[tile_off(tr00, tc00) + i];
                        old[1] = data[tile_off(tr01, tc00) + i];
                        old[2] = data[tile_off(tr10, tc00) + i];
                        old[3] = data[tile_off(tr11, tc00) + i];
                        old[4] = data[tile_off(tr00, tc01) + i];
                        old[5] = data[tile_off(tr01, tc01) + i];
                        old[6] = data[tile_off(tr10, tc01) + i];
                        old[7] = data[tile_off(tr11, tc01) + i];
                        old[8] = data[tile_off(tr00, tc10) + i];
                        old[9] = data[tile_off(tr01, tc10) + i];
                        old[10] = data[tile_off(tr10, tc10) + i];
                        old[11] = data[tile_off(tr11, tc10) + i];
                        old[12] = data[tile_off(tr00, tc11) + i];
                        old[13] = data[tile_off(tr01, tc11) + i];
                        old[14] = data[tile_off(tr10, tc11) + i];
                        old[15] = data[tile_off(tr11, tc11) + i];

                        for (unsigned char col = 0; col < 4; ++col) {
                            for (unsigned char row = 0; row < 4; ++row) {
                                for (unsigned char k = 0; k < 4; ++k) {
                                    for (unsigned char l = 0; l < 4; ++l) {
                                        new[row + col * 4] += mat[row + k * 4] * old[k + l * 4] \
                                                              * conj(mat[col + l * 4]);
                                    }
                                }
                            }
                        }

                        data[tile_off(tr00, tc00) + i] = new[0];
                        data[tile_off(tr01, tc00) + i] = new[1];
                        data[tile_off(tr10, tc00) + i] = new[2];
                        data[tile_off(tr11, tc00) + i] = new[3];
                        data[tile_off(tr00, tc01) + i] = new[4];
                        data[tile_off(tr01, tc01) + i] = new[5];
                        data[tile_off(tr10, tc01) + i] = new[6];
                        data[tile_off(tr11, tc01) + i] = new[7];
                        data[tile_off(tr00, tc10) + i] = new[8];
                        data[tile_off(tr01, tc10) + i] = new[9];
                        data[tile_off(tr10, tc10) + i] = new[10];
                        data[tile_off(tr11, tc10) + i] = new[11];
                        data[tile_off(tr00, tc11) + i] = new[12];
                        data[tile_off(tr01, tc11) + i] = new[13];
                        data[tile_off(tr10, tc11) + i] = new[14];
                        data[tile_off(tr11, tc11) + i] = new[15];
                    }
                }

                else if (tr00 > tc10) {

                    for (gate_idx_t i = 0; i < TILE_DIM; ++i) {
                        for (gate_idx_t j = 0; j < TILE_DIM; ++j) {
                            cplx_t old[16], new[16] = {0};

                            old[0] = data[tile_off(tr00, tc00) + elem_off(i, j)];
                            old[1] = data[tile_off(tr01, tc00) + elem_off(i, j)];
                            old[2] = data[tile_off(tr10, tc00) + elem_off(i, j)];
                            old[3] = data[tile_off(tr11, tc00) + elem_off(i, j)];
                            old[4] = data[tile_off(tr00, tc01) + elem_off(i, j)];
                            old[5] = data[tile_off(tr01, tc01) + elem_off(i, j)];
                            old[6] = data[tile_off(tr10, tc01) + elem_off(i, j)];
                            old[7] = data[tile_off(tr11, tc01) + elem_off(i, j)];
                            old[8] = data[tile_off(tr00, tc10) + elem_off(i, j)];
                            old[9] = data[tile_off(tr01, tc10) + elem_off(i, j)];
                            old[10] = data[tile_off(tr10, tc10) + elem_off(i, j)];
                            old[11] = data[tile_off(tr11, tc10) + elem_off(i, j)];
                            old[12] = conj(data[tile_off(tc11, tr00) + elem_off(j, i)]);
                            old[13] = data[tile_off(tr01, tc11) + elem_off(i, j)];
                            old[14] = data[tile_off(tr10, tc11) + elem_off(i, j)];
                            old[15] = data[tile_off(tr11, tc11) + elem_off(i, j)];

                            for (unsigned char col = 0; col < 4; ++col) {
                                for (unsigned char row = 0; row < 4; ++row) {
                                    for (unsigned char k = 0; k < 4; ++k) {
                                        for (unsigned char l = 0; l < 4; ++l) {
                                            new[row + col * 4] += mat[row + k * 4] * old[k + l * 4] \
                                                                  * conj(mat[col + l * 4]);
                                        }
                                    }
                                }
                            }

                            data[tile_off(tr00, tc00) + elem_off(i, j)] = new[0];
                            data[tile_off(tr01, tc00) + elem_off(i, j)] = new[1];
                            data[tile_off(tr10, tc00) + elem_off(i, j)] = new[2];
                            data[tile_off(tr11, tc00) + elem_off(i, j)] = new[3];
                            data[tile_off(tr00, tc01) + elem_off(i, j)] = new[4];
                            data[tile_off(tr01, tc01) + elem_off(i, j)] = new[5];
                            data[tile_off(tr10, tc01) + elem_off(i, j)] = new[6];
                            data[tile_off(tr11, tc01) + elem_off(i, j)] = new[7];
                            data[tile_off(tr00, tc10) + elem_off(i, j)] = new[8];
                            data[tile_off(tr01, tc10) + elem_off(i, j)] = new[9];
                            data[tile_off(tr10, tc10) + elem_off(i, j)] = new[10];
                            data[tile_off(tr11, tc10) + elem_off(i, j)] = new[11];
                            data[tile_off(tc11, tr00) + elem_off(j, i)] = conj(new[12]);
                            data[tile_off(tr01, tc11) + elem_off(i, j)] = new[13];
                            data[tile_off(tr10, tc11) + elem_off(i, j)] = new[14];
                            data[tile_off(tr11, tc11) + elem_off(i, j)] = new[15];
                        }
                    }
                }

                else if (tr00 > tc01) {

                    if(tr01 > tc10) {

                        for (gate_idx_t i = 0; i < TILE_DIM; ++i) {
                            for (gate_idx_t j = 0; j < TILE_DIM; ++j) {
                                cplx_t old[16], new[16] = {0};

                                old[0] = data[tile_off(tr00, tc00) + elem_off(i, j)];
                                old[1] = data[tile_off(tr01, tc00) + elem_off(i, j)];
                                old[2] = data[tile_off(tr10, tc00) + elem_off(i, j)];
                                old[3] = data[tile_off(tr11, tc00) + elem_off(i, j)];
                                old[4] = data[tile_off(tr00, tc01) + elem_off(i, j)];
                                old[5] = data[tile_off(tr01, tc01) + elem_off(i, j)];
                                old[6] = data[tile_off(tr10, tc01) + elem_off(i, j)];
                                old[7] = data[tile_off(tr11, tc01) + elem_off(i, j)];
                                old[8] = conj(data[tile_off(tc10, tr00) + elem_off(j, i)]);
                                old[9] = data[tile_off(tr01, tc10) + elem_off(i, j)];
                                old[10] = data[tile_off(tr10, tc10) + elem_off(i, j)];
                                old[11] = data[tile_off(tr11, tc10) + elem_off(i, j)];
                                old[12] = conj(data[tile_off(tc11, tr00) + elem_off(j, i)]);
                                old[13] = conj(data[tile_off(tc11, tr01) + elem_off(j, i)]);
                                old[14] = data[tile_off(tr10, tc11) + elem_off(i, j)];
                                old[15] = data[tile_off(tr11, tc11) + elem_off(i, j)];

                                for (unsigned char col = 0; col < 4; ++col) {
                                    for (unsigned char row = 0; row < 4; ++row) {
                                        for (unsigned char k = 0; k < 4; ++k) {
                                            for (unsigned char l = 0; l < 4; ++l) {
                                                new[row + col * 4] += mat[row + k * 4] * old[k + l * 4] \
                                                                      * conj(mat[col + l * 4]);
                                            }
                                        }
                                    }
                                }

                                data[tile_off(tr00, tc00) + elem_off(i, j)] = new[0];
                                data[tile_off(tr01, tc00) + elem_off(i, j)] = new[1];
                                data[tile_off(tr10, tc00) + elem_off(i, j)] = new[2];
                                data[tile_off(tr11, tc00) + elem_off(i, j)] = new[3];
                                data[tile_off(tr00, tc01) + elem_off(i, j)] = new[4];
                                data[tile_off(tr01, tc01) + elem_off(i, j)] = new[5];
                                data[tile_off(tr10, tc01) + elem_off(i, j)] = new[6];
                                data[tile_off(tr11, tc01) + elem_off(i, j)] = new[7];
                                data[tile_off(tc10, tr00) + elem_off(j, i)] = conj(new[8]);
                                data[tile_off(tr01, tc10) + elem_off(i, j)] = new[9];
                                data[tile_off(tr10, tc10) + elem_off(i, j)] = new[10];
                                data[tile_off(tr11, tc10) + elem_off(i, j)] = new[11];
                                data[tile_off(tc11, tr00) + elem_off(j, i)] = conj(new[12]);
                                data[tile_off(tc11, tr01) + elem_off(j, i)] = conj(new[13]);
                                data[tile_off(tr10, tc11) + elem_off(i, j)] = new[14];
                                data[tile_off(tr11, tc11) + elem_off(i, j)] = new[15];
                            }
                        }
                    }

                    else {

                        for (gate_idx_t i = 0; i < TILE_DIM; ++i) {
                            for (gate_idx_t j = 0; j < TILE_DIM; ++j) {
                                cplx_t old[16], new[16] = {0};

                                old[0] = data[tile_off(tr00, tc00) + elem_off(i, j)];
                                old[1] = data[tile_off(tr01, tc00) + elem_off(i, j)];
                                old[2] = data[tile_off(tr10, tc00) + elem_off(i, j)];
                                old[3] = data[tile_off(tr11, tc00) + elem_off(i, j)];
                                old[4] = data[tile_off(tr00, tc01) + elem_off(i, j)];
                                old[5] = data[tile_off(tr01, tc01) + elem_off(i, j)];
                                old[6] = data[tile_off(tr10, tc01) + elem_off(i, j)];
                                old[7] = data[tile_off(tr11, tc01) + elem_off(i, j)];
                                old[8] = conj(data[tile_off(tc10, tr00) + elem_off(j, i)]);
                                old[9] = conj(data[tile_off(tc10, tr01) + elem_off(j, i)]);
                                old[10] = data[tile_off(tr10, tc10) + elem_off(i, j)];
                                old[11] = data[tile_off(tr11, tc10) + elem_off(i, j)];
                                old[12] = conj(data[tile_off(tc11, tr00) + elem_off(j, i)]);
                                old[13] = conj(data[tile_off(tc11, tr01) + elem_off(j, i)]);
                                old[14] = data[tile_off(tr10, tc11) + elem_off(i, j)];
                                old[15] = data[tile_off(tr11, tc11) + elem_off(i, j)];

                                for (unsigned char col = 0; col < 4; ++col) {
                                    for (unsigned char row = 0; row < 4; ++row) {
                                        for (unsigned char k = 0; k < 4; ++k) {
                                            for (unsigned char l = 0; l < 4; ++l) {
                                                new[row + col * 4] += mat[row + k * 4] * old[k + l * 4] \
                                                                      * conj(mat[col + l * 4]);
                                            }
                                        }
                                    }
                                }

                                data[tile_off(tr00, tc00) + elem_off(i, j)] = new[0];
                                data[tile_off(tr01, tc00) + elem_off(i, j)] = new[1];
                                data[tile_off(tr10, tc00) + elem_off(i, j)] = new[2];
                                data[tile_off(tr11, tc00) + elem_off(i, j)] = new[3];
                                data[tile_off(tr00, tc01) + elem_off(i, j)] = new[4];
                                data[tile_off(tr01, tc01) + elem_off(i, j)] = new[5];
                                data[tile_off(tr10, tc01) + elem_off(i, j)] = new[6];
                                data[tile_off(tr11, tc01) + elem_off(i, j)] = new[7];
                                data[tile_off(tc10, tr00) + elem_off(j, i)] = conj(new[8]);
                                data[tile_off(tc10, tr01) + elem_off(j, i)] = conj(new[9]);
                                data[tile_off(tr10, tc10) + elem_off(i, j)] = new[10];
                                data[tile_off(tr11, tc10) + elem_off(i, j)] = new[11];
                                data[tile_off(tc11, tr00) + elem_off(j, i)] = conj(new[12]);
                                data[tile_off(tc11, tr01) + elem_off(j, i)] = conj(new[13]);
                                data[tile_off(tr10, tc11) + elem_off(i, j)] = new[14];
                                data[tile_off(tr11, tc11) + elem_off(i, j)] = new[15];
                            }
                        }
                    }
                }

                else if (btr != btc) {

                    if(tr01 > tc10) {

                        for (gate_idx_t i = 0; i < TILE_DIM; ++i) {
                            for (gate_idx_t j = 0; j < TILE_DIM; ++j) {
                                cplx_t old[16], new[16] = {0};

                                old[0] = data[tile_off(tr00, tc00) + elem_off(i, j)];
                                old[1] = data[tile_off(tr01, tc00) + elem_off(i, j)];
                                old[2] = data[tile_off(tr10, tc00) + elem_off(i, j)];
                                old[3] = data[tile_off(tr11, tc00) + elem_off(i, j)];
                                old[4] = conj(data[tile_off(tc01, tr00) + elem_off(j, i)]);
                                old[5] = data[tile_off(tr01, tc01) + elem_off(i, j)];
                                old[6] = data[tile_off(tr10, tc01) + elem_off(i, j)];
                                old[7] = data[tile_off(tr11, tc01) + elem_off(i, j)];
                                old[8] = conj(data[tile_off(tc10, tr00) + elem_off(j, i)]);
                                old[9] = data[tile_off(tr01, tc10) + elem_off(i, j)];
                                old[10] = data[tile_off(tr10, tc10) + elem_off(i, j)];
                                old[11] = data[tile_off(tr11, tc10) + elem_off(i, j)];
                                old[12] = conj(data[tile_off(tc11, tr00) + elem_off(j, i)]);
                                old[13] = conj(data[tile_off(tc11, tr01) + elem_off(j, i)]);
                                old[14] = conj(data[tile_off(tc11, tr10) + elem_off(j, i)]);
                                old[15] = data[tile_off(tr11, tc11) + elem_off(i, j)];

                                for (unsigned char col = 0; col < 4; ++col) {
                                    for (unsigned char row = 0; row < 4; ++row) {
                                        for (unsigned char k = 0; k < 4; ++k) {
                                            for (unsigned char l = 0; l < 4; ++l) {
                                                new[row + col * 4] += mat[row + k * 4] * old[k + l * 4] \
                                                                      * conj(mat[col + l * 4]);
                                            }
                                        }
                                    }
                                }

                                data[tile_off(tr00, tc00) + elem_off(i, j)] = new[0];
                                data[tile_off(tr01, tc00) + elem_off(i, j)] = new[1];
                                data[tile_off(tr10, tc00) + elem_off(i, j)] = new[2];
                                data[tile_off(tr11, tc00) + elem_off(i, j)] = new[3];
                                data[tile_off(tc01, tr00) + elem_off(j, i)] = conj(new[4]);
                                data[tile_off(tr01, tc01) + elem_off(i, j)] = new[5];
                                data[tile_off(tr10, tc01) + elem_off(i, j)] = new[6];
                                data[tile_off(tr11, tc01) + elem_off(i, j)] = new[7];
                                data[tile_off(tc10, tr00) + elem_off(j, i)] = conj(new[8]);
                                data[tile_off(tr01, tc10) + elem_off(i, j)] = new[9];
                                data[tile_off(tr10, tc10) + elem_off(i, j)] = new[10];
                                data[tile_off(tr11, tc10) + elem_off(i, j)] = new[11];
                                data[tile_off(tc11, tr00) + elem_off(j, i)] = conj(new[12]);
                                data[tile_off(tc11, tr01) + elem_off(j, i)] = conj(new[13]);
                                data[tile_off(tc11, tr10) + elem_off(j, i)] = conj(new[14]);
                                data[tile_off(tr11, tc11) + elem_off(i, j)] = new[15];
                            }
                        }
                    }

                    else {

                        for (gate_idx_t i = 0; i < TILE_DIM; ++i) {
                            for (gate_idx_t j = 0; j < TILE_DIM; ++j) {
                                cplx_t old[16], new[16] = {0};

                                old[0] = data[tile_off(tr00, tc00) + elem_off(i, j)];
                                old[1] = data[tile_off(tr01, tc00) + elem_off(i, j)];
                                old[2] = data[tile_off(tr10, tc00) + elem_off(i, j)];
                                old[3] = data[tile_off(tr11, tc00) + elem_off(i, j)];
                                old[4] = conj(data[tile_off(tc01, tr00) + elem_off(j, i)]);
                                old[5] = data[tile_off(tr01, tc01) + elem_off(i, j)];
                                old[6] = data[tile_off(tr10, tc01) + elem_off(i, j)];
                                old[7] = data[tile_off(tr11, tc01) + elem_off(i, j)];
                                old[8] = conj(data[tile_off(tc10, tr00) + elem_off(j, i)]);
                                old[9] = conj(data[tile_off(tc10, tr01) + elem_off(j, i)]);
                                old[10] = data[tile_off(tr10, tc10) + elem_off(i, j)];
                                old[11] = data[tile_off(tr11, tc10) + elem_off(i, j)];
                                old[12] = conj(data[tile_off(tc11, tr00) + elem_off(j, i)]);
                                old[13] = conj(data[tile_off(tc11, tr01) + elem_off(j, i)]);
                                old[14] = conj(data[tile_off(tc11, tr10) + elem_off(j, i)]);
                                old[15] = data[tile_off(tr11, tc11) + elem_off(i, j)];

                                for (unsigned char col = 0; col < 4; ++col) {
                                    for (unsigned char row = 0; row < 4; ++row) {
                                        for (unsigned char k = 0; k < 4; ++k) {
                                            for (unsigned char l = 0; l < 4; ++l) {
                                                new[row + col * 4] += mat[row + k * 4] * old[k + l * 4] \
                                                                      * conj(mat[col + l * 4]);
                                            }
                                        }
                                    }
                                }

                                data[tile_off(tr00, tc00) + elem_off(i, j)] = new[0];
                                data[tile_off(tr01, tc00) + elem_off(i, j)] = new[1];
                                data[tile_off(tr10, tc00) + elem_off(i, j)] = new[2];
                                data[tile_off(tr11, tc00) + elem_off(i, j)] = new[3];
                                data[tile_off(tc01, tr00) + elem_off(j, i)] = conj(new[4]);
                                data[tile_off(tr01, tc01) + elem_off(i, j)] = new[5];
                                data[tile_off(tr10, tc01) + elem_off(i, j)] = new[6];
                                data[tile_off(tr11, tc01) + elem_off(i, j)] = new[7];
                                data[tile_off(tc10, tr00) + elem_off(j, i)] = conj(new[8]);
                                data[tile_off(tc10, tr01) + elem_off(j, i)] = conj(new[9]);
                                data[tile_off(tr10, tc10) + elem_off(i, j)] = new[10];
                                data[tile_off(tr11, tc10) + elem_off(i, j)] = new[11];
                                data[tile_off(tc11, tr00) + elem_off(j, i)] = conj(new[12]);
                                data[tile_off(tc11, tr01) + elem_off(j, i)] = conj(new[13]);
                                data[tile_off(tc11, tr10) + elem_off(j, i)] = conj(new[14]);
                                data[tile_off(tr11, tc11) + elem_off(i, j)] = new[15];
                            }
                        }
                    }
                }

                /* If btr==btc, adjoints of the tiles in the upper triangle are exactly the lower triangles involved */
                else {
                    
                    for (gate_idx_t i = 0; i < TILE_DIM; ++i) {
                        for(gate_idx_t j = 0; j <= i; ++j) {
                            cplx_t old[16], new[16] = {0};

                            old[0] = data[tile_off(tr00, tc00) + elem_off(i, j)];
                            old[1] = data[tile_off(tr01, tc00) + elem_off(i, j)];
                            old[2] = data[tile_off(tr10, tc00) + elem_off(i, j)];
                            old[3] = data[tile_off(tr11, tc00) + elem_off(i, j)];
                            old[4] = conj(data[tile_off(tc01, tr00) + elem_off(j, i)]);
                            old[5] = data[tile_off(tr01, tc01) + elem_off(i, j)];
                            old[6] = data[tile_off(tr10, tc01) + elem_off(i, j)];
                            old[7] = data[tile_off(tr11, tc01) + elem_off(i, j)];
                            old[8] = conj(data[tile_off(tc10, tr00) + elem_off(j, i)]);
                            old[9] = conj(data[tile_off(tc10, tr01) + elem_off(j, i)]);
                            old[10] = data[tile_off(tr10, tc10) + elem_off(i, j)];
                            old[11] = data[tile_off(tr11, tc10) + elem_off(i, j)];
                            old[12] = conj(data[tile_off(tc11, tr00) + elem_off(j, i)]);
                            old[13] = conj(data[tile_off(tc11, tr01) + elem_off(j, i)]);
                            old[14] = conj(data[tile_off(tc11, tr10) + elem_off(j, i)]);
                            old[15] = data[tile_off(tr11, tc11) + elem_off(i, j)];

                            for (unsigned char col = 0; col < 4; ++col) {
                                for (unsigned char row = 0; row < 4; ++row) {
                                    for (unsigned char k = 0; k < 4; ++k) {
                                        for (unsigned char l = 0; l < 4; ++l) {
                                            new[row + col * 4] += mat[row + k * 4] * old[k + l * 4] \
                                                                  * conj(mat[col + l * 4]);
                                        }
                                    }
                                }
                            }

                            data[tile_off(tr00, tc00) + elem_off(i, j)] = new[0];
                            data[tile_off(tr00, tc00) + elem_off(j, i)] = conj(new[0]);
                            data[tile_off(tr01, tc00) + elem_off(i, j)] = new[1];
                            data[tile_off(tr01, tc00) + elem_off(j, i)] = conj(new[4]);
                            data[tile_off(tr10, tc00) + elem_off(i, j)] = new[2];
                            data[tile_off(tr10, tc00) + elem_off(j, i)] = conj(new[8]);
                            data[tile_off(tr11, tc00) + elem_off(i, j)] = new[3];
                            data[tile_off(tr11, tc00) + elem_off(j, i)] = conj(new[12]);
                            data[tile_off(tr01, tc01) + elem_off(i, j)] = new[5];
                            data[tile_off(tr01, tc01) + elem_off(j, i)] = conj(new[5]);
                            data[tile_off(tr10, tc01) + elem_off(i, j)] = new[6];
                            data[tile_off(tr10, tc01) + elem_off(j, i)] = conj(new[9]);
                            data[tile_off(tr11, tc01) + elem_off(i, j)] = new[7];
                            data[tile_off(tr11, tc01) + elem_off(j, i)] = conj(new[13]);
                            data[tile_off(tr10, tc10) + elem_off(i, j)] = new[10];
                            data[tile_off(tr10, tc10) + elem_off(j, i)] = conj(new[10]);
                            data[tile_off(tr11, tc10) + elem_off(i, j)] = new[11];
                            data[tile_off(tr11, tc10) + elem_off(j, i)] = conj(new[14]);
                            data[tile_off(tr11, tc11) + elem_off(i, j)] = new[15];
                            data[tile_off(tr11, tc11) + elem_off(j, i)] = conj(new[15]);
                        }
                    }
                }
            }
        }
    }
}
