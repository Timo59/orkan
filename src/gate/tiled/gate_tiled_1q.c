#include "gate_tiled.h"

void single_from_mat(state_t *state, const qubit_t target, const cplx_t *mat) {
    const gate_idx_t dim = (gate_idx_t)1 << state->qubits;
    const gate_idx_t dim_tile = (dim + TILE_DIM - 1) >> LOG_TILE_DIM;
    cplx_t * restrict data = state->data;

    const cplx_t mat00_2 = pow(cabs(mat[0]), 2);
    const cplx_t mat01_2 = pow(cabs(mat[2]), 2);
    const cplx_t mat10_2 = pow(cabs(mat[1]), 2);
    const cplx_t mat11_2 = pow(cabs(mat[3]), 2);

    /* Case 1:
     * All constituents of the same weighted combination are within the same tile
     */
    if (target < LOG_TILE_DIM) {
        const gate_idx_t n_tiles = dim_tile * (dim_tile + 1) / 2;
        const gate_idx_t nbl = (dim < TILE_DIM) ? dim >> 1 : TILE_DIM >> 1;
        const gate_idx_t incr = (gate_idx_t)1 << target;

        #pragma omp parallel for if(dim >= OMP_THRESHOLD)
        for (gate_idx_t i = 0; i < n_tiles; ++i) {
            cplx_t * restrict tile = data + i * TILE_SIZE;

            for (gate_idx_t blr = 0; blr < nbl; ++blr) {
                const gate_idx_t offset_lr0 = insertBit0(blr, target) * TILE_DIM;
                const gate_idx_t offset_lr1 = offset_lr0 + incr * TILE_DIM;

                for (gate_idx_t blc = 0; blc < nbl; ++blc) {
                    const gate_idx_t idx00 = insertBit0(blc, target) + offset_lr0;
                    const gate_idx_t idx01 = idx00 + incr;
                    const gate_idx_t idx10 = insertBit0(blc, target) + offset_lr1;
                    const gate_idx_t idx11 = idx10 + incr;

                    cplx_t new00 = tile[idx00] * mat00_2 \
                                   + tile[idx01] * mat[0] * conj(mat[2]) \
                                   + tile[idx10] * mat[2] * conj(mat[0]) \
                                   + tile[idx11] * mat01_2;

                    cplx_t new01 = tile[idx00] * mat[0] * conj(mat[1]) \
                                   + tile[idx01] * mat[0] * conj(mat[3]) \
                                   + tile[idx10] * mat[2] * conj(mat[1]) \
                                   + tile[idx11] * mat[2] * conj(mat[3]);

                    cplx_t new10 = tile[idx00] * mat[1] * conj(mat[0]) \
                                   + tile[idx01] * mat[1] * conj(mat[2]) \
                                   + tile[idx10] * mat[3] * conj(mat[0]) \
                                   + tile[idx11] * mat[3] * conj(mat[2]);

                    cplx_t new11 = tile[idx00] * mat10_2 \
                                   + tile[idx01] * mat[1] * conj(mat[3]) \
                                   + tile[idx10] * mat[3] * conj(mat[1]) \
                                   + tile[idx11] * mat11_2;

                    tile[idx00] = new00;
                    tile[idx01] = new01;
                    tile[idx10] = new10;
                    tile[idx11] = new11;
                }
            }
        }
    }

    /* Case 2:
     * All constituents are spread over different tiles
     */
    else{
        const gate_idx_t nbt = dim_tile >> 1;
        const qubit_t target_tile = target - LOG_TILE_DIM;
        const gate_idx_t incr = (gate_idx_t)1 << target_tile;

        #pragma omp parallel for schedule(static, 1) if(dim >= OMP_THRESHOLD)
        for (gate_idx_t btr = 0; btr < nbt; ++btr) {
            const gate_idx_t tr0 = insertBit0(btr, target_tile);
            const gate_idx_t tr1 = tr0 | incr;

            for (gate_idx_t btc = 0; btc <= btr; ++btc) {
                const gate_idx_t tc0 = insertBit0(btc, target_tile);
                const gate_idx_t tc1 = tc0 | incr;

                cplx_t * restrict t00 = data + tile_off(tr0, tc0);
                cplx_t * restrict t10 = data + tile_off(tr1, tc0);
                cplx_t * restrict t11 = data + tile_off(tr1, tc1);

                if (tr0 > tc1) {
                    cplx_t * restrict t01 = data + tile_off(tr0, tc1); 

                    for (gate_idx_t i = 0; i < TILE_SIZE; ++i) {
                        cplx_t new00 = t00[i] * mat00_2 \
                                       + t01[i] * mat[0] * conj(mat[2]) \
                                       + t10[i] * mat[2] * conj(mat[0]) \
                                       + t11[i] * mat01_2;

                        cplx_t new01 = t00[i] * mat[0] * conj(mat[1]) \
                                       + t01[i] * mat[0] * conj(mat[3]) \
                                       + t10[i] * mat[2] * conj(mat[1]) \
                                       + t11[i] * mat[2] * conj(mat[3]);

                        cplx_t new10 = t00[i] * mat[1] * conj(mat[0]) \
                                       + t01[i] * mat[1] * conj(mat[2]) \
                                       + t10[i] * mat[3] * conj(mat[0]) \
                                       + t11[i] * mat[3] * conj(mat[2]);

                        cplx_t new11 = t00[i] * mat10_2 \
                                       + t01[i] * mat[1] * conj(mat[3]) \
                                       + t10[i] * mat[3] * conj(mat[1]) \
                                       + t11[i] * mat11_2;

                        t00[i] = new00;
                        t01[i] = new01;
                        t10[i] = new10;
                        t11[i] = new11;
                    }
                }

                else if (btc != btr) {
                    cplx_t * restrict t01 = data + tile_off(tc1, tr0);

                    for (gate_idx_t i = 0; i < TILE_DIM; ++i) {
                        for (gate_idx_t j = 0; j < TILE_DIM; ++j) {
                            cplx_t new00 = t00[elem_off(i, j)] * mat00_2 \
                                           + conj(t01[elem_off(j, i)]) * mat[0] * conj(mat[2]) \
                                           + t10[elem_off(i, j)] * mat[2] * conj(mat[0]) \
                                           + t11[elem_off(i, j)] * mat01_2;

                            cplx_t new01 = t00[elem_off(i, j)] * mat[0] * conj(mat[1]) \
                                           + conj(t01[elem_off(j, i)]) * mat[0] * conj(mat[3]) \
                                           + t10[elem_off(i, j)] * mat[2] * conj(mat[1]) \
                                           + t11[elem_off(i, j)] * mat[2] * conj(mat[3]);

                            cplx_t new10 = t00[elem_off(i, j)] * mat[1] * conj(mat[0]) \
                                           + conj(t01[elem_off(j, i)]) * mat[1] * conj(mat[2]) \
                                           + t10[elem_off(i, j)] * mat[3] * conj(mat[0]) \
                                           + t11[elem_off(i, j)] * mat[3] * conj(mat[2]);

                            cplx_t new11 = t00[elem_off(i, j)] * mat10_2 \
                                           + conj(t01[elem_off(j, i)]) * mat[1] * conj(mat[3]) \
                                           + t10[elem_off(i, j)] * mat[3] * conj(mat[1]) \
                                           + t11[elem_off(i, j)] * mat11_2;

                            t00[elem_off(i, j)] = new00;
                            t01[elem_off(j, i)] = conj(new01);
                            t10[elem_off(i, j)] = new10;
                            t11[elem_off(i, j)] = new11;
                        }
                    }
                }

                else {
                    for (gate_idx_t i = 0; i < TILE_DIM; ++i) {
                        for (gate_idx_t j = 0; j <= i; ++j) {
                            cplx_t new00 = t00[elem_off(i, j)] * mat00_2 \
                                           + conj(t10[elem_off(j, i)]) * mat[0] * conj(mat[2]) \
                                           + t10[elem_off(i, j)] * mat[2] * conj(mat[0]) \
                                           + t11[elem_off(i, j)] * mat01_2;

                            cplx_t new01 = t00[elem_off(i, j)] * mat[0] * conj(mat[1]) \
                                           + conj(t10[elem_off(j, i)]) * mat[0] * conj(mat[3]) \
                                           + t10[elem_off(i, j)] * mat[2] * conj(mat[1]) \
                                           + t11[elem_off(i, j)] * mat[2] * conj(mat[3]);

                            cplx_t new10 = t00[elem_off(i, j)] * mat[1] * conj(mat[0]) \
                                           + conj(t10[elem_off(j, i)]) * mat[1] * conj(mat[2]) \
                                           + t10[elem_off(i, j)] * mat[3] * conj(mat[0]) \
                                           + t11[elem_off(i, j)] * mat[3] * conj(mat[2]);

                            cplx_t new11 = t00[elem_off(i, j)] * mat10_2 \
                                           + conj(t10[elem_off(j, i)]) * mat[1] * conj(mat[3]) \
                                           + t10[elem_off(i, j)] * mat[3] * conj(mat[1]) \
                                           + t11[elem_off(i, j)] * mat11_2;

                            t00[elem_off(i, j)] = new00;
                            t00[elem_off(j, i)] = conj(new00);
                            t10[elem_off(i, j)] = new10;
                            t10[elem_off(j, i)] = conj(new01);
                            t11[elem_off(i, j)] = new11;
                            t11[elem_off(j, i)] = conj(new11);
                        }
                    }
                }
            }
        }
    }
}
