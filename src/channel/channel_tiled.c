#include "channel.h"

void channel_tiled_1q(state_t *state, const cplx_t * restrict sop, const qubit_t target) {
  const idx_t dim = (idx_t)1 << state->qubits;
  const idx_t dim_tile = (dim + TILE_DIM - 1) >> LOG_TILE_DIM;
  cplx_t * restrict data = state->data;

  /* Case 1:
   * All constituents of the same weighted combination are within the same tile
   */
  if (target < LOG_TILE_DIM) {
    const idx_t n_tiles = dim_tile * (dim_tile + 1) / 2;
    const idx_t nbl = (dim < TILE_DIM) ? dim >> 1 : TILE_DIM >> 1;
    const idx_t incr = (idx_t)1 << target;

    #pragma omp parallel for if(dim >= OMP_THRESHOLD)
    for (idx_t i = 0; i < n_tiles; ++i) {
      cplx_t * restrict tile = data + i * TILE_SIZE;

      for (idx_t blr = 0; blr < nbl; ++blr) {
        const idx_t offset_lr0 = insertBit0(blr, target) * TILE_DIM;
        const idx_t offset_lr1 = offset_lr0 + incr * TILE_DIM;

        for (idx_t blc = 0; blc < nbl; ++blc) {
          const idx_t idx00 = insertBit0(blc, target) + offset_lr0;
          const idx_t idx01 = idx00 + incr;
          const idx_t idx10 = insertBit0(blc, target) + offset_lr1;
          const idx_t idx11 = idx10 + incr;

          cplx_t cur[4] = {tile[idx00], tile[idx01], tile[idx10], tile[idx11]};
          cplx_t upd[4] = {0};

          for (unsigned char row = 0; row < 4; ++row) {
            for (unsigned char col = 0; col < 4; ++col) {
              upd[row] += sop[col + row * 4] * cur[col];
            }
          }

          tile[idx00] = upd[0];
          tile[idx01] = upd[1];
          tile[idx10] = upd[2];
          tile[idx11] = upd[3];
        }
      }
    }
  }

  /* Case 2:
   * All constituents are spread over different tiles
   */
  else {
    const idx_t nbt = dim_tile >> 1;
    const qubit_t target_tile = target - LOG_TILE_DIM;
    const idx_t incr = (idx_t)1 << target_tile;

    #pragma omp parallel for schedule(static, 1) if(dim >= OMP_THRESHOLD)
    for (idx_t btr = 0; btr < nbt; ++btr) {
      const idx_t tr0 = insertBit0(btr, target_tile);
      const idx_t tr1 = tr0 | incr;

      for (idx_t btc = 0; btc <= btr; ++btc) {
        const idx_t tc0 = insertBit0(btc, target_tile);
        const idx_t tc1 = tc0 | incr;

        cplx_t * restrict t00 = data + tile_off(tr0, tc0);
        cplx_t * restrict t10 = data + tile_off(tr1, tc0);
        cplx_t * restrict t11 = data + tile_off(tr1, tc1);

        if (tr0 > tc1) {
          cplx_t * restrict t01 = data + tile_off(tr0, tc1); 

          for (idx_t i = 0; i < TILE_SIZE; ++i) {
            cplx_t cur[4] = {t00[i], t01[i], t10[i], t11[i]};
            cplx_t upd[4] = {0};

            for (unsigned char row = 0; row < 4; ++row) {
              for (unsigned char col = 0; col < 4; ++col) {
                upd[row] += sop[col + row * 4] * cur[col];
              }
            }

            t00[i] = upd[0];
            t01[i] = upd[1];
            t10[i] = upd[2];
            t11[i] = upd[3];
          }
        }
        else if (btc != btr) {
          cplx_t * restrict t01 = data + tile_off(tc1, tr0);

          for (idx_t i = 0; i < TILE_DIM; ++i) {
            for (idx_t j = 0; j < TILE_DIM; ++j) {
              cplx_t cur[4] = {t00[elem_off(i,j)], conj(t01[elem_off(j,i)]), t10[elem_off(i,j)], t11[elem_off(i,j)]};
              cplx_t upd[4] = {0};

              for (unsigned char row = 0; row < 4; ++row) {
                for (unsigned char col = 0; col < 4; ++col) {
                  upd[row] += sop[col + row * 4] * cur[col];
                }
              }

              t00[elem_off(i,j)] = upd[0];
              t01[elem_off(j,i)] = conj(upd[1]);
              t10[elem_off(i,j)] = upd[2];
              t11[elem_off(i,j)] = upd[3];
            }
          }
        }
        else {
          for (idx_t i = 0; i < TILE_DIM; ++i) {
            for (idx_t j = 0; j <= i; ++j) {
              cplx_t cur[4] = {t00[elem_off(i,j)], conj(t10[elem_off(j,i)]), t10[elem_off(i,j)], t11[elem_off(i,j)]};
              cplx_t upd[4] = {0};

              for (unsigned char row = 0; row < 4; ++row) {
                for (unsigned char col = 0; col < 4; ++col) {
                  upd[row] += sop[col + row * 4] * cur[col];
                }
              }

              t00[elem_off(i, j)] = upd[0];
              t00[elem_off(j, i)] = conj(upd[0]);
              t10[elem_off(i, j)] = upd[2];
              t10[elem_off(j, i)] = conj(upd[1]);
              t11[elem_off(i, j)] = upd[3];
              t11[elem_off(j, i)] = conj(upd[3]);
            }
          }
        }
      }
    }
  }
}
