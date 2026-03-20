#include "channel.h"

void channel_packed_1q(state_t *state, const cplx_t * restrict sop, const qubit_t target) {
  const idx_t dim = (idx_t)1 << state->qubits;
  const idx_t incr = (idx_t)1 << target;
  const idx_t n_base = dim >> 1;
  cplx_t * restrict data = state->data;

  #pragma omp parallel for schedule(static, 1) if(dim >= OMP_THRESHOLD)
  for (idx_t bc = 0; bc < n_base; ++bc) {
    const idx_t c0 = insertBit0(bc, target);
    const idx_t c1 = c0 | incr;

    const idx_t col0_off = col_off(dim, c0);
    const idx_t col1_off = col_off(dim, c1);

    for (idx_t br = bc; br < n_base; ++br) {
      const idx_t r0 = insertBit0(br, target);
      const idx_t r1 = r0 | incr;

      /* Packed indices for the 2x2 block */
      const idx_t idx00 = col0_off + (r0 - c0);      /* rho(r0, c0) */
      const idx_t idx10 = idx00 + incr;               /* rho(r1, c0) */
      const idx_t idx11 = col1_off + (r1 - c1);      /* rho(r1, c1) */

      cplx_t cur[4], upd[4] = {0};

      /* Three elements are always in the lower triangle */
      cur[0] = data[idx00]; cur[2] = data[idx10]; cur[3] = data[idx11];

      if (r0 > c1) {
        /* rho(r0, c1) is in the lower triangle */
        const idx_t idx01 = col1_off + (r0 - c1);
        cur[1] = data[idx01];

        for (unsigned char i = 0; i < 4; ++i) {
          for (unsigned char j = 0; j < 4; ++j) {
            upd[i] += sop[j + i * 4] * cur[j];
          }
        }

        data[idx00] = upd[0]; data[idx01] = upd[1]; data[idx10] = upd[2]; data[idx11] = upd[3];
      }
      else if (br != bc) {
        /* rho(r0, c1) is in the upper triangle: read conj(rho(c1, r0)) */
        const idx_t idx01_conj = col_off(dim, r0) + (c1 - r0);
        cur[1] = conj(data[idx01_conj]);

        for (unsigned char i = 0; i < 4; ++i) {
          for (unsigned char j = 0; j < 4; ++j) {
            upd[i] += sop[j + i * 4] * cur[j];
          }
        }

        data[idx00] = upd[0]; data[idx01_conj] = conj(upd[1]); data[idx10] = upd[2]; data[idx11] = upd[3];
      }
      else {
        /* Diagonal block: rho(r0, c1) = conj(rho(r1, c0)) */
        cur[1] = conj(cur[2]);

        for (unsigned char i = 0; i < 4; ++i) {
          for (unsigned char j = 0; j < 4; ++j) {
            upd[i] += sop[j + i * 4] * cur[j];
          }
        }

        data[idx00] = upd[0]; data[idx10] = upd[2]; data[idx11] = upd[3];
      }
    }
  }
}
