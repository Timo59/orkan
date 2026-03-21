#include "channel.h"

/*
 * Preload 4x4 superoperator into 32 real-valued locals.
 * Captures from sop pointer; produces _s0r.._s15r, _s0i.._s15i.
 */
#define SOP_PRELOAD(sop)                                                    \
    const double _s0r  = creal((sop)[ 0]), _s0i  = cimag((sop)[ 0]);       \
    const double _s1r  = creal((sop)[ 1]), _s1i  = cimag((sop)[ 1]);       \
    const double _s2r  = creal((sop)[ 2]), _s2i  = cimag((sop)[ 2]);       \
    const double _s3r  = creal((sop)[ 3]), _s3i  = cimag((sop)[ 3]);       \
    const double _s4r  = creal((sop)[ 4]), _s4i  = cimag((sop)[ 4]);       \
    const double _s5r  = creal((sop)[ 5]), _s5i  = cimag((sop)[ 5]);       \
    const double _s6r  = creal((sop)[ 6]), _s6i  = cimag((sop)[ 6]);       \
    const double _s7r  = creal((sop)[ 7]), _s7i  = cimag((sop)[ 7]);       \
    const double _s8r  = creal((sop)[ 8]), _s8i  = cimag((sop)[ 8]);       \
    const double _s9r  = creal((sop)[ 9]), _s9i  = cimag((sop)[ 9]);       \
    const double _s10r = creal((sop)[10]), _s10i = cimag((sop)[10]);        \
    const double _s11r = creal((sop)[11]), _s11i = cimag((sop)[11]);        \
    const double _s12r = creal((sop)[12]), _s12i = cimag((sop)[12]);        \
    const double _s13r = creal((sop)[13]), _s13i = cimag((sop)[13]);        \
    const double _s14r = creal((sop)[14]), _s14i = cimag((sop)[14]);        \
    const double _s15r = creal((sop)[15]), _s15i = cimag((sop)[15])

/*
 * Apply 4x4 superoperator in real arithmetic.
 * Inputs:  8 doubles (r0,i0, r1,i1, r2,i2, r3,i3)
 * Outputs: 8 doubles (o0r,o0i, o1r,o1i, o2r,o2i, o3r,o3i)
 * Requires SOP_PRELOAD in enclosing scope.
 */
#define SOP_APPLY(r0, i0, r1, i1, r2, i2, r3, i3,                          \
                  o0r, o0i, o1r, o1i, o2r, o2i, o3r, o3i)                   \
    (o0r) = _s0r*(r0)  - _s0i*(i0)  + _s1r*(r1)  - _s1i*(i1)               \
          + _s2r*(r2)  - _s2i*(i2)  + _s3r*(r3)  - _s3i*(i3);              \
    (o0i) = _s0r*(i0)  + _s0i*(r0)  + _s1r*(i1)  + _s1i*(r1)               \
          + _s2r*(i2)  + _s2i*(r2)  + _s3r*(i3)  + _s3i*(r3);              \
    (o1r) = _s4r*(r0)  - _s4i*(i0)  + _s5r*(r1)  - _s5i*(i1)               \
          + _s6r*(r2)  - _s6i*(i2)  + _s7r*(r3)  - _s7i*(i3);              \
    (o1i) = _s4r*(i0)  + _s4i*(r0)  + _s5r*(i1)  + _s5i*(r1)               \
          + _s6r*(i2)  + _s6i*(r2)  + _s7r*(i3)  + _s7i*(r3);              \
    (o2r) = _s8r*(r0)  - _s8i*(i0)  + _s9r*(r1)  - _s9i*(i1)               \
          + _s10r*(r2) - _s10i*(i2) + _s11r*(r3) - _s11i*(i3);             \
    (o2i) = _s8r*(i0)  + _s8i*(r0)  + _s9r*(i1)  + _s9i*(r1)               \
          + _s10r*(i2) + _s10i*(r2) + _s11r*(i3) + _s11i*(r3);             \
    (o3r) = _s12r*(r0) - _s12i*(i0) + _s13r*(r1) - _s13i*(i1)              \
          + _s14r*(r2) - _s14i*(i2) + _s15r*(r3) - _s15i*(i3);             \
    (o3i) = _s12r*(i0) + _s12i*(r0) + _s13r*(i1) + _s13i*(r1)              \
          + _s14r*(i2) + _s14i*(r2) + _s15r*(i3) + _s15i*(r3)

void channel_tiled_1q(state_t *state, const cplx_t * restrict sop, const qubit_t target) {
  const idx_t dim = (idx_t)1 << state->qubits;
  const idx_t dim_tile = (dim + TILE_DIM - 1) >> LOG_TILE_DIM;
  cplx_t * restrict data = state->data;

  SOP_PRELOAD(sop);

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

          const double r0 = creal(tile[idx00]), i0 = cimag(tile[idx00]);
          const double r1 = creal(tile[idx01]), i1 = cimag(tile[idx01]);
          const double r2 = creal(tile[idx10]), i2 = cimag(tile[idx10]);
          const double r3 = creal(tile[idx11]), i3 = cimag(tile[idx11]);
          double o0r, o0i, o1r, o1i, o2r, o2i, o3r, o3i;

          SOP_APPLY(r0, i0, r1, i1, r2, i2, r3, i3,
                    o0r, o0i, o1r, o1i, o2r, o2i, o3r, o3i);

          tile[idx00] = CMPLX(o0r, o0i);
          tile[idx01] = CMPLX(o1r, o1i);
          tile[idx10] = CMPLX(o2r, o2i);
          tile[idx11] = CMPLX(o3r, o3i);
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

    #pragma omp parallel for schedule(dynamic, 1) if(dim >= OMP_THRESHOLD)
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

          #pragma clang loop vectorize(enable)
          for (idx_t i = 0; i < TILE_SIZE; ++i) {
            const double r0 = creal(t00[i]), i0 = cimag(t00[i]);
            const double r1 = creal(t01[i]), i1 = cimag(t01[i]);
            const double r2 = creal(t10[i]), i2 = cimag(t10[i]);
            const double r3 = creal(t11[i]), i3 = cimag(t11[i]);
            double o0r, o0i, o1r, o1i, o2r, o2i, o3r, o3i;

            SOP_APPLY(r0, i0, r1, i1, r2, i2, r3, i3,
                      o0r, o0i, o1r, o1i, o2r, o2i, o3r, o3i);

            t00[i] = CMPLX(o0r, o0i);
            t01[i] = CMPLX(o1r, o1i);
            t10[i] = CMPLX(o2r, o2i);
            t11[i] = CMPLX(o3r, o3i);
          }
        }
        else if (btc != btr) {
          cplx_t * restrict t01 = data + tile_off(tc1, tr0);

          for (idx_t i = 0; i < TILE_DIM; ++i) {
            for (idx_t j = 0; j < TILE_DIM; ++j) {
              const double r0 = creal(t00[elem_off(i,j)]), i0 = cimag(t00[elem_off(i,j)]);
              const double r1 = creal(t01[elem_off(j,i)]), i1 = -cimag(t01[elem_off(j,i)]);
              const double r2 = creal(t10[elem_off(i,j)]), i2 = cimag(t10[elem_off(i,j)]);
              const double r3 = creal(t11[elem_off(i,j)]), i3 = cimag(t11[elem_off(i,j)]);
              double o0r, o0i, o1r, o1i, o2r, o2i, o3r, o3i;

              SOP_APPLY(r0, i0, r1, i1, r2, i2, r3, i3,
                        o0r, o0i, o1r, o1i, o2r, o2i, o3r, o3i);

              t00[elem_off(i,j)] = CMPLX(o0r, o0i);
              t01[elem_off(j,i)] = CMPLX(o1r, -o1i);
              t10[elem_off(i,j)] = CMPLX(o2r, o2i);
              t11[elem_off(i,j)] = CMPLX(o3r, o3i);
            }
          }
        }
        else {
          for (idx_t i = 0; i < TILE_DIM; ++i) {
            for (idx_t j = 0; j <= i; ++j) {
              const double r0 = creal(t00[elem_off(i,j)]), i0 = cimag(t00[elem_off(i,j)]);
              const double r1 = creal(t10[elem_off(j,i)]), i1 = -cimag(t10[elem_off(j,i)]);
              const double r2 = creal(t10[elem_off(i,j)]), i2 = cimag(t10[elem_off(i,j)]);
              const double r3 = creal(t11[elem_off(i,j)]), i3 = cimag(t11[elem_off(i,j)]);
              double o0r, o0i, o1r, o1i, o2r, o2i, o3r, o3i;

              SOP_APPLY(r0, i0, r1, i1, r2, i2, r3, i3,
                        o0r, o0i, o1r, o1i, o2r, o2i, o3r, o3i);

              t00[elem_off(i, j)] = CMPLX(o0r, o0i);
              t00[elem_off(j, i)] = CMPLX(o0r, -o0i);
              t10[elem_off(i, j)] = CMPLX(o2r, o2i);
              t10[elem_off(j, i)] = CMPLX(o1r, -o1i);
              t11[elem_off(i, j)] = CMPLX(o3r, o3i);
              t11[elem_off(j, i)] = CMPLX(o3r, -o3i);
            }
          }
        }
      }
    }
  }
}
