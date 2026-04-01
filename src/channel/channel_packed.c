#include "channel.h"
#include "index.h"

#ifndef OMP_THRESHOLD
#define OMP_THRESHOLD 512
#endif

/*
 * Preload 4x4 superoperator into 32 real-valued locals.
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

void channel_packed_1q(state_t *state, const cplx_t * restrict sop, const qubit_t target) {
  const idx_t dim = (idx_t)1 << state->qubits;
  const idx_t incr = (idx_t)1 << target;
  const idx_t n_base = dim >> 1;
  cplx_t * restrict data = state->data;

  SOP_PRELOAD(sop);

  #pragma omp parallel for schedule(dynamic, 1) if(dim >= OMP_THRESHOLD)
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

      /* Three elements are always in the lower triangle */
      const double cr0 = creal(data[idx00]), ci0 = cimag(data[idx00]);
      const double cr2 = creal(data[idx10]), ci2 = cimag(data[idx10]);
      const double cr3 = creal(data[idx11]), ci3 = cimag(data[idx11]);

      if (r0 > c1) {
        /* rho(r0, c1) is in the lower triangle */
        const idx_t idx01 = col1_off + (r0 - c1);
        const double cr1 = creal(data[idx01]), ci1 = cimag(data[idx01]);
        double o0r, o0i, o1r, o1i, o2r, o2i, o3r, o3i;

        SOP_APPLY(cr0, ci0, cr1, ci1, cr2, ci2, cr3, ci3,
                  o0r, o0i, o1r, o1i, o2r, o2i, o3r, o3i);

        data[idx00] = CMPLX(o0r, o0i);
        data[idx01] = CMPLX(o1r, o1i);
        data[idx10] = CMPLX(o2r, o2i);
        data[idx11] = CMPLX(o3r, o3i);
      }
      else if (br != bc) {
        /* rho(r0, c1) is in the upper triangle: read conj(rho(c1, r0)) */
        const idx_t idx01_conj = col_off(dim, r0) + (c1 - r0);
        const double cr1 = creal(data[idx01_conj]), ci1 = -cimag(data[idx01_conj]);
        double o0r, o0i, o1r, o1i, o2r, o2i, o3r, o3i;

        SOP_APPLY(cr0, ci0, cr1, ci1, cr2, ci2, cr3, ci3,
                  o0r, o0i, o1r, o1i, o2r, o2i, o3r, o3i);

        data[idx00]      = CMPLX(o0r, o0i);
        data[idx01_conj] = CMPLX(o1r, -o1i);
        data[idx10]      = CMPLX(o2r, o2i);
        data[idx11]      = CMPLX(o3r, o3i);
      }
      else {
        /* Diagonal block: rho(r0, c1) = conj(rho(r1, c0)) */
        const double cr1 = cr2, ci1 = -ci2;
        double o0r, o0i, o1r, o1i, o2r, o2i, o3r, o3i;

        SOP_APPLY(cr0, ci0, cr1, ci1, cr2, ci2, cr3, ci3,
                  o0r, o0i, o1r, o1i, o2r, o2i, o3r, o3i);

        data[idx00] = CMPLX(o0r, o0i);
        data[idx10] = CMPLX(o2r, o2i);
        data[idx11] = CMPLX(o3r, o3i);
      }
    }
  }
}
