// gate.h   - Functions representing the action of local quantum gates to a general quantum state
//
// It includes <stdio.h> and <stdlib.h> because GATE_VALIDATE needs fprintf() and exit().

#ifndef GATE_H
#define GATE_H

/*
 * =====================================================================================================================
 * Includes
 * =====================================================================================================================
 */
#include "state.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

/*
 * =====================================================================================================================
 * C++ check
 * =====================================================================================================================
 */
#ifdef __cplusplus
extern "C" {
#endif

/*
 * =====================================================================================================================
 * Fortran BLAS function zrot_
 * =====================================================================================================================
 */
#ifdef __linux__
void zrot_(
    const blasint *n,
    openblas_complex_double *x, const blasint *incx,
    openblas_complex_double *y, const blasint *incy,
    const double *c, const openblas_complex_double *s);
#endif

/*
 * =====================================================================================================================
 * Index type for gate internals
 * =====================================================================================================================
 *
 * gate_idx_t is an unsigned 64-bit type used for ALL index computations within the gate
 * module (loop variables, dim, stride, step, pack_idx results). This avoids signed 32-bit
 * overflow that would occur with dim_t (LAPACK int) for 16+ qubit mixed states where
 * dim*(dim+1)/2 exceeds 2^31. dim_t is retained only at LAPACK/BLAS call boundaries.
 */
typedef uint64_t gate_idx_t;

/*
 * =====================================================================================================================
 * Bit-insertion helper (shared across gate backends)
 * =====================================================================================================================
 */

/**
 * @brief Insert a 0 bit at position `pos` in value `val`
 *
 * For k from 0 to dim/2-1, insertBit0(k, target) produces all indices with bit
 * `target` equal to 0. This enables direct enumeration without wasted iterations.
 *
 * Example: pos=2, val=0b101 -> 0b1001 (insert 0 at bit 2)
 */
static inline gate_idx_t insertBit0(gate_idx_t val, qubit_t pos) {
    gate_idx_t mask = ((gate_idx_t)1 << pos) - 1;
    return (val & mask) | ((val & ~mask) << 1);
}

/**
 * @brief Insert 0 bits at two positions in value `val`
 *
 * For two-qubit gates, we iterate over indices where both qubit bits are 0.
 * This helper inserts 0 bits at positions lo and hi (lo < hi).
 */
static inline gate_idx_t insertBits2_0(gate_idx_t val, qubit_t lo, qubit_t hi) {
    gate_idx_t mask_lo = ((gate_idx_t)1 << lo) - 1;
    val = (val & mask_lo) | ((val & ~mask_lo) << 1);
    gate_idx_t mask_hi = ((gate_idx_t)1 << hi) - 1;
    return (val & mask_hi) | ((val & ~mask_hi) << 1);
}

/**
 * @brief Compute packed array index for element (r, c) where r >= c
 *
 * LAPACK packed lower-triangle, column-major: index = c*(2*dim - c + 1)/2 + (r - c)
 */
static inline gate_idx_t pack_idx(gate_idx_t dim, gate_idx_t r, gate_idx_t c) {
    return c * (2 * dim - c + 1) / 2 + (r - c);
}

/*
 * =====================================================================================================================
 * Error handling
 * =====================================================================================================================
 *
 * GATE_VALIDATE calls exit() on failure. This is a deliberate design choice for a
 * scientific computing library where correctness is paramount: passing an out-of-range
 * qubit index or null pointer indicates a programming error (not a recoverable runtime
 * condition), and continuing would silently produce wrong simulation results. The
 * fail-fast behavior ensures such bugs surface immediately during development and testing.
 */

#define GATE_VALIDATE(cond, msg) do {                           \
    if (!(cond)) {                                              \
        fprintf(stderr, "qlib: gate: %s\n", (msg));            \
        exit(EXIT_FAILURE);                                     \
    }                                                           \
} while(0)

/*
 * =====================================================================================================================
 * Pauli gates
 * =====================================================================================================================
 */

void x(state_t* state, qubit_t target);
void y(state_t* state, qubit_t target);
void z(state_t* state, qubit_t target);

/*
 * =====================================================================================================================
 * Clifford gates
 * =====================================================================================================================
 */

void h(state_t* state, qubit_t target);
void s(state_t* state, qubit_t target);
void sdg(state_t* state, qubit_t target);
void t(state_t* state, qubit_t target);
void tdg(state_t* state, qubit_t target);
void hy(state_t* state, qubit_t target);

/*
 * =====================================================================================================================
 * Rotation gates
 * =====================================================================================================================
 */

void rx(state_t* state, qubit_t target, double theta);
void ry(state_t* state, qubit_t target, double theta);
void rz(state_t* state, qubit_t target, double theta);
void p(state_t* state, qubit_t target, double theta);

/*
 * =====================================================================================================================
 * Two-qubit gates
 * =====================================================================================================================
 */

void cx(state_t* state, qubit_t control, qubit_t target);
void cy(state_t* state, qubit_t control, qubit_t target);
void cz(state_t* state, qubit_t control, qubit_t target);
void swap_gate(state_t* state, qubit_t q1, qubit_t q2);

/*
 * =====================================================================================================================
 * Three-qubit gates
 * =====================================================================================================================
 */

void ccx(state_t* state, qubit_t ctrl1, qubit_t ctrl2, qubit_t target);

#ifdef __cplusplus
}
#endif

#endif // GATE_H
