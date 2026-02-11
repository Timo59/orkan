/**
 * @file gate_packed_3q.c
 * @brief Three-qubit gates for mixed states in packed lower-triangular storage
 *
 * Storage: LAPACK packed lower-triangle, column-major, N(N+1)/2 complex elements.
 *
 * All two-qubit gates iterate over quarter-dim blocks indexed by (br, bc) with
 * br >= bc. Each block addresses 4 rows (r00, r01, r10, r11) × 4 columns,
 * producing 6 element-swap pairs. Three pairs always land in the lower triangle;
 * three may cross the diagonal and require conjugation on read/write.
 *
 * Key pitfalls:
 *   - Upper-triangle access: read/write as conj(data[pack_idx(dim, c, r)])
 *   - Diagonal blocks (bc == br): Hermitian pairs share storage; skip one side
 *   - insertBits2_0() requires lo < hi
 */

#include "gate.h"
#include <complex.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/* Minimum dimension to enable OpenMP parallelization (avoid thread overhead for small systems) */
#define OMP_THRESHOLD 64


/**
 * @brief Compute packed array index for element (r, c) where r >= c
 */
static inline gate_idx_t pack_idx(gate_idx_t dim, gate_idx_t r, gate_idx_t c) {
    return c * (2 * dim - c + 1) / 2 + (r - c);
}

void ccx_packed(state_t *state, const qubit_t ctrl1, const qubit_t ctrl2, const qubit_t target) {
    const gate_idx_t dim = (gate_idx_t)1 << state->qubits;

    cplx_t * restrict data = state->data;
}
