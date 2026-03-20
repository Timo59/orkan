#ifndef CHANNEL_H
#define CHANNEL_H

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
 * Index type for channel internals
 * =====================================================================================================================
 *
 * idx_t is an unsigned 64-bit type used for ALL index computations within the channel module (loop variables, dim,
 * stride, step, pack_idx results).
 */
typedef uint64_t idx_t;

/*
 * =====================================================================================================================
 * Kraus representation type
 * =====================================================================================================================
 *
 * kraus_t is a struct storing the Kraus representation of a quantum channel acting on n_qubits qubits as a
 * concatenation of the n_terms Kraus matrices in row-major format
 */
typedef struct {
  qubit_t n_qubits;
  uint64_t n_terms;
  cplx_t* data;
} kraus_t;

/*
 * =====================================================================================================================
 * Superoperator matrix representation type
 * =====================================================================================================================
 *
 * superop_t is a struct storing the superoperator matrix representation of a quantum channel acting on n_qubits qubits
 * in row-major format
 */
typedef struct {
  qubit_t n_qubits;
  cplx_t* data;
} superop_t;

/*
 * =====================================================================================================================
 * Bit-insertion helper (shared across channel backends)
 * =====================================================================================================================
 */

/*
 * Minimum dimension to enable OpenMP parallelization (avoid thread overhead for small systems)
 */
#ifndef OMP_THRESHOLD
#define OMP_THRESHOLD 512
#endif

/**
 * @brief Insert a 0 bit at position `pos` in value `val`
 *
 * For k from 0 to dim/2-1, insertBit0(k, target) produces all indices with bit
 * `target` equal to 0. This enables direct enumeration without wasted iterations.
 *
 * Example: pos=2, val=0b101 -> 0b1001 (insert 0 at bit 2)
 */
static inline idx_t insertBit0(idx_t val, qubit_t pos) {
    idx_t mask = ((idx_t)1 << pos) - 1;
    return (val & mask) | ((val & ~mask) << 1);
}

/*
 * =====================================================================================================================
 * Index helpers
 * =====================================================================================================================
 */

/*
 * @brief Returns the index of the first element of column col in packed column-major format
 */
static inline idx_t col_off(idx_t dim, idx_t col) {
  return col * (2 * dim - col + 1) / 2;
}

/*
 * @brief Returns the index of the first element in tile(tr,tc) in packed row-major format
 */
static inline idx_t tile_off(idx_t tr, idx_t tc) {
    return (tr * (tr + 1) / 2 + tc) * TILE_SIZE;
}

/*
 * @brief Returns the index of element(lr,lc) in a tile of size TILE_DIM x TILE_DIM
 */
static inline idx_t elem_off(idx_t lr, idx_t lc) {
    return lr * TILE_DIM + lc;
}

/*
 * =====================================================================================================================
 * Compute superoperator from Kraus representation
 * =====================================================================================================================
 */

/*
 * @brief Returns the superoperator from a Kraus representation
 */
superop_t kraus_to_superop(const kraus_t *kraus);

/*
 * =====================================================================================================================
 * Quantum channels
 * =====================================================================================================================
 */

/*
 * @brief Apply the 1-local quantum channel given in terms of its superoperator to the quantum state
 */
void channel_1q(state_t *state, const superop_t *sop, qubit_t target);

#ifdef __cplusplus
}
#endif

#endif // CHANNEL_H

