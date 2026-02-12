// test_mixed_utils.h - Utility functions for mixed state test harnesses (packed and tiled)

#ifndef TEST_MIXED_UTILS_H
#define TEST_MIXED_UTILS_H

/*
 * =====================================================================================================================
 * includes
 * =====================================================================================================================
 */

#ifndef Q_TEST_H
#include "test.h"
#endif

#ifndef STATE_H
#include "state.h"
#endif

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
 * Function declarations
 * =====================================================================================================================
 */

/*
 * @brief   Unpacks a packed lower-triangular Hermitian matrix to a full n×n matrix
 *
 * @param[in]   n       Dimension of the matrix (n×n)
 * @param[in]   packed  Packed lower-triangular storage (column-major, length n(n+1)/2)
 *
 * @returns Newly allocated full n×n matrix (column-major), or NULL on failure.
 *          Upper triangle is filled via Hermitian conjugation.
 *          Caller is responsible for freeing the returned pointer.
 */
cplx_t* density_unpack(unsigned n, const cplx_t *packed);


/*
 * @brief   Packs the lower triangle of a full n×n matrix into packed storage
 *
 * @param[in]   n       Dimension of the matrix (n×n)
 * @param[in]   full    Full n×n matrix in column-major layout
 *
 * @returns Newly allocated packed lower-triangular storage (length n(n+1)/2), or NULL on failure.
 *          Caller is responsible for freeing the returned pointer.
 */
cplx_t* density_pack(unsigned n, const cplx_t *full);

/*
 * Maximum qubits for tiled tests. Set to 6 (dim=64) to exercise cross-tile
 * code paths, which require dim > TILE_DIM (32 with LOG_TILE_DIM=5).
 * For 6-qubit mixed states the full density matrix is 64x64 = 4096 elements,
 * and the reference computation (U rho U^dag) runs in O(dim^3) which is
 * fast enough for testing.
 */
#define MAXQUBITS_TILED 6

/*
 * @brief   Create a MIXED_TILED state from a full density matrix
 *
 * Allocates a tiled state and populates it element by element using state_set().
 */
void tiled_state_from_full(state_t *state, unsigned nqubits, const cplx_t *full);

/*
 * @brief   Compare a tiled state against a full reference matrix element by element
 */
void assert_tiled_equals_full(const state_t *state, const cplx_t *ref_full, unsigned dim);

#ifdef __cplusplus
}
#endif

#endif // TEST_MIXED_UTILS_H
