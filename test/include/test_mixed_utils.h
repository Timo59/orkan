// test_mixed_utils.h - Utility functions for mixed state test harnesses (packed and tiled)

#ifndef TEST_MIXED_UTILS_H
#define TEST_MIXED_UTILS_H

/*
 * =====================================================================================================================
 * includes
 * =====================================================================================================================
 */

#include <assert.h>

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
 * Helper functions
 * =====================================================================================================================
 */

/*
 * @brief   Calculate the index in packed lower-triangular storage for element (row, col)
 *
 * @param[in]   n       Matrix dimension (order)
 * @param[in]   row     Row index (0-indexed, must be >= col)
 * @param[in]   col     Column index (0-indexed)
 *
 * @returns Index in the packed array (column-major order, LAPACK convention)
 *
 * @note    For Hermitian matrix stored in lower-triangular packed format (UPLO='L')
 *          Elements are stored column by column:
 *          Column 0: A(0,0), A(1,0), ..., A(n-1,0)
 *          Column 1: A(1,1), A(2,1), ..., A(n-1,1)
 *          etc.
 *          Element A(i,j) with i >= j is at position: i - j + j*(2n - j + 1)/2
 */
static inline unsigned packed_index(unsigned n, unsigned row, unsigned col) {
    assert(row >= col && "packed_index(): row must be >= col for lower-triangular storage");
    assert(row < n && col < n && "packed_index(): indices out of bounds");
    return row - col + (col * (2 * n - col + 1)) / 2;
}

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
