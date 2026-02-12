// test_mixed_utils.h - Utility functions for converting between packed and full density matrix formats

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

#ifdef __cplusplus
}
#endif

#endif // TEST_MIXED_UTILS_H
