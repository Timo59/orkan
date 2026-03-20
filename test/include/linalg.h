// linalg.h - Basic linear algebra functions

#ifndef LINALG_H
#define LINALG_H

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
 * Function definitions
 * =====================================================================================================================
 */

/*
 * @brief   Returns the matrix-vector product M*v
 *
 * @param[in]   n   Order of the square matrix and the vector
 * @param[in]   m   Double complex array of size n*n; the square matrix in column major format
 * @param[in]   v   Double complex array of size n; the vector
 *
 * @returns Double complex array of size n; the matrix-vector product w=M*v
 */
cplx_t* zmv(unsigned n, const cplx_t *m, const cplx_t *v);

/*
 * @brief   Returns the similarity transformation U*M*U†
 *
 * @param[in]   n   Order of the square matrices
 * @param[in]   u   Double complex array of size n*n; the matrix U in column major format
 * @param[in]   m   Double complex array of size n*n; the matrix M in column major format
 *
 * @returns Double complex array of size n*n; the result R=U*M*U† in column major format
 */
cplx_t* zumu(unsigned n, const cplx_t *u, const cplx_t *m);

cplx_t* zsumu(unsigned n, unsigned r, const cplx_t **u, const cplx_t *m);

/*
 * @brief   Returns the Kronecker product of the k-by-l matrix A with the m-by-n matrix B.
 *
 * @param[in]   k   Number of rows in A
 * @param[in]   l   Number of columns in A
 * @param[in]   A   Double complex array of size k*l; the k-by-l matrix A in column major format
 * @param[in]   m   Number of rows in B
 * @param[in]   n   Number of columns in B
 * @param[in]   B   Double complex array of size m*n; the m-by-n matrix B in column major format
 *
 * @returns Double complex array of size k*l*m*n; the (k*m)-by-(l*n) Kronecker product matrix of A and B in column major
 *          format.
 */
cplx_t* zkron(unsigned k, unsigned l, const cplx_t A[], unsigned m, unsigned n, const cplx_t B[]);

#ifdef __cplusplus
}
#endif

#endif // LINALG_H
