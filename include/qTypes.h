#pragma once

#ifndef QTYPES_H
#define QTYPES_H

/*
 * =====================================================================================================================
 *                                                      includes
 * =====================================================================================================================
 */
#ifndef __COMPLEX_H__
#include <complex.h>
#endif

/*
 * =====================================================================================================================
 *                                                      C++ check
 * =====================================================================================================================
 */

#ifdef __cplusplus
extern "C" {
#endif

/*
 * =====================================================================================================================
 *                                                  type definitions
 * =====================================================================================================================
 */
typedef unsigned long long          complength_t;
typedef double complex              cplx_t;
typedef unsigned short              depth_t;
#ifdef MACOS
#define ACCELERATE_NEW_LAPACK                                   // Required to use cblas_new
#define ACCELERATE_LAPACK_ILP64                                 // __LAPACK_int is a 64-bit integer
#include <vecLib/lapack_types.h>
    typedef __LAPACK_double_complex cplx_t;
    typedef __LAPACK_int            dim_t;
#include <vecLib/blas_new.h>
#else
#define OPENBLAS_USE64BITINT
#include <openblas-pthread/openblas_config.h>
    typedef double complex          cplx_t;
    typedef blasint                 dim_t;
#include <openblas-pthread/f77blas.h>
#endif
typedef enum pauli {
    ID,
    X,
    Y,
    Z
} pauli_t;
typedef unsigned char           qubit_t;

/*
 * =====================================================================================================================
 *                                                      Macros
 * =====================================================================================================================
 */

#define INVSQRT2                0.7071067811865475
#define PI                      3.1415926535897932

#ifndef _SYS_PARAM_H_
#define MIN(a,b)                ((a) < (b) ? (a) : (b))
#define MAX(a,b)                ((a) > (b) ? (a) : (b))
#endif

#define SWAP(a, b, T)           do { register T q; q = *(a); *(a) = *(b); *(b) = q; } while(0)
#define POW2(a, T)              (1 << (T) (a))

#define SETREAL(a, b)           (*((double*) &(a)) = (b))
#define SETIMAG(a, b)           (*(((double*) &(a)) + 1) = (b))

/*
 * =====================================================================================================================
 *                                                      Wrapper
 * =====================================================================================================================
 */
inline cplx_t zdotc(const dim_t *N, cplx_t *X, const dim_t *INCX, cplx_t *Y, const dim_t *INCY) {
#ifdef MACOS
    cplx_t out;
    zdotc_(&out, N, X, INCX, Y, INCY);
    return out;
#else
    return zdotc_(N, X, INCX, Y, INCY);
#endif
}

#ifdef __cplusplus
}
#endif

#endif /* QTYPES_H */
