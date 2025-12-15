#pragma once

#ifndef QTYPES_H
#define QTYPES_H

/*
 * =====================================================================================================================
 * Include
 * =====================================================================================================================
 */

#if defined(__APPLE__)
    #include <vecLib/lapack_types.h>
#elif defined(__linux__)
    #include <openblas_config.h>
#endif

/*
 * =====================================================================================================================
 *  C++ check
 * =====================================================================================================================
 */

#ifdef __cplusplus
extern "C" {
#endif

/*
 * =====================================================================================================================
 *  Type definitions
 * =====================================================================================================================
 */

typedef unsigned char               qubit_t;    // Qubit identifier

#if defined(__APPLE__)
    typedef __LAPACK_double_complex cplx_t;
    typedef __LAPACK_int            dim_t;
#elif defined(__linux__)
    typedef openblas_complex_double cplx_t;
    typedef blasint                 dim_t;
#endif

/*
 * =====================================================================================================================
 * Macros
 * =====================================================================================================================
 */

#ifndef _SYS_PARAM_H_
    #define MIN(a,b)                ((a) < (b) ? (a) : (b))
    #define MAX(a,b)                ((a) > (b) ? (a) : (b))
#endif

#define SWAP(a, b, T)           do { register T q; q = *(a); *(a) = *(b); *(b) = q; } while(0)
#define POW2(a, T)              (1 << (T) (a))

#define SETREAL(a, b)           (*((double*) &(a)) = (b))
#define SETIMAG(a, b)           (*(((double*) &(a)) + 1) = (b))

#ifdef __cplusplus
}
#endif

#endif /* QTYPES_H */
