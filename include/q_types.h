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
 *  Error codes
 * =====================================================================================================================
 */

typedef enum {
    QS_OK           =  0,   // Success
    QS_ERR_NULL     = -1,   // Null pointer argument
    QS_ERR_OOM      = -2,   // Out of memory
    QS_ERR_QUBIT    = -3,   // Invalid qubit index
    QS_ERR_TYPE     = -4,   // Invalid state type for operation
    QS_ERR_FILE     = -5,   // File I/O error
    QS_ERR_FORMAT   = -6,   // Invalid file format
    QS_ERR_PARAM    = -7,   // Invalid parameter value
} qs_error_t;

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

#ifndef MIN
    #define MIN(a,b)                ((a) < (b) ? (a) : (b))
#endif
#ifndef MAX
    #define MAX(a,b)                ((a) > (b) ? (a) : (b))
#endif

#define SWAP(a, b, T)           do { T q; q = *(a); *(a) = *(b); *(b) = q; } while(0)
#define POW2(a, T)              (((T)1) << (a))

#define SETREAL(a, b)           (*((double*) &(a)) = (b))
#define SETIMAG(a, b)           (*(((double*) &(a)) + 1) = (b))

#ifdef __cplusplus
}
#endif

#endif /* QTYPES_H */
