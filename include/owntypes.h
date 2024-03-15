#ifndef QOWNTYPES_H
#define QOWNTYPES_H

/*
 * =================================================================================================
 *                                              includes
 * =================================================================================================
 */

#include <complex.h>

/*
 * =================================================================================================
 *                                              C++ check
 * =================================================================================================
 */

#ifdef __cplusplus
extern "C" {
#endif

/*
 * =================================================================================================
 *                                         macros
 * =================================================================================================
 */

#define INVSQRT2                0.7071067811865475

#define MIN(a,b)                ((a) < (b) ? (a) : (b))
#define MAX(a,b)                ((a) > (b) ? (a) : (b))
#define SWAP(a, b, T)           do { register T q; q = *(a); *(a) = *(b); *(b) = q; } while(0)
#define POW2(a, T)              (1 << (T) (a))

#define SETREAL(a, b)           (*((double*) &(a)) = (b))
#define SETIMAG(a, b)           (*(((double*) &(a)) + 1) = (b))

/*
 * =================================================================================================
 *                                              type definitions
 * =================================================================================================
 */

typedef unsigned        complength_t;
typedef unsigned        qubit_t;
typedef unsigned        dim_t;
typedef unsigned        mat_t;
typedef unsigned        depth_t;
typedef unsigned        iter_t;
typedef double _Complex cplx_t;

typedef enum pauli {
    ID,
    X,
    Y,
    Z
} pauli_t;

#ifdef __cplusplus
}
#endif

#endif
