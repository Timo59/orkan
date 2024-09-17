#ifndef CPLXUTIL_H
#define CPLXUTIL_H

/*
 * =====================================================================================================================
 *                                                      Includes
 * =====================================================================================================================
 */

#ifndef QTYPES_H
#include "qTypes.h"
#endif

#include <math.h>

#ifndef __STDBOOL_H
#include <stdbool.h>
#endif

#ifndef _STDIO_H_
#include <stdio.h>
#endif

#ifndef _STDLIB_H_
#include <stdlib.h>
#endif

/*
 * =====================================================================================================================
 *                                                      Macros
 * =====================================================================================================================
 */

#define numberPrint(X) _Generic((X), \
    double: rnumberPrint, \
    cplx_t: cnumberPrint \
)

#define vectorPrint(X, Y) _Generic((X), \
    double*: rvectorPrint, \
    cplx_t*: cvectorPrint \
) ((X), (Y))

#define matrixPrint(X, Y) _Generic((X), \
    double*: rmatrixPrint, \
    cplx_t*: cmatrixPrint \
) (X, Y)

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
 *                                              Imaginary power multiplier
 * =====================================================================================================================
 */

void multiplyByIPower(cplx_t* number, int8_t power);

/*
 * =====================================================================================================================
 *                                              Real almost equalities
 * =====================================================================================================================
 */

bool ralmostEqual(double a, double b, double precision);

bool rvectorAlmostEqual(double v[], double w[], dim_t dim, double precision);

/*
 * =====================================================================================================================
 *                                              Complex almost equalities
 * =====================================================================================================================
 */

bool calmostEqual(cplx_t a, cplx_t b, double precision);

bool cvectorAlmostEqual(cplx_t v[], cplx_t w[], dim_t dim, double precision);

/*
 * =====================================================================================================================
 *                                              Double prints
 * =====================================================================================================================
 */

void rnumberPrint(double number);

void rvectorPrint(double vector[], dim_t dim);

void rmatrixPrint(double matrix[], dim_t dim);

/*
 * =====================================================================================================================
 *                                                  Complex prints
 * =====================================================================================================================
 */

void cnumberPrint(cplx_t number);

void cvectorPrint(cplx_t vector[], dim_t dim);

void cmatrixPrint(cplx_t matrix[], dim_t dim);

#ifdef __cplusplus
}
#endif

#endif
