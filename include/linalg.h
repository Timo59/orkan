#ifndef LINALG_H
#define LINALG_H

/*
 * =================================================================================================
 *                                              includes
 * =================================================================================================
 */

#include "owntypes.h"
#include <complex.h>
#include <math.h>
#include <stdlib.h>

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
 *                                              purely real matrix/vector operations
 * =================================================================================================
 */

double* matAdd(const double a[], const double b[], dim_t dim);

void matAddInPlace(double a[], const double b[], dim_t dim);

double* vecAdd(const double a[], const double b[], dim_t dim);

void vecAddInPlace(double a[], const double b[], dim_t dim);

double vecProduct(const double a[], const double b[], dim_t dim);

double* scalarMatMul(double a, const double b[], dim_t dim);

void scalarMatMulInPlace(double a, double b[], dim_t dim);

double* scalarVecMul(double a, const double b[], dim_t dim);

void scalarVecMulInPlace(double a, double b[], dim_t dim);

double* matMul(const double a[], const double b[], dim_t dim);

void matMulInPlace(double a[], const double b[], dim_t dim);

double* matVecMul(const double a[], const double b[], dim_t dim);

double* kronecker(const double a[], const double b[], dim_t dimA, dim_t dimB);

/* 
 * =================================================================================================
 *                                              purely complex matrix/vector operations
 * =================================================================================================
 */

cplx_t* cmatAdd(const cplx_t a[], const cplx_t b[], dim_t dim);

void cmatAddInPlace(cplx_t a[], const cplx_t b[], dim_t dim);

cplx_t* cvecAdd(const cplx_t a[], const cplx_t b[], dim_t dim);

void cvecAddInPlace(cplx_t a[], const cplx_t b[], dim_t dim);

cplx_t cvecProduct(const cplx_t a[], const cplx_t b[], dim_t dim);

cplx_t cinnerProduct(const cplx_t a[], const cplx_t b[], dim_t dim);

cplx_t* cscalarMatMul(cplx_t a, const cplx_t b[], dim_t dim);

void cscalarMatMulInPlace(cplx_t a, cplx_t b[], dim_t dim);

cplx_t* cscalarVecMul(cplx_t a, const cplx_t b[], dim_t dim);

void cscalarVecMulInPlace(cplx_t a, cplx_t b[], dim_t dim);

cplx_t* cmatMul(const cplx_t a[], const cplx_t b[], dim_t dim);

void cmatMulInPlace(cplx_t a[], const cplx_t b[], dim_t dim);

cplx_t* cmatVecMul(const cplx_t a[], const cplx_t b[], dim_t dim);

void cmatVecMulInPlace(const cplx_t a[], cplx_t v[], dim_t dim);

cplx_t* ckronecker(const cplx_t a[], const cplx_t b[], dim_t dimA, dim_t dimB);

#ifdef __cplusplus
}
#endif

#endif
