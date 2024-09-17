#ifndef LINALG_H
#define LINALG_H

/*
 * =====================================================================================================================
 *                                                      includes
 * =====================================================================================================================
 */

#ifndef __COMPLEX_H__
#include <complex.h>
#endif

#ifndef _STDLIB_H_
#include <stdlib.h>
#endif

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

double* matAdd(const double a[], const double b[], unsigned long long dim);

void matAddInPlace(double a[], const double b[], unsigned long long dim);

double* vecAdd(const double a[], const double b[], unsigned long long dim);

void vecAddInPlace(double a[], const double b[], unsigned long long dim);

double vecProduct(const double a[], const double b[], unsigned long long dim);

double* scalarMatMul(double a, const double b[], unsigned long long dim);

void scalarMatMulInPlace(double a, double b[], unsigned long long dim);

double* scalarVecMul(double a, const double b[], unsigned long long dim);

void scalarVecMulInPlace(double a, double b[], unsigned long long dim);

double* matMul(const double a[], const double b[], unsigned long long dim);

void matMulInPlace(double a[], const double b[], unsigned long long dim);

double* matVecMul(const double a[], const double b[], unsigned long long dim);

double* kronecker(const double a[], const double b[], unsigned long long dimA, unsigned long long dimB);

/* 
 * =================================================================================================
 *                                              purely complex matrix/vector operations
 * =================================================================================================
 */

complex double* cmatAdd(const complex double a[], const complex double b[], unsigned long long dim);

void cmatAddInPlace(complex double a[], const complex double b[], unsigned long long dim);

complex double * cvecAdd(const complex double a[], const complex double b[], unsigned long long dim);

void cvecAddInPlace(complex double a[], const complex double b[], unsigned long long dim);

complex double cvecProduct(const complex double a[], const complex double b[], unsigned long long dim);

complex double cinnerProduct(const complex double a[], const complex double b[], unsigned long long dim);

complex double* cscalarMatMul(complex double a, const complex double b[], unsigned long long dim);

void cscalarMatMulInPlace(complex double a, complex double b[], unsigned long long dim);

complex double* cscalarVecMul(complex double a, const complex double b[], unsigned long long dim);

void cscalarVecMulInPlace(complex double a, complex double b[], unsigned long long dim);

complex double* cmatMul(const complex double a[], const complex double b[], unsigned long long dim);

void cmatMulInPlace(complex double a[], const complex double b[], unsigned long long dim);

complex double* cmatVecMul(const complex double a[], const complex double b[], unsigned long long dim);

void cmatVecMulInPlace(const complex double a[], complex double v[], unsigned long long dim);

complex double* ckronecker(const complex double a[],
                           const complex double b[],
                           unsigned long long dimA,
                           unsigned long long dimB);

#ifdef __cplusplus
}
#endif

#endif
