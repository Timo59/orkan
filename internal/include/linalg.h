#ifndef LINALG_H
#define LINALG_H

/*
 * =====================================================================================================================
 *                                                      Includes
 * =====================================================================================================================
 */

#ifndef __COMPLEX_H__
#include <complex.h>
#endif

#ifdef MACOS
#define ACCELERATE_NEW_LAPACK                                   // Required to use cblas_new
#define ACCELERATE_LAPACK_ILP64                                 // __LAPACK_int is a 64-bit integer
#include <vecLib/lapack_types.h>
typedef __LAPACK_double_complex cplx_t;
typedef __LAPACK_int            dim_t;
#include <vecLib/cblas_new.h>
#include <vecLib/lapack.h>
#else
#define OPENBLAS_USE64BITINT
#include <openblas_config.h>
typedef double complex          cplx_t;
typedef blasint                 dim_t;
#include <cblas.h>
#include <lapacke.h>
#endif

/*
 * =====================================================================================================================
 *                                              C++ check
 * =====================================================================================================================
 */
#ifdef __cplusplus
extern "C" {
#endif

/*
 * =====================================================================================================================
 *                                          Fortran LAPACK function zheev
 * =====================================================================================================================
 */
#ifdef LINUX
void zheev_(
    const char *jobz, const char *uplo,
    const blasint *n,
    openblas_complex_double *a, const blasint *lda,
    const double *w,
    openblas_complex_double *work, const blasint *lwork, const double *rwork,
    const blasint *info
    );
#endif


/*
 * =====================================================================================================================
 *                                              Real vector operations
 * =====================================================================================================================
 */
double* vecAdd(const double a[], const double b[], unsigned long long dim);

void vecAddInPlace(double a[], const double b[], unsigned long long dim);

double vecProduct(const double a[], const double b[], unsigned long long dim);

double* scalarVecMul(double a, const double b[], unsigned long long dim);

void scalarVecMulInPlace(double a, double b[], unsigned long long dim);

/*
 * =====================================================================================================================
 *                                              Complex vector operations
 * =====================================================================================================================
 */
complex double* cscalarVecMul(complex double a, const complex double b[], unsigned long long dim);

void cscalarVecMulInPlace(complex double a, complex double b[], unsigned long long dim);

complex double* cvecAdd(const complex double a[], const complex double b[], unsigned long long dim);

void cvecAddInPlace(complex double a[], const complex double b[], unsigned long long dim);

complex double cvecProduct(const complex double a[], const complex double b[], unsigned long long dim);

complex double cInner(const complex double a[], const complex double b[], unsigned long long dim);

/*
 * =====================================================================================================================
 *                                              Real matrix operations
 * =====================================================================================================================
 */
void scalarMatMulInPlace(double a, double b[], unsigned long long dim);

double* scalarMatMul(double a, const double b[], unsigned long long dim);

double* matAdd(const double a[], const double b[], unsigned long long dim);

void matAddInPlace(double a[], const double b[], unsigned long long dim);

double* matVecMul(const double a[], const double b[], unsigned long long dim);

double* matMul(const double a[], const double b[], unsigned long long dim);

void matMulInPlace(double a[], const double b[], unsigned long long dim);

double* kronecker(const double a[], const double b[], unsigned long long dimA, unsigned long long dimB);

/*
 * =====================================================================================================================
 *                                              Complex matrix operations
 * =====================================================================================================================
 */
complex double* cscalarMatMul(complex double a, const complex double b[], unsigned long long dim);

void cscalarMatMulInPlace(complex double a, complex double b[], unsigned long long dim);

complex double* cmatAdd(const complex double a[], const complex double b[], unsigned long long dim);

void cmatAddInPlace(complex double a[], const complex double b[], unsigned long long dim);

complex double* cmatVecMul(const complex double a[], const complex double b[], unsigned long long dim);

void cmatVecMulInPlace(const complex double a[], complex double v[], unsigned long long dim);

complex double* cmatMul(const complex double a[], const complex double b[], unsigned long long dim);

void cmatMulInPlace(complex double a[], const complex double b[], unsigned long long dim);

complex double* ckronecker(const complex double a[],
                           const complex double b[],
                           unsigned long long dimA,
                           unsigned long long dimB);

double complex* zexpm(double complex* m, const double complex a, const dim_t dim);

#ifdef __cplusplus
}
#endif

#endif
