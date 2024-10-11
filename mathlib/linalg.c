/*
 * =====================================================================================================================
 *                                                      Includes
 * =====================================================================================================================
 */

#ifndef LINALG_H
#include "linalg.h"
#endif

/*
 * =====================================================================================================================
 *                                              Real vector operations
 * =====================================================================================================================
 */

double* vecAdd(const double a[], const double b[], unsigned long long dim) {
    double* result = malloc(dim * sizeof(double));
    for (unsigned long long i = 0; i < dim; ++i) {
        result[i] = a[i] + b[i];
    }
    return result;
}

void vecAddInPlace(double a[], const double b[], unsigned long long dim) {
    for (unsigned long long i = 0; i < dim; ++i) {
        a[i] += b[i];
    }
}

double vecProduct(const double a[], const double b[], unsigned long long dim) {
    double result = 0;
    for (unsigned long long i = 0; i < dim; ++i) {
        result += a[i] * b[i];
    }
    return result;
}

double* scalarVecMul(double a, const double b[], unsigned long long dim) {
    double* result = malloc(dim * sizeof(double));
    for (unsigned long long i = 0; i < dim; ++i) {
        result[i] = a * b[i];
    }
    return result;
}

void scalarVecMulInPlace(double a, double b[], unsigned long long dim) {
    for (unsigned long long i = 0; i < dim; ++i) {
        b[i] *= a;
    }
}

/*
 * =====================================================================================================================
 *                                              Complex vector operations
 * =====================================================================================================================
 */

complex double* cscalarVecMul(complex double a, const complex double b[], unsigned long long dim) {
    complex double* result = malloc(dim * sizeof(complex double));
    for (unsigned long long i = 0; i < dim; ++i) {
        result[i] = a * b[i];
    }
    return result;
}

void cscalarVecMulInPlace(complex double a, complex double b[], unsigned long long dim) {
    for (unsigned long long i = 0; i < dim; ++i) {
        b[i] *= a;
    }
}

complex double* cvecAdd(const complex double a[], const complex double b[], unsigned long long dim) {
    complex double* result = malloc(dim * sizeof(complex double));
    for (unsigned long long i = 0; i < dim; ++i) {
        result[i] = a[i] + b[i];
    }
    return result;
}

void cvecAddInPlace(complex double a[], const complex double b[], unsigned long long dim) {
    for (unsigned long long i = 0; i < dim; ++i) {
        a[i] += b[i];
    }
}

complex double cvecProduct(const complex double a[], const complex double b[], unsigned long long dim) {
    complex double result = 0;
    for (unsigned long long int i = 0; i < dim; ++i) {
        result += a[i] * b[i];
    }
    return result;
}

complex double cInner(const complex double a[], const complex double b[], unsigned long long dim) {
    complex double result = 0;
    for (unsigned long long i = 0; i < dim; ++i) {
        result += conj(a[i]) * b[i];
    }
    return result;
}

/*
 * =====================================================================================================================
 *                                              Real matrix operations
 * =====================================================================================================================
 */

void scalarMatMulInPlace(double a, double b[], unsigned long long dim) {
    for (unsigned long long i = 0; i < dim * dim; ++i) {
        b[i] *= a;
    }
}

double* scalarMatMul(double a, const double b[], unsigned long long dim) {
    double* result = malloc(dim * dim * sizeof(double));
    for (unsigned long long i = 0; i < dim * dim; ++i) {
        result[i] = a * b[i];
    }
    return result;
}

double* matAdd(const double a[], const double b[], unsigned long long dim) {
    double* result = malloc(dim * dim * sizeof(double));
    for (unsigned long long i = 0; i < dim * dim; ++i) {
        result[i] = a[i] + b[i];
    }
    return result;
}

void matAddInPlace(double a[], const double b[], unsigned long long dim) {
    for (unsigned long long i = 0; i < dim * dim; ++i) {
        a[i] += b[i];
    }
}

double* matVecMul(const double a[], const double b[], unsigned long long dim) {
    double* result = (double*) malloc(dim * sizeof(double));
    for (unsigned long long i = 0; i < dim; ++i) {
        double tmp = 0;
        for (unsigned long long j = 0; j < dim; ++j) {
            tmp += a[i * dim + j] * b[j];
        }
        result[i] = tmp;
    }
    return result;
}

double* matMul(const double a[], const double b[], unsigned long long dim) {
    double* result = malloc(dim * dim * sizeof(double));
    for (unsigned long long i = 0; i < dim; ++i) {
        for (unsigned long long j = 0; j < dim; ++j) {
            double tmp = 0;
            for (unsigned long long k = 0; k < dim; ++k) {
                tmp += a[i * dim + k] * b[k * dim + j];
            }
            result[i * dim + j] = tmp;
        }
    }
    return result;
}

void matMulInPlace(double a[], const double b[], unsigned long long dim) {
    double* tmpRow = malloc(dim * sizeof(double));
    for (unsigned long long i = 0; i < dim; ++i) {
        for (unsigned long long j = 0; j < dim; ++j) {
            double tmp = 0;
            for (unsigned long long k = 0; k < dim; ++k){
                tmp += a[i * dim + k] * b[k * dim + j];
            }
            tmpRow[j] = tmp;
        }
        for (unsigned long long l = 0; l < dim; ++l){
            a[i * dim + l] = tmpRow[l];
        }
    }
    free(tmpRow);
}

/*
 * This function computes the kronecker product of two square matrices A,B of dimension dimA and dimA, respectively.
 * Input:
 *      const cplx_t a[]: read only 1d array holding the matrix A as a concatenation of its rows
 *      const cplx_t b[]: read only 1d array holding the matrix B as a concatenation of its rows
 *      dim_t dimA: number of rows and columns of matrix A
 *      dim_t dimB: number of rows and columns of matrix
 * Output:
 *      Complex double pointer to a 1d-array holding the kronecker product of A and B as a concatenation
 *      of its rows, i.e.,
 *      [a_{11}b_{11}, a_{11}b_{12}, ..., a_{11}b_{1n}, a_{12}b_{11}, a_{12}b_{12} ..., a_{1m}b_{1n},
 *      a_{11}b_{21}, a_{11}b_{22}, ..., a_{1m}b_{nn}, ..., a_{21}b_{11}, ..., a_{mm}b_{nn}].
 */
double* kronecker(const double a[], const double b[], unsigned long long dimA, unsigned long long dimB) {
    unsigned long long rowIndex, colIndex, entry;
    double* result = malloc(dimA * dimA * dimB * dimB * sizeof(double));
    for (unsigned long long i = 0; i < dimA; ++i) {
        for (unsigned long long j = 0; j < dimA; ++j) {
            for (unsigned long long k = 0; k < dimB; ++k) {
                for (unsigned long long l = 0; l < dimB; ++l)
                {
                    rowIndex = dimB * i + k;
                    colIndex = dimB * j + l;
                    entry = rowIndex * dimA * dimB + colIndex;
                    result[entry] = a[i * dimA + j] * b[k * dimB + l];
                }
            }
        }
    }
    return result;
}

/*
 * =====================================================================================================================
 *                                              Complex matrix operations
 * =====================================================================================================================
 */

complex double* cscalarMatMul(complex double a, const complex double b[], unsigned long long dim) {
    complex double* result = malloc(dim * dim * sizeof(complex double));
    for (unsigned long long i = 0; i < dim * dim; ++i) {
        result[i] = a * b[i];
    }
    return result;
}

void cscalarMatMulInPlace(complex double a, complex double b[], unsigned long long dim) {
    for (unsigned long long i = 0; i < dim * dim; ++i) {
        b[i] *= a;
    }
}

complex double* cmatAdd(const complex double a[], const complex double b[], unsigned long long dim) {
    complex double* result = malloc(dim * dim * sizeof(complex double));
    for (unsigned long long i = 0; i < dim * dim; ++i) {
        result[i] = a[i] + b[i];
    }
    return result;
}

void cmatAddInPlace(complex double a[], const complex double b[], unsigned long long dim) {
    for (unsigned long long i = 0; i < dim * dim; ++i) {
        a[i] += b[i];
    }
}

complex double* cmatVecMul(const complex double a[], const complex double b[], unsigned long long dim) {
    complex double* result = malloc(dim * sizeof(complex double));
    for (unsigned long long i = 0; i < dim; ++i) {
        complex double tmp = 0;
        for (unsigned long long j = 0; j < dim; ++j) {
            tmp += a[i * dim + j] * b[j];
        }
        result[i] = tmp;
    }
    return result;
}

void cmatVecMulInPlace(const complex double a[], complex double v[], unsigned long long dim) {
    complex double* tmp = malloc(dim * sizeof(complex double));
    for (unsigned long long j = 0; j < dim; ++j) {
        tmp[j] = v[j];
    }
    for (unsigned long long i = 0; i < dim; ++i) {
        v[i] = (complex double) 0.;
        for (unsigned long long j = 0; j < dim; ++j) {
            v[i] += a[i * dim + j] * tmp[j];
        }
    }
    free(tmp);
}

complex double* cmatMul(const complex double a[], const complex double b[], unsigned long long dim) {
    complex double* result = malloc(dim * dim * sizeof(complex double));
    for (unsigned long long i = 0; i < dim; ++i) {
        for (unsigned long long j = 0; j < dim; ++j) {
            complex double tmp = 0;
            for (unsigned long long k = 0; k < dim; ++k) {
                tmp += a[i * dim + k] * b[k * dim + j];
            }
            result[i * dim + j] = tmp;
        }
    }
    return result;
}

void cmatMulInPlace(complex double a[], const complex double b[], unsigned long long dim) {
    complex double* tmpRow = malloc(dim * sizeof(complex double));
    for (unsigned long long i = 0; i < dim; ++i) {
        for (unsigned long long j = 0; j < dim; ++j) {
            complex double tmp = 0;
            for (unsigned long long k = 0; k < dim; ++k){
                tmp += a[i * dim + k] * b[k * dim + j];
            }
            tmpRow[j] = tmp;
        }
        for (unsigned long long l = 0; l < dim; ++l){
            a[i * dim + l] = tmpRow[l];
        }
    }
    free(tmpRow);
}

/*
 * This function computes the kronecker product of two square matrices A,B of dimension dimA and dimA, respectively.
 * Input:
 *      const cplx_t a[]: read only 1d array holding the matrix A as a concatenation of its rows
 *      const cplx_t b[]: read only 1d array holding the matrix B as a concatenation of its rows
 *      dim_t dimA: number of rows and columns of matrix A
 *      dim_t dimB: number of rows and columns of matrix B
 * Output:
 *      Complex double pointer to a 1d-array holding the kronecker product of A and B as a concatenation
 *      of its rows, i.e.,
 *      [a_{11}b_{11}, a_{11}b_{12}, ..., a_{11}b_{1n}, a_{12}b_{11}, a_{12}b_{12} ..., a_{1m}b_{1n},
 *      a_{11}b_{21}, a_{11}b_{22}, ..., a_{1m}b_{nn}, ..., a_{21}b_{11}, ..., a_{mm}b_{nn}].
 */
complex double* ckronecker(const complex double a[],
                           const complex double b[],
                           unsigned long long dimA,
                           unsigned long long dimB)
{
    complex double* result = malloc(dimA * dimA * dimB * dimB * sizeof(complex double));
    unsigned long long dim = dimA * dimB;
    for (unsigned long long i = 0; i < dimA; ++i) {
        for (unsigned long long j = 0; j < dimA; ++j) {
            for (unsigned long long k = 0; k < dimB; ++k) {
                for (unsigned long long l = 0; l < dimB; ++l) {
                    result[i * dimB * dim + k * dim + j * dimB + l] = a[i * dimA + j] * b[k * dimB + l];
                }
            }
        }
    }
    return result;
}