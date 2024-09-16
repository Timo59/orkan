/*
 * =================================================================================================
 *                                              includes
 * =================================================================================================
 */

#include "linalg.h"

/*
 * =================================================================================================
 *                                              purely real matrix/vector operations
 * =================================================================================================
 */

double* matAdd(const double a[], const double b[], dim_t dim) {
    double* result = (double*) malloc(((mat_t) dim) * dim * sizeof(double));
    for (mat_t i = 0; i < ((mat_t) dim) * dim; ++i) {
        result[i] = a[i] + b[i];
    }
    return result;
}

void matAddInPlace(double a[], const double b[], dim_t dim) {
    for (mat_t i = 0; i < ((mat_t) dim) * dim; ++i) {
        a[i] += b[i];
    }
}

double* vecAdd(const double a[], const double b[], dim_t dim) {
    double* result = (double*) malloc(dim * sizeof(double));
    for (dim_t i = 0; i < dim; ++i) {
        result[i] = a[i] + b[i];
    }
    return result;
}

void vecAddInPlace(double a[], const double b[], dim_t dim) {
    for (dim_t i = 0; i < dim; ++i) {
        a[i] += b[i];
    }
}

double vecProduct(const double a[], const double b[], dim_t dim) {
    double result = 0;
    for (dim_t i = 0; i < dim; ++i) {
        result += a[i] * b[i];
    }
    return result;
}

double* scalarMatMul(double a, const double b[], dim_t dim) {
    double* result = (double*) malloc(((mat_t) dim) * dim * sizeof(double));
    for (mat_t i = 0; i < ((mat_t) dim) * dim; ++i) {
        result[i] = a * b[i];
    }
    return result;
}

void scalarMatMulInPlace(double a, double b[], dim_t dim) {
    for (mat_t i = 0; i < ((mat_t) dim) * dim; ++i) {
        b[i] *= a;
    }
}

double* scalarVecMul(double a, const double b[], dim_t dim) {
    double* result = (double*) malloc(dim * sizeof(double));
    for (dim_t i = 0; i < dim; ++i) {
        result[i] = a * b[i];
    }
    return result;
}

void scalarVecMulInPlace(double a, double b[], dim_t dim) {
    for (dim_t i = 0; i < dim; ++i) {
        b[i] *= a;
    }
}

double* matMul(const double a[], const double b[], dim_t dim) {
    double tmp;
    double* result = (double*) malloc(((mat_t) dim) * dim * sizeof(double));
    for (dim_t i = 0; i < dim; ++i) {
        for (dim_t j = 0; j < dim; ++j) {
            tmp = 0;
            for (dim_t k = 0; k < dim; ++k) {
                tmp += a[((mat_t) i) * dim + k] * b[((mat_t) k) * dim + j];
            }
            result[((mat_t) i) * dim + j] = tmp;
        }
    }
    return result;
}

void matMulInPlace(double a[], const double b[], dim_t dim) {
    double tmp;
    double* tmpRow = (double*) malloc(dim * sizeof(double)); 
    for (dim_t i = 0; i < dim; ++i) {               
        for (dim_t j = 0; j < dim; ++j) {
            tmp = 0;
            for (dim_t k = 0; k < dim; ++k){
                tmp += a[i * dim + k] * b[k * dim + j];
            }
            tmpRow[j] = tmp;
        }
        for (dim_t l = 0; l < dim; ++l){
            a[i * dim + l] = tmpRow[l];
        }
    }
    free(tmpRow);
}

double* matVecMul(const double a[], const double b[], dim_t dim) {
    double tmp;
    double* result = (double*) malloc(dim * sizeof(double));
    for (dim_t i = 0; i < dim; ++i) {
        tmp = 0;
        for (dim_t j = 0; j < dim; ++j) {
            tmp += a[((mat_t) i) * dim + j] * b[j];
        }
        result[i] = tmp;
    }
    return result;
} 

/*
This function computes the kronecker product of two square matrices A,B of dimension dimA and dimA,
respectively.

Input:
    const cplx_t a[]: read only 1d array holding the matrix A as a concatenation of its rows
    const cplx_t b[]: read only 1d array holding the matrix B as a concatenation of its rows
    dim_t dimA: number of rows and columns of matrix A
    dim_t dimB: number of rows and columns of matrix B
Output:
    Complex double pointer to a 1d-array holding the kronecker product of A and B as a concatenation
    of its rows, i.e.,
    [a_{11}b_{11}, a_{11}b_{12}, ..., a_{11}b_{1n}, a_{12}b_{11}, a_{12}b_{12} ..., a_{1m}b_{1n},
    a_{11}b_{21}, a_{11}b_{22}, ..., a_{1m}b_{nn}, ..., a_{21}b_{11}, ..., a_{mm}b_{nn}].
*/
double* kronecker(const double a[], const double b[], dim_t dimA, dim_t dimB) {
    mat_t rowIndex, colIndex, entry;                    // unsigned integers holding the indices
                                                        // of row and column of the kronecker
                                                        // product's entry

    double* result = (double*) malloc(((mat_t) dimA) * dimA * dimB * dimB * sizeof(double));
                                                        // Complex double pointer to a 1d array
                                                        // of the size dimA * dimB * dimA * dimB
                                                        // holding the entries of the kronecker
                                                        // product as a concatenation of its rows

    /*
    Starting in the first row, this loop iterates through A's rows.
    */
    for (mat_t i = 0; i < dimA; ++i) {
        /*
        Starting in the first column, this loop iterates through A's columns.
        */
        for (mat_t j = 0; j < dimA; ++j) {
            /*
            Starting in the first row, this loop iterates through B's rows.
            */
            for (mat_t k = 0; k < dimB; ++k) {
                /*
                Starting in the first column, this loop iterates through B's columns.
                */
                for (mat_t l = 0; l < dimB; ++l) 
                {
                    rowIndex = ((mat_t) dimB) * i + k;          // The kronecker product's current
                                                                // row is A's current row times the
                                                                // number of rows in B plus B's
                                                                // current row

                    colIndex = ((mat_t) dimB) * j + l;          // The kronecker product's current
                                                                // column is is A's current column
                                                                // times the number of columns in B
                                                                // plus B's current column

                    entry = rowIndex * dimA * dimB + colIndex;  // Each row of the kronecker product
                                                                // has dimA * dimB entries, so its
                                                                // concatenated form's entry in row
                                                                // x and column y is x * dimA * dimB
                                                                // blocks away from the first block
                                                                // and y cells away from the
                                                                // beginning of a block

                    result[entry] = a[((mat_t) i) * dimA + j] * b[((mat_t) k) * dimB + l];
                                                                // The kronecker product's entry at
                                                                // at position at the current row
                                                                // and column is the product of the
                                                                // the A's and B's entries
                }
            }
        }
    }
    return result;
}

/*
 * =================================================================================================
 *                                              purely complex matrix/vector operations
 * =================================================================================================
 */

cplx_t* cmatAdd(const cplx_t a[], const cplx_t b[], dim_t dim) {
    cplx_t* result = (cplx_t*) malloc(((mat_t) dim) * dim * sizeof(cplx_t));
    for (mat_t i = 0; i < ((mat_t) dim) * dim; ++i) {
        result[i] = a[i] + b[i];
    }
    return result;
}

void cmatAddInPlace(cplx_t a[], const cplx_t b[], dim_t dim) {
    for (mat_t i = 0; i < ((mat_t) dim) * dim; ++i) {
        a[i] += b[i];
    }
}

cplx_t* cvecAdd(const cplx_t a[], const cplx_t b[], dim_t dim) {
    cplx_t* result = (cplx_t*) malloc(dim * sizeof(cplx_t));
    for (dim_t i = 0; i < dim; ++i) {
        result[i] = a[i] + b[i];
    }
    return result;
}

void cvecAddInPlace(cplx_t a[], const cplx_t b[], dim_t dim) {
    for (dim_t i = 0; i < dim; ++i) {
        a[i] += b[i];
    }
}

cplx_t cvecProduct(const cplx_t a[], const cplx_t b[], dim_t dim) {
    cplx_t result = 0;
    for (dim_t i = 0; i < dim; ++i) {
        result += a[i] * b[i];
    }
    return result;
}

cplx_t cinnerProduct(const cplx_t a[], const cplx_t b[], dim_t dim) {
    cplx_t result = 0;
    for (dim_t i = 0; i < dim; ++i) {
        result += conj(a[i]) * b[i];
    }
    return result;
}

cplx_t* cscalarMatMul(cplx_t a, const cplx_t b[], dim_t dim) {
    cplx_t* result = (cplx_t*) malloc(((mat_t) dim) * dim * sizeof(cplx_t));
    for (mat_t i = 0; i < ((mat_t) dim) * dim; ++i) {
        result[i] = a * b[i];
    }
    return result;
}

void cscalarMatMulInPlace(cplx_t a, cplx_t b[], dim_t dim) {
    for (mat_t i = 0; i < ((mat_t) dim) * dim; ++i) {
        b[i] *= a;
    }
}

cplx_t* cscalarVecMul(cplx_t a, const cplx_t b[], dim_t dim) {
    cplx_t* result = (cplx_t*) malloc(dim * sizeof(cplx_t));
    for (dim_t i = 0; i < dim; ++i) {
        result[i] = a * b[i];
    }
    return result;
}

void cscalarVecMulInPlace(cplx_t a, cplx_t b[], dim_t dim) {
    for (dim_t i = 0; i < dim; ++i) {
        b[i] *= a;
    }
}

cplx_t* cmatMul(const cplx_t a[], const cplx_t b[], dim_t dim) {
    cplx_t tmp;
    cplx_t* result = (cplx_t*) malloc(((mat_t) dim) * dim * sizeof(cplx_t));
    for (dim_t i = 0; i < dim; ++i) {
        for (dim_t j = 0; j < dim; ++j) {
            tmp = 0;
            for (dim_t k = 0; k < dim; ++k) {
                tmp += a[((mat_t) i) * dim + k] * b[((mat_t) k) * dim + j];
            }
            result[((mat_t) i) * dim + j] = tmp;
        }
    }
    return result;
}

void cmatMulInPlace(cplx_t a[], const cplx_t b[], dim_t dim) {
    cplx_t tmp;
    cplx_t* tmpRow = (cplx_t*) malloc(dim * sizeof(cplx_t)); 
    for (dim_t i = 0; i < dim; ++i) {               
        for (dim_t j = 0; j < dim; ++j) {
            tmp = 0;
            for (dim_t k = 0; k < dim; ++k){
                tmp += a[i * dim + k] * b[k * dim + j];
            }
            tmpRow[j] = tmp;
        }
        for (dim_t l = 0; l < dim; ++l){
            a[i * dim + l] = tmpRow[l];
        }
    }
    free(tmpRow);
}

cplx_t* cmatVecMul(const cplx_t a[], const cplx_t b[], dim_t dim) {
    cplx_t tmp;
    cplx_t* result = (cplx_t*) malloc(dim * sizeof(cplx_t));
    for (dim_t i = 0; i < dim; ++i) {
        tmp = 0;
        for (dim_t j = 0; j < dim; ++j) {
            tmp += a[((mat_t) i) * dim + j] * b[j];
        }
        result[i] = tmp;
    }
    return result;
}

void cmatVecMulInPlace(const cplx_t a[], cplx_t v[], dim_t dim) {
    cplx_t* tmp = (cplx_t*) malloc(dim * sizeof(cplx_t));
    for (dim_t j = 0; j < dim; ++j) {
        tmp[j] = v[j];
    }
    for (dim_t i = 0; i < dim; ++i) {
        v[i] = (cplx_t) 0.;
        for (dim_t j = 0; j < dim; ++j) {
            v[i] += a[i * dim + j] * tmp[j];
        }
    }
    free(tmp);
}

/*
This function computes the kronecker product of two square matrices A,B of dimension dimA and dimA,
respectively.

Input:
    const cplx_t a[]: read only 1d array holding the matrix A as a concatenation of its rows
    const cplx_t b[]: read only 1d array holding the matrix B as a concatenation of its rows
    dim_t dimA: number of rows and columns of matrix A
    dim_t dimB: number of rows and columns of matrix B
Output:
    Complex double pointer to a 1d-array holding the kronecker product of A and B as a concatenation
    of its rows, i.e.,
    [a_{11}b_{11}, a_{11}b_{12}, ..., a_{11}b_{1n}, a_{12}b_{11}, a_{12}b_{12} ..., a_{1m}b_{1n},
    a_{11}b_{21}, a_{11}b_{22}, ..., a_{1m}b_{nn}, ..., a_{21}b_{11}, ..., a_{mm}b_{nn}].
*/
cplx_t* ckronecker(const cplx_t a[], const cplx_t b[], dim_t dimA, dim_t dimB)
{
    cplx_t* result = malloc(dimA * dimA * dimB * dimB * sizeof(cplx_t ));
    mat_t dim = dimA * dimB;
    for (mat_t i = 0; i < dimA; ++i) {
        for (mat_t j = 0; j < dimA; ++j) {
            for (mat_t k = 0; k < dimB; ++k) {
                for (mat_t l = 0; l < dimB; ++l) {
                    result[i * dimB * dim + k * dim + j * dimB + l] = a[i * dimA + j] * b[k * dimB + l];
                }
            }
        }
    }
    return result;
}
