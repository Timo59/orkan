// linalg.c - Basic linear algebra functions

/*
 * =====================================================================================================================
 * includes
 * =====================================================================================================================
 */
#ifndef LINALG_H
#include "linalg.h"
#endif

#include <stdio.h>
#include <stdlib.h>

#if defined(__APPLE__)
#include <vecLib/cblas_new.h>
#elif defined(__linux__)
#include <cblas.h>
#endif

/*
 * =====================================================================================================================
 * Function definitions
 * =====================================================================================================================
 */

cplx_t* zmv(const unsigned n, const cplx_t *m, const cplx_t *v) {
    cplx_t *out = NULL;

    // Allocate memory for the output
    if (!((out = malloc(n * sizeof (*out))))) {
        fprintf(stderr, "mv(): out allocation failed\n");
        return out;
    }

    // Define multipliers: w = ALPHA*M.v + BETA*w
    const dim_t N = (dim_t) n;
    const cplx_t ALPHA = 1.0 + I*0.0;
    const cplx_t BETA = 0.0 + I*0.0;
    cblas_zgemv(CblasColMajor, CblasNoTrans, N, N, &ALPHA, m, N, v, 1, &BETA, out, 1);

    return out;
}


cplx_t* zkron(const unsigned k, const unsigned l, const cplx_t A[],
             const unsigned m, const unsigned n, const cplx_t B[]) {
    // Allocate memory for the Kronecker product matrix
    cplx_t *out = malloc(k * l * m * n * sizeof(*out));
    if (!out) {
        fprintf(stderr, "kron(): out allocation failed\n");
        return out;
    }

    // Iterate A's columns
    for (unsigned a_col = 0; a_col < l; ++a_col) {
        // Iterate A's rows
        for (unsigned a_row = 0; a_row < k; ++a_row) {
            const cplx_t a = A[a_row + a_col * k];  // (a_row,a_col)-th element of A

            // Iterate B's columns
            for (unsigned b_col = 0; b_col < n; ++b_col) {
                // Iterate B's rows
                for (unsigned b_row = 0; b_row < m; ++b_row) {
                    // b_row is the row in the a_row-th block of rows with each block comprising m rows
                    // b_col is the column in the a_col-th block of columns with each block comprising n columns
                    const unsigned col = b_col + a_col * n;
                    const unsigned row = b_row + a_row * m;
                    out[row + col * k * m] = a * B[b_row + b_col * m];
                }   // b_row
            }   // b_col
        }   // a_row
    }   // a_col

    return out;
}
