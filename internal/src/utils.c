//
// Created by Timo Ziegler on 18.09.24.
//

#include "../../test/utils.h"

/*
 * =====================================================================================================================
 *                                                  Test functions
 * =====================================================================================================================
 */

bool ralmostEqual(double a, double b, double precision) {
    return (fabs(a - b) < precision);
}

bool rvectorAlmostEqual(double v[], double w[], dim_t dim, double precision) {
    for (dim_t i = 0; i < dim; ++i) {
        if (!ralmostEqual(v[i], w[i], precision)) {
            return false;
        }
    }
    return true;
}

bool calmostEqual(cplx_t a, cplx_t b, double precision) {
    return (fabs(creal(a) - creal(b)) < precision) && (fabs(cimag(a) - cimag(b)) < precision);
}

bool cvectorAlmostEqual(cplx_t v[], cplx_t w[], dim_t dim, double precision) {
    for (dim_t i = 0; i < dim; ++i) {
        if (!calmostEqual(v[i], w[i], precision)) {
            return false;
        }
    }
    return true;
}

/*
 * =====================================================================================================================
 *                                                      Prints
 * =====================================================================================================================
 */

void rnumberPrint(double number) {
    printf("%.17g", number);
}

void cnumberPrint(cplx_t number) {
    printf("%.17g ", creal(number));
    if (cimag(number) > 0) {
        printf("+ i%.17g", cimag(number));
    }
    else if (cimag(number) < 0) {
        printf("- i%.17g", -cimag(number));
    }
    else {
        printf("+ i0.");
    }
}

void rvectorPrint(double vector[], dim_t dim) {
    printf("[");
    rnumberPrint(vector[0]);
    for (dim_t i = 1; i < dim; ++i) {
        printf(", ");
        rnumberPrint(vector[i]);
    }
    printf("]\n");
}

void cvectorPrint(cplx_t vector[], dim_t dim) {
    printf("[");
    cnumberPrint(vector[0]);
    for (dim_t i = 1; i < dim; ++i) {
        printf(", ");
        cnumberPrint(vector[i]);
    }
    printf("]\n");
}

void rmatrixPrint(double matrix[], dim_t dim) {
    printf("[");
    rvectorPrint(matrix, dim);
    for (dim_t i = 1; i < dim; ++i) {
        rvectorPrint(matrix + (i * dim), dim);
    }
    printf("]\n");
}

void cmatrixPrint(cplx_t matrix[], dim_t dim) {
    printf("[");
    cvectorPrint(matrix, dim);
    for (dim_t i = 1; i < dim; ++i) {
        cvectorPrint(matrix + (i * dim), dim);
    }
    printf("]\n");
}

/*
 * =====================================================================================================================
 *                                                      Test vectors
 * =====================================================================================================================
 */

/*
 * This function generates a set of vectors used to unittest the elementary operations on pure states.
 * Input:
 *      qubit_t qubits:     The number of qubits
 * Output:
 *      An array of pointers to all computational basis states plus one with alternating 0.5 and 0.5I in its entries
*/
cplx_t** generateTestVectors(qubit_t qubits) {
    dim_t dim = POW2(qubits, dim_t);
    cplx_t** vectors = calloc(dim + 1, sizeof(cplx_t*));
    vectors[dim] = calloc(dim, sizeof(*(vectors[dim])));

    for (dim_t i = 0; i < dim; ++i) {
        vectors[i] = calloc(dim, sizeof(*(vectors[i])));
        vectors[i][i] = 1.0 + 0.0 * I;

        vectors[dim][i] = (i % 2) ? 0.5 + 0.0 * I : 0.0 + 0.5 * I;
    }
    return vectors;
}

/*
 * This function frees the memory allocated for the test vectors as well as for the array holding them
 *
 * Input:
 *      cplx_t** vectors:   Array of pointers to the test vectors
 *      qubit_t qubits:     The number of qubits
 * Output:
 *      There is no output value, but the previuosly allocated memory for test vectors is freed
*/
void freeTestVectors(cplx_t** vectors, qubit_t qubits) {
    dim_t dim = POW2(qubits, dim_t);
    for (dim_t i = 0; i < dim + 1; ++i) {
        free(vectors[i]);
    }
    free(vectors);
}