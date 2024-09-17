/*
 * =====================================================================================================================
 *                                                      Includes
 * =====================================================================================================================
 */

#ifndef CPLXUTIL_H
#include "cplxutil.h"
#endif

/*
 * =====================================================================================================================
 *                                              Imaginary power multiplier
 * =====================================================================================================================
 */

void multiplyByIPower(cplx_t* number, int8_t power) {
    power %= 4;
    int8_t prefactor = (power != 0) ? -((power + 1) % 2) : ((power + 1) % 2);
    int8_t rotation = (power != 3) ? -(power % 2) : (power % 2);

    SETREAL(*number, prefactor * creal(*number) + rotation * cimag(*number));
    SETIMAG(*number, -rotation * creal(*number) + prefactor * cimag(*number));
}

/*
 * =====================================================================================================================
 *                                              Real almost equalities
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

/*
 * =====================================================================================================================
 *                                              Complex almost equalities
 * =====================================================================================================================
 */

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
 *                                              Double prints
 * =====================================================================================================================
 */

void rnumberPrint(double number) {
    printf("%.17g", number);
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

void rmatrixPrint(double matrix[], dim_t dim) {
    printf("[");
    vectorPrint(matrix, dim);
    for (dim_t i = 1; i < dim; ++i) {
        vectorPrint(matrix + (i * dim), dim);
    }
    printf("]\n");
}

/*
 * =====================================================================================================================
 *                                                  Complex prints
 * =====================================================================================================================
 */

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

void cvectorPrint(cplx_t vector[], dim_t dim) {
    printf("[");
    cnumberPrint(vector[0]);
    for (dim_t i = 1; i < dim; ++i) {
        printf(", ");
        cnumberPrint(vector[i]);
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