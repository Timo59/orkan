#ifndef UTILS_H
#define UTILS_H

/*
 * =====================================================================================================================
 * includes
 * =====================================================================================================================
 */

#include <complex.h>
#include <stdio.h>

/*
 * =====================================================================================================================
 * C++ check
 * =====================================================================================================================
 */

#ifdef __cplusplus
extern "C" {
#endif

/*
 * =====================================================================================================================
 * Print functions
 * =====================================================================================================================
 */

// Print double number with 7 digits after the decimal point
inline void dnprint(const double x) {
    printf("%.7f", x);
}

// Print double complex number with 7 digits after the decimal point
inline void znprint(const double complex x) {
    printf("%.7f", creal(x));
    cimag(x) >= 0 ? printf("+i%.7f", cimag(x)) : printf("-i%.7f", fabs(cimag(x)));
}

#define nprint(_x) \
    _Generic((_x), \
        double: dnprint, \
        double complex: znprint \
    )(_x)

#define vprint(_v, _n) \
    do { \
        printf("["); \
        for (unsigned _i = 0; _i < (_n); ++_i) { \
            nprint((_v)[_i]); \
            if (_i < (_n) - 1) printf(", "); \
        } \
        printf("]\n"); \
    } while(0)

#define mprint(_m, _k, _n) \
    do { \
        printf("["); \
        for (unsigned _i = 0; _i < (_k); ++_i) { \
            if (_i != 0) printf(" "); \
            printf("["); \
            for (unsigned _j = 0; _j < (_n); ++_j) { \
                nprint((_m)[_j * (_k) + _i]); \
                if (_j < (_n) - 1) printf(", "); \
            } \
            printf("]"); \
            if (_i < (_k) - 1) printf("\n"); \
        } \
        printf("]\n"); \
    } while(0)

#ifdef __cplusplus
}
#endif

#endif //UTILS_H
