#ifndef UTILS_H
#define UTILS_H

/*
 * =====================================================================================================================
 * includes
 * =====================================================================================================================
 */

#include <complex.h>
#include <math.h>
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
        const unsigned _len = (_n); \
        printf("["); \
        for (unsigned _i = 0; _i < _len; ++_i) { \
            nprint((_v)[_i]); \
            if (_i < _len - 1) printf(", "); \
        } \
        printf("]\n"); \
    } while(0)

#define mprint(_m, _k, _n) \
    do { \
        const unsigned _rows = (_k); \
        const unsigned _cols = (_n); \
        printf("["); \
        for (unsigned _i = 0; _i < _rows; ++_i) { \
            if (_i != 0) printf(" "); \
            printf("["); \
            for (unsigned _j = 0; _j < _cols; ++_j) { \
                nprint((_m)[_j * _rows + _i]); \
                if (_j < _cols - 1) printf(", "); \
            } \
            printf("]"); \
            if (_i < _rows - 1) printf("\n"); \
        } \
        printf("]\n"); \
    } while(0)

#define mprint_packed(_packed, _n) \
    do { \
        const unsigned _dim = (_n); \
        printf("["); \
        for (unsigned _row = 0; _row < _dim; ++_row) { \
            if (_row != 0) printf(" "); \
            printf("["); \
            for (unsigned _col = 0; _col < _dim; ++_col) { \
                if (_col <= _row) { \
                    unsigned _idx = _col * _dim - _col * (_col + 1) / 2 + _row; \
                    nprint((_packed)[_idx]); \
                } else { \
                    unsigned _idx = _row * _dim - _row * (_row + 1) / 2 + _col; \
                    nprint(conj((_packed)[_idx])); \
                } \
                if (_col < _dim - 1) printf(", "); \
            } \
            printf("]"); \
            if (_row < _dim - 1) printf("\n"); \
        } \
        printf("]\n"); \
    } while(0)

#ifdef __cplusplus
}
#endif

#endif //UTILS_H
