//
// Created by Timo Ziegler on 09.10.24.
//
/*
 * =====================================================================================================================
 *                                                      includes
 * =====================================================================================================================
 */
#ifndef GATELIB_H
#include "gatelib.h"
#endif

/*
 * =====================================================================================================================
 *                                                  Pauli gates
 * =====================================================================================================================
 */

void applyX(state_t* state, qubit_t qubit) {
    __LAPACK_int flipDistance = POW2(qubit, __LAPACK_int);
    __LAPACK_int blockDistance = POW2(qubit + 1, __LAPACK_int);

#pragma omp parallel for default(none) shared(blockDistance, flipDistance, qubit, state)
    for (__LAPACK_int i = 0; i < state->dim; i += blockDistance) {
        cblas_zswap(flipDistance, state->vec + i, 1, state->vec + (i + flipDistance), 1);
    }
}

void applyY(state_t* state, qubit_t qubit) {
    __LAPACK_int flipDistance = POW2(qubit, __LAPACK_int);
    __LAPACK_int blockDistance = POW2(qubit + 1, __LAPACK_int);
    __LAPACK_double_complex plus = I;
    __LAPACK_double_complex minus = -I;

#pragma omp parallel for default(none) shared(blockDistance, flipDistance, minus, plus, qubit, state)
    for (__LAPACK_int i = 0; i < state->dim; i += blockDistance) {
        cblas_zswap(flipDistance, state->vec + i, 1, state->vec + (i + flipDistance), 1);
        cblas_zscal(flipDistance, &minus, state->vec + i, 1);
        cblas_zscal(flipDistance, &plus, state->vec + (i + flipDistance), 1);
    }
}

void applyZ(state_t* state, qubit_t qubit) {
    __LAPACK_int flipDistance = POW2(qubit, __LAPACK_int);
    __LAPACK_int blockDistance = POW2(qubit + 1, __LAPACK_int);
    __LAPACK_double_complex minus = -1;

#pragma omp parallel for default(none) shared(blockDistance, flipDistance, minus, qubit, state)
    for (__LAPACK_int i = 0; i < state->dim; i += blockDistance) {
        cblas_zscal(flipDistance, &minus, state->vec + (i + flipDistance), 1);
    }
}