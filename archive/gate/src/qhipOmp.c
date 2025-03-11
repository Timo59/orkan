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

#ifndef PSTATE_H
#include "pstate.h"
#endif
/*
 * =====================================================================================================================
 *                                                  Pauli gates
 * =====================================================================================================================
 */

void applyX(state_t* state, const qubit_t qubit) {
    dim_t flipDistance = POW2(qubit, dim_t);
    dim_t indices = POW2(state->qubits - 1, dim_t);

#pragma omp parallel for default(none) shared(flipDistance, indices, qubit, state)
    for (dim_t i = 0; i < indices; ++i) {
        dim_t left = (i & (~0U << qubit)) << 1;
        dim_t right = i & ~(~0U << qubit);
        dim_t j = (left | right);
        SWAP(state->vec + j, state->vec + (j + flipDistance), cplx_t);
    }
}

void applyY(state_t* state, const qubit_t qubit) {
    dim_t flipDistance = POW2(qubit, dim_t);
    dim_t indices = POW2(state->qubits - 1, dim_t);

#pragma omp parallel for default(none) shared(flipDistance, indices, qubit, state)
    for (dim_t i = 0; i < indices; ++i) {
        dim_t left = (i & (~0U << qubit)) << 1;
        dim_t right = i & ~(~0U << qubit);
        dim_t j = (left | right);
        SWAP(state->vec + j, state->vec + (j + flipDistance), cplx_t);
        state->vec[j] *= -I;
        state->vec[j + flipDistance] *= I;
    }
}

void applyZ(state_t* state, const qubit_t qubit) {
    dim_t flipDistance = POW2(qubit, dim_t);
    dim_t indices = POW2(state->qubits - 1, dim_t);

#pragma omp parallel for default(none) shared(flipDistance, indices, qubit, state)
    for (dim_t i = 0; i < indices; ++i) {
        dim_t left = (i & (~0U << qubit)) << 1;
        left |= (1 << qubit);
        dim_t right = i & ~(~0U << qubit);
        dim_t j = (left | right);
        state->vec[j] *= -1;
    }
}