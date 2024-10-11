//
// Created by Timo Ziegler on 09.10.24.
//

/*
 * =====================================================================================================================
 *                                                      includes
 * =====================================================================================================================
 */

#ifndef QHIPGCD_H
#include "qhipGcd.h"
#endif

/*
 * =====================================================================================================================
 *                                                  Pauli gates
 * =====================================================================================================================
 */

void applyX_gcd(state_t* state, qubit_t qubit) {
    dim_t flipDistance = POW2(qubit, dim_t);
    dim_t indices = POW2(state->qubits - 1, dim_t);
    dim_t stride = (state->qubits > 6) ? POW2(state->qubits - 5, dim_t) : indices;

    /* Get a concurrent dispatch queue */
    dispatch_queue_t queue = dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0);

    dispatch_apply(indices / stride, queue, ^(size_t i) {
        size_t j = i * stride;
        size_t j_stop = j + stride;

        do {
            dim_t left = (j & (~0U << qubit)) << 1;
            dim_t right = j & ~(~0U << qubit);
            dim_t k = (left | right);
            SWAP(state->vec + k, state->vec + (k + flipDistance), cplx_t);
            ++j;
        }while (j < j_stop);
    });
}

void applyY_gcd(state_t* state, qubit_t qubit) {
    dim_t flipDistance = POW2(qubit, dim_t);
    dim_t indices = POW2(state->qubits - 1, dim_t);
    dim_t stride = (state->qubits > 6) ? POW2(state->qubits - 5, dim_t) : indices;

    /* Get a concurrent dispatch queue */
    dispatch_queue_t queue = dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0);

    dispatch_apply(indices / stride, queue, ^(size_t i) {
        size_t j = i * stride;
        size_t j_stop = j + stride;

        do {
            dim_t left = (j & (~0U << qubit)) << 1;
            dim_t right = j & ~(~0U << qubit);
            dim_t k = (left | right);
            SWAP(state->vec + k, state->vec + (k + flipDistance), cplx_t);
            state->vec[k] *= -I;
            state->vec[k + flipDistance] *= I;
            ++j;
        }while (j < j_stop);
    });
}

void applyZ_gcd(state_t* state, qubit_t qubit) {
    dim_t flipDistance = POW2(qubit, dim_t);
    dim_t indices = POW2(state->qubits - 1, dim_t);
    dim_t stride = (state->qubits > 6) ? POW2(state->qubits - 5, dim_t) : indices;

    /* Get a concurrent dispatch queue */
    dispatch_queue_t queue = dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0);

    dispatch_apply(indices / stride, queue, ^(size_t i) {
        size_t j = i * stride;
        size_t j_stop = j + stride;

        do {
            dim_t left = (j & (~0U << qubit)) << 1;
            left |= (1 << qubit);
            dim_t right = j & ~(~0U << qubit);
            dim_t k = (left | right);
            state->vec[k] *= -1;
            ++j;
        }while (j < j_stop);
    });
}