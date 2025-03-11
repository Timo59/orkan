//
// Created by Timo Ziegler on 17.10.24.
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
 *                                                  Pauli X-gate
 * =====================================================================================================================
 */
void applyX(state_t* state, qubit_t qubit) {
    dim_t blockDistance = POW2(qubit + 1, dim_t);       // Range of indices whose state vector entries after the
                                                        // execution of the single qubit gate will be linear
                                                        // combinations of entries with indices only within the range
    dim_t flipDistance = POW2(qubit, dim_t);            // Distance of indices whose state vector entries after the
                                                        // execution of the single qubit gate will be linear
                                                        // combinations of their entries prior to the execution
    for (dim_t i = 0; i < state->dim; i += blockDistance) {     // Iterate the computational basis states with all zeros
                                                                // to the right and including the target qubit
        for (dim_t j = i; j < i + flipDistance; ++j) {              // Iterate the computational basis states whose
                                                                    // qubits left to and including the target qubit
                                                                    // are constant
            SWAP(state->vec + j, state->vec + (j + flipDistance), cplx_t);      // Swap entry with entry at flipDistance
        }
    }
}

void applyX_blas(state_t* state, qubit_t qubit) {
    __LAPACK_int flipDistance = POW2(qubit, __LAPACK_int);
    __LAPACK_int blockDistance = POW2(qubit + 1, __LAPACK_int);

    for (__LAPACK_int i = 0; i < state->dim; i += blockDistance) {
        cblas_zswap(flipDistance, state->vec + i, 1, state->vec + (i + flipDistance), 1);
    }
}

void applyX_blas_omp(state_t* state, qubit_t qubit) {
    __LAPACK_int flipDistance = POW2(qubit, __LAPACK_int);
    __LAPACK_int blockDistance = POW2(qubit + 1, __LAPACK_int);

#pragma omp parallel for default(none) shared(blockDistance, flipDistance, qubit, state)
    for (__LAPACK_int i = 0; i < state->dim; i += blockDistance) {
        cblas_zswap(flipDistance, state->vec + i, 1, state->vec + (i + flipDistance), 1);
    }
}

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

void applyX_omp(state_t* state, qubit_t qubit) {
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

/*
 * =====================================================================================================================
 *                                                  Pauli Y-gate
 * =====================================================================================================================
 */
void applyY(state_t* state, qubit_t qubit) {
    dim_t blockDistance = POW2(qubit + 1, dim_t);
    dim_t flipDistance = POW2(qubit, dim_t);
    for (dim_t i = 0; i < state->dim; i += blockDistance) {
        for (dim_t j = i; j < i + flipDistance; ++j) {
            SWAP(state->vec + j, state->vec + (j + flipDistance), cplx_t);
            state->vec[j] *= -I;
            state->vec[j + flipDistance] *= I;
        }
    }
}

void applyY_blas(state_t* state, qubit_t qubit) {
    __LAPACK_int flipDistance = POW2(qubit, __LAPACK_int);
    __LAPACK_int blockDistance = POW2(qubit + 1, __LAPACK_int);
    __LAPACK_double_complex plus = I;
    __LAPACK_double_complex minus = -I;

    for (__LAPACK_int i = 0; i < state->dim; i += blockDistance) {
        cblas_zswap(flipDistance, state->vec + i, 1, state->vec + (i + flipDistance), 1);
        cblas_zscal(flipDistance, &minus, state->vec + i, 1);
        cblas_zscal(flipDistance, &plus, state->vec + (i + flipDistance), 1);
    }
}

void applyY_blas_omp(state_t* state, qubit_t qubit) {
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

void applyY_omp(state_t* state, qubit_t qubit) {
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

/*
 * =====================================================================================================================
 *                                                  Pauli Z-gate
 * =====================================================================================================================
 */
void applyZ(state_t* state, qubit_t qubit) {
    dim_t blockDistance = POW2(qubit + 1, dim_t);
    dim_t flipDistance = POW2(qubit, dim_t);
    for (dim_t i = 0; i < state->dim; i += blockDistance) {
        for (dim_t j = i; j < i + flipDistance; ++j) {
            state->vec[j + flipDistance] *= -1;
        }
    }
}

void applyZ_blas(state_t* state, qubit_t qubit) {
    __LAPACK_int flipDistance = POW2(qubit, __LAPACK_int);
    __LAPACK_int blockDistance = POW2(qubit + 1, __LAPACK_int);
    __LAPACK_double_complex minus = -1;

    for (__LAPACK_int i = 0; i < state->dim; i += blockDistance) {
        cblas_zscal(flipDistance, &minus, state->vec + (i + flipDistance), 1);
    }
}

void applyZ_blas_omp(state_t* state, qubit_t qubit) {
    __LAPACK_int flipDistance = POW2(qubit, __LAPACK_int);
    __LAPACK_int blockDistance = POW2(qubit + 1, __LAPACK_int);
    __LAPACK_double_complex minus = -1;

#pragma omp parallel for default(none) shared(blockDistance, flipDistance, minus, qubit, state)
    for (__LAPACK_int i = 0; i < state->dim; i += blockDistance) {
        cblas_zscal(flipDistance, &minus, state->vec + (i + flipDistance), 1);
    }
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

void applyZ_omp(state_t* state, qubit_t qubit) {
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