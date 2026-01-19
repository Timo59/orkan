// qhipster.c - Functions representing fundamental quantum gates applied to pure states, i.e., Hilbert space vectors

#ifndef GATE_H
#include "gate.h"
#endif

#include <math.h>

#include <stdio.h>

#if defined(__APPLE__)
    #include <vecLib/cblas_new.h>
    #include <vecLib/lapack.h>
#elif defined(__linux__)
    #include <cblas.h>
#endif

/*
 * =====================================================================================================================
 * Pauli gates
 * =====================================================================================================================
 */

void x_pure(state_t* state, const qubit_t target) {
    // Swap entries with target = 0 and 1 but all other qubits fixed
    const dim_t dim = POW2(state->qubits, dim_t);   // Number of entries in the state vector
    const dim_t stride = POW2(target, dim_t);    // Distance between elements only differing in the targeted qubit
    const dim_t step = POW2(target + 1, dim_t);  // Size of the mutually independent blocks

    // Iterate blocks with all qubits left to the addressed qubit fixed
    for (dim_t i = 0; i < dim; i += step) {
        cblas_zswap(stride, state->data + i, 1, state->data + i + stride, 1);
    }
}


void y_pure(state_t* state, const qubit_t target) {
    const dim_t dim = POW2(state->qubits, dim_t);   // Number of entries in the state vector
    const dim_t stride = POW2(target, dim_t);    // Distance between elements only differing in the targeted qubit
    const dim_t step = POW2(target + 1, dim_t);  // Size of the mutually independent blocks
    const cplx_t iplus = 0.0 + I*1.0;   // Positive imaginary unit
    const cplx_t iminus = 0.0 - I*1.0;  // Negative imaginary unit

    // Iterate blocks with all qubits left to the addressed qubit fixed
    for (dim_t i = 0; i < dim; i += step) {
        // Swap entries with target bit = 0 and target bit = 1
        cblas_zswap(stride, state->data + i, 1, state->data + i + stride, 1);

        // Multiply the new entries with target bit = 0 by -i
        cblas_zscal(stride, &iminus, state->data + i, 1);

        // Multiply the new entries with target bit = 1 by i
        cblas_zscal(stride, &iplus, state->data + stride + i, 1);
    }
}


void z_pure(state_t* state, const qubit_t target) {
    const dim_t dim = POW2(state->qubits, dim_t);   // Number of entries in the state vector
    const dim_t stride = POW2(target, dim_t);    // Distance between elements only differing in the targeted qubit
    const dim_t step = POW2(target + 1, dim_t);  // Size of the mutually independent blocks
    const cplx_t alpha = -1.0;

    // Iterate blocks with all qubits left to the addressed qubit fixed
    for (dim_t i = 0; i < dim; i += step) {
        // Multiply entries with target bit = 1 by -1
        cblas_zscal(stride, &alpha, state->data + stride + i, 1);
    }
}