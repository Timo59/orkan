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

void applyX(state_t* state, const qubit_t target) {
    // Check that target is in the system's scope
    if (target >= state->qubits) {
        fprintf(stderr, "applyX(): Target is out of scope; Expected %u, Was %u", state->qubits, target);
        state_free(state);
    }

    // Swap entries with target = 0 and 1 but all other qubits fixed
    const dim_t dim = POW2(state->qubits, dim_t);   // Number of entries in the state vector
    const dim_t stride = POW2(target, dim_t);    // Stride between elements only differing in the targeted qubit
    const dim_t step = POW2(target + 1, dim_t);  // Size of the mutually independent blocks

    // Iterate blocks with all qubits left to the addressed qubit fixed
    for (dim_t i = 0; i < dim; i += step) {
        cblas_zswap(stride, state->data + i, 1, state->data + i + stride, 1);
    }
}


void applyY(state_t* state, const qubit_t target) {
    // Check that target is in the system's scope
    if (target >= state->qubits) {
        fprintf(stderr, "applyY(): Target is out of scope; Expected %u, Was %u", state->qubits, target);
        state_free(state);
    }

    const dim_t dim = POW2(state->qubits, dim_t);
    const dim_t stride = POW2(target, dim_t);
    const dim_t step = POW2(target + 1, dim_t);

    const dim_t incr = 1;   // Increment for qubits right to the target
    const cplx_t minus = -1;
    const double c = 0;
    const cplx_t s = I;

    // |...0...> --> -i |...1...> and |...1...> --> i |...0...>
    for (dim_t i = 0; i < dim; i += step) {
        zrot_(&stride, state->data + i, &incr, state->data + i + stride, &incr, &c, &s);
        cblas_zscal(stride, &minus, state->data + i, 1);
    }
}


void applyZ(state_t* state, const qubit_t target) {
    // Check that target is in the system's scope
    if (target >= state->qubits) {
        fprintf(stderr, "applyX(): Target is out of scope; Expected %u, Was %u", state->qubits, target);
        state_free(state);
    }

    const dim_t dim = POW2(state->qubits, dim_t);
    const dim_t stride = POW2(target, dim_t);
    const dim_t step = POW2(target + 1, dim_t);
    const cplx_t minus = -1;

    // |...0...> --> |...0...> and |...1...> --> - |...1...>
    for (dim_t i = 0; i < dim; i += step) {
        cblas_zscal(stride, &minus, state->data + i + stride, 1);
    }
}
