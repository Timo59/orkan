// qhipster.c - Functions representing fundamental quantum gates applied to pure states, i.e., Hilbert space vectors

#ifndef GATE_H
#include "gate.h"
#endif

#include <math.h>
#include <stdio.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#if defined(__APPLE__)
    #include <vecLib/cblas_new.h>
    #include <vecLib/lapack.h>
#elif defined(__linux__)
    #include <cblas.h>
#endif

/* Minimum dimension to enable OpenMP parallelization (avoid thread overhead for small systems)
 * Pure states: O(dim) work per gate, need dim >= 2048 (~11 qubits) to offset thread overhead
 */
#define OMP_THRESHOLD 2048

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

    if (step >= dim && stride >= OMP_THRESHOLD) {
        // Leftmost qubit: single iteration, parallelize element-wise
        #pragma omp parallel for
        for (dim_t j = 0; j < stride; ++j) {
            cplx_t tmp = state->data[j];
            state->data[j] = state->data[j + stride];
            state->data[j + stride] = tmp;
        }
    } else {
        // Normal case: parallelize over independent blocks
        #pragma omp parallel for if(dim >= OMP_THRESHOLD)
        for (dim_t i = 0; i < dim; i += step) {
            cblas_zswap(stride, state->data + i, 1, state->data + i + stride, 1);
        }
    }
}


void y_pure(state_t* state, const qubit_t target) {
    const dim_t dim = POW2(state->qubits, dim_t);   // Number of entries in the state vector
    const dim_t stride = POW2(target, dim_t);    // Distance between elements only differing in the targeted qubit
    const dim_t step = POW2(target + 1, dim_t);  // Size of the mutually independent blocks
    const cplx_t iplus = 0.0 + I*1.0;   // Positive imaginary unit
    const cplx_t iminus = 0.0 - I*1.0;  // Negative imaginary unit

    if (step >= dim && stride >= OMP_THRESHOLD) {
        // Leftmost qubit: single iteration, parallelize element-wise
        #pragma omp parallel for
        for (dim_t j = 0; j < stride; ++j) {
            cplx_t tmp = state->data[j];
            state->data[j] = state->data[j + stride] * iminus;
            state->data[j + stride] = tmp * iplus;
        }
    } else {
        // Normal case: parallelize over independent blocks
        #pragma omp parallel for if(dim >= OMP_THRESHOLD)
        for (dim_t i = 0; i < dim; i += step) {
            // Swap entries with target bit = 0 and target bit = 1
            cblas_zswap(stride, state->data + i, 1, state->data + i + stride, 1);

            // Multiply the new entries with target bit = 0 by -i
            cblas_zscal(stride, &iminus, state->data + i, 1);

            // Multiply the new entries with target bit = 1 by i
            cblas_zscal(stride, &iplus, state->data + stride + i, 1);
        }
    }
}


void z_pure(state_t* state, const qubit_t target) {
    const dim_t dim = POW2(state->qubits, dim_t);   // Number of entries in the state vector
    const dim_t stride = POW2(target, dim_t);    // Distance between elements only differing in the targeted qubit
    const dim_t step = POW2(target + 1, dim_t);  // Size of the mutually independent blocks

    if (step >= dim && stride >= OMP_THRESHOLD) {
        // Leftmost qubit: single iteration, parallelize element-wise
        #pragma omp parallel for
        for (dim_t j = 0; j < stride; ++j) {
            state->data[j + stride] = -state->data[j + stride];
        }
    } else {
        // Normal case: parallelize over independent blocks
        const cplx_t alpha = -1.0;
        #pragma omp parallel for if(dim >= OMP_THRESHOLD)
        for (dim_t i = 0; i < dim; i += step) {
            // Multiply entries with target bit = 1 by -1
            cblas_zscal(stride, &alpha, state->data + stride + i, 1);
        }
    }
}


/*
 * =====================================================================================================================
 * Clifford gates
 * =====================================================================================================================
 */

void h_pure(state_t* state, const qubit_t target) {
    const dim_t dim = POW2(state->qubits, dim_t);
    const dim_t stride = POW2(target, dim_t);
    const dim_t step = POW2(target + 1, dim_t);
    const double inv_sqrt2 = M_SQRT1_2;  // 1/√2

    // Hadamard: |0⟩ → (|0⟩+|1⟩)/√2, |1⟩ → (|0⟩-|1⟩)/√2
    // For each pair (a,b): new_a = (a+b)/√2, new_b = (a-b)/√2
    if (step >= dim && stride >= OMP_THRESHOLD) {
        // Leftmost qubit: single outer iteration, parallelize over pairs
        #pragma omp parallel for
        for (dim_t j = 0; j < stride; ++j) {
            cplx_t a = state->data[j];
            cplx_t b = state->data[j + stride];
            state->data[j] = (a + b) * inv_sqrt2;
            state->data[j + stride] = (a - b) * inv_sqrt2;
        }
    } else {
        // Normal case: parallelize over independent blocks
        #pragma omp parallel for if(dim >= OMP_THRESHOLD)
        for (dim_t i = 0; i < dim; i += step) {
            for (dim_t j = 0; j < stride; ++j) {
                cplx_t a = state->data[i + j];
                cplx_t b = state->data[i + j + stride];
                state->data[i + j] = (a + b) * inv_sqrt2;
                state->data[i + j + stride] = (a - b) * inv_sqrt2;
            }
        }
    }
}


void s_pure(state_t* state, const qubit_t target) {
    const dim_t dim = POW2(state->qubits, dim_t);
    const dim_t stride = POW2(target, dim_t);
    const dim_t step = POW2(target + 1, dim_t);
    const cplx_t alpha = I;  // i

    // S gate: |0⟩ → |0⟩, |1⟩ → i|1⟩
    if (step >= dim && stride >= OMP_THRESHOLD) {
        // Leftmost qubit: single iteration, parallelize element-wise
        #pragma omp parallel for
        for (dim_t j = 0; j < stride; ++j) {
            state->data[j + stride] *= alpha;
        }
    } else {
        // Normal case: parallelize over independent blocks
        #pragma omp parallel for if(dim >= OMP_THRESHOLD)
        for (dim_t i = 0; i < dim; i += step) {
            cblas_zscal(stride, &alpha, state->data + stride + i, 1);
        }
    }
}


void sdg_pure(state_t* state, const qubit_t target) {
    const dim_t dim = POW2(state->qubits, dim_t);
    const dim_t stride = POW2(target, dim_t);
    const dim_t step = POW2(target + 1, dim_t);
    const cplx_t alpha = -I;  // -i

    // S† gate: |0⟩ → |0⟩, |1⟩ → -i|1⟩
    if (step >= dim && stride >= OMP_THRESHOLD) {
        // Leftmost qubit: single iteration, parallelize element-wise
        #pragma omp parallel for
        for (dim_t j = 0; j < stride; ++j) {
            state->data[j + stride] *= alpha;
        }
    } else {
        // Normal case: parallelize over independent blocks
        #pragma omp parallel for if(dim >= OMP_THRESHOLD)
        for (dim_t i = 0; i < dim; i += step) {
            cblas_zscal(stride, &alpha, state->data + stride + i, 1);
        }
    }
}


void t_pure(state_t* state, const qubit_t target) {
    const dim_t dim = POW2(state->qubits, dim_t);
    const dim_t stride = POW2(target, dim_t);
    const dim_t step = POW2(target + 1, dim_t);
    const cplx_t alpha = M_SQRT1_2 + M_SQRT1_2 * I;  // e^(iπ/4) = (1+i)/√2

    // T gate: |0⟩ → |0⟩, |1⟩ → e^(iπ/4)|1⟩
    if (step >= dim && stride >= OMP_THRESHOLD) {
        // Leftmost qubit: single iteration, parallelize element-wise
        #pragma omp parallel for
        for (dim_t j = 0; j < stride; ++j) {
            state->data[j + stride] *= alpha;
        }
    } else {
        // Normal case: parallelize over independent blocks
        #pragma omp parallel for if(dim >= OMP_THRESHOLD)
        for (dim_t i = 0; i < dim; i += step) {
            cblas_zscal(stride, &alpha, state->data + stride + i, 1);
        }
    }
}


void tdg_pure(state_t* state, const qubit_t target) {
    const dim_t dim = POW2(state->qubits, dim_t);
    const dim_t stride = POW2(target, dim_t);
    const dim_t step = POW2(target + 1, dim_t);
    const cplx_t alpha = M_SQRT1_2 - M_SQRT1_2 * I;  // e^(-iπ/4) = (1-i)/√2

    // T† gate: |0⟩ → |0⟩, |1⟩ → e^(-iπ/4)|1⟩
    if (step >= dim && stride >= OMP_THRESHOLD) {
        // Leftmost qubit: single iteration, parallelize element-wise
        #pragma omp parallel for
        for (dim_t j = 0; j < stride; ++j) {
            state->data[j + stride] *= alpha;
        }
    } else {
        // Normal case: parallelize over independent blocks
        #pragma omp parallel for if(dim >= OMP_THRESHOLD)
        for (dim_t i = 0; i < dim; i += step) {
            cblas_zscal(stride, &alpha, state->data + stride + i, 1);
        }
    }
}