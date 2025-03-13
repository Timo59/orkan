//
// Created by Timo Ziegler on 12.03.25.
//
/*
 * =====================================================================================================================
 *                                                      Includes
 * =====================================================================================================================
 */
#ifndef EXACT_H
#include "exact.h"
#endif

#ifndef __MATH__
#include <math.h>
#endif

#ifndef __OMP_H
#include <omp.h>
#endif

#ifndef _STDLIB_H_
#include <stdlib.h>
#endif

#ifndef _STDIO_H_
#include <stdio.h>
#endif

#ifdef MACOS
#include <vecLib/blas_new.h>
#include <vecLib/lapack.h>
#else
#include <openblas-pthread/f77blas.h>
#endif
#include "../../internal/include/utils.h"

/*
 * =====================================================================================================================
 *                                                  Function definitions
 * =====================================================================================================================
 */
/*
 * This function returns the mean value of an observable with respect to a quantum state
 *
 * @param[in,out]   state   Quantum state; untouched by the function
 * @param[in]       obs     Diagonal entries of an observanle (is required to be diagonalizable in computational basis)
 *
 * @return  The mean value of the observable with respect to the quantum state with DOUBLE PRECISION
 */
double meanObs(const state_t* state, const double obs[]) {
    double out = 0;
#pragma omp parallel for reduction(+:out) default(none) shared(state, obs)
    for (dim_t i = 0; i < state->dim; ++i) {
        out += pow(cabs(state->vec[i]), 2) * obs[i];
    }
    return out;
}

/*
 * This function returns the gradient of an observable's mean value with respect to a quantum state after the
 * application of a parametrized quantum circuit
 *
 * @param[in,out]   state   Input quantum state to the PQC. On exit, the output of the PQC with the CURRENT parameter
 *                          setting
 * @param[in]       d       Number of Parametrized Quantum Blocks in the PQC
 * @param[in]       pqbs[]  Functions resembling the application of a PQB
 * @param[in]       par[]   Parameters for the PQC
 * @param[in]       qbs[]   Functions resembling the application of the evolution operators (Order must be in accordance
 *                          with pqbs
 * @param[in]       obs[]   Observable
 * @param[in,out]   grad*   On exit, gradient of the mean value
 */
void gradPQC(state_t* state, const depth_t d, const applyPQB pqbs[], const double par[], const applyQB qbs[],
    const double obs[], double** grad) {
    const dim_t incr = 1;
    state_t bra;                                    // state, initialized to the input state, evolved with all evolution
    stateInitEmpty(&bra, state->qubits);            // operators and acted on with the observable
    stateInitVector(&bra, state->vec);

    double* tmp = malloc(d * sizeof (double));
    if (!tmp) {
        fprintf(stderr, "gradPQC: grad reallocation failed\n");
        return;
    }
    *grad = tmp;

    applyPQC(&bra, d, pqbs, par);
    for (dim_t i = 0; i < bra.dim; ++i) {                           // Apply the diagonal observable to bra
        bra.vec[i] *= obs[i];
    }
#pragma omp parallel for default(none) shared(state, d, pqbs, par, qbs, grad, incr, bra) num_threads(d)
    for (depth_t j = 0; j < d; ++j) {                               // Iterate the components of the PQC
        state_t ket;                                                // Initialized to the input state; acted on by the
        stateInitEmpty(&ket, state->qubits);                        // derivative of the PQC wrt to the j-th parameter
        stateInitVector(&ket, state->vec);
        applyPQC(&ket, j, pqbs, par);                               // Apply parametrized blocks up to j
        qbs[j](&ket);                                               // and apply the evolution operator to ket

        applyPQC(&ket, d - j, pqbs + j, par + j);                   // Evolve ket with rest of the evolution operators

        *(*grad + j) = 2 * cimag(stateOverlap(bra, ket));
        stateFreeVector(&ket);
    }
    stateFreeVector(&bra);
}

// void gradPQC(state_t* state, const depth_t d, const applyPQB pqbs[], const double par[], const applyQB qbs[],
//     const double obs[], double** grad) {
//     const dim_t incr = 1;
//     state_t bra;                                    // state, initialized to the input state, evolved with all evolution
//     stateInitEmpty(&bra, state->qubits);            // operators and acted on with the observable
//     stateInitVector(&bra, state->vec);
//
//     state_t ket;                            // state, inheriting each intermediate evolution from tmp, acted on with the
//     stateInitEmpty(&ket, state->qubits);    // respective evolution operator and finally evolved with the remaining
//     // evolution operators
//
//     double* tmp = malloc(d * sizeof (double));
//     if (!tmp) {
//         fprintf(stderr, "gradPQC: grad reallocation failed\n");
//         return;
//     }
//     *grad = tmp;
//
//     applyPQC(&bra, d, pqbs, par);
//     for (dim_t i = 0; i < bra.dim; ++i) {                           // Apply the diagonal observable to bra
//         bra.vec[i] *= obs[i];
//     }
//
//     for (depth_t k = 0; k < d; ++k) {                               // Iterate the components of the PQC
//         pqbs[k](state, par[k]);                                     // Evolve tmp with the current component,
//         stateInitVector(&ket, state->vec);                          // copy its vector to ket
//         qbs[k](state);                                              // and apply the observable to ket
//
//         applyPQC(&ket, d - k - 1, pqbs + k + 1, par + k + 1);
//
//         *(*grad + k) = 2 * cimag(stateOverlap(bra, ket));           // The k-th entry is the overlap of bra and ket
//     }
//     stateFreeVector(&bra);
//     stateFreeVector(&ket);
// }