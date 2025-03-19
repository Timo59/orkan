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
 * @param[in]       pqc[]   Parametrized Quantum Blocks; i.e., quantum blocks internally depending on parameter
 * @param[in]       qbs[]   Quantum blocks defining the partial derivative of the PQC
 * @param[in]       obs[]   Observable
 * @param[in,out]   grad*   Array of d DOUBLE PRECISION values. On exit, gradient of the mean value
 */
void gradPQC(state_t* state, const depth_t d, const applyPQB pqc[], const double par[], const applyQB qbs[],
    const double obs[], double* grad) {
    state_t bra;                                    // state, initialized to the input state, evolved with all evolution
    stateInitEmpty(&bra, state->qubits);            // operators and acted on with the observable
    stateInitVector(&bra, state->vec);

    for (depth_t i = 0; i < d; ++i){                                // Apply the PQC to bra
        pqc[i](&bra, par[i]);
    }
    for (dim_t i = 0; i < bra.dim; ++i) {                           // Apply the diagonal observable to bra
        bra.vec[i] *= obs[i];
    }

#pragma omp parallel for default(none) shared(state, d, pqc, par, qbs, grad, bra) num_threads(d)
    for (depth_t i = 0; i < d; ++i) {                               // Iterate the components of the PQC
        state_t ket;                                                // Initialized to the input state; acted on by the
        stateInitEmpty(&ket, state->qubits);                        // derivative of the PQC wrt to the i-th parameter
        stateInitVector(&ket, state->vec);

        for (depth_t j = 0; j < i; ++j) {                           // Apply parametrized blocks up to i-1
            pqc[j](&ket, par[j]);
        }
        qbs[i](&ket);                                               // and apply the evolution operator to ket
        for (depth_t j = i; j < d; ++j) {                           // Evolve ket with rest of the evolution operators
            pqc[j](&ket, par[j]);
        }

        grad[i] = 2 * cimag(stateOverlap(bra, ket));
        stateFreeVector(&ket);
    }
    stateFreeVector(&bra);
}