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

#ifndef UTILS_H
#include "utils.h"
#endif

#include <unistd.h>

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
void gradPQC(state_t* state, const depth_t d, const applyPQB pqc[], const double par[], const applyQB qb[],
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

#pragma omp parallel for default(none) shared(state, d, pqc, par, qb, grad, bra) num_threads(d)
    for (depth_t i = 0; i < d; ++i) {                               // Iterate the components of the PQC
        state_t ket;                                                // Initialized to the input state; acted on by the
        stateInitEmpty(&ket, state->qubits);                        // derivative of the PQC wrt to the i-th parameter
        stateInitVector(&ket, state->vec);

        for (depth_t j = 0; j < i; ++j) {                           // Apply parametrized blocks up to i-1
            pqc[j](&ket, par[j]);
        }
        qb[i](&ket);                                               // and apply the evolution operator to ket
        for (depth_t j = i; j < d; ++j) {                           // Evolve ket with rest of the evolution operators
            pqc[j](&ket, par[j]);
        }

        grad[i] = 2 * cimag(stateOverlap(bra, ket));
        stateFreeVector(&ket);
    }
    stateFreeVector(&bra);
}

/*
 * This function calculates the moment matrices of a set of observables wrt to a sequence of LCU.
 *
 *  @param[in, out] state:      Input quantum state to the PQC. On exit, the output of the PQC up to link with the
 *                               CURRENT parameter setting
 *  @param[in]      obsc:       Number of moment matrices
 *  @param[in, out] obs[]:      Set of observables constituting the moment matrices
 *  @param[in]      circdepth:  Number of LCU channels
 *  @param[in, out] c[]:        Coefficients for LCU; concatenation for all channels
 *  @param[in]      uc:         Number of unitary operators in each LCU TODO: Allow for different numbers across channel
 *  @param[in, out] u[]:        Set of unitary operators; concatenation for all LCU
 *  @param[in, out] link:       Element of the LCU sequence constituting the basis for the moment matrices
 *  @param[in, out] momMat[]:   Array of moment matrices in packed column major form
 *
 * Output:
 *  An array of moment matrices, whose entries are the overlaps of two quantum states with the corresponding observable
 *  applied. In particular, all LCU preceding the one specified by link are applied with their coefficients setting on
 *  input. Then, the set of unitary operators at link are applied one by one to the quantum state to give the states for
 *  each column and row, respectively. The remaining LCU are also applied with their coefficients determined on inpt.
 */

void mmseq(state_t* state,
           const depth_t obsc,
           const applyQB obs[],
           const depth_t circdepth,
           const cplx_t c[],
           const depth_t uc,
           const applyQB u[],
           const depth_t link,
           cplx_t* momMat[])
{
    for (depth_t k = 0; k < link; ++k) {                            // Apply all LCU channels prior to 'link'
        lcQB(state, uc, u + k * uc, c + k * uc);
    }

//#pragma omp parallel for default(none) shared(state, obsc, obs, circdepth, c, uc, u, link, momMat)
    for (depth_t j = 0; j < uc; ++j) {                              // Iterate the moment matrices' columns
        state_t bra, ket;
        stateInitEmpty(&ket, state->qubits);
        stateInitVector(&ket, state->vec);
        const unsigned short offset = j * uc - j * (j - 1) / 2;     // Elements prior to the j-th column

        u[j + link * uc](&ket);                                     // Apply the search unitary according to the link

        for (depth_t k = link + 1; k < circdepth; ++k) {            // and the column
            lcQB(&ket, uc, u + k * uc, c + k * uc);             // Apply all remaining LCU channels
        }
        printf("test[%d]\n", j);

        for (depth_t k = 0; k < obsc; ++k) {                        // Calculate the j-th diagonal element for all
            stateInitEmpty(&bra, state->qubits);                    // observables
            stateInitVector(&bra, ket.vec);
            obs[k](&bra);

            printf("\tobs[%d] = ", k);
            vectorPrint(bra.vec, bra.dim);

            momMat[k][offset] = stateOverlap(bra, ket);
            stateFreeVector(&bra);
        }

        for (depth_t i = 0; i < uc - j; ++i) {
            for (depth_t k = 0; k < obsc; ++k) {
                stateInitEmpty(&bra, state->qubits);
                stateInitVector(&bra, state->vec);
                u[i + j + link * uc](&bra);
                for (depth_t l = link + 1; l < circdepth; ++l) {
                    lcQB(&bra, uc, u + l * uc, c + l * uc);
                }
                obs[k](&bra);
                momMat[k][i + offset] = stateOverlap(bra, ket);
                stateFreeVector(&bra);
            }
        }
        stateFreeVector(&ket);
    }
}