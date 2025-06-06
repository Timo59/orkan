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
 *                                              Observable mean value
 * =====================================================================================================================
 */

double meanDiagObs(const state_t* state, const double obs[]) {
    double out = 0;
#pragma omp parallel for reduction(+:out) default(none) shared(state, obs)
    for (dim_t i = 0; i < state->dim; ++i) {
        out += pow(cabs(state->vec[i]), 2) * obs[i];
    }
    return out;
}

double meanCG(const state_t* state, const applyCG obs) {
    state_t tmp;                                                    // Temporary state the observable is applied to
    stateInitEmpty(&tmp, state->qubits);
    stateInitVector(&tmp, state->vec);
    obs(&tmp);
    const cplx_t overlap = stateOverlap(tmp, *state);
    if (fabs(cimag(overlap)) < 1e-9) {
        stateFreeVector(&tmp);
        return creal(overlap);
    }
    fprintf(stderr, "meanCG(): Mean value of observable has non-negligible imaginary component\n");
    stateFreeVector(&tmp);
    exit(EXIT_FAILURE);
}

double meanObsHerm(const state_t* state, const herm_t* obs) {
    state_t tmp;                                                    // Temporary state the observable is applied to
    stateInitEmpty(&tmp, state->qubits);
    stateInitVector(&tmp, state->vec);
    applyHerm(&tmp, obs);
    const cplx_t overlap = stateOverlap(tmp, *state);
    if (fabs(cimag(overlap)) < 1e-9) {
        stateFreeVector(&tmp);
        return creal(overlap);
    }
    fprintf(stderr, "meanObsHerm(): Mean value of observable has non-negligible imaginary component\n");
    stateFreeVector(&tmp);
    exit(EXIT_FAILURE);
}

/*
 * =====================================================================================================================
 *                                              Gradient PQC
 * =====================================================================================================================
 */

double* gradPQCDiag(state_t* state, const double obs[], const depth_t d, const applyPCG pqc[], const double par[],
    const applyCG gen[]) {
    double* out = malloc(d * sizeof(double));
    if (out == NULL) {
        fprintf(stderr, "gradPQCDiag(): out allocation failed\n");
        exit(EXIT_FAILURE);
    }
    state_t bra;                                    // Initialized to the input state, evolved with all evolution
    stateInitEmpty(&bra, state->qubits);            // operators and acted on with the observable
    stateInitVector(&bra, state->vec);

    for (depth_t i = 0; i < d; ++i){                                // Apply the PQC to bra
        pqc[i](&bra, par[i]);
    }
    applyDiag(&bra, obs);                                           // Apply the diagonal observable to bra

#pragma omp parallel for default(none) shared(state, d, pqc, par, gen, bra, out) num_threads(d)
    for (depth_t i = 0; i < d; ++i) {                               // Iterate the components of the PQC
        state_t ket;                                                // Initialized to the input state; acted on by the
        stateInitEmpty(&ket, state->qubits);                        // derivative of the PQC wrt to the i-th parameter
        stateInitVector(&ket, state->vec);

        for (depth_t j = 0; j < i; ++j) {                           // Apply parametrized blocks up to i-1
            pqc[j](&ket, par[j]);
        }
        gen[i](&ket);                                               // and apply the evolution operator to ket
        for (depth_t j = i; j < d; ++j) {                           // Evolve ket with rest of the evolution operators
            pqc[j](&ket, par[j]);
        }

        out[i] = 2 * cimag(stateOverlap(bra, ket));
        stateFreeVector(&ket);
    }
    stateFreeVector(&bra);
    return out;
}

double* gradPQCHerm(state_t* state, const herm_t* obs, const depth_t d, const applyPCG pqc[], const double par[],
    const applyCG gen[]) {
    double* out = malloc(d * sizeof(double));
    if (out == NULL) {
        fprintf(stderr, "gradPQCHerm(): out allocation failed\n");
        exit(EXIT_FAILURE);
    }

    state_t bra;                                    // Initialized to the input state, evolved with all evolution
    stateInitEmpty(&bra, state->qubits);            // operators and acted on with the observable
    stateInitVector(&bra, state->vec);

    for (depth_t i = 0; i < d; ++i){                                // Apply the PQC to bra
        pqc[i](&bra, par[i]);
    }
    applyHerm(&bra, obs);                                           // Apply the diagonal observable to bra

#pragma omp parallel for default(none) shared(state, d, pqc, par, gen, bra, out) num_threads(d)
    for (depth_t i = 0; i < d; ++i) {                               // Iterate the components of the PQC
        state_t ket;                                                // Initialized to the input state; acted on by the
        stateInitEmpty(&ket, state->qubits);                        // derivative of the PQC wrt to the i-th parameter
        stateInitVector(&ket, state->vec);

        for (depth_t j = 0; j < i; ++j) {                           // Apply parametrized blocks up to i-1
            pqc[j](&ket, par[j]);
        }
        gen[i](&ket);                                               // and apply the evolution operator to ket
        for (depth_t j = i; j < d; ++j) {                           // Evolve ket with rest of the evolution operators
            pqc[j](&ket, par[j]);
        }

        out[i] = 2 * cimag(stateOverlap(bra, ket));
        stateFreeVector(&ket);
    }
    stateFreeVector(&bra);
    return out;
}

/*
 * This function calculates the moment matrices of a set of observables wrt to a sequence of LCU.
 *
 *  @param[in, out] state       : Input quantum state to the PQC. On exit, the output of the PQC up to link with the
 *                                CURRENT parameter setting
 *  @param[in]      obsc        : Number of moment matrices
 *  @param[in, out] obs[]       : Set of observables constituting the moment matrices
 *  @param[in]      circdepth   : Number of LCU channels
 *  @param[in, out] lccg[]      : Array of LCU channels
 *  @param[in, out] link        : Element of the LCU sequence constituting the basis for the moment matrices
 *  @param[in, out] momMat[]    : Array of moment matrices in column major form
 *
 * Output:
 *  An array of moment matrices, whose entries are the overlaps of two quantum states with the corresponding observable
 *  applied. In particular, all LCU preceding the one specified by link are applied with their coefficients setting on
 *  input. Then, the set of unitary operators at link are applied one by one to the quantum state to give the states for
 *  each column and row, respectively. The remaining LCU are also applied with their coefficients determined on input.
 */

void mmseq(state_t* state,
           const depth_t obsc,
           const applyCG obs[],
           const depth_t circdepth,
           const lccg_t lcu[],
           const depth_t link,
           cplx_t* momMat[])
{
    // Apply all LCU channels prior to 'link'
    for (depth_t k = 0; k < link; ++k) {
        applyLCCG(state, lcu + k);
    }

#pragma omp parallel for default(none) shared(state, obsc, obs, circdepth, lcu, link, momMat)
    // Iterate the moment matrices' columns
    for (depth_t j = 0; j < lcu[link].len; ++j) {
        state_t bra, ket;
        stateInitEmpty(&ket, state->qubits);
        stateInitVector(&ket, state->vec);

        // Apply the search unitary according to link and column
        lcu[link].comp[j](&ket);

        // Apply the remaining LCU channels
        for (depth_t k = link + 1; k < circdepth; ++k) {
            applyLCCG(state, lcu + k);
        }

        // Calculate the j-th diagonal element for all observables
        for (depth_t k = 0; k < obsc; ++k) {
            stateInitEmpty(&bra, state->qubits);
            stateInitVector(&bra, ket.vec);
            obs[k](&bra);

            momMat[k][j + j * lcu[link].len] = stateOverlap(bra, ket);
            stateFreeVector(&bra);
        }

        // Iterate the moment matrices' rows
        for (depth_t i = j + 1; i < lcu[link].len; ++i) {
            // Iterate the moment matrices themselves
            for (depth_t k = 0; k < obsc; ++k) {
                stateInitEmpty(&bra, state->qubits);
                stateInitVector(&bra, state->vec);
                // Apply the search unitary according to link and row
                lcu[link].comp[i](&bra);
                // Apply the remaining LCU channels
                for (depth_t l = link + 1; l < circdepth; ++l) {
                    applyLCCG(state, lcu + l);
                }
                obs[k](&bra);
                momMat[k][i + j * lcu[link].len] = stateOverlap(bra, ket);
                momMat[k][j + i * lcu[link].len] = conj(momMat[k][i + j * lcu[link].len]);
                stateFreeVector(&bra);
            }
        }
        stateFreeVector(&ket);
    }
}