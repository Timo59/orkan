//
// Created by Timo Ziegler on 11.10.24.
//
/*
 * =====================================================================================================================
 *                                                      includes
 * =====================================================================================================================
 */

#ifndef EXACT_H
#include "exact.h"
#endif

/*
 * =====================================================================================================================
 *                                                  Mean value operator
 * =====================================================================================================================
 */
cplx_t expValPauli(const state_t* state, const pauliOp_t* op) {
    cplx_t result = 0;                          // Output variable

#pragma omp parallel default(none) shared(op, result, state)
    {
        state_t tmp;                            // Temporary state, initialized with the input state
        stateInitEmpty(&tmp, state->qubits);    // vector, the obperator's terms are applied to

#pragma omp for
        for (complength_t i = 0; i < op->length; ++i) {                 // Iterate the operator's terms
            stateCopyVector(&tmp, state->vec);                          // Copy the input state's vector to tmp
            applyPauliStr(&tmp, op->comps + (i * op->qubits));          // Apply the current Pauli string to tmp
            cplx_t term = cInner(state->vec, tmp.vec, state->dim);      // Store inner product of state and ket to tmp

#pragma omp critical(sum)
            {
                result += (op->coeffs[i] * creal(term));                // Add tmp times the current coefficient to
            }                                                           // result
        }
        stateFreeVector(&tmp);
    }
    return result;
}

/*
 * =====================================================================================================================
 *                                                  Mean value observable
 * =====================================================================================================================
 */
double expValObsPauli(const state_t* state, const pauliObs_t* obs) {
    double result = 0;                          // Output variable

#pragma omp parallel default(none) shared(obs, result, state)
    {
        state_t tmp;                            // Temporary state, initialized with the input state
        stateInitEmpty(&tmp, state->qubits);    // vector, the observable's terms are applied to

#pragma omp for
        for (complength_t i = 0; i < obs->length; ++i) {                // Iterate the observable's terms
            stateCopyVector(&tmp, state->vec);                          // Copy the input state's vector to tmp
            applyPauliStr(&tmp, obs->comps + (i * obs->qubits));        // Apply the current Pauli string to tmp
            cplx_t term = cInner(state->vec, tmp.vec, state->dim);      // Store inner product of state and ket to tmp

#pragma omp critical(sum)
            {
                result += (obs->coeffs[i] * creal(term));               // Add tmp times the current coefficient to
            }                                                           // result
        }
        stateFreeVector(&tmp);
    }
    return result;
}

/*
 * =====================================================================================================================
 *                                                      gradient PQC
 * =====================================================================================================================
 */
double* gradPQC(const state_t* state,
                const double par[],
                const obs_t* obs,
                const obs_t *evoOps[],
                depth_t circdepth)
{
    state_t bra;                            // state, initialized to the input state, evolved with all evolution
    stateInitEmpty(&bra, state->qubits);    // operators and acted on with the observable
    stateCopyVector(&bra, state->vec);

    double* result = malloc(circdepth * sizeof(double));

    for (depth_t i = 0; i < circdepth; ++i) {           // Evolve bra with all evolutional operators
        evolveObsTrotter(&bra, evoOps[i], par[i]);
    }
    applyObs(&bra, obs);                                // And apply the observable to it

#pragma omp parallel for default(none) shared(bra, circdepth, evoOps, par, result, state) num_threads(circdepth)
    for (depth_t i = 0; i < circdepth; ++i) {   // Iterate the components of the PQC
        state_t ket;                            // state, inheriting vector from input state, evolved with operators up
        stateInitEmpty(&ket, state->qubits);    // to and including the current, acted on with the current evolution
        stateCopyVector(&ket, state->vec);      // operator and finally evolved with the remaining operators


        for (depth_t j = 0; j <= i; ++j) {
            evolveObsTrotter(&ket, evoOps[j], par[j]);  // Evolve ket with operators up to and including the current
        }                                               // evolution operator
        applyObs(&ket, evoOps[i]);                      // and apply the current evolution operator to it

        for (depth_t j = i + 1; j < circdepth; ++j) {       // Evolve with the remaining operators
            evolveObsTrotter(&ket, evoOps[j], par[j]);
        }

        result[i] = 2 * cimag(stateOverlap(&bra, &ket));    // The k-th entry is the overlap of bra and ket
        stateFreeVector(&ket);
    }
    stateFreeVector(&bra);
    return result;
}