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
    state_t tmp;                                                // Temporary state, initialized with the input state
                                                                // vector, the operator's terms are applied to
    register cplx_t result = 0;                                 // Output variable

    stateInitEmpty(&tmp, state->qubits);                        // Initialize the temporary state
    for (complength_t i = 0; i < op->length; ++i) {             // Iterate the operator's terms
        stateCopyVector(&tmp, state->vec);                      // Copy the input state's vector to the temporary state
        applyPauliStr(&tmp, op->comps + (i * op->qubits));      // Apply the current Pauli string to the temporary state
        cplx_t term = cInner(state->vec, tmp.vec, state->dim);  // Store inner product of state and ket to tmp
        result += (op->coeffs[i] * term);                       // Add tmp times the current coefficient to result
    }
    stateFreeVector(&tmp);

    return result;
}

/*
 * =====================================================================================================================
 *                                                  Mean value observable
 * =====================================================================================================================
 */
double expValObsPauli(const state_t* state, const pauliObs_t* obs) {
    state_t tmp;                                                // Temporary state, initialized with the input state
                                                                // vector, the observable's terms are applied to
    double result = 0;                                          // Output variable

    stateInitEmpty(&tmp, state->qubits);                        // Initialize the temporary state

    for (complength_t i = 0; i < obs->length; ++i) {            // Iterate the observable's terms
        stateCopyVector(&tmp, state->vec);                      // Copy the input state's vector to the temporary state
        applyPauliStr(&tmp, obs->comps + (i * obs->qubits));    // Apply the current Pauli string to the temporary state
        cplx_t term = cInner(state->vec, tmp.vec, state->dim);  // Store inner product of state and ket to tmp
        result += (obs->coeffs[i] * creal(term));               // Add tmp times the current coefficient to result
    }
    stateFreeVector(&tmp);

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

    state_t tmp;                            // state, initialized to the input state, holding each intermediate
    stateInitEmpty(&tmp, state->qubits);    // evolution
    stateCopyVector(&tmp, state->vec);

    state_t ket;                            // state, inheriting each intermediate evolution from tmp, acted on with the
    stateInitEmpty(&ket, state->qubits);    // respective evolution operator and finally evolved with the remaining
    // evolution operators
    double* result = malloc(circdepth * sizeof(double));

    for (depth_t i = 0; i < circdepth; ++i) {           // Evolve bra with all evolutional operators
        evolveObsTrotter(&bra, evoOps[i], par[i]);
    }
    applyObs(&bra, obs);                                // And apply the observable to it

    for (depth_t k = 0; k < circdepth; ++k) {               // Iterate the components of the PQC
        evolveObsTrotter(&tmp, evoOps[k], par[k]);          // Evolve tmp with the current component,
        stateCopyVector(&ket, tmp.vec);                     // copy its vector to ket
        applyObs(&ket, evoOps[k]);                          // and apply the observable to ket

        for (depth_t l = k + 1; l < circdepth; ++l) {       // Evolve ket with the rest of the evolution operators
            evolveObsTrotter(&ket, evoOps[l], par[l]);
        }

        result[k] = 2 * cimag(stateOverlap(&bra, &ket));    // The k-th entry is the overlap of bra and ket
    }
    stateFreeVector(&bra);
    stateFreeVector(&ket);
    stateFreeVector(&tmp);
    return result;
}
