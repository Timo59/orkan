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
    cplx_t term;                                                // Temporary var holding the mean value of each term of
                                                                // the operator
    register double result = 0;                                 // Output variable
    __LAPACK_int N = (__LAPACK_int) state->dim;

    stateInitEmpty(&tmp, state->qubits);                        // Initialize the temporary state
    for (complength_t i = 0; i < op->length; ++i) {             // Iterate the obperator's terms
        cblas_zcopy(N, state->vec, 1, tmp.vec, 1);              // Copy the input state's vector to the temporary state
        applyPauliStr(&tmp, op->comps + (i * op->qubits));      // Apply the current Pauli string to the temporary state
        cblas_zdotc_sub(N, state->vec, 1, tmp.vec, 1, &term);   // Store inner product of state and ket to tmp
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
    cplx_t term;                                                // Temporary var holding the mean value of each term of
                                                                // the observable
    register double result = 0;                                 // Output variable
    __LAPACK_int N = (__LAPACK_int) state->dim;
    stateInitEmpty(&tmp, state->qubits);                        // Initialize the temporary state
    for (complength_t i = 0; i < obs->length; ++i) {            // Iterate the observable's terms
        cblas_zcopy(N, state->vec, 1, tmp.vec, 1);              // Copy the input state's vector to the temporary state
        applyPauliStr(&tmp, obs->comps + (i * obs->qubits));    // Apply the current Pauli string to the temporary state
        cblas_zdotc_sub(N, state->vec, 1, tmp.vec, 1, &term);   // Store inner product of state and ket to tmp
        result += (obs->coeffs[i] * creal(term));                // Add tmp times the current coefficient to result
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
    __LAPACK_int N = (__LAPACK_int) state->dim;

    state_t bra;                                    // state, initialized to the input state, evolved with all evolution
    stateInitEmpty(&bra, state->qubits);            // operators and acted on with the observable
    cblas_zcopy(N, state->vec, 1, bra.vec, 1);

    state_t tmp;                                    // state, initialized to the input state, holding each intermediate
    stateInitEmpty(&tmp, state->qubits);            // evolution
    cblas_zcopy(N, state->vec, 1, tmp.vec, 1);

    state_t ket;                                    // state, inheriting each intermediate evolution from tmp, acted on
    stateInitEmpty(&ket, state->qubits);            // with the respective evolution operator and finally evolved with
                                                    // the remaining evolution operators
    double* result = malloc(circdepth * sizeof(double));

    for (depth_t i = 0; i < circdepth; ++i) {       // Evolve bra with all evolutional operators
        evolveObsTrotter(&bra, evoOps[i], par[i]);
    }
    applyObs(&bra, obs);                            // And apply the observable to it

    for (depth_t k = 0; k < circdepth; ++k) {           // Iterate the components of the PQC
        evolveObsTrotter(&tmp, evoOps[k], par[k]);      // Evolve tmp with the current component,
        cblas_zcopy(N, tmp.vec, 1, ket.vec, 1);         // copy its vector to ket
        applyObs(&ket, evoOps[k]);                      // and apply the observable to ket

        for (depth_t l = k + 1; l < circdepth; ++l) {   // Evolve ket with the rest of the evolution operators
            evolveObsTrotter(&ket, evoOps[l], par[l]);
        }

        cplx_t dotc;                                        // The k-th entry is the overlap of bra and ket
        cblas_zdotc_sub(N, bra.vec, 1, ket.vec, 1, &dotc);
        result[k] = 2 * cimag(dotc);
    }
    stateFreeVector(&bra);
    stateFreeVector(&ket);
    stateFreeVector(&tmp);
    return result;
}
