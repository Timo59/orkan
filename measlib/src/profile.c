//
// Created by Timo Ziegler on 17.10.24.
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

cplx_t expValPauli_blas(const state_t* state, const pauliOp_t* op) {
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

cplx_t expValPauli_blas_omp(const state_t* state, const pauliOp_t* op) {
    cplx_t result = 0;                              // Output variable
    __LAPACK_int N = (__LAPACK_int) state->dim;
    __LAPACK_int one = (__LAPACK_int) 1;

#pragma omp parallel default(none) shared(N, op, one, result, state)
    {
        state_t tmp;                            // Temporary state, initialized with the input state
        stateInitEmpty(&tmp, state->qubits);    // vector, the operator's terms are applied to
        cplx_t term;

#pragma omp for
        for (complength_t i = 0; i < op->length; ++i) {             // Iterate the operator's terms
            cblas_zcopy(N, state->vec, 1, tmp.vec, 1);              // Copy the input state's vector to tmp
            applyPauliStr(&tmp, op->comps + (i * op->qubits));      // Apply the current Pauli string to tmp
            cblas_zdotc_sub(N, state->vec, 1, tmp.vec, 1, &term);   // Store inner product of state and ket to tmp

#pragma omp critical(sum)
            {
                result += (op->coeffs[i] * term);                   // Add tmp times the current coefficient to
            }                                                       // result
        }
        stateFreeVector(&tmp);
    }
    return result;
}

cplx_t expValPauli_omp(const state_t* state, const pauliOp_t* op) {
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

double expValObsPauli_blas(const state_t* state, const pauliObs_t* obs) {
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

double expValObsPauli_blas_omp(const state_t* state, const pauliObs_t* obs) {
    double result = 0;                              // Output variable
    __LAPACK_int N = (__LAPACK_int) state->dim;
    __LAPACK_int one = (__LAPACK_int) 1;

#pragma omp parallel default(none) shared(N, obs, one, result, state)
    {
        state_t tmp;                            // Temporary state, initialized with the input state
        stateInitEmpty(&tmp, state->qubits);    // vector, the observable's terms are applied to
        cplx_t term;

#pragma omp for
        for (complength_t i = 0; i < obs->length; ++i) {            // Iterate the observable's terms
            cblas_zcopy(N, state->vec, 1, tmp.vec, 1);              // Copy the input state's vector to tmp
            applyPauliStr(&tmp, obs->comps + (i * obs->qubits));    // Apply the current Pauli string to tmp
            cblas_zdotc_sub(N, state->vec, 1, tmp.vec, 1, &term);   // Store inner product of state and ket to tmp

#pragma omp critical(sum)
            {
                result += (obs->coeffs[i] * creal(term));           // Add tmp times the current coefficient to
            }                                                       // result
        }
        stateFreeVector(&tmp);
    }
    return result;
}
double expValObsPauli_omp(const state_t* state, const pauliObs_t* obs) {
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

double* gradPQC_blas(const state_t* state,
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

double* gradPQC_blas_omp(const state_t* state,
                         const double par[],
                         const obs_t* obs,
                         const obs_t *evoOps[],
                         depth_t circdepth)
{
    __LAPACK_int N = (__LAPACK_int) state->dim;
    state_t bra;                                // state, initialized to the input state, evolved with all evolution
    stateInitEmpty(&bra, state->qubits);        // operators and acted on with the observable
    cblas_zcopy(N, state->vec, 1, bra.vec, 1);

    double* result = malloc(circdepth * sizeof(double));

    for (depth_t i = 0; i < circdepth; ++i) {           // Evolve bra with all evolutional operators
        evolveObsTrotter(&bra, evoOps[i], par[i]);
    }
    applyObs(&bra, obs);                                // And apply the observable to it

#pragma omp parallel for default(none) shared(bra, circdepth, evoOps, N, par, result, state) num_threads(circdepth)
    for (depth_t i = 0; i < circdepth; ++i) {   // Iterate the components of the PQC
        state_t ket;                                    // state, inheriting vector from input state, evolved with
        stateInitEmpty(&ket, state->qubits);            // operators up to and including the current, acted on with the
        cblas_zcopy(N, state->vec, 1, ket.vec, 1);      // current evolution operator and finally evolved with the
                                                        // remaining operators
        for (depth_t j = 0; j <= i; ++j) {
            evolveObsTrotter(&ket, evoOps[j], par[j]);  // Evolve ket with operators up to and including the current
        }                                               // evolution operator
        applyObs(&ket, evoOps[i]);                      // and apply the current evolution operator to it

        for (depth_t j = i + 1; j < circdepth; ++j) {   // Evolve with the remaining operators
            evolveObsTrotter(&ket, evoOps[j], par[j]);
        }

        cplx_t dotc;                                        // The k-th entry is the overlap of bra and ket
        cblas_zdotc_sub(N, bra.vec, 1, ket.vec, 1, &dotc);
        result[i] = 2 * cimag(dotc);
        stateFreeVector(&ket);
    }
    stateFreeVector(&bra);
    return result;
}

double* gradPQC_omp(const state_t* state,
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
