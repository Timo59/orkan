//
// Created by Timo Ziegler on 17.10.24.
//
/*
 * =====================================================================================================================
 *                                                      includes
 * =====================================================================================================================
 */
#ifndef PAULI_H
#include "pauli.h"
#endif

/*
 * =====================================================================================================================
 *                                                  Apply Pauli string
 * =====================================================================================================================
 */
void applyPauliStr(state_t* state, const pauli_t paulistr[]) {
    for (qubit_t qubit = 0; qubit < state->qubits; ++qubit) {
        switch (paulistr[qubit]) {
            case ID: {
                break;
            }
            case Z: {
#if defined(BLAS_GATES)
                applyZ_blas(state, qubit);
#elif defined(BLAS_OMP_GATES)
                applyZ_blas_omp(state, qubit);
#elif defined(OMP_GATES)
                applyZ_omp(state, qubit);
#else
                applyZ(state, qubit);
#endif
                break;
            }
            case X: {
#if defined(BLAS_GATES)
                applyX_blas(state, qubit);
#elif defined(BLAS_OMP_GATES)
                applyX_blas_omp(state, qubit);
#elif defined(OMP_GATES)
                applyX_omp(state, qubit);
#else
                applyX(state, qubit);
#endif
                break;
            }
            case Y: {
#if defined(BLAS_GATES)
                applyY_blas(state, qubit);
#elif defined(BLAS_OMP_GATES)
                applyY_blas_omp(state, qubit);
#elif defined(OMP_GATES)
                applyY_omp(state, qubit);
#else
                applyY(state, qubit);
#endif
                break;
            }
            default: {
                printf("Error in applyPauliStr: Invalid value for Pauli.\n");
                exit(1);
            }
        }
    }
}

/*
 * =====================================================================================================================
 *                                                  Apply Pauli operator
 * =====================================================================================================================
 */
void applyOpPauli(state_t* state, const pauliOp_t* op) {
    state_t tmp;                                                        // Temporary state to apply each Pauli string to
    stateInitEmpty(&tmp, state->qubits);
    cplx_t* tmpSum = calloc(state->dim, sizeof(cplx_t));     // Temporary vector holding the intermediate sum
                                                                        // of state vectors
    for (complength_t i = 0; i < op->length; ++i) {                     // Iterate the operator's terms
        stateCopyVector(&tmp, state->vec);                              // Copy input's state vector to tmp
        applyPauliStr(&tmp, op->comps + (i * op->qubits));              // Apply the Pauli string to it
        cscalarVecMulInPlace(op->coeffs[i], tmp.vec, state->dim);       // Multiply the terms' coefficient
        cvecAddInPlace(tmpSum, tmp.vec, state->dim);                            // and add the result to the
                                                                                        // intermediate sum
    }
    stateFreeVector(&tmp);
    stateFreeVector(state);
    state->vec = tmpSum;
}

void applyOpPauli_blas(state_t* state, const pauliOp_t* op) {
    __LAPACK_int N = (__LAPACK_int) state->dim;
    state_t tmp;                                                        // Temporary state to apply each Pauli string to
    stateInitEmpty(&tmp, state->qubits);
    cplx_t* tmpSum = calloc(state->dim, sizeof(cplx_t));                // Temporary vector holding the intermediate sum
                                                                        // of state vectors
    for (complength_t i = 0; i < op->length; ++i) {                     // Iterate the observable's terms
        cblas_zcopy(N, state->vec, 1, tmp.vec, 1);                      // Copy input's state vector to tmp
        applyPauliStr(&tmp, op->comps + (i * op->qubits));              // Apply the Pauli string to it
        __LAPACK_double_complex alpha = op->coeffs[i];                  // Multiply the terms' coefficient
        cblas_zaxpy(N, &alpha, tmp.vec, 1, tmpSum, 1);                  // and add the result to the
                                                                        // intermediate sum
    }
    stateFreeVector(&tmp);
    stateFreeVector(state);                                             // Free the input's state vector
    state->vec = tmpSum;                                                // and replace it by tmp
}

void applyOpPauli_blas_omp(state_t* state, const pauliOp_t* op) {
    __LAPACK_int N = (__LAPACK_int) state->dim;
    cplx_t* tmpSum = calloc(state->dim, sizeof(cplx_t));                // Temporary vector holding the intermediate sum
# pragma omp parallel default(none) shared(N, op, state, tmpSum)
    {
        state_t tmp;
        stateInitEmpty(&tmp, state->qubits);

#pragma omp for
        for (complength_t i = 0; i < op->length; ++i) {                             // Iterate the observable's terms
            cblas_zcopy(N, state->vec, 1, tmp.vec, 1);                              // Copy input's state vector to tmp
            applyPauliStr(&tmp, op->comps + (i * op->qubits));                      // Apply the Pauli string to it
            __LAPACK_double_complex alpha = op->coeffs[i];

#pragma omp critical(sum)
            {
                cblas_zaxpy(N, &alpha, tmp.vec, 1, tmpSum, 1);                      // and add the result to the
            }                                                                       // intermediate sum
        }
        stateFreeVector(&tmp);
    }
    stateFreeVector(state);                                                         // Free the input's state vector
    state->vec = tmpSum;                                                            // and replace it by tmpSum
}

void applyOpPauli_omp(state_t* state, const pauliOp_t* op) {
    cplx_t *tmpSum = calloc(state->dim, sizeof(cplx_t));                // Temporary vector holding the intermediate sum

# pragma omp parallel default(none) shared(op, state, tmpSum)
    {
        state_t tmp;                                                    // Temporary state to apply each Pauli string to
        stateInitEmpty(&tmp, state->qubits);

#pragma omp for
        for (complength_t i = 0; i < op->length; ++i) {                             // Iterate the observable's terms
            stateCopyVector(&tmp, state->vec);                                      // Copy input's state vector to tmp
            applyPauliStr(&tmp, op->comps + (i * op->qubits));                      // Apply the Pauli string to it
            cscalarVecMulInPlace(op->coeffs[i], tmp.vec, state->dim);               // Multiply the terms' coefficient

#pragma omp critical(sum)
            {
                cvecAddInPlace(tmpSum, tmp.vec, state->dim);                        // and add the result to the
            }                                                                       // intermediate sum
        }
        stateFreeVector(&tmp);
    }
    stateFreeVector(state);                                                         // Free the input's state vector
    state->vec = tmpSum;                                                            // and replace it by tmp
}

/*
 * =====================================================================================================================
 *                                                  Apply Pauli observable
 * =====================================================================================================================
 */
void applyObsPauli(state_t* state, const pauliObs_t* obs) {
    state_t tmp;                                                        // Temporary state to apply each Pauli string to
    stateInitEmpty(&tmp, state->qubits);
    cplx_t* tmpSum = calloc(state->dim, sizeof(cplx_t));                // Temporary vector holding the intermediate sum
                                                                        // of state vectors
    for (complength_t i = 0; i < obs->length; ++i) {                    // Iterate the observable's terms
        stateCopyVector(&tmp, state->vec);                              // Copy input's state vector to tmp
        applyPauliStr(&tmp, obs->comps + (i * obs->qubits));            // Apply the Pauli string to it
        cscalarVecMulInPlace((cplx_t) obs->coeffs[i], tmp.vec, state->dim);     // Multiply the terms' coefficient
        cvecAddInPlace(tmpSum, tmp.vec, state->dim);                            // and add the result to the
                                                                                        // intermediate sum
    }
    stateFreeVector(&tmp);
    stateFreeVector(state);
    state->vec = tmpSum;
}

void applyObsPauli_blas(state_t* state, const pauliObs_t* obs) {
    __LAPACK_int N = (__LAPACK_int) state->dim;
    state_t tmp;                                                        // Temporary state to apply each Pauli string to
    stateInitEmpty(&tmp, state->qubits);
    cplx_t* tmpSum = calloc(state->dim, sizeof(cplx_t));                // Temporary vector holding the intermediate sum
                                                                        // of state vectors
    for (complength_t i = 0; i < obs->length; ++i) {                    // Iterate the observable's terms
        cblas_zcopy(N, state->vec, 1, tmp.vec, 1);                      // Copy input's state vector to tmp
        applyPauliStr(&tmp, obs->comps + (i * obs->qubits));            // Apply the Pauli string to it
        __LAPACK_double_complex alpha = obs->coeffs[i];                 // Multiply the terms' coefficient
        cblas_zaxpy(N, &alpha, tmp.vec, 1, tmpSum, 1);                  // and add the result to the
                                                                        // intermediate sum
    }
    stateFreeVector(&tmp);
    stateFreeVector(state);                                             // Free the input's state vector
    state->vec = tmpSum;                                                // and replace it by tmp
}

void applyObsPauli_blas_omp(state_t* state, const pauliObs_t* obs) {
    __LAPACK_int N = (__LAPACK_int) state->dim;
    cplx_t* tmpSum = calloc(state->dim, sizeof(cplx_t));                // Temporary vector holding the intermediate sum

# pragma omp parallel default(none) shared(N, obs, state, tmpSum)
    {
        state_t tmp;
        stateInitEmpty(&tmp, state->qubits);

#pragma omp for
        for (complength_t i = 0; i < obs->length; ++i) {                            // Iterate the observable's terms
            cblas_zcopy(N, state->vec, 1, tmp.vec, 1);                              // Copy input's state vector to tmp
            applyPauliStr(&tmp, obs->comps + (i * obs->qubits));                    // Apply the Pauli string to it
            __LAPACK_double_complex alpha = obs->coeffs[i];

#pragma omp critical(sum)
            {
                cblas_zaxpy(N, &alpha, tmp.vec, 1, tmpSum, 1);                      // and add the result to the
            }                                                                       // intermediate sum
        }
        stateFreeVector(&tmp);
    }
    stateFreeVector(state);                                                         // Free the input's state vector
    state->vec = tmpSum;                                                            // and replace it by tmpSum
}

void applyObsPauli_omp(state_t* state, const pauliObs_t* obs) {
    cplx_t* tmpSum = calloc(state->dim, sizeof(cplx_t));                // Temporary vector holding the intermediate sum

# pragma omp parallel default(none) shared(obs, state, tmpSum)
    {
        state_t tmp;                                                    // Temporary state to apply each Pauli string to
        stateInitEmpty(&tmp, state->qubits);

#pragma omp for
        for (complength_t i = 0; i < obs->length; ++i) {                            // Iterate the observable's terms
            stateCopyVector(&tmp, state->vec);                                      // Copy input's state vector to tmp
            applyPauliStr(&tmp, obs->comps + (i * obs->qubits));                    // Apply the Pauli string to it
            cscalarVecMulInPlace((cplx_t) obs->coeffs[i], tmp.vec, state->dim);     // Multiply the terms' coefficient

#pragma omp critical(sum)
            {
                cvecAddInPlace(tmpSum, tmp.vec, state->dim);                        // and add the result to the
            }                                                                       // intermediate sum
        }
        stateFreeVector(&tmp);
    }
    stateFreeVector(state);                                                         // Free the input's state vector
    state->vec = tmpSum;                                                            // and replace it by tmp
}

/*
 * =====================================================================================================================
 *                                              Evolve with Pauli string
 * =====================================================================================================================
 */
void evolvePauliStr(state_t* state, const pauli_t paulistr[], double angle) {
    cplx_t* tmp = malloc(sizeof(cplx_t) * state->dim);          // Temporary vector holding state vector multiplied
                                                                    // with cos(angle)
    for (dim_t entry = 0; entry < state->dim; ++entry) {        // Iterate entries of the input's state vector
        tmp[entry] = cos(angle) * state->vec[entry];            // Store it times cos(angle) in tmp
        state->vec[entry] *= -I * sin(angle);                   // Multiply it with -i * sin(angle)
    }
    applyPauliStr(state, paulistr);                             // Apply the Pauli string to the input state
    for (dim_t entry = 0; entry < state->dim; ++entry) {        // Iterate the entries of the input's state vector
        state->vec[entry] += tmp[entry];                        // Add the tmp's entry to it
    }
    free(tmp);
}

void evolvePauliStr_blas(state_t* state, const pauli_t paulistr[], double angle) {
    __LAPACK_int N = (__LAPACK_int) state->dim;
    __LAPACK_double_complex alpha = cos(angle);
    cplx_t* tmp = calloc(sizeof(cplx_t), state->dim);               // Temporary vector holding state vector multiplied
                                                                    // by cos(angle)
    cblas_zaxpy(N, &alpha, state->vec, 1, tmp, 1);                  // tmp is input state's vector time cos(alpha)
    applyPauliStr(state, paulistr);                                 // Apply the Pauli string to the input state
    alpha = -I * sin(angle);
    cblas_zaxpy(N, &alpha, state->vec, 1, tmp, 1);                  // Multiply the state's vector by -I * sin(alpha)
                                                                    // and add it to tmp
    free(state->vec);                                               // Free the input's state vector
    state->vec = tmp;                                               // and replace it by tmp
}

void evolvePauliStr_blas_omp(state_t* state, const pauli_t paulistr[], double angle) {
    cplx_t* tmp = malloc(sizeof(cplx_t) * state->dim);              // Temporary vector holding state vector multiplied
    // with cos(angle)
#pragma omp parallel for default(none) shared(angle, state, tmp)
    for (dim_t entry = 0; entry < state->dim; ++entry) {            // Iterate entries of the input's state vector
        tmp[entry] = state->vec[entry];                             // Store it time cos(angle) in tmp
        tmp[entry] *= cos(angle);                                   // Store it times cos(angle) in tmp
        state->vec[entry] *= -I * sin(angle);                       // Multiply it with -i * sin(angle)
    }
    applyPauliStr(state, paulistr);                                 // Apply the Pauli string to the input state

#pragma omp parallel for default(none) shared(state, tmp)
    for (dim_t entry = 0; entry < state->dim; ++entry) {            // Iterate the entries of the input's state vector
        state->vec[entry] += tmp[entry];                            // Add the tmp's entry to it
    }

    free(tmp);
}

void evolvePauliStr_omp(state_t* state, const pauli_t paulistr[], double angle) {
    cplx_t* tmp = malloc(sizeof(cplx_t) * state->dim);              // Temporary vector holding state vector multiplied
    // with cos(angle)
#pragma omp parallel for default(none) shared(angle, state, tmp)
    for (dim_t entry = 0; entry < state->dim; ++entry) {            // Iterate entries of the input's state vector
        tmp[entry] = state->vec[entry];                             // Store it time cos(angle) in tmp
        tmp[entry] *= cos(angle);                                   // Store it times cos(angle) in tmp
        state->vec[entry] *= -I * sin(angle);                       // Multiply it with -i * sin(angle)
    }
    applyPauliStr(state, paulistr);                                 // Apply the Pauli string to the input state

#pragma omp parallel for default(none) shared(state, tmp)
    for (dim_t entry = 0; entry < state->dim; ++entry) {            // Iterate the entries of the input's state vector
        state->vec[entry] += tmp[entry];                            // Add the tmp's entry to it
    }

    free(tmp);
}