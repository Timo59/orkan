//
// Created by Timo Ziegler on 10.10.24.
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
 *                                              Apply Pauli structs
 * =====================================================================================================================
 */
void applyPauliStr(state_t* state, const pauli_t paulistr[]) {
    for (qubit_t qubit = 0; qubit < state->qubits; ++qubit) {
        switch (paulistr[qubit]) {
            case ID: {
                break;
            }
            case Z: {
                applyZ(state, qubit);
                break;
            }
            case X: {
                applyX(state, qubit);
                break;
            }
            case Y: {
                applyY(state, qubit);
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
 * This function applies a Pauli operator to a quantum state.
 *
 * Input:
 *      state_t* state:         State of the qubit system
 *      pauliOp_t* op:          Pauli operator, i.e., a complex linear combination of Pauli strings
 *
 * Output:
 *      This function has no return value, but sums the outcome of single Pauli string applications weighted by the
 *      coefficients in place.
 */

void applyOpPauli(state_t* state, const pauliOp_t* op) {
    __LAPACK_int N = (__LAPACK_int) state->dim;
    state_t tmp;                                                        // Temporary state to apply each Pauli string to
//    stateInitEmpty(&tmp, state->qubits);
    tmp.qubits = state->qubits;
    tmp.dim = state->dim;
    cplx_t tmpVec[tmp.dim];
    tmp.vec = tmpVec;

    cplx_t* tmpSum;                                                     // Temporary vector holding the intermediate sum
    if ((tmpSum = calloc(state->dim, sizeof(cplx_t))) == NULL) {                // of state vectors
        fprintf(stderr, "tmpSum allocation in applyObsPauli failed\n");
        return;
    }
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

/*
 * This function applies a Pauli observable to a quantum state.
 *
 * Input:
 *      state_t* state:         State of the qubit system
 *      pauliObs_t* obs:        Pauli observable, i.e., a real linear combination of Pauli strings
 *
Output:
 *      This function has no return value, but sums the outcome of single Pauli string applications weighted by the
 *      coefficients in place.
 */
void applyObsPauli(state_t* state, const pauliObs_t* obs) {
    __LAPACK_int N = (__LAPACK_int) state->dim;
    state_t tmp;                                                        // Temporary state to apply each Pauli string to
    stateInitEmpty(&tmp, state->qubits);

    cplx_t* tmpSum;                                                     // Temporary vector holding the intermediate sum
    if ((tmpSum = calloc(state->dim, sizeof(cplx_t))) == NULL) {                // of state vectors
        fprintf(stderr, "tmpSum allocation in applyObsPauli failed\n");
        return;
    }
    for (complength_t i = 0; i < obs->length; ++i) {                    // Iterate the observable's terms
        cblas_zcopy(N, state->vec, 1, tmp.vec, 1);                      // Copy input's state vector to tmp
        applyPauliStr(&tmp, obs->comps + (i * obs->qubits));            // Apply the Pauli string to it
        __LAPACK_double_complex alpha = obs->coeffs[i];                 // Multiply the terms' coefficient
        cblas_zaxpy(N, &alpha, tmp.vec, 1, tmpSum, 1);                  // and add the result to the
    }                                                                    // intermediate sum
    stateFreeVector(&tmp);
    stateFreeVector(state);                                             // Free the input's state vector
    state->vec = tmpSum;                                                // and replace it by tmp
}

/*
 * =====================================================================================================================
 *                                              Evolve with Pauli string
 * =====================================================================================================================
 */

/*
 * This function evolves a state with the with a Pauli string.
 *
 * Input:
 *      state_t* state:         State of the qubit system
 *      pauli_t paulistr[]:     Pauli string, i.e., array of Pauli operators
 * Output:
 *      This function has no return value, but changes the state vector in place.
 */

void evolvePauliStr(state_t* state, const pauli_t paulistr[], double angle) {
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