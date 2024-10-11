//
// Created by Timo Ziegler on 10.10.24.
//
/*
 * =====================================================================================================================
 *                                                      includes
 * =====================================================================================================================
 */

#ifndef PAULIOLD_H
#include "pauliOld.h"
#endif

/*
 * =====================================================================================================================
 *                                              Apply Pauli structs
 * =====================================================================================================================
 */

void applyPauliStr_old(state_t* state, const pauli_t paulistr[]) {
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
void applyOpPauli_old(state_t* state, const pauliOp_t* op) {
    state_t tmp;                                                        // Temporary state to apply each Pauli string to
    stateInitEmpty(&tmp, state->qubits);
    cplx_t* tmpSum = calloc(state->dim, sizeof(cplx_t));     // Temporary vector holding the intermediate sum
                                                                        // of state vectors
    for (complength_t i = 0; i < op->length; ++i) {                     // Iterate the operator's terms
        stateCopyVector(&tmp, state->vec);                              // Copy input's state vector to tmp
        applyPauliStr_old(&tmp, op->comps + (i * op->qubits));              // Apply the Pauli string to it
        cscalarVecMulInPlace(op->coeffs[i], tmp.vec, state->dim);       // Multiply the terms' coefficient
        cvecAddInPlace(tmpSum, tmp.vec, state->dim);                            // and add the result to the
                                                                                // intermediate sum
    }
    stateFreeVector(&tmp);
    stateFreeVector(state);
    state->vec = tmpSum;
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
void applyObsPauli_old(state_t* state, const pauliObs_t* obs) {
    state_t tmp;                                                        // Temporary state to apply each Pauli string to
    stateInitEmpty(&tmp, state->qubits);
    cplx_t* tmpSum = calloc(state->dim, sizeof(cplx_t));     // Temporary vector holding the intermediate sum
                                                                        // of state vectors
    printf("Hello from Pauli OLD\n");
    for (complength_t i = 0; i < obs->length; ++i) {                    // Iterate the observable's terms
        stateCopyVector(&tmp, state->vec);                              // Copy input's state vector to tmp
        applyPauliStr_old(&tmp, obs->comps + (i * obs->qubits));            // Apply the Pauli string to it
        cscalarVecMulInPlace((cplx_t) obs->coeffs[i], tmp.vec, state->dim);     // Multiply the terms' coefficient
        cvecAddInPlace(tmpSum, tmp.vec, state->dim);                            // and add the result to the
                                                                                // intermediate sum
    }
    stateFreeVector(&tmp);
    stateFreeVector(state);
    state->vec = tmpSum;
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
void evolvePauliStr_old(state_t* state, const pauli_t paulistr[], double angle) {
    cplx_t* tmp = malloc(sizeof(cplx_t) * state->dim);          // Temporary vector holding state vector multiplied
                                                                // with cos(angle)
    for (dim_t entry = 0; entry < state->dim; ++entry) {        // Iterate entries of the input's state vector
        tmp[entry] = cos(angle) * state->vec[entry];            // Store it times cos(angle) in tmp
        state->vec[entry] *= -I * sin(angle);                   // Multiply it with -i * sin(angle)
    }
    applyPauliStr_old(state, paulistr);                             // Apply the Pauli string to the input state
    for (dim_t entry = 0; entry < state->dim; ++entry) {        // Iterate the entries of the input's state vector
        state->vec[entry] += tmp[entry];                        // Add the tmp's entry to it
    }
    free(tmp);
}