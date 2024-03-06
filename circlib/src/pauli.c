/*
 * =================================================================================================
 *                                              includes
 * =================================================================================================
 */

#include "pauli.h"

/*
 * =================================================================================================
 *                                              compound gate constructor
 * =================================================================================================
 */

// void initCompgate(compgate_t* compgate, complength_t length, qubit_t qubits) {
//     compgate->components = (pauli_t*) calloc(length * qubits, sizeof(pauli_t));
//     compgate->coefficients = (cplx_t*) calloc(length, sizeof(cplx_t));
//     compgate->length = length;
//     compgate->qubits = qubits;
// }

// void freeCompgate(compgate_t compgate) {
//     free(compgate.components);
//     free(compgate.coefficients);
// }

/* 
 * =================================================================================================
 *                                              compound gate fusion
 * =================================================================================================
 */

// void fuseCompgates(compgate_t* result, compgate_t compgate1, compgate_t compgate2) {
//     complength_t entry;
//     int8_t imagCounter;
//     for (complength_t i = 0; i < compgate1.length; ++i) {
//         for (complength_t j = 0; j < compgate2.length; ++j) {
//             entry = i * compgate2.length + j;
//             result->coefficients[entry] = compgate1.coefficients[i] * compgate2.coefficients[j];
//             imagCounter = 0;
//             for (qubit_t k = 0; k < compgate1.components->qubits; ++k) {
//                 switch (compgate1.components[i].paulis[k]) {
//                     case ID: {
//                         result->components[entry].paulis[k] = compgate2.components[j].paulis[k];
//                         break;
//                     }
//                     case Z: {
//                         switch (compgate2.components[j].paulis[k]) {
//                             case ID: {
//                                 result->components[entry].paulis[k] = Z;
//                                 break;
//                             }
//                             case Z: {
//                                 result->components[entry].paulis[k] = ID;
//                                 break;
//                             }
//                             case X: {
//                                 result->components[entry].paulis[k] = Y;
//                                 imagCounter += 1;
//                                 break;
//                             }
//                             case Y: {
//                                 result->components[entry].paulis[k] = X;
//                                 imagCounter -= 1;
//                                 break;
//                             }
//                             default: {
//                                 printf("Error in fuseCompgates: Invalid value for Pauli.\n");
//                                 exit(1);
//                             }  
//                         }
//                         break;
//                     }
//                     case X: {
//                         switch (compgate2.components[j].paulis[k]) {
//                             case ID: {
//                                 result->components[entry].paulis[k] = X;
//                                 break;
//                             }
//                             case Z: {
//                                 result->components[entry].paulis[k] = Y;
//                                 imagCounter -= 1;
//                                 break;
//                             }
//                             case X: {
//                                 result->components[entry].paulis[k] = ID;
//                                 break;
//                             }
//                             case Y: {
//                                 result->components[entry].paulis[k] = Z;
//                                 imagCounter += 1;
//                                 break;
//                             }
//                             default: {
//                                 printf("Error in fuseCompgates: Invalid value for Pauli.\n");
//                                 exit(1);
//                             } 
//                         }
//                         break;
//                     }
//                     case Y: {
//                         switch (compgate2.components[j].paulis[k]) {
//                             case ID: {
//                                 result->components[entry].paulis[k] = Y;
//                                 break;
//                             }
//                             case Z: {
//                                 result->components[entry].paulis[k] = X;
//                                 imagCounter += 1;
//                                 break;
//                             }
//                             case X: {
//                                 result->components[entry].paulis[k] = Z;
//                                 imagCounter -= 1;
//                                 break;
//                             }
//                             case Y: {
//                                 result->components[entry].paulis[k] = ID;
//                                 break;
//                             }
//                             default: {
//                                 printf("Error in fuseCompgates: Invalid value for Pauli.\n");
//                                 exit(1);
//                             } 
//                         }
//                         break;                        
//                     }
//                     default: {
//                         printf("Error in fuseCompgates: Invalid value for Pauli.\n");
//                         exit(1);
//                     }
//                 }
//             }
//             multiplyByIPower(&(result->coefficients[entry]), imagCounter);
//         }
//     }
// }

/*
 * =================================================================================================
 *                                              component applications
 * =================================================================================================
 */

void applyPauliString(state_t* state, const pauli_t paulistr[]) {
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
                printf("Error in applyPauliString: Invalid value for Pauli.\n");
                exit(1);
            }
        }
    }
}

void applyOperatorPauli(state_t* state, const pauliOp_t* op) {
    state_t tmp;
    state_t tmpSum;

    stateInitEmpty(&tmp, state->qubits);
    stateInitEmpty(&tmpSum, state->qubits);

    for (complength_t i = 0; i < op->length; ++i) {
        stateCopyVector(&tmp, state->vector);
        applyPauliString(&tmp, op->components + (i * op->qubits));
        cscalarVecMulInPlace(op->coefficients[i], tmp.vector, state->dimension);
        cvecAddInPlace(tmpSum.vector, tmp.vector, state->dimension);
    }

    stateFreeVector(&tmp);
    stateCopyVector(state, tmpSum.vector);
    stateFreeVector(&tmpSum);
}

/*
This function applies an observable, which is a sum of Pauli strings with real coefficients, to a 
state.

Input:
    state_t state:                  state that the observable is to be applied on
    const pauliObs_t observable:    observable that is applied to the state    
*/
void applyObservablePauli(state_t* state, const pauliObs_t* observable) {
    state_t tmp;                                // temporary state each Pauli string of the
                                                // observable is applied to

    state_t tmpSum;                             // temporary state holding each intermediate sum
                                                // of the observable's summands application to the
                                                // input state

    stateInitEmpty(&tmp, state->qubits);         // initialize tmp's statevector to all zeros

    stateInitEmpty(&tmpSum, state->qubits);      // initialize tmpSum's statevector to all zeros

    /*
    Starting with the first Pauli string, corresponding to the first qubits entries in the
    observable's components, apply the Pauli string, multiply by the corresponding coefficient and
    add the resulting state vector to the temporary sum state for each Pauli string.
    Note that the Pauli strings are in a continuous array, which means that the applyPauliString
    function takes blocks of length observable.qubits and then moves forward by the same amount.
    */
    for (complength_t i = 0; i < observable->length; ++i) {
        stateCopyVector(&tmp, state->vector);   // copy the input state's vector to tmp's
                                                // statevector

        applyPauliString(&tmp, observable->components + (i * observable->qubits));
                                                // apply the Pauli string according to the current
                                                // block in observable.components to the tmp state

        cscalarVecMulInPlace((cplx_t) observable->coefficients[i], tmp.vector, state->dimension);
                                                // multiply the tmp state by the coefficient
                                                // to observable.coefficient according
        cvecAddInPlace(tmpSum.vector, tmp.vector, state->dimension);
                                                // add the state vector in tmp state to the
                                                // intermediate state vector in tmpSum
    }
    stateFreeVector(&tmp);
    stateCopyVector(state, tmpSum.vector);      // copy the final state in tmpSum to the input
                                                // state's vector
    stateFreeVector(&tmpSum);
}

void evolveWithPauliString(state_t* state, const pauli_t paulistr[], double angle) {
    cplx_t* tmpVector = (cplx_t*) malloc(sizeof(cplx_t) * state->dimension);

    for (dim_t entry = 0; entry < state->dimension; ++entry) {
        tmpVector[entry] = state->vector[entry];
        tmpVector[entry] *= cos(angle);
        state->vector[entry] *= -I * sin(angle);
    }

    applyPauliString(state, paulistr);
    for (dim_t entry = 0; entry < state->dimension; ++entry) {
        state->vector[entry] += tmpVector[entry];
    }

    free(tmpVector);
}

void evolveWithTrotterizedObservablePauli(state_t* state, const pauliObs_t* observable, \
                                          double angle) {
    for (complength_t i = 0; i < observable->length; ++i) {
        evolveWithPauliString(state, observable->components + (i * observable->qubits), \
                              angle * observable->coefficients[i]);
    }
}
