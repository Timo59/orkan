/*
 * =====================================================================================================================
 *                                                      includes
 * =====================================================================================================================
 */

#include "pauli.h"

/*
 * =====================================================================================================================
 *                                              compound gate constructor
 * =====================================================================================================================
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
 * =====================================================================================================================
 *                                              compound gate fusion
 * =====================================================================================================================
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
 * =====================================================================================================================
 *                                              component applications
 * =====================================================================================================================
 */

/*
 * This function applies a Pauli string to a quantum state.
 *
 * Input:
 *      state_t* state:         State of the qubit system
 *      pauli_t paulistr[]:     Pauli string, i.e., array of Pauli operators
 *
 * Output:
 *      This function has no return value, but changes the statevector in place.
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


/*
 * This function applies a Pauli operator to a quantum state.
 *
 * Input:
 *      state_t* state:         State of the qubit system
 *      pauliOp_t* op:          Pauli operator, i.e., a complex linear combination of Pauli strings
 *
 * Output:
 *      This function has no return value, but changes the statevector in place.
 */
void applyOperatorPauli(state_t* state, const pauliOp_t* op) {
    state_t tmp;        // temporary vector to apply each Pauli string to
    stateInitEmpty(&tmp, state->qubits);
    state_t tmpSum;     // temporary vector holding each intermediate sum of altered input vectors
    stateInitEmpty(&tmpSum, state->qubits);

    for (complength_t i = 0; i < op->length; ++i) {
        /* For each term, copy the input statevector, apply the Pauli string and multiply the resulting vector with the
         * Pauli string's coefficient. Add the result to the temporary sum vector.
         */
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
 * This function applies a Pauli observable to a quantum state.
 *
 * Input:
 *      state_t* state:         State of the qubit system
 *      pauliOp_t* op:          Pauli operator, i.e., a real linear combination of Pauli strings
 *
 * Output:
 *      This function has no return value, but changes the statevector in place.
 */
void applyObservablePauli(state_t* state, const pauliObs_t* observable) {
    state_t tmp;                                        // temporary vector to apply each Pauli string to
    stateInitEmpty(&tmp, state->qubits);        // initialize tmp's statevector to all zeros
    state_t tmpSum;                                     // temporary vector holding each intermediate sum of altered input
                                                        // vectors
    stateInitEmpty(&tmpSum, state->qubits);     // initialize tmpSum's statevector to all zeros

    for (complength_t i = 0; i < observable->length; ++i) {
        /* For each term, copy the input statevector, apply the Pauli string and multiply the resulting vector with the
         * Pauli string's coefficient. Add the result to the temporary sum vector.
         */
        stateCopyVector(&tmp, state->vector);
        applyPauliString(&tmp, observable->components + (i * observable->qubits));
        cscalarVecMulInPlace((cplx_t) observable->coefficients[i], tmp.vector, state->dimension);
        cvecAddInPlace(tmpSum.vector, tmp.vector, state->dimension);
    }
    stateFreeVector(&tmp);
    stateFreeVector(state);
    state->vector = tmpSum.vector;
}

void applyObservablePauliBlas(state_t* state, const pauliObs_t* observable) {
    __LAPACK_int N = (__LAPACK_int) state->dimension;
    state_t tmp;                                        // temporary vector to apply each Pauli string to
    stateInitEmpty(&tmp, state->qubits);        // initialize tmp's statevector to all zeros
    state_t tmpSum;                                     // temporary vector holding each intermediate sum of altered input
                                                        // vectors
    stateInitEmpty(&tmpSum, state->qubits);     // initialize tmpSum's statevector to all zeros

    for (complength_t i = 0; i < observable->length; ++i) {
        /* For each term, copy the input statevector, apply the Pauli string and multiply the resulting vector with the
         * Pauli string's coefficient. Add the result to the temporary sum vector.
         */
        cblas_zcopy(N, state->vector, 1, tmp.vector, 1);
        applyPauliString(&tmp, observable->components + (i * observable->qubits));
        __LAPACK_double_complex alpha = observable->coefficients[i];
        cblas_zaxpy(N, &alpha, tmp.vector, 1, tmpSum.vector, 1);
    }
    stateFreeVector(&tmp);
    stateFreeVector(state);
    state->vector = tmpSum.vector;
}

void applyObservablePauliOmp(state_t* state, const pauliObs_t* observable) {
    state_t tmp;                                        // temporary vector to apply each Pauli string to
    stateInitEmpty(&tmp, state->qubits);        // initialize tmp's statevector to all zeros
    state_t tmpSum;                                     // temporary vector holding each intermediate sum of altered input
    // vectors
    stateInitEmpty(&tmpSum, state->qubits);     // initialize tmpSum's statevector to all zeros

# pragma omp parallel for default(none) private(tmp) shared(observable, state, tmpSum)
    for (complength_t i = 0; i < observable->length; ++i) {
        /* For each term, copy the input statevector, apply the Pauli string and multiply the resulting vector with the
         * Pauli string's coefficient. Add the result to the temporary sum vector.
         */
        stateCopyVector(&tmp, state->vector);
        applyPauliString(&tmp, observable->components + (i * observable->qubits));
        cscalarVecMulInPlace((cplx_t) observable->coefficients[i], tmp.vector, state->dimension);
        cvecAddInPlace(tmpSum.vector, tmp.vector, state->dimension);
    }
    stateFreeVector(&tmp);
    stateFreeVector(state);
    state->vector = tmpSum.vector;
}

/*
 * This function evolves a state with the with a Pauli string.
 *
 * Input:
 *      state_t* state: State to be evolved
 *      pauli_t paulistr[]:
 */
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
