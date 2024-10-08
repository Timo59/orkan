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
    stateInitEmpty(&tmp, state->qubits);
    cplx_t* tmpSum = calloc(state->dim, sizeof(cplx_t));    // Temporary vector holding the intermediate sum
                                                                        // of state vectors
    for (complength_t i = 0; i < op->length; ++i) {                         // Iterate the operator's terms
        stateCopyVector(&tmp, state->vec);                      // Copy input's state vector to tmp
        applyPauliStr(&tmp, op->comps + (i * op->qubits));     // Apply the Pauli string to it
        __LAPACK_double_complex alpha = op->coeffs[i];
        cblas_zaxpy(N, &alpha, tmp.vec, 1, tmpSum, 1);  // Multiply the terms' coefficient and add the
    }                                                                       // result to the intermediate sum

    stateFreeVector(&tmp);
    stateFreeVector(state);
    state->vec = tmpSum;
}

/*
 * =====================================================================================================================
 *                                              Apply Pauli observable
 * =====================================================================================================================
 */

/*
 * This function applies a Pauli observable to a quantum state.
 *
 * Input:
 *      state_t* state:         State of the qubit system
 *      pauliOp_t* op:          Pauli observable, i.e., a real linear combination of Pauli strings
 *
Output:
 *      This function has no return value, but sums the outcome of single Pauli string applications weighted by the
 *      coefficients in place.
 */
void applyObsPauli(state_t* state, const pauliObs_t* obs) {
    state_t tmp;                                                        // Temporary state to apply each Pauli string to
    stateInitEmpty(&tmp, state->qubits);
    cplx_t* tmpSum = calloc(state->dim, sizeof(cplx_t));    // Temporary vector holding the intermediate sum
                                                                        // of state vectors
    for (complength_t i = 0; i < obs->length; ++i) {                                // Iterate the observable's terms
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

void applyObsPauliBlas(state_t* state, const pauliObs_t* observable) {
    __LAPACK_int N = (__LAPACK_int) state->dim;
    state_t tmp;                                        // temporary vector to apply each Pauli string to
    stateInitEmpty(&tmp, state->qubits);        // initialize tmp's statevector to all zeros
    /* temporary vector holding each intermediate sum of altered input vectors */
    cplx_t* tmpSum = calloc(state->dim, sizeof(cplx_t));


    for (complength_t i = 0; i < observable->length; ++i) {
        /* For each term, copy the input statevector, apply the Pauli string and multiply the resulting vector with the
         * Pauli string's coefficient. Add the result to the temporary sum vector.
         */
        cblas_zcopy(N, state->vec, 1, tmp.vec, 1);
        applyPauliStr(&tmp, observable->comps + (i * observable->qubits));
        __LAPACK_double_complex alpha = observable->coeffs[i];
        cblas_zaxpy(N, &alpha, tmp.vec, 1, tmpSum, 1);
    }
    stateFreeVector(&tmp);
    stateFreeVector(state);
    state->vec = tmpSum;
}

void applyObsPauliOmp(state_t* state, const pauliObs_t* observable) {
    cplx_t* tmpSum = calloc(state->dim, sizeof(cplx_t));

# pragma omp parallel default(none) shared(observable, state, tmpSum)
    {
        state_t tmp;
        stateInitEmpty(&tmp, state->qubits);

#pragma omp for
        for (complength_t i = 0; i < observable->length; ++i) {
            stateCopyVector(&tmp, state->vec);
            applyPauliStr(&tmp, observable->comps + (i * observable->qubits));
            cscalarVecMulInPlace((cplx_t) observable->coeffs[i], tmp.vec, state->dim);

#pragma omp critical(sum)
            {
                cvecAddInPlace(tmpSum, tmp.vec, state->dim);
            }
        }
        stateFreeVector(&tmp);
    }
    stateFreeVector(state);
    state->vec = tmpSum;
}

void applyObsPauliBlas_omp(state_t* state, const pauliObs_t* observable) {
    __LAPACK_int N = (__LAPACK_int) state->dim;
    cplx_t* tmpSum = calloc(state->dim, sizeof(cplx_t));

# pragma omp parallel default(none) shared(N, observable, state, tmpSum)
    {
        state_t tmp;
        stateInitEmpty(&tmp, state->qubits);

#pragma omp for
        for (complength_t i = 0; i < observable->length; ++i) {
            cblas_zcopy(N, state->vec, 1, tmp.vec, 1);
            applyPauliStr(&tmp, observable->comps + (i * observable->qubits));
            __LAPACK_double_complex alpha = observable->coeffs[i];

#pragma omp critical(sum)
            {
                cblas_zaxpy(N, &alpha, tmp.vec, 1, tmpSum, 1);
            }
        }
        stateFreeVector(&tmp);

    }
    cblas_zcopy(N, tmpSum, 1, state->vec, 1);
}

void applyObsPauliBlas_omp_new(state_t* state, const pauliObs_t* observable) {
    __LAPACK_int N = (__LAPACK_int) state->dim;
    __LAPACK_double_complex one = 1;
    cplx_t* result = calloc(state->dim, sizeof(cplx_t));

# pragma omp parallel default(none) shared(N, observable, one, state, result)
    {
        state_t tmp;
        stateInitEmpty(&tmp, state->qubits);
        cplx_t* tmpSum = calloc(state->dim, sizeof(cplx_t));

#pragma omp for
        for (complength_t i = 0; i < observable->length; ++i) {
            cblas_zcopy(N, state->vec, 1, tmp.vec, 1);
            applyPauliStr(&tmp, observable->comps + (i * observable->qubits));
            __LAPACK_double_complex alpha = observable->coeffs[i];
            cblas_zaxpy(N, &alpha, tmp.vec, 1, tmpSum, 1);
        }
        stateFreeVector(&tmp);

#pragma omp critical(sum)
        {
            cblas_zaxpy(N, &one, tmpSum, 1, result, 1);
        }
        free(tmpSum);
    }
    stateFreeVector(state);
    stateInitVector(state, result, state->qubits);
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
    cplx_t* tmp = malloc(sizeof(cplx_t) * state->dim);      // Temporary vector holding state vector multiplicated
                                                                // with cos(angle)
    for (dim_t entry = 0; entry < state->dim; ++entry) {        // Iterate entries of the input's state vector
        tmp[entry] = cos(angle) * state->vec[entry];            // Store it time cos(angle) in tmp
        state->vec[entry] *= -I * sin(angle);                   // Multiply it with -i * sin(angle)
    }
    applyPauliStr(state, paulistr);                             // Apply the Pauli string to the input state
    for (dim_t entry = 0; entry < state->dim; ++entry) {        // Iterate the entries of the input's state vector
        state->vec[entry] += tmp[entry];                        // Add the tmp's entry to it
    }
    free(tmp);
}

void evolvePauliStrBlas(state_t* state,
                        const pauli_t paulistr[],
                        double angle)
{
    __LAPACK_int N = (__LAPACK_int) state->dim;
    __LAPACK_double_complex alpha = cos(angle);
    cplx_t* tmp = calloc(sizeof(cplx_t), state->dim);

    cblas_zaxpy(N, &alpha, state->vec, 1, tmp, 1);
    applyPauliStr(state, paulistr);
    alpha = -I * sin(angle);
    cblas_zaxpy(N, &alpha, state->vec, 1, tmp, 1);

    free(state->vec);
    state->vec = tmp;
}

void evolvePauliStrOmp(state_t* state,
                       const pauli_t paulistr[],
                       double angle)
{
    cplx_t* tmp = malloc(sizeof(cplx_t) * state->dim);

#pragma omp parallel for default(none) shared(angle, state, tmp)
    for (dim_t entry = 0; entry < state->dim; ++entry) {
        tmp[entry] = state->vec[entry];
        tmp[entry] *= cos(angle);
        state->vec[entry] *= -I * sin(angle);
    }
    applyPauliStr(state, paulistr);

#pragma omp parallel for default(none) shared(state, tmp)
    for (dim_t entry = 0; entry < state->dim; ++entry) {
        state->vec[entry] += tmp[entry];
    }

    free(tmp);
}

/*
 * =====================================================================================================================
 *                      Evolve with Pauli observable with Trotter approximation to first order
 * =====================================================================================================================
 */

/*
 * This function evolves a state with an observable in Trotter approximation to first order
 *
 * Input:
 *      state_t* state:         State of the qubit system
 *      pauliObs_t*:            Pauli observable, i.e., a real linear combination of Pauli strings
 *      double angle:           Time of the evolution
 * Output:
 *      This function has no return value, but sequentially evolves the input state with Pauli strings as if the
 *      observable's components commute.
 */

void evolveObsPauliTrotter(state_t* state, const pauliObs_t* obs, double angle) {
    for (complength_t i = 0; i < obs->length; ++i) {
        evolvePauliStr(state, obs->comps + (i * obs->qubits), angle * obs->coeffs[i]);
    }
}

void evolveObsPauliTrotterDagger(state_t* state, const pauliObs_t* obs, double angle) {
    for (complength_t i = 0; i < obs->length; ++i) {
        evolvePauliStr(state,
                       obs->comps + ((obs->length - i - 1) * obs->qubits),
                       -angle * obs->coeffs[obs->length - i - 1]);
    }
}
