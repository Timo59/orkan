/*
 * =================================================================================================
 *                                              includes
 * =================================================================================================
 */

#include "exact.h"

/*
 * =================================================================================================
 *                                              expectation value
 * =================================================================================================
 */

/*
This function computes the complex expectation value of an operator.

Input:
    state_t state: a user defined structure holding the pointer to a double complex valued array
    corresponding to the state vector, the number of qubits as an unsigned integer of type qubits_t
    and the Hilbert space dimension aka the number of entries in the state vector as an unsigned
    integer of type dim_t
    compgate_t operator: a user defined structure holding the pointer to an array of components 
    comprising gates, a pointer to an array of coefficients for each sequence and the number of
    sequences in the compound gate

Output:
    A double complex carrying the expectation value of a possibly non-hermitian operator.
*/
cplx_t expValPauli(const state_t* state, const pauliOp_t* op) {
    register cplx_t tmp;                                        // complex double holding the
                                                                // expectation value of a single
                                                                // component on the input state

    state_t ket;                                                // a state struct holding the vector
                                                                // that the components are applied
                                                                // to

    register cplx_t result = 0;                                 // complex double holding the sum of 
                                                                // expectation values of components
                                                                // on the input state times the 
                                                                // coefficients

    /*
    Starting at the first sequence of gates, the integer i goes through all components and their
    respective coefficents. The result of the component applied to the input's state vector times
    the coefficient is stored in tmp which is used to gradually add up the resulting state vector
    after applying the compgate in ket.
    */

    stateInitEmpty(&ket, state->qubits);
    for (complength_t i = 0; i < op->length; ++i) {
        stateCopyVector(&ket, state->vector);

        applyPauliString(&ket, op->components + (i * op->qubits));

        tmp = cinnerProduct(state->vector, ket.vector, state->dimension);
                                                                // store the inner product of the
                                                                // the input state with ket, aka the
                                                                // expectation value, in tmp               

        result += (op->coefficients[i] * tmp);            // add the single component's
                                                                // expectation value times its
                                                                // coefficient to the existing
                                                                // result
    }
    stateFreeVector(&ket);

    return result;
}

/*
This function computes the real expectation value of an observable.

Input:
    state_t state: a user defined structure holding the pointer to a double complex valued array
    corresponding to the state vector, the number of qubits as an unisgned integer of type qubits_t
    and the Hilbert space dimension aka the number of entries in the state vector as an unsigned
    integer of type dim_t
    obs_t observable: a user defined structure holding the pointer to an array of observables
    comprising hermitian gates, a pointer to an array of real coefficients for each component and
    the number of components in the compound observable

Output:
    A double carrying the expectation value of a hermitian operator (observable).
*/
double expValObsPauli(const state_t* state, const pauliObs_t* observable) {
    register cplx_t tmp;                                        // double complex holding the
                                                                // expectation value of a single
                                                                // component on the input state

    state_t ket;                                                // a state struct holding the vector
                                                                // that the components are applied
                                                                // to

    register double result = 0;                                 // double holding the sum of
                                                                // expectation values of components
                                                                // on the input state times the 
                                                                // coefficients

    /*
    Starting at the first sequence of gates, the integer i goes through all components and their
    respective coefficents. The result of the component applied to the input's state vector times
    the coefficient is stored in tmp which is used to gradually add up the resulting state vector
    after applying the compgate in ket.
    */

    stateInitEmpty(&ket, state->qubits);
    for (complength_t i = 0; i < observable->length; ++i) {
        stateCopyVector(&ket, state->vector);
        applyPauliString(&ket, observable->components + (i * observable->qubits));
                                                                // applies the i-th component to the
                                                                // state vector stored in tmp

        tmp = cinnerProduct(state->vector, ket.vector, state->dimension);
                                                                // store the inner product of the
                                                                // the input state with ket, aka the
                                                                // expectation value, in tmp            

        result += (observable->coefficients[i] * creal(tmp));   // add the single component's
                                                                // expectation value times its
                                                                // coefficient to the existing
                                                                // result

    }
    stateFreeVector(&ket);
    
    return result;
}

cplx_t expValDiag(const state_t* state, const cplx_t op[]) {
    register cplx_t result = 0;

    for (dim_t i = 0; i < state->dimension; ++i) {        
        result += conj(state->vector[i]) * state->vector[i] * op[i];
    }
    
    return result;
}

double expValObsDiag(const state_t* state, const double observable[]) {
    register double result = 0;

    for (dim_t i = 0; i < state->dimension; ++i) {        
        result += pow(cabs(state->vector[i]), 2) * observable[i];
    }
    
    return result;
}

double expValObs(const state_t* state, const obs_t* observable) {
    switch(observable->type) {
        case DIAG: {
            return expValObsDiag(state, observable->diagObs);
        }
        case PAULI: {
            return expValObsPauli(state, observable->pauliObs);
        }
        default: {
            perror("Error in applyOperator: Invalid operator type");
            exit(1);
        }
    }
}

/*
 * =================================================================================================
 *                                              apply PQC
 * =================================================================================================
 */

void applyPQC(state_t* state, const double params[], const obs_t* evoOps[], depth_t circdepth) {
    for (depth_t i = 0; i < circdepth; ++i) {
        evolveWithTrotterizedObservable(state, evoOps[i], params[i]);
    }
}

/*
 * =================================================================================================
 *                                              expectation value PQC
 * =================================================================================================
 */

double expValPQC(const state_t* state, const double params[], const obs_t* observable, \
                 const obs_t* evoOps[], depth_t circdepth) {
    state_t copy;
    stateInitEmpty(&copy, state->qubits);
    stateCopyVector(&copy, state->vector);

    for (depth_t i = 0; i < circdepth; ++i) {
        evolveWithTrotterizedObservable(&copy, evoOps[i], params[i]);
    }

    double result = expValObs(&copy, observable);
    stateFreeVector(&copy);
    return result;
}

/*
 * =================================================================================================
 *                                              gradient PQC
 * =================================================================================================
 */

/*
This function computes the gradient's specified component of a parametrised quantum circuit at the
current point in parameter space.
We have to assume that each parametrised quantum gate consists of a sum of commuting Pauli strings.

Input:
    state_t* state:         user defined structure holding the vector of the initial state inserted
                            into the parametrised quantum circuit, the number of qubits it is
                            defined on and the corresponding Hilbert space dimension
    double point:           pointer to an array of length circuit depth holding the current values
                            of all parameters
    obs_t observable:       observable whose expectation value's descent direction in parameter
                            space one is interested in, i.e., whose expectation value is to be 
                            optimised
    compgate_t* evoOps:     array of length circuit depth holding the evolution operators (for now
                            only sums of Pauli strings) of the PQC
    depth_t direction:      unsigned integer ranging from 0 to circuit depth, determining with
                            respect to which parameter component the derivative should be computed
    depth_t circdepth:      unsigned integer holding the number of evolution operators in the PQC
                            aka its circuit depth

Output:
    An array representing the PQC's gradient at the specified point in parameter space.
*/
double* gradientPQC(const state_t* state, const double params[], const obs_t* observable, \
                    const obs_t* evoOps[], depth_t circdepth) {
    state_t bra;
                                                                        
    state_t ket;

    state_t tmp;

    register double* result = (double*) malloc(circdepth * sizeof(double));

//    printf("Input state to gradientPQC:\n");
//    cvectorPrint(state->vector, state->dimension);
//
//    printf("\nInput parameters to gradientPQC:\n");
//    rvectorPrint(params, circdepth);
//
//    printf("\nInput observable to gradient PQC:\n");
//    double* obsInput = observable->diagObs;
//    rvectorPrint(obsInput, state->dimension);
//
//    printf("\nInput evolution operators to gradient PQC:\n");
//    for (depth_t i = 0; i < circdepth; ++i) {
//        switch (evoOps[i]->type) {
//            case DIAG: {
//                double* evoOpsInput = evoOps[i]->diagObs;
//                rvectorPrint(evoOpsInput, state->dimension);
//                break;
//            }
//            case PAULI: {
//                pauliObs_t *evoOpsInput = evoOps[i]->pauliObs;
//                for (complength_t j = 0; j < evoOpsInput->length; ++j) {
//                    for (qubit_t k = 0; k < evoOpsInput->qubits; ++k) {
//                        printf("%d", evoOpsInput->components[k + (j * evoOpsInput->qubits)]);
//                    }
//                    printf(", %f\n", evoOpsInput->coefficients[j]);
//                }
//                break;
//            }
//            default: {
//                printf("Got %d instead of %d or %d", evoOps[i]->type, DIAG, PAULI);
//                exit(1);
//            }
//        }
//    }
//
//    printf("\nInput circdepth to gradient PQC:\n %d\n", circdepth);

    stateInitEmpty(&bra, state->qubits);        // Initialize bra with the input state's vector
    stateCopyVector(&bra, state->vector);

    stateInitEmpty(&ket, state->qubits);        // Initialize ket to the all zero state

    stateInitEmpty(&tmp, state->qubits);        // Initialize tmp with the input state's vector
    stateCopyVector(&tmp, state->vector);

    /*
    Evolve bra with all evolutional operators
    */
    for (depth_t i = 0; i < circdepth; ++i) {
        evolveWithTrotterizedObservable(&bra, evoOps[i], params[i]);
    }

    /*
    Apply the observable to bra
    */
    applyObservable(&bra, observable);

    /*
    Starting at k=0, for every integer k upt to the circuit depth, evolve tmp with the k-th 
    evolutional operator. Then copy its state vector to ket, apply the observable to ket and evolve
    ket with the remaining evolutional operators and add two times the imaginary part of bra's and
    ket's overlap to the gradient. 
    */
    for (depth_t k = 0; k < circdepth; ++k) {
        evolveWithTrotterizedObservable(&tmp, evoOps[k], params[k]);

        stateCopyVector(&ket, tmp.vector);
        applyObservable(&ket, evoOps[k]);

        for (depth_t l = k + 1; l < circdepth; ++l) {
            evolveWithTrotterizedObservable(&ket, evoOps[l], params[l]);
        }

        result[k] = 2 * cimag(stateOverlap(&bra, &ket));
    }
    stateFreeVector(&bra);
    stateFreeVector(&ket);
    stateFreeVector(&tmp);
    return result;
}

/*
 * =================================================================================================
 *                                              hessian PQC
 * =================================================================================================
 */

/*
This function returns the hessian according to a PQC at the current point the parameter space.
*/

double* hessianPQC(const state_t* state, const double params[], const obs_t* observable, \
                   const obs_t* evoOps[], depth_t circdepth) {

    state_t tmpBra;     // temporary state holding each intermediate evolution according to the
                        // hessian's column index both bra vectors have in common
    
    state_t tmpKet;     // temporary state holding each intermediate evolution according to the
                        // hessian's row index both ket vectors have in common
    
    state_t tmpTmpKet;  // temporary state holding each intermediate evolution and applying the  
                        // current evolutional observable according to the hessian's row as well as
                        // the evolution according to the difference of the hessian's row and column
                        // index

    state_t bra1;       // state inheriting the intermediate evolution from tmpBra, then altered by
                        // the evolution operator according to the hessian's column index and
                        // evolved with the remaining evolution operators and used for the first
                        // stateOverlap

    state_t bra2;       // state evolved with all evolution operators and used for the second
                        // stateOverlap

    state_t ket1;       // state inheriting the intermediate evolution from tmpTmpKet, then evolved
                        // by the remaining evolution operators according to the difference of the
                        // hessian's row and column index and used for the first stateOverlap

    state_t ket2;       // state inheriting the intermediate evolution form tmpTmpKet, then altered
                        // by the evolution operator according to the hessian's column index and 
                        // evolved by the remaining evolution operators

    register double* result = (double*) malloc((circdepth * circdepth) * sizeof(double));
                        // array holding the upper triangle of the hessian

    stateInitEmpty(&tmpBra, state->qubits);
    stateInitEmpty(&tmpKet, state->qubits);
    stateInitEmpty(&tmpTmpKet, state->qubits);
    stateInitEmpty(&bra1, state->qubits);
    stateInitEmpty(&bra2, state->qubits);
    stateInitEmpty(&ket1, state->qubits);
    stateInitEmpty(&ket2, state->qubits);

    stateCopyVector(&tmpBra, state->vector);
    stateCopyVector(&tmpKet, state->vector);

    stateCopyVector(&bra2, state->vector);

    /*
    Evolve the bra vector of the second stateOverlap with all evolution operators and corresponding
    parameters
    */
    for (depth_t i = 0; i < circdepth; ++i) {
        evolveWithTrotterizedObservable(&bra2, evoOps[i], params[i]);
    }

    /*
    Starting in the hessian's first row, counted by k, evolve the tmpKet with all operators up to k,
    copy its state vector to tmpTmpKet and apply the k-th evolution operator to tmpTmpKet. Then,
    compute the diagonal element of the hessian. After allcolumns in the hessian's k-th row are
    computed, reset tmpTmpKet as well as tmpBra and increase the index by entries in the upper
    triangle's row k, i.e., circdepth - k.
    */
    for (depth_t k = 0; k < circdepth; ++k) {
        evolveWithTrotterizedObservable(&tmpKet, evoOps[k], params[k]);
        stateCopyVector(&tmpTmpKet, tmpKet.vector);
        applyObservable(&tmpTmpKet, evoOps[k]);

        /*
        The diagonal elements, i.e., k=l are computed seperately at each new iteration of the first
        index because they only appear once in the matrix.
        */
        stateCopyVector(&bra1, tmpTmpKet.vector);
        stateCopyVector(&ket1, tmpTmpKet.vector);
        stateCopyVector(&ket2, tmpTmpKet.vector);
        applyObservable(&ket2, evoOps[k]);

        for (depth_t kk = k + 1; kk < circdepth; ++kk) {
            evolveWithTrotterizedObservable(&bra1, evoOps[kk], params[kk]);
            evolveWithTrotterizedObservable(&ket1, evoOps[kk], params[kk]);
            evolveWithTrotterizedObservable(&ket2, evoOps[kk], params[kk]);
        }
        result[k * circdepth + k] = 2 * creal(stateOverlap(&bra1, &ket1));
        result[k * circdepth + k] -= 2 * creal(stateOverlap(&bra2, &ket2));
        stateFreeVector(&bra1);
        stateFreeVector(&ket1);
        stateFreeVector(&ket2);

        /*
        Starting at the k-th column in the hessian's k-th row, counted by l, evolve tmpBra with all
        operators up to l, copy its state vector to bra1 and apply the l-th evolution operator to 
        bra1 before evolving bra1 with all the remaining operators.
        Moreover, evolve tmpTmpKet with all operators whose index is within [k,l], copy its state
        vector to both, ket1 and ket2 and evolve ket1 with the remaining operators before applying
        the observable to it. Then apply the operator corresponding to l to ket2, evolve it with the
        remaining operators and apply the observable to it. The sum of the two stateOverlaps' real parts
        are appended to the result and state vectors of bra1, ket1 and ket2 are reset.
        */
        for (depth_t l = k + 1; l < circdepth; ++l) {
            evolveWithTrotterizedObservable(&tmpBra, evoOps[l], params[l]);
            stateCopyVector(&bra1, tmpBra.vector);
            applyObservable(&bra1, evoOps[l]);

            for (depth_t ll = l + 1; ll < circdepth; ++ll) {
                evolveWithTrotterizedObservable(&bra1, evoOps[ll], params[ll]);
            }

            for (depth_t kk = k + 1; kk <= l; ++kk) {
                evolveWithTrotterizedObservable(&tmpTmpKet, evoOps[kk], params[kk]);
            }
            stateCopyVector(&ket1, tmpTmpKet.vector);
            stateCopyVector(&ket2, tmpTmpKet.vector);

            for (depth_t kk = l + 1; kk < circdepth; ++kk) {
                evolveWithTrotterizedObservable(&ket1, evoOps[kk], params[kk]);
            }
            applyObservable(&ket1, observable);

            applyObservable(&ket2, evoOps[l]);
            for (depth_t ll = l + 1; ll < circdepth; ++ll) {
                evolveWithTrotterizedObservable(&ket2, evoOps[ll], params[ll]);
            }
            applyObservable(&ket2, observable);

            result[k * circdepth + l] = 2 * creal(stateOverlap(&bra1, &ket1));
            result[k * circdepth + l] -= 2 * creal(stateOverlap(&bra2, &ket2));
            result[l * circdepth + k] = 2 * creal(stateOverlap(&bra1, &ket1));
            result[l * circdepth + k] -= 2 * creal(stateOverlap(&bra2, &ket2));

            stateFreeVector(&bra1);
            stateFreeVector(&ket1);
            stateFreeVector(&ket2);
        }
        stateFreeVector(&tmpTmpKet);
        stateFreeVector(&tmpBra);
    }
    stateFreeVector(&bra2);
    stateFreeVector(&tmpKet);

    return result;
}

/*
 * =================================================================================================
 *                                    Moment matrices
 * =================================================================================================
 */

/*
This function returns the moment matrices E_{ij} and H_{ij} w.r.t. the identity, the QAOA mixer and
the phase separator as search unitaries.

 Input:
    state_t* state:         state of the system
    obs_t hamiltonian:      the objective function's hamiltonian
    double angleB:          mixer's evolution parameter
    double angleH:          phase separator's evolution parameter

Output:
    Array holding the matrix E in column major order concatenated with the matrix H in column major
    order.
*/
cplx_t* momentMatQAOA(state_t* state, const obs_t* hamiltonian, double angleB, double angleH) {
    state_t tmpB;       // temporary state holding input state's evolution by angleB w.r.t. B

    state_t tmpH;       // temporary state holding input state's evolution by angleB w.r.t. H

    register cplx_t* result = (cplx_t*) malloc((18) * sizeof(cplx_t));

    /*
     * Initialize tmpB
     */
    stateInitEmpty(&tmpB, state->qubits);
    stateCopyVector(&tmpB, state->vector);
    for (qubit_t i = 0; i < tmpB.qubits; ++i) {
        applyRX(&tmpB, i, 2*angleB);
    }

    /*
     * Initialize tmpH
     */
    stateInitEmpty(&tmpH, state->qubits);
    stateCopyVector(&tmpH, state->vector);
    evolveWithTrotterizedObservable(&tmpH, hamiltonian, angleH);

    /*
     * Save E in column major form to the first 9 entries of result
     */
    result[0] = 1. + 0. * I;
    result[1] = stateOverlap(&tmpB, state);
    result[2] = stateOverlap(&tmpH, state);
    result[3] = conj(result[1]);
    result[4] = 1. + 0. * I;
    result[5] = stateOverlap(&tmpH, &tmpB);
    result[6] = conj(result[2]);
    result[7] = conj(result[5]);
    result[8] = 1. + 0. * I;

    /*
     * Save H in column major form to the other 9 entries of result starting with the diagonal
     * elements
     */
    result[9] = expValObs(state, hamiltonian);
    result[13] = expValObs(&tmpB, hamiltonian);
    result[17] = expValObs(&tmpH, hamiltonian);
    /*
     * Apply the hamiltonian to tmpB and compute all expectation values involving it
     */
    applyObservable(&tmpB, hamiltonian);
    result[10] = stateOverlap(&tmpB, state);
    result[12] = conj(result[10]);
    result[14] = stateOverlap(&tmpH, &tmpB);
    result[16] = conj(result[14]);
    /*
     * Apply the hamiltonian to tmpH and compute the remaining expectation value
     */
    applyObservable(&tmpH, hamiltonian);
    result[11] = stateOverlap(&tmpH, state);
    result[15] = conj(result[11]);

    /*
     * Free temporary state vectors
     */
    stateFreeVector(&tmpB);
    stateFreeVector(&tmpH);
    return result;
}

/*
This function computes the moment matrices of the set of unitaries generated by exponentiating Pauli
strings.

Input:
    state_t* state          Current state of the system, number of its qubits and the resulting
                            Hilbert space dimension
    obs_t* searchHerms      Observables to be exponentiated to give the search unitaries
    double* angles          Angles defining the search unitaries
    depth_t searchDim       Number of search unitaries
    obs_t* observable       Observable whose expectation value is the objective

Output: 
    The moment matrix <state|U_i^* U_j|state> appended by <state|U_i^* H U_j|state> in column major
    form.

*/

cplx_t* momentMatrices(const state_t* state, const obs_t searchHerms[], double* angles, \
                       depth_t searchDim, obs_t* observable) {

    state_t bra;        // state holding the bra, i.e. the complex conjugate of the input states 
                        // times search unitary

    state_t ket;        // state holding the ket, i.e. the input states times search unitary

    register cplx_t* result = malloc((searchDim*searchDim) * sizeof(double));

    stateInitEmpty(&bra, state->qubits);
    stateInitEmpty(&ket, state->qubits);

    /*
    Starting in the first column, the integer i counts the column of the moment matrices.
    */
    for (depth_t i = 0; i < searchDim; ++i) {
        stateCopyVector(&ket, state->vector);
        evolveWithTrotterizedObservable(&ket, searchHerms + i, angles[i]);
                                                                        // Evolve ket with the i-th
                                                                        // observable by the i-th
                                                                        // angle

        result[i*searchDim + i] = 1.;                                   // The diagonal entries
                                                                        // <state|U_i^* U_i|state>
                                                                        // are always equal to one

        result[(searchDim*searchDim) + i*searchDim + i] = (cplx_t) expValObs(&ket, observable);
                                                                        // The diagonal entries of
                                                                        // the objective matrix are
                                                                        // expectation values of the
                                                                        // observable   

        /*
        Starting at the i-th entry of the i-th column, j counts through all entries of the moment
        matrices' lower triangle.
        */
        for (depth_t j = i+1; j < searchDim; ++j) {
            stateCopyVector(&bra, state->vector);
            evolveWithTrotterizedObservable(&bra, searchHerms + j, angles[j]);
                                                                        // Evolve bra with the j-th
                                                                        // observable by the j-th
                                                                        // angle

            result[i*searchDim + j] = stateOverlap(&bra, &ket);         // To the j-th entry of the
                                                                        // i-th column add overlap
                                                                        // <state|U_j^* U_i|state>

            result[j*searchDim + i] = conj(result[i*searchDim + j]);    // By hermiticity of E, its 
                                                                        // i-th entry of the j-th
                                                                        // column is the complex
                                                                        // conjugate of the j-th
                                                                        // entry in the i-th column

            applyObservable(&bra, observable);                          // Apply the observable to 
                                                                        // bra
            result[(searchDim*searchDim) + i*searchDim + j] = stateOverlap(&bra, &ket);
                                                                        // and store its overlap
                                                                        // with ket in the j-th
                                                                        // entry in the i-th column
                                                                        // of the objective matrix

            result[(searchDim*searchDim) + j*searchDim + i] = conj(result[(searchDim*searchDim) \
                                                                      + i*searchDim + j]);
                                                                        // By hermiticity of H, its
                                                                        // i-th entry of the j-th
                                                                        // column is the complex
                                                                        // conjugate of the j-th
                                                                        // entry in the i-th column
            stateFreeVector(&bra);

        }
        stateFreeVector(&ket);
    }
    return result;
}