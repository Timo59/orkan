/*
 * =====================================================================================================================
 *                                                      includes
 * =====================================================================================================================
 */

#include "exact.h"

/*
 * =====================================================================================================================
 *                                                  Mean value operator
 * =====================================================================================================================
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
cplx_t expValPauli(const state_t* state, const pauliOp_t* operator) {
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
    for (complength_t i = 0; i < operator->length; ++i) {
        stateCopyVector(&ket, state->vector);

        applyPauliString(&ket, operator->components + (i * operator->qubits));

        tmp = cinnerProduct(state->vector, ket.vector, state->dimension);
                                                                // store the inner product of the
                                                                // the input state with ket, aka the
                                                                // expectation value, in tmp               

        result += (operator->coefficients[i] * tmp);            // add the single component's
                                                                // expectation value times its
                                                                // coefficient to the existing
                                                                // result
    }
    stateFreeVector(&ket);

    return result;
}

cplx_t expValDiag(const state_t* state, const cplx_t operator[]) {
    register cplx_t result = 0;

    for (dim_t i = 0; i < state->dimension; ++i) {
        result += conj(state->vector[i]) * state->vector[i] * operator[i];
    }

    return result;
}

cplx_t expVal(const state_t* state, const op_t* operator) {
    switch(operator->type) {
        case DIAG: {
            return expValDiag(state, operator->diagOp);
        }
        case PAULI: {
            return expValPauli(state, operator->pauliOp);
        }
        default: {
            perror("Error in applyOperator: Invalid operator type");
            exit(1);
        }
    }
}

/*
 * =====================================================================================================================
 *                                                  Mean value observable
 * =====================================================================================================================
 */

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
 * =====================================================================================================================
 *                                                  Mean value PQC
 * =====================================================================================================================
 */

double expValObsPauliPQC(const state_t* state, const double params[], const pauliObs_t* observable, \
                         const obs_t* evoOps[], depth_t circdepth) {
    state_t bra;                    // temporary state holding the evolved input state
    state_t ket;                    // temporary state holding the evolved input state altered with each summand of obs

    cplx_t tmp;                     // cplx_t holding each summand of the observable's mean value

    register double result = 0;     // double holding the result

    /*
     * 1. Initialize tmp states
     */
    stateInitEmpty(&bra, state->qubits);
    stateInitEmpty(&ket, state->qubits);
    stateCopyVector(&bra, state->vector);

    /*
     * 2. Evolve tmp states with PQC
     */
    applyPQC(&bra, params, evoOps, circdepth);

    /*
     * 3. Calculate the mean of each weighted Pauli string in the observable on the input state
     */
    for (complength_t i = 0; i < observable->length; ++i) {
        stateCopyVector(&ket, bra.vector);
        applyPauliString(&ket, observable->components + (i * observable->qubits));

        tmp = cinnerProduct(bra.vector, ket.vector, bra.dimension);

        result += (observable->coefficients[i] * creal(tmp));

    }

    return result;
}

double expValObsDiagPQC(const state_t* state, const double params[], const double observable[], \
                        const obs_t* evoOps[], depth_t circdepth) {
    state_t tmp;                    // temporary state holding the evolved input state

    register double result = 0;     // double holding the result

    /*
     * 1. Initialize tmp state
     */
    stateInitEmpty(&tmp, state->qubits);
    stateCopyVector(&tmp, state->vector);

    /*
     * 2. Evolve tmp state with PQC
     */
    applyPQC(&tmp, params, evoOps, circdepth);

    /*
     * 3. Calculate the mean of each weighted Pauli string in the observable on the input state
     */
    for (dim_t i = 0; i < tmp.dimension; ++i) {
        result += pow(cabs(tmp.vector[i]), 2) * observable[i];
    }

    return result;
}

double expValObsPQC(const state_t* state, const double params[], const obs_t* observable, \
                    const obs_t* evoOps[], depth_t circdepth) {
    switch(observable->type) {
        case DIAG: {
            return expValObsDiagPQC(state, params, observable->diagObs, evoOps, circdepth);
        }
        case PAULI: {
            return expValObsPauliPQC(state, params, observable->pauliObs, evoOps, circdepth);
        }
        default: {
            perror("Error in applyOperator: Invalid operator type");
            exit(1);
        }
    }
}

/*
 * =====================================================================================================================
 *                                                      gradient PQC
 * =====================================================================================================================
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

    stateInitEmpty(&bra, state->qubits);        // Initialize bra with the input state's vector
    stateCopyVector(&bra, state->vector);

    stateInitEmpty(&ket, state->qubits);        // Initialize ket

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


double* approxGradientPQC(const state_t* state, const double params[], const obs_t* observable, \
                               const obs_t* evoOps[], depth_t circdepth, double epsilon) {

    register double* result = (double*) malloc(circdepth * sizeof(double));

    double* tmp_params = malloc(circdepth * sizeof(double ));       // temporary copy of the parameters
    for (depth_t i = 0; i < circdepth; ++i) {
        tmp_params[i] = params[i];
    }


    /*
     * 1. Subtract the mean value at the current parameter setting from all entries of the latter gradient
     */
    double mean = expValObsPQC(state, params, observable, evoOps, circdepth);

    /*
     * 2. For each direction shift the corresponding parameter calculate the slope by the difference method
     */
    for (depth_t i = 0; i < circdepth; ++i) {
        tmp_params[i] += epsilon;
        result[i] = (1./epsilon) * (expValObsPQC(state, tmp_params, observable, evoOps, circdepth) - mean);
        tmp_params[i] -= epsilon;
    }

    free(tmp_params);

    return result;
}

/*
 * =====================================================================================================================
 *                                                      hessian PQC
 * =====================================================================================================================
 */

/*
 * This function returns the hessian according to a PQC at the current point the parameter space.
 *
 * Input:
 *  state_t* state:         initial state of the parameterized quantum circuit
 *  const double params:    array of length circuit depth holding the current parameter values
 *  obs_t observable:       observable whose expectation value's hessian with respect to the parameterized quantum
 *                           circuit is the objective
 *  compgate_t* evoOps:     array of length circuit depth holding pointers to the evolution operators of the PQC
 *  depth_t circdepth:      unsigned integer holding the number of evolution operators in the PQC aka its circuit depth
 *
 * Output:
 *  An array representing the PQC's gradient at the specified point in parameter space.
 *
*/

double* hessianPQC(const state_t* state, const double params[], const obs_t* observable, const obs_t* evoOps[], \
                   depth_t circdepth) {

    state_t tmp;        // state holding each intermediate evolution of the state vector according to the outer loop aka
                        // the hessian's row index

    state_t tmpBra;     // state holding each intermediate evolution of the state vector according to the inner loop aka
                        // the hessian's column index
    
    state_t tmpKet;     // state descendant from tmp, altered by the evolution operator according to the hessian's row
                        // index and evolved with the remaining evolution operators up to the one according to the
                        // hessian's column index

    state_t bra1;       // state descendant from tmpBra, altered by the evolution operator according to the hessian's
                        // column index and evolved with the remaining evolution operators; used for the first overlap

    state_t bra2;       // state evolved with all evolution operators; used for the second overlap

    state_t ket1;       // state descendant from tmpKet, evolved with the remaining evolution operators; used for the
                        // first overlap

    state_t ket2;       // state descendant from tmpKet, altered by the evolution operator according to the hessian's
                        // column index and evolved with the remaining evolution operators; used for the second overlap

    register double* result = (double*) malloc((circdepth * circdepth) * sizeof(double));

    stateInitEmpty(&tmp, state->qubits);
    stateInitEmpty(&tmpBra, state->qubits);
    stateInitEmpty(&tmpKet, state->qubits);
    stateInitEmpty(&bra1, state->qubits);
    stateInitEmpty(&bra2, state->qubits);
    stateInitEmpty(&ket1, state->qubits);
    stateInitEmpty(&ket2, state->qubits);

    stateCopyVector(&tmp, state->vector);

    /*
     * 1a) Evolve the bra vector of the second stateOverlap with all evolution operators and corresponding parameters.
    */
    stateCopyVector(&bra2, state->vector);
    for (depth_t i = 0; i < circdepth; ++i) {
        evolveWithTrotterizedObservable(&bra2, evoOps[i], params[i]);
    }

    /*
     * 1b) The unsigned integer k iterates through all row indices of the hessian's upper triangle. In each iteration
     * tmp is evolved with the k-th evolution operator.
    */
    for (depth_t k = 0; k < circdepth; ++k) {
        evolveWithTrotterizedObservable(&tmp, evoOps[k], params[k]);

        /*
         * 2a) tmpKet inherits its state vector from tmp and applies the k-th evolution operator to it.
         */
        stateCopyVector(&tmpKet, tmp.vector);
        applyObservable(&tmpKet, evoOps[k]);

        /*
         * 2b) tmpBra inherits its state vector from tmp.
         */
        stateCopyVector(&tmpBra, tmp.vector);

        /*
         * 3a) The k-th diagonal element is computed separately because it only appears once in the matrix. bra1, ket1
         * and ket2 are all descendants of tmpKet. The k-th evolution operator is applied to ket2 a second time.
        */
        stateCopyVector(&bra1, tmpKet.vector);
        stateCopyVector(&ket1, tmpKet.vector);
        stateCopyVector(&ket2, tmpKet.vector);
        applyObservable(&ket2, evoOps[k]);

        /*
         * 3b) bra1, ket1 and ket2 are all evolved with the remaining evolution operators.
         */
        for (depth_t kk = k + 1; kk < circdepth; ++kk) {
            evolveWithTrotterizedObservable(&bra1, evoOps[kk], params[kk]);
            evolveWithTrotterizedObservable(&ket1, evoOps[kk], params[kk]);
            evolveWithTrotterizedObservable(&ket2, evoOps[kk], params[kk]);
        }

        /*
         * 3d) The observable is applied to ket1 and ket2.
         */
        applyObservable(&ket1, observable);
        applyObservable(&ket2, observable);

        /*
         * 3e) The diagonal entries are added once to the array.
         */
        result[k * circdepth + k] = 2 * creal(stateOverlap(&bra1, &ket1));
        result[k * circdepth + k] -= 2 * creal(stateOverlap(&bra2, &ket2));

        /*
         * 4) The unsigned integer l iterates through all column indices of the hessian's upper triangle without the
         * diagonal, i.e., l > k. In each iteration tmpBra and tmpKet are evolved with the l-th evolution operator
        */
        for (depth_t l = k + 1; l < circdepth; ++l) {
            evolveWithTrotterizedObservable(&tmpBra, evoOps[l], params[l]);
            evolveWithTrotterizedObservable(&tmpKet, evoOps[l], params[l]);

            /*
             * 4a) bra1 is a descendant of tmpBra, it is applied the l-th evolution operator and evolved with the
             * remaining evolution operators.
             */
            stateCopyVector(&bra1, tmpBra.vector);
            applyObservable(&bra1, evoOps[l]);
            for (depth_t ll = l + 1; ll < circdepth; ++ll) {
                evolveWithTrotterizedObservable(&bra1, evoOps[ll], params[ll]);
            }

            /*
             * 4b) ket1 is a descendant of tmpKet, is evolved with the remaining evolution operators and the observable
             * is applied to it
             */
            stateCopyVector(&ket1, tmpKet.vector);
            for (depth_t kk = l + 1; kk < circdepth; ++kk) {
                evolveWithTrotterizedObservable(&ket1, evoOps[kk], params[kk]);
            }

            applyObservable(&ket1, observable);


            /*
             * 4c) ket2 is a descendant of tmpKet and the l-th evolution operator is applied to it. It is evolved with
             * the remaining evolution operators and the observable is applied to it.
             */
            stateCopyVector(&ket2, tmpKet.vector);
            applyObservable(&ket2, evoOps[l]);
            for (depth_t ll = l + 1; ll < circdepth; ++ll) {
                evolveWithTrotterizedObservable(&ket2, evoOps[ll], params[ll]);
            }
            applyObservable(&ket2, observable);

            /*
             * 4d) The difference of the two overlap's real part are added to the positions corresponding to the l-th
             * entry of the k-th row and vice versa of the array.
             */
            result[k * circdepth + l] = 2 * creal(stateOverlap(&bra1, &ket1));
            result[k * circdepth + l] -= 2 * creal(stateOverlap(&bra2, &ket2));
            result[l * circdepth + k] = 2 * creal(stateOverlap(&bra1, &ket1));
            result[l * circdepth + k] -= 2 * creal(stateOverlap(&bra2, &ket2));
        }
    }

    /*
     * The state vectors of bra2 and tmp are freed. No more state vector is attached
     */
    stateFreeVector(&tmp);
    stateFreeVector(&tmpBra);
    stateFreeVector(&tmpKet);
    stateFreeVector(&bra1);
    stateFreeVector(&bra2);
    stateFreeVector(&ket1);
    stateFreeVector(&ket2);

    return result;
}


double* approxHessianPQC(const state_t* state, const double params[], const obs_t* observable, \
                         const obs_t* evoOps[], depth_t circdepth, gradPQC grad, double epsilon) {

    register double* result = (double*) malloc(circdepth*circdepth * sizeof(double));

    double* tmp_params = malloc(circdepth * sizeof(double ));       // temporary copy of the parameters
    for (depth_t i = 0; i < circdepth; ++i) {
        tmp_params[i] = params[i];
    }

    /*
     * 1. Calculate the gradient at the current parameter setting
     */
    double* refGradient = grad(state, tmp_params, observable, evoOps, circdepth);

    /*
     * 2. For each direction shift the corresponding parameter; calculate the gradient's finite difference. This will be
     * the i-th column of the hessian.
     */
    for (depth_t i = 0; i < circdepth; ++i) {
        tmp_params[i] += epsilon;
        double* gradient = grad(state, tmp_params, observable, evoOps, circdepth);

        for(depth_t j = 0; j < circdepth; ++j) {
            result[i*circdepth + j] = (1./epsilon) * (gradient[j] - refGradient[j]);
        }

        tmp_params[i] -= epsilon;
    }

    free(tmp_params);
    return(result);
}

/*
 * =====================================================================================================================
 *                                                  Moment matrices
 * =====================================================================================================================
 */

/*
 * This function calculates the moment matrices of a set of observables corresponding to the application of a set of
 * unitaries.
 *
 * Input:
 *  state_t* state          System's state, including statevector and Hilbert space dimension
 *  obs_t* obs[]            Set of observables
 *  depth_t obsc            Number of observables to calculate the moment matrix of
 *  obs_t* srchOps[]        Set of observables making up for the unitaries
 *  depth_t matDim          Number of operators making up for the unitaries, i.e., dimension of the resulting matrices
 *
 * Output:
 *  Array of moment matrices in column major form corresponding to the input array of observables.
 */

cplx_t** momMat(const state_t* state, \
                const obs_t* obs[], \
                depth_t obsc, \
                const obs_t* srchOps[], \
                depth_t matDim, \
                srchU applyU) {

    /*
     * 1. Set up the array of moment matrices
     */
    register cplx_t** result = (cplx_t**) malloc(obsc * sizeof(cplx_t*));

    for (depth_t i = 0; i < obsc; ++i) {
        result[i] = (cplx_t*) malloc(matDim * matDim * sizeof(cplx_t));
    }

    /*
     * 2. Initialize the temporary states altered by the unitaries
     */
    state_t tmpBra;
    stateInitEmpty(&tmpBra, state->qubits);
    state_t tmpTmpBra;
    stateInitEmpty(&tmpTmpBra, state->qubits);
    state_t tmpKet;
    stateInitEmpty(&tmpKet, state->qubits);

    /*
     * 3a. Iterate the columns of all moment matrices, applying the respective unitary to the statevector
     */
    for (depth_t i = 0; i < matDim; ++i) {
        stateCopyVector(&tmpKet, state->vector);
        applyU(&tmpKet, srchOps[i]);

        /*
         * 3b. Iterate the columns' entries of all moment matrices starting at the current column's index, applying the
         * respective unitary to the statevector
         */
        for (depth_t j = i; j < matDim; ++j) {
            stateCopyVector(&tmpKet, state->vector);
            applyU(&tmpBra, srchOps[j]);

            /*
             * 3c. Iterate the moment matrices by applying the respective observable to tmpBra's statevector and compute
             * the overlap with tmpKet
             */
            for (depth_t k = 0; k < obsc; ++k) {
                stateCopyVector(&tmpTmpBra, tmpBra.vector);
                applyObservable(&tmpTmpBra, obs[k]);
                result[k][i * matDim + j] = stateOverlap(&tmpTmpBra, &tmpKet);
                result[k][j * matDim + i] = result[k][i * matDim + j];
            }
        }
    }

    stateFree(&tmpBra);
    stateFree(&tmpTmpBra);
    stateFree(&tmpKet);
    return result;
}

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

    register cplx_t *result = (cplx_t *) malloc((18) * sizeof(cplx_t));

    /*
     * Initialize tmpB
     */
    stateInitEmpty(&tmpB, state->qubits);
    stateCopyVector(&tmpB, state->vector);
    for (qubit_t i = 0; i < tmpB.qubits; ++i) {
        applyRX(&tmpB, i, 2 * angleB);
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
