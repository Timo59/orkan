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

/*
 * This function computes the complex mean value of an operator wrt to the state
 * Input:
 *      state_t* state:         State of the qubit system
 *      pauliOp_t* operator:    Pauli operator; i.e., complex linear combination of Pauli strings
 * Output:
 *      The mean value of a Pauli operator without affecting the state itself.
*/
//cplx_t expValPauli(const state_t* state, const pauliOp_t* op) {
//#if defined(EXACT_BLAS)
//    return expValPauli_blas(state, op);
//#elif defined(EXACT_OMP)
//    return expValPauli_omp(state, op);
//#elif defined(EXACT_BLAS_OMP)
//    return expValPauli_blas_omp(state, op);
//#else
//    return expValPauli_old(state, op);
//#endif
//}

cplx_t expValDiag(const state_t* state, const cplx_t op[]) {
    register cplx_t result = 0;

    for (dim_t i = 0; i < state->dim; ++i) {
        result += conj(state->vec[i]) * state->vec[i] * op[i];
    }

    return result;
}

cplx_t expVal(const state_t* state, const op_t* op) {
    switch(op->type) {
        case DIAG: {
            return expValDiag(state, op->diagOp);
        }
        case PAULI: {
            return expValPauli(state, op->pauliOp);
        }
        default: {
            perror("Error in applyOp: Invalid operator type");
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
 * This function returns the expected value of a Pauli observable on a quantum state.
 *
 * Input:
 *      state_t* state:             State of a qubit system
 *      pauliObs_t* observable:     Observable as a sum of weighted Pauli strings
 *
 * Output:
 *      Expected value of the observable on the state.
 */
//double expValObsPauli(const state_t* state, const pauliObs_t* obs) {
//#if defined(EXACT_BLAS)
//    return expValObsPauli_blas(state, obs);
//#elif defined(EXACT_OMP)
//    return expValObsPauli_omp(state, obs);
//#elif defined(EXACT_BLAS_OMP)
//    return expValObsPauli_blas_omp(state, obs);
//#else
//    return expValObsPauli_old(state, obs);
//#endif
//}

/*
 * This function returns the expected value of a diagonal observable on a quantum state.
 *
 * Input:
 *      state_t* state:         State of a qubit system
 *      double observable[]:    Observable as its diagonal entries in the computational basis
 *
 * Output:
 *      Expected value of the observable on the state.
 */
double expValObsDiag(const state_t* state, const double observable[]) {
    double result = 0;      // double holding the sum of the weighted state vector's absolute squares

    /* For each of the state vector's entries: Multiply its absolute square with the corresponding observable's diagonal
     * entry and add this to the result. */
    for (dim_t i = 0; i < state->dim; ++i) {
        result += pow(cabs(state->vec[i]), 2) * observable[i];
    }
    
    return result;
}

/*
 * This function returns the expected value of an observable on a quantum state, using the function determined
 * by the type of the observable.
 *
 * Input:
 *      state_t* state:     State of a qubit system
 *      obs_t* observable:  Observable with a type union
 *
 * Output:
 *      Expected value of the observable on the state.
 */
double expValObs(const state_t* state, const obs_t* obs) {
    /* Depending on the type of the observable, calls either function returning the expected value */
    switch(obs->type) {
        case DIAG: {
            return expValObsDiag(state, obs->diagObs);
        }
        case PAULI: {
            return expValObsPauli(state, obs->pauliObs);
        }
        default: {
            perror("Error in applyOp: Invalid operator type");
            exit(1);
        }
    }
}

/*
 * =====================================================================================================================
 *                                                  Mean value PQC
 * =====================================================================================================================
 */

/*
 * This function returns the expected value of a Pauli observable on a state evolved with a PQC.
 *
 * Input:
 *      state_t* state:             Initial state of the qubit system
 *      double params[]:            Array of evolution parameters for the PQC
 *      pauliObs_t* observable:     Observable as a sum of Paulis strings
 *      obs_t* evoOps[]:            Array of hermitian evolution operators
 *      depth_t circdepth:          Number of evolution operators in the PQC
 *
 * Output:
 *      Expected value of the observable on the parametrized quantum circuit.
 */
double expValObsPqcPauli(const state_t* state,
                         const double par[],
                         const pauliObs_t* obs,
                         const obs_t *evoOps[],
                         depth_t circdepth)
{
    state_t tmp;                                    // Temporary state, initialized to the input state and evolved with
                                                    // the PQC
    stateInitEmpty(&tmp, state->qubits);            // Initialize temporary state with the input state's vector
    stateCopyVector(&tmp, state->vec);

    applyPQC(&tmp, par, evoOps, circdepth);         // Apply the PQC to the temporary state
    double result = expValObsPauli(&tmp, obs);
    stateFreeVector(&tmp);

    return result;
}

/*
 * This function returns the expected value of a diagonal observable on a state evolved with a PQC.
 *
 * Input:
 *      state_t* state:         Initial state of a qubit system
 *      double params[]:        Array of evolution parameters for the PQC
 *      double observable[]:    Observable as its diagonal entries in the computational basis
 *      obs_t* evoOps[]:        Array of hermitian evolution operators
 *      depth_t circdepth:      Number of evolution operators in the PQC
 *
 * Output:
 *      Expected value of the observable on the parametrized quantum circuit.
 */
double expValObsPqcDiag(const state_t* state,
                        const double par[],
                        const double obs[],
                        const obs_t *evoOps[],
                        depth_t circdepth)
{
    state_t tmp;                    // temporary state holding the evolved input state

    register double result = 0;     // double holding the result

    /* Initialize the temporary state and copy the input state's vector to it */
    stateInitEmpty(&tmp, state->qubits);
    stateCopyVector(&tmp, state->vec);

    /* Evolve temporary state with the PQC */
    applyPQC(&tmp, par, evoOps, circdepth);

    /* For each of the evolved state vector's entries: Multiply its absolute square with the corresponding observable's
     * diagonal entry and add this to the result. */
    for (dim_t i = 0; i < tmp.dim; ++i) {
        result += pow(cabs(tmp.vec[i]), 2) * obs[i];
    }

    /* Free the memory allocated for tmp's state vector */
    stateFreeVector(&tmp);

    return result;
}

/*
 * This function returns the expected value of an observable on a state evolved with a PQC, using the function determined
 * by the type of the observable.
 *
 * Input:
 *      state_t* state:         Initial state of a qubit system
 *      double params[]:        Array of evolution parameters for the PQC
 *      obs_t* observable:      Observable with a type union
 *      obs_t* evoOps[]:        Array of hermitian evolution operators
 *      depth_t circdepth:      Number of evolution operators in the PQC
 *
 * Output:
 *      Expected value of the observable on the parametrized quantum circuit.
 */
double expValObsPQC(const state_t* state,
                    const double par[],
                    const obs_t* obs,
                    const obs_t* evoOps[],
                    depth_t circdepth)
{
    /* Depending on the type of the observable, calls either function returning the expected value */
    switch(obs->type) {
        case DIAG: {
            return expValObsPqcDiag(state, par, obs->diagObs, evoOps, circdepth);
        }
        case PAULI: {
            return expValObsPqcPauli(state, par, obs->pauliObs, evoOps, circdepth);
        }
        default: {
            perror("Error in applyOp: Invalid operator type");
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
 * This function returns the gradient of an observable's expected value on a parametrized quantum circuit
 *
 * Input:
 *      state_t* state:         Initial state of a qubit system
 *      double par[]:           Array of evolution parameters for the PQC
 *      obs_t* obs:             Observable with a type union
 *      obs_t* evoOps[]:        Array of hermitian evolution operators
 *      depth_t circdepth:      Number of evolution operators in the PQC
 *
 * Output:
 *      Gradient of the observable's expected value on the parametrized quantum circuit wrt the evolution parameters at
 *      the current point in parameter space.
 */
//double* gradPQC(const state_t* state,
//                const double par[],
//                const obs_t* obs,
//                const obs_t *evoOps[],
//                depth_t circdepth)
//{
//#if defined(EXACT_BLAS)
//    return gradPQC_blas(state, par, obs, evoOps, circdepth);
//#elif defined(EXACT_OMP)
//    return gradPQC_omp(state, par, obs, evoOps, circdepth);
//#elif defined(EXACT_BLAS_OMP)
//    return gradPQC_blas_omp(state, par, obs, evoOps, circdepth);
//#else
//    return gradPQC_old(state, par, obs, evoOps, circdepth);
//#endif
//}

/*
 * This function returns the approximate gradient of an observable's expected value on a parametrized quantum circuit
 * using a finite difference method.
 *
 * Input:
 *      state_t* state:         Initial state of a qubit system
 *      double par[]:           Array of evolution parameters for the PQC
 *      obs_t* observable:      Observable with a type union
 *      obs_t* evoOps[]:        Array of hermitian evolution operators
 *      depth_t circdepth:      Number of evolution operators in the PQC
 *      double epsilon:         Absolute difference the parameter values are shifted in finite difference method
 *
 * Output:
 *      Gradient of the observable's expected value on the parametrized quantum circuit wrt the evolution parameters at
 *      the current point in parameter space using a finite difference method.
 */
double* approxGradientPQC(const state_t* state,
                          const double par[],
                          const obs_t* observable,
                          const obs_t* evoOps[],
                          depth_t circdepth,
                          double epsilon)
{
    double* result = (double*) malloc(circdepth * sizeof(double));

    /* Create a copy of the array holding the parameter values */
    double* tmpPar = malloc(circdepth * sizeof(double));
    for (depth_t i = 0; i < circdepth; ++i) {
        tmpPar[i] = par[i];
    }

    /* Calculate the mean value of the current parameter setting from all entries of the gradient */
    double mean = expValObsPQC(state, par, observable, evoOps, circdepth);

    /* For each direction in parameter space, shift the parameter value by epsilon and calculate the slope by the
     * difference method */
    for (depth_t i = 0; i < circdepth; ++i) {
        tmpPar[i] += epsilon;
        result[i] = (1./epsilon) * (expValObsPQC(state, tmpPar, observable, evoOps, circdepth) - mean);
        tmpPar[i] -= epsilon;
    }

    /* Free the memory allocated to the copy of the parameter array */
    free(tmpPar);

    return result;
}

/*
 * =====================================================================================================================
 *                                                      hessian PQC
 * =====================================================================================================================
 */

/*
 * This function returns the hessian of an observable's expected value on a parametrized quantum circuit
 *
 * Input:
 *      state_t* state:         Initial state of a qubit system
 *      double params[]:        Array of evolution parameters for the PQC
 *      obs_t* observable:      Observable with a type union
 *      obs_t* evoOps[]:        Array of hermitian evolution operators
 *      depth_t circdepth:      Number of evolution operators in the PQC
 *
 * Output:
 *      Hessian of the observable's expected value on the parametrized quantum circuit wrt the evolution parameters at
 *      the current point in parameter space.
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

    double* result = malloc(circdepth * circdepth * sizeof(double));

    /* Initialize all temporary states and copy the input state's vector to tmp */
    stateInitEmpty(&tmp, state->qubits);
    stateInitEmpty(&tmpBra, state->qubits);
    stateInitEmpty(&tmpKet, state->qubits);
    stateInitEmpty(&bra1, state->qubits);
    stateInitEmpty(&bra2, state->qubits);
    stateInitEmpty(&ket1, state->qubits);
    stateInitEmpty(&ket2, state->qubits);

    stateCopyVector(&tmp, state->vec);

    /* 1a) Evolve the bra vector of the second stateOverlap with all evolution operators and corresponding parameters */
    stateCopyVector(&bra2, state->vec);
    for (depth_t i = 0; i < circdepth; ++i) {
        evolveObsTrotter(&bra2, evoOps[i], params[i]);
    }

    /* 1b) The unsigned integer k iterates through all row indices of the hessian's upper triangle. In each iteration
     * tmp is evolved with the k-th evolution operator */
    for (depth_t k = 0; k < circdepth; ++k) {
        evolveObsTrotter(&tmp, evoOps[k], params[k]);

        /* 2a) tmpKet inherits its state vector from tmp and applies the k-th evolution operator to it */
        stateCopyVector(&tmpKet, tmp.vec);
        applyObs(&tmpKet, evoOps[k]);

        /* 2b) tmpBra inherits its state vector from tmp */
        stateCopyVector(&tmpBra, tmp.vec);

        /* 3a) The k-th diagonal element is computed separately because it only appears once in the matrix. bra1, ket1
         * and ket2 are all descendants of tmpKet. The k-th evolution operator is applied to ket2 a second time */
        stateCopyVector(&bra1, tmpKet.vec);
        stateCopyVector(&ket1, tmpKet.vec);
        stateCopyVector(&ket2, tmpKet.vec);
        applyObs(&ket2, evoOps[k]);

        /* 3b) bra1, ket1 and ket2 are all evolved with the remaining evolution operators */
        for (depth_t kk = k + 1; kk < circdepth; ++kk) {
            evolveObsTrotter(&bra1, evoOps[kk], params[kk]);
            evolveObsTrotter(&ket1, evoOps[kk], params[kk]);
            evolveObsTrotter(&ket2, evoOps[kk], params[kk]);
        }

        /* 3d) The observable is applied to ket1 and ket2 */
        applyObs(&ket1, observable);
        applyObs(&ket2, observable);

        /* 3e) The diagonal entries are added once to the array */
        result[k * circdepth + k] = 2 * creal(stateOverlap(&bra1, &ket1));
        result[k * circdepth + k] -= 2 * creal(stateOverlap(&bra2, &ket2));

        /* 4) The unsigned integer l iterates through all column indices of the hessian's upper triangle without the
         * diagonal, i.e., l > k. In each iteration tmpBra and tmpKet are evolved with the l-th evolution operator */
        for (depth_t l = k + 1; l < circdepth; ++l) {
            evolveObsTrotter(&tmpBra, evoOps[l], params[l]);
            evolveObsTrotter(&tmpKet, evoOps[l], params[l]);

            /* 4a) bra1 is a descendant of tmpBra, it is applied the l-th evolution operator and evolved with the
             * remaining evolution operators */
            stateCopyVector(&bra1, tmpBra.vec);
            applyObs(&bra1, evoOps[l]);
            for (depth_t ll = l + 1; ll < circdepth; ++ll) {
                evolveObsTrotter(&bra1, evoOps[ll], params[ll]);
            }

            /* 4b) ket1 is a descendant of tmpKet, is evolved with the remaining evolution operators and the observable
             * is applied to it */
            stateCopyVector(&ket1, tmpKet.vec);
            for (depth_t kk = l + 1; kk < circdepth; ++kk) {
                evolveObsTrotter(&ket1, evoOps[kk], params[kk]);
            }
            applyObs(&ket1, observable);

            /* 4c) ket2 is a descendant of tmpKet and the l-th evolution operator is applied to it. It is evolved with
             * the remaining evolution operators and the observable is applied to it */
            stateCopyVector(&ket2, tmpKet.vec);
            applyObs(&ket2, evoOps[l]);
            for (depth_t ll = l + 1; ll < circdepth; ++ll) {
                evolveObsTrotter(&ket2, evoOps[ll], params[ll]);
            }
            applyObs(&ket2, observable);

            /* 4d) The difference of the two overlap's real part are added to the positions corresponding to the l-th
             * entry of the k-th row and vice versa of the array */
            result[k * circdepth + l] = 2 * creal(stateOverlap(&bra1, &ket1));
            result[k * circdepth + l] -= 2 * creal(stateOverlap(&bra2, &ket2));
            result[l * circdepth + k] = 2 * creal(stateOverlap(&bra1, &ket1));
            result[l * circdepth + k] -= 2 * creal(stateOverlap(&bra2, &ket2));
        }
    }

    /* Free all the memory allocated for the temporary state's vectors */
    stateFreeVector(&tmp);
    stateFreeVector(&tmpBra);
    stateFreeVector(&tmpKet);
    stateFreeVector(&bra1);
    stateFreeVector(&bra2);
    stateFreeVector(&ket1);
    stateFreeVector(&ket2);

    return result;
}

/*
 * This function returns the approximate hessian of an observable's expected value on a parametrized quantum circuit
 * using a finite difference method
 *
 * Input:
 *      state_t* state:         Initial state of a qubit system
 *      double params[]:        Array of evolution parameters for the PQC
 *      obs_t* observable:      Observable with a type union
 *      obs_t* evoOps[]:        Array of hermitian evolution operators
 *      depth_t circdepth:      Number of evolution operators in the PQC
 *      gradPQC grad:           Function that computes the gradient of a PQC
 *      double epsilon:         Absolute difference the parameter values are shifted in finite difference method
 *
 * Output:
 *      Approximate hessian of the observable's expected value on the parametrized quantum circuit wrt the evolution
 *      parameters at the current point in parameter space using a finite difference method on the gradient.
 */
double* approxHessianPQC(const state_t* state, const double params[], const obs_t* observable, \
                         const obs_t* evoOps[], depth_t circdepth, gradientPQC grad, double epsilon) {

    double* result = malloc(circdepth * circdepth * sizeof(double));

    /* Create a copy of the array holding the parameter values */
    double* tmpParams = malloc(circdepth * sizeof(double));
    for (depth_t i = 0; i < circdepth; ++i) {
        tmpParams[i] = params[i];
    }


    /* Compute the gradient at the current parameter setting */
    double* refGradient = grad(state, tmpParams, observable, evoOps, circdepth);

    /* For each direction in parameter space, shift the parameter value by epsilon and calculate the gradient's finite
     * difference. This is the i-th column of the hessian */
    for (depth_t i = 0; i < circdepth; ++i) {
        tmpParams[i] += epsilon;
        double* gradient = grad(state, tmpParams, observable, evoOps, circdepth);

        for(depth_t j = 0; j < circdepth; ++j) {
            result[i * circdepth + j] = (1. / epsilon) * (gradient[j] - refGradient[j]);
        }
        free(gradient);
        tmpParams[i] -= epsilon;
    }

    /* Free the memory allocated to the copy of the parameter array */
    free(refGradient);
    free(tmpParams);
    return(result);
}

/*
 * =====================================================================================================================
 *                                                  Moment matrices
 * =====================================================================================================================
 */

/*
 * This function calculates the moment matrices of a set of observables.
 *
 * Input:
 *      state_t* state:     Initial state of a qubit system
 *      obs_t* obs[]:       Array of observables constituting the moment matrices
 *      depth_t obsc:       Number of observables to calculate the moment matrix of
 *      obs_t* srchOps[]:   Set of observables generating the search unitaries
 *      depth_t matDim:     Dimension of the moment matrices; i.e., number of operators generating the search unitaries
 *      srchU applyU:       Rule to transform search operators to search unitaries
 *
 * Output:
 *  Array of moment matrices with respect to search unitaries generated from a set of hermitian operators in column
 *  major form in order corresponding to the input array of observables.
 */

cplx_t** momMat(const state_t* state,
                const obs_t* obs[],
                depth_t obsc,
                const obs_t* srchOps[],
                depth_t matDim,
                srchU applyU)
{
    cplx_t** result = (cplx_t**) malloc(obsc * sizeof(cplx_t*));            // Array holding as many moment
    for (depth_t i = 0; i < obsc; ++i) {                                        // matrices as there are observables
        result[i] = (cplx_t*) malloc(matDim * matDim * sizeof(cplx_t));
    }

    state_t bra;
    stateInitEmpty(&bra, state->qubits);
    state_t ket;
    stateInitEmpty(&ket, state->qubits);
    state_t tmpKet;
    stateInitEmpty(&tmpKet, state->qubits);

    /* Iterate the columns of all moment matrices; copy the input state's vector to tmpKet and apply the column's search
     * unitary to it */
    for (depth_t i = 0; i < matDim; ++i) {
        stateCopyVector(&tmpKet, state->vec);
        applyU(&tmpKet, srchOps[i]);

        /* Iterate the moment matrices; copy tmpKet's vector to ket, apply the moment matrix's observable to it and
         * store the real part of the overlap with tmpKet to the moment matrix's diagonal entry */
        for (depth_t j = 0; j < obsc; ++j) {
            expValObs(&tmpKet, obs[j]);
        }

        /* Iterate the columns' entries of all moment matrices starting at the current column's index; copy the input
         * state's vector to tmpBra and apply the row's search unitary to it */
        for (depth_t j = i; j < matDim; ++j) {
            stateCopyVector(&bra, state->vec);
            applyU(&bra, srchOps[j]);

            /* Iterate the moment matrices; copy tmpBra's vector to tmpTmpBra and apply the moment matrix's observable
             * to it. */
            for (depth_t k = 0; k < obsc; ++k) {
                stateCopyVector(&ket, tmpKet.vec);
                applyObs(&ket, obs[k]);
                result[k][i * matDim + j] = stateOverlap(&bra, &ket);
                result[k][j * matDim + i] = conj(result[k][i * matDim + j]);
            }
        }
    }

    /* Free all the memory allocated for the temporary state's vectors */
    stateFreeVector(&bra);
    stateFreeVector(&ket);
    stateFreeVector(&tmpKet);
    return result;
}

/*
 * This function calculates the moment matrices with respect to the PQC's evolution operators as search unitaries.
 *
 * Input:
 *      state_t* state:         Initial state of a qubit system
 *      obs_t* obs[]:           Array of observables constituting the moment matrices
 *      depth_t obsc:           Number of observables to calculate the moment matrix of
 *      obs_t* evoOps[]:        Set of observables generating the search unitaries
 *      depth_t circdepth:      Dimension of the moment matrices; i.e., number of evolution operators in the PQC
 *      double angles[]:        Array of fixed angles in the evolution
 *
 * Output:
 *  Array of moment matrices with respect to fixed angle evolution unitaries and preceding identity operator in column
 *  major form in order corresponding to the input array of observables.
 */
cplx_t** momMatPQC(const state_t* state,
                   const obs_t* obs[],
                   depth_t obsc,
                   const obs_t* evoOps[],
                   depth_t circdepth,
                   double angles[]) {

    /* Set up the array of moment matrices */
    register cplx_t** result = (cplx_t**) malloc(obsc * sizeof(cplx_t*));
    for (depth_t i = 0; i < obsc; ++i) {
        result[i] = malloc((circdepth + 1) * (circdepth + 1) * sizeof(*(result[i])));
    }

    /* Initialize the temporary states altered by the unitaries */
    state_t bra;
    stateInitEmpty(&bra, state->qubits);
    state_t ket;
    stateInitEmpty(&ket, state->qubits);
    state_t tmpKet;
    stateInitEmpty(&tmpKet, state->qubits);

    /* The topmost diagonal entry of each moment matrix is the corresponding observable's mean value wrt to the input
     * state */
    for (depth_t i = 0; i < obsc; ++i) {
        result[i][0] = expValObs(state, obs[i]);
    }

    /* Iterate the columns of all moment matrices; copy the input state's vector to tmpKet and apply the column's search
     * unitary to it */
    for (depth_t i = 0; i < circdepth; ++i) {
        stateCopyVector(&tmpKet, state->vec);
        evolveObsTrotter(&tmpKet, evoOps[i], angles[i]);

        /* Iterate the moment matrices; copy tmpKet to bra and apply the corresponding observable to compute the first
         * column/row of each moment matrix */
        for (depth_t j = 0; j < obsc; ++j) {
            stateCopyVector(&bra, tmpKet.vec);
            applyObs(&bra, obs[j]);
            result[j][i + 1] = stateOverlap(&bra, state);
            result[j][(i + 1) * (circdepth + 1)] = conj(result[j][i + 1]);
        }

        /* Iterate the moment matrices; copy tmpKet's vector to ket, apply the moment matrix's observable to it and
         * store the real part of the overlap with tmpKet to the moment matrix's diagonal entry */
        for (depth_t j = 0; j < obsc; ++j) {
            result[j][(i + 1) * (circdepth + 1) + (i + 1)] = expValObs(&tmpKet, obs[j]);
        }

        /* Iterate the columns' entries of all moment matrices starting at the current column's index; copy the input
         * state's vector to tmpBra and apply the row's search unitary to it */
        for (depth_t j = i + 1; j < circdepth; ++j) {
            stateCopyVector(&bra, state->vec);
            evolveObsTrotter(&bra, evoOps[j], angles[j]);

            /* Iterate the moment matrices; copy tmpBra's vector to tmpTmpBra and apply the moment matrix's observable
             * to it. */
            for (depth_t k = 0; k < obsc; ++k) {
                stateCopyVector(&ket, tmpKet.vec);
                applyObs(&ket, obs[k]);
                result[k][(i + 1) * (circdepth + 1) + (j + 1)] = stateOverlap(&bra, &ket);
                result[k][(j + 1) * (circdepth + 1) + (i + 1)] = conj(result[k][(i + 1) * (circdepth + 1) + (j + 1)]);
            }
        }
    }

    /* Free all the memory allocated for the temporary state's vectors */
    stateFreeVector(&bra);
    stateFreeVector(&ket);
    stateFreeVector(&tmpKet);
    return result;
}