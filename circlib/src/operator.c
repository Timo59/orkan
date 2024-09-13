/*
 * =====================================================================================================================
 *                                                      includes
 * =====================================================================================================================
 */

#include "circlib.h"

/*
 * =====================================================================================================================
 *                                                  apply operator
 * =====================================================================================================================
 */

void applyOperator(state_t* state, const op_t* op) {
    switch(op->type) {
        case DIAG: {
            applyOpDiag(state, op->diagOp);
            break;
        }
        case PAULI: {
            applyOpPauli(state, op->pauliOp);
            break;
        }
        default: {
            perror("Error in applyOperator: Invalid operator type");
            exit(1);
        }
    }
}

/*
 * =====================================================================================================================
 *                                                  apply observable
 * =====================================================================================================================
 */

void applyObservable(state_t* state, const obs_t* observable) {
    switch(observable->type) {
        case DIAG: {
            applyObsDiag(state, observable->diagObs);
            break;
        }
        case PAULI: {
            applyObsPauli(state, observable->pauliObs);
            break;
        }
        default: {
            printf("Error in applyObservable: Invalid observable type.\n");
            exit(1);
        }
    }
}

/*
 * =====================================================================================================================
 *                                      evolve with trotterized observable
 * =====================================================================================================================
 */

void evolveWithTrotterizedObservable(state_t* state, const obs_t* observable, double angle) {
    switch(observable->type) {
        case DIAG: {
            evolveObsDiag(state, observable->diagObs, angle);
            break;
        }
        case PAULI: {
            evolveObsPauliTrotter(state, observable->pauliObs, angle);
            break;
        }
        default: {
            printf("Error in evolveWithTrotterizedObservable: Invalid observable type.\n");
            exit(1);
        }
    }
}

/*
 * =====================================================================================================================
 *                                                      apply PQC
 * =====================================================================================================================
 */

void applyPQC(state_t* state, const double par[], const obs_t* evoOps[], depth_t circdepth) {
    for (depth_t i = 0; i < circdepth; ++i) {
        evolveWithTrotterizedObservable(state, evoOps[i], par[i]);
    }
}

/*
 * =====================================================================================================================
 *                                      Linear combination of parametrized quantum gates
 * =====================================================================================================================
 */

/*
 * This function performs a linear combination of parametrized quantum gates on a quantum state returning the quantum
 * state conditioned on the outcome of an ancillary measurement (normalized).
 *
 * Input:
 *      state_t* state:         Input state of a quantum system
 *      cplx_t coeff[]:         Coefficients of the evolution operators in the linear combination
 *      double par[]:           Angles of the evolution operators
 *      obs_t* evoOps[]:        Array of observables for the evolution operators
 *      depth_t parc:           Number of evolution operators
 *
 * Output:
 *      The function has no output value, instead it performs the linear combination of unitaries on the state vector of
 *      the quantum state.
 */

void lcupqg(state_t* state, const cplx_t coeff[], const double angles[], const obs_t* evoOps[], depth_t anglesc) {
    cplx_t* stateVec = (cplx_t*) calloc(state->dim, sizeof(cplx_t));
    state_t tmp;
    stateInitEmpty(&tmp, state->qubits);

    cblas_zaxpy((__LAPACK_int) state->dim,
                (__LAPACK_double_complex*) coeff,
                (__LAPACK_double_complex*) state->vec,
                (__LAPACK_int) 1,
                (__LAPACK_double_complex*) stateVec,
                (__LAPACK_int) 1);

    /* For each evolution operator, evolve the input state and store the resulting state vector multiplied with the
     * coefficient of the LCU to stateVec */
    for (depth_t i = 0; i < anglesc; ++i) {
        stateCopyVector(&tmp, state->vec);
        evolveWithTrotterizedObservable(&tmp, evoOps[i], angles[i]);

        cblas_zaxpy((__LAPACK_int) state->dim,
                    (__LAPACK_double_complex*) coeff + (i + 1),
                    (__LAPACK_double_complex*) tmp.vec,
                    (__LAPACK_int) 1,
                    (__LAPACK_double_complex*) stateVec,
                    (__LAPACK_int) 1);
    }
    stateFreeVector(&tmp);

    state->vec = stateVec;
}