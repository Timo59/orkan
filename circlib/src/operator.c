/*
 * =================================================================================================
 *                                              includes
 * =================================================================================================
 */

#include "circlib.h"

/*
 * =================================================================================================
 *                                              apply operator
 * =================================================================================================
 */

void applyOperator(state_t* state, const op_t* op) {
    switch(op->type) {
        case DIAG: {
            applyOperatorDiag(state, op->diagOp);
            break;
        }
        case PAULI: {
            applyOperatorPauli(state, op->pauliOp);
            break;
        }
        default: {
            perror("Error in applyOperator: Invalid operator type");
            exit(1);
        }
    }
}

/*
 * =================================================================================================
 *                                              apply observable
 * =================================================================================================
 */

void applyObservable(state_t* state, const obs_t* observable) {
    switch(observable->type) {
        case DIAG: {
            applyObservableDiag(state, observable->diagObs);
            break;
        }
        case PAULI: {
            applyObservablePauli(state, observable->pauliObs);
            break;
        }
        default: {
            printf("Error in applyObservable: Invalid observable type.\n");
            exit(1);
        }
    }
}

/*
 * =================================================================================================
 *                                evolve with trotterized observable
 * =================================================================================================
 */

void evolveWithTrotterizedObservable(state_t* state, const obs_t* observable, double angle) {
    switch(observable->type) {
        case DIAG: {
            evolveWithObservableDiag(state, observable->diagObs, angle);
            break;
        }
        case PAULI: {
            evolveWithTrotterizedObservablePauli(state, observable->pauliObs, angle);
            break;
        }
        default: {
            printf("Error in evolveWithTrotterizedObservable: Invalid observable type.\n");
            exit(1);
        }
    }
}
