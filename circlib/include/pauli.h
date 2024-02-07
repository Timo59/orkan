#ifndef PAULI_H
#define PAULI_H

/*
 * =================================================================================================
 *                                              includes
 * =================================================================================================
 */

#include "gatelib.h"
#include "cplxutil.h"
#include "linalg.h"
#include <stdio.h>

/*
 * =================================================================================================
 *                                              C++ check
 * =================================================================================================
 */

#ifdef __cplusplus
extern "C" {
#endif
    
/*
 * =================================================================================================
 *                                              type definitions
 * =================================================================================================
 */

typedef struct pauliOp {
    pauli_t* components;
    cplx_t* coefficients;
    complength_t length;
    qubit_t qubits;
} pauliOp_t;

typedef struct pauliObs {
    pauli_t* components;
    double* coefficients;
    complength_t length;
    qubit_t qubits;
} pauliObs_t;

/*
 * =================================================================================================
 *                                              compound gate constructor
 * =================================================================================================
 */

// void initOperatorPauli(pauliOp_t* operator, complength_t length, qubit_t qubits);

// void freeOperatorPauli(pauliOp_t operator);

/*
 * =================================================================================================
 *                                              compound gate fusion
 * =================================================================================================
 */

// void fuseCompgates(compgate_t* result, compgate_t compgate1, compgate_t compgate2);

/*
 * =================================================================================================
 *                                              component applications
 * =================================================================================================
 */

void applyPauliString(state_t* state, const pauli_t paulistr[]);

void applyOperatorPauli(state_t* state, const pauliOp_t* op);

void applyObservablePauli(state_t* state, const pauliObs_t* observable);

void evolveWithPauliString(state_t* state, const pauli_t paulistr[], double angle);

void evolveWithTrotterizedObservablePauli(state_t* state, const pauliObs_t* observable, double angle);

#ifdef __cplusplus
}
#endif

#endif
