/*
 * =====================================================================================================================
 *                                                      includes
 * =====================================================================================================================
 */
#ifndef __MATH__
#include <math.h>
#endif

#ifndef PSTATE_H
#include "pstate.h"
#endif

#ifndef _STDIO_H_
#include <stdio.h>
#endif

#ifndef _STDLIB_H_
#include <stdlib.h>
#endif
/*
 * =====================================================================================================================
 *                                                  State preparation
 * =====================================================================================================================
 */
void stateInitEmpty(state_t* state, qubit_t qubits) {
	state->qubits = qubits;
	state->dim = POW2(qubits, dim_t);
	if((state->vec = calloc(state->dim, sizeof(cplx_t))) == NULL) {
        fprintf(stderr, "Statevector allocation failed\n");
        return;
    }
}

void stateInitZero(state_t* state, qubit_t qubits) {
	stateInitEmpty(state, qubits);
	state->vec[0] = 1.0;
}

void stateInitPlus(state_t* state, qubit_t qubits) {
	stateInitEmpty(state, qubits);
	double prefactor = 1. / sqrt((double) state->dim);
	for (dim_t i = 0; i < state->dim; ++i) {
		state->vec[i] = prefactor;
	}
}

void stateInitVector(state_t* state, cplx_t vector[], qubit_t qubits) {
	state->qubits = qubits;
	state->dim = POW2(qubits, dim_t);
	state->vec = vector;
}

void stateCopyVector(state_t* state, const cplx_t vector[]) {
	for (dim_t i = 0; i < state->dim; ++i) {
		state->vec[i] = vector[i];
	}
}

/*
 * =====================================================================================================================
 *                                  		        Free state
 * =====================================================================================================================
 */

void stateFreeVector(state_t* state) {
	free(state->vec);
}

void stateFree(state_t* state) {
	stateFreeVector(state);
	free(state);
}

/*
 * =====================================================================================================================
 *                                              Binary operations
 * =====================================================================================================================
 */

cplx_t stateOverlap(const state_t* state1, const state_t* state2) {
	cplx_t result = 0;
	for (dim_t i = 0; i < state1->dim; ++i) {
		result += conj(state1->vec[i]) * state2->vec[i];
	}
	return result;
}
