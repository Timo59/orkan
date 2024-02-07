/*
 * =================================================================================================
 *                                              includes
 * =================================================================================================
 */

#include "statelib.h"
#include <stdio.h>

/*
 * =================================================================================================
 *                                              state construction
 * =================================================================================================
 */

void stateInitEmpty(state_t* state, qubit_t qubits) {
	state->qubits = qubits;
	state->dimension = POW2(qubits, dim_t);
	state->vector = (cplx_t*) calloc(state->dimension, sizeof(cplx_t));
}

void stateInitZero(state_t* state, qubit_t qubits) {
	stateInitEmpty(state, qubits);
	state->vector[0] = 1.0; 
}

void stateInitPlus(state_t* state, qubit_t qubits) {
	stateInitEmpty(state, qubits);
	double prefactor = 1. / sqrt(state->dimension);
	for (dim_t i = 0; i < state->dimension; ++i) {
		state->vector[i] = prefactor;
	}
}

void stateInitVector(state_t* state, cplx_t vector[], qubit_t qubits) {
	state->qubits = qubits;
	state->dimension = POW2(qubits, dim_t);
	state->vector = vector;
}

void stateCopyVector(state_t* state, const cplx_t vector[]) {
	for (dim_t i = 0; i < state->dimension; ++i) {
		state->vector[i] = vector[i];
	}
}

/*
 * =================================================================================================
 *                                              free state
 * =================================================================================================
 */

void stateFreeVector(state_t* state) {
	free(state->vector);
}

void stateFree(state_t* state) {
	stateFreeVector(state);
	free(state);
}

/*
 * =================================================================================================
 *                                              binary operations
 * =================================================================================================
 */

cplx_t stateOverlap(const state_t* state1, const state_t* state2) {
	cplx_t result = 0;
	for (dim_t i = 0; i < state1->dimension; ++i) {
		result += conj(state1->vector[i]) * state2->vector[i];
	}
	return result;
}
