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

#ifdef MACOS
#include <vecLib/cblas_new.h>
#else
#include <cblas.h>
#endif

/*
 * =====================================================================================================================
 *                                                  State preparation
 * =====================================================================================================================
 */
/*
 * @brief	Defines a pure quantum state of a specified number of qubits and allocates the required memory for the
 *			statevector
 *
 * @param[in,out]	state	Address of the state
 * @param[in]		qubits	Number of qubits
 *
 * @returns	This function has no return value; alters the state in place
 */
void stateInitEmpty(state_t* state, const qubit_t qubits) {
	state->qubits = qubits;
	state->dim = POW2(qubits, dim_t);
	if((state->vec = calloc(state->dim, sizeof(cplx_t))) == NULL) {
		fprintf(stderr, "stateInitEmpty(): state->vec allocation failed\n");
        exit(EXIT_FAILURE);
    }
}

/*
 * @brief	Initialize a pure quantum state to the plus state (uniform superposition of all computational basis states)
 *
 * @param[in,out]	state		Address of the state
 * @param[in]		qubits		Number of qubits
 *
 * @returns	This function has no return value; alters the state in place
 */
void stateInitPlus(state_t* state, const qubit_t qubits) {
	stateInitEmpty(state, qubits);
	const double prefactor = 1. / sqrt((double) state->dim);
	for (dim_t i = 0; i < state->dim; ++i) {
		state->vec[i] = prefactor;
	}
}

/*
 * @brief	Initialize a pure quantum state to the state given by a complex vector
 *
 * @param[in,out]	state	Address of the state
 * @param[in]		vector	Complex vector with DOUBLE PRECISION that is copied to the state's vector representation
 *
 * @returns	This function has no return value; alters the state in place
 */
void stateInitVector(state_t* state, const cplx_t vector[]) {
	const dim_t inc = 1;
	cblas_zcopy(state->dim, vector, inc, state->vec, inc);
}

/*
 * =====================================================================================================================
 *                                  		        Free state
 * =====================================================================================================================
 */
void stateFreeVector(state_t* state) {
	free(state->vec);
	state->vec = NULL;
}

/*
 * =====================================================================================================================
 *                                              Binary operations
 * =====================================================================================================================
 */

/*
 * @brief	Calculate the overlap of two pure quantum states, i.e., <psi|phi>
 *
 * param[in,out]	state1	bra
 * param[in,out]	state2	ket
 */
cplx_t stateOverlap(const state_t state1, const state_t state2) {
	cplx_t out;
	cblas_zdotc_sub(state1.dim, state1.vec, 1, state2.vec, 1, &out);
	return out;
}
