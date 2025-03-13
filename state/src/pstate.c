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
#include <vecLib/blas_new.h>
#else
#include <f77blas.h>
#endif

/*
 * =====================================================================================================================
 *                                                  State preparation
 * =====================================================================================================================
 */
/*
 * This function defines a pure quantum state of a specified number of qubits and allocates the required memory for the
 * statevector.
 *
 * @param[in,out] state		Address of the state_t struct to be initialized
 * @param[in] qubits		Number of qubits
 * @return					This function has no return value; alters the state in place
 */
void stateInitEmpty(state_t* state, const qubit_t qubits) {
	state->qubits = qubits;
	state->dim = POW2(qubits, dim_t);
	if((state->vec = calloc(state->dim, sizeof(cplx_t))) == NULL) {
        fprintf(stderr, "Statevector allocation failed\n");
    }
}

/*
 * This function initializes a pure quantum state of a specified number of qubits to the plus state (uniform
 * superposition of all computational basis states).
 *
 * @param[in,out] state		Address of the state_t struct to be initialized
 * @param[in] qubits		Number of qubits
 * @return					This function has no return value; alters the state in place
 */
void stateInitPlus(state_t* state, const qubit_t qubits) {
	stateInitEmpty(state, qubits);
	const double prefactor = 1. / sqrt((double) state->dim);
	for (dim_t i = 0; i < state->dim; ++i) {
		state->vec[i] = prefactor;
	}
}

/*
 * This function initializes a pure quantum state of a specified number of qubits to the state given by a complex vector
 *
 * @param[in,out] state		Address of the state_t struct to be initialized
 * @param[in] vector		Complex vector with DOUBLE PRECISION that is copied to the state's vector representation
 * @return					This function has no return value; alters the state in place
 */
void stateInitVector(state_t* state, const cplx_t vector[]) {
	const dim_t inc = 1;
	zcopy_(&state->dim, vector, &inc, state->vec, &inc);
}

/*
 * =====================================================================================================================
 *                                  		        Free state
 * =====================================================================================================================
 */
void stateFreeVector(state_t* state) {
	free(state->vec);
}

/*
 * =====================================================================================================================
 *                                              Binary operations
 * =====================================================================================================================
 */

cplx_t stateOverlap(const state_t state1, const state_t state2) {
	const dim_t incr = 1;
	return zdotc(&state1.dim, state1.vec, &incr, state2.vec, &incr);
}
