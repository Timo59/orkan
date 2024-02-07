#ifndef STATELIB_H
#define STATELIB_H

/*
 * =================================================================================================
 *                                              includes
 * =================================================================================================
 */

#include "owntypes.h"
#include <complex.h>
#include <math.h>
#include <stdlib.h>

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

/*
 * Struct:      state_t
 * --------------------
 * Description: This struct represents a multi-qubit pure state.
 * Contents:    
 *      vector:             Vector representation of the state.
 *      qubits:             Number of qubits of the underlying system.
 *      length:             Length of the state's vector representation.
 */
typedef struct state {
    cplx_t* vector;
    qubit_t qubits;
    dim_t dimension;
} state_t;

/*
 * =================================================================================================
 *                                              state construction
 * =================================================================================================
 */

void stateInitEmpty(state_t* state, qubit_t qubits);

void stateInitZero(state_t* state, qubit_t qubits);

void stateInitPlus(state_t* state, qubit_t qubits);

void stateInitVector(state_t* state, cplx_t vector[], qubit_t qubits);

void stateCopyVector(state_t* state, const cplx_t vector[]);

/*
 * =================================================================================================
 *                                  		free state
 * =================================================================================================
 */

void stateFreeVector(state_t* state);

void stateFree(state_t* state);

/*
 * =================================================================================================
 *                                          binary operations
 * =================================================================================================
 */

cplx_t stateOverlap(const state_t* state1, const state_t* state2);

#ifdef __cplusplus
}
#endif

#endif
