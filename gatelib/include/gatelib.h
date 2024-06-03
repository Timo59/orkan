#ifndef GATELIB_H
#define GATELIB_H

/*
 * =================================================================================================
 *                                              includes
 * =================================================================================================
 */

#include "owntypes.h"
#include "statelib.h"
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
 *                                              Pauli gates
 * =================================================================================================
 */

void applyX(state_t* state, qubit_t qubit);

void applyY(state_t* state, qubit_t qubit);

void applyZ(state_t* state, qubit_t qubit);

/*
 * =================================================================================================
 *                                              Clifford gates
 * =================================================================================================
 */

void applyS(state_t* state, qubit_t qubit);

void applySdagger(state_t* state, qubit_t qubit);

void applyH(state_t* state, qubit_t qubit);

/*
 * =================================================================================================
 *                                              Hadamard-Y gate
 * =================================================================================================
 */

void applyHy(state_t* state, qubit_t qubit);

/*
 * =================================================================================================
 *                                              T gate
 * =================================================================================================
 */

void applyT(state_t* state, qubit_t qubit);

void applyTdagger(state_t* state, qubit_t qubit);

/*
 * =================================================================================================
 *                                              P gate
 * =================================================================================================
 */

void applyP(state_t* state, qubit_t qubit, double angle);

void applyPdagger(state_t* state, qubit_t qubit, double angle);

/*
 * =================================================================================================
 *                                              rotation gates
 * =================================================================================================
 */

void applyRX(state_t* state, qubit_t qubit, double angle);

void applyRXdagger(state_t* state, qubit_t qubit, double angle);

void applyRY(state_t* state, qubit_t qubit, double angle);

void applyRYdagger(state_t* state, qubit_t qubit, double angle);

void applyRZ(state_t* state, qubit_t qubit, double angle);

void applyRZdagger(state_t* state, qubit_t qubit, double angle);

/*
 * =================================================================================================
 *                                              controlled Pauli gates
 * =================================================================================================
 */

void applyCX(state_t* state, qubit_t control, qubit_t target);

void applyCY(state_t* state, qubit_t control, qubit_t target);

void applyCZ(state_t* state, qubit_t control, qubit_t target);

/*
 * =================================================================================================
 *                                              controlled Clifford gates
 * =================================================================================================
 */

void applyCS(state_t* state, qubit_t control, qubit_t target);

void applyCSdagger(state_t* state, qubit_t control, qubit_t target);

void applyCH(state_t* state, qubit_t control, qubit_t target);

/*
 * =================================================================================================
 *                                              controlled Hadamard-Y gate
 * =================================================================================================
 */

void applyCHy(state_t* state, qubit_t control, qubit_t target);

/*
 * =================================================================================================
 *                                              controlled T gate
 * =================================================================================================
 */

void applyCT(state_t* state, qubit_t control, qubit_t target);

void applyCTdagger(state_t* state, qubit_t control, qubit_t target);

/*
 * =================================================================================================
 *                                              controlled P gate
 * =================================================================================================
 */

void applyCP(state_t* state, qubit_t control, qubit_t target, double angle);

void applyCPdagger(state_t* state, qubit_t control, qubit_t target, double angle);

/*
 * =================================================================================================
 *                                              SWAP gate
 * =================================================================================================
 */

void applySWAP(state_t* state, qubit_t qubit1, qubit_t qubit2);

/*
 * =================================================================================================
 *                                              Toffoli gate
 * =================================================================================================
 */

void applyToffoli(state_t* state, qubit_t control1, qubit_t control2, qubit_t target);

#ifdef __cplusplus
}
#endif

#endif
