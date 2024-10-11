//
// Created by Timo Ziegler on 09.10.24.
//

#ifndef QHIPSTER_H
#define QHIPSTER_H

/*
 * =====================================================================================================================
 *                                                      Includes
 * =====================================================================================================================
 */

#ifndef STATELIB_H
#include "statelib.h"
#endif

#ifndef _STDIO_H_
#include <stdio.h>
#endif


/*
 * =====================================================================================================================
 *                                                      C++ check
 * =====================================================================================================================
 */

#ifdef __cplusplus
extern "C" {
#endif

/*
 * =====================================================================================================================
 *                                                  Pauli gates
 * =====================================================================================================================
 */

void applyX_old(state_t* state, qubit_t qubit);

void applyY_old(state_t* state, qubit_t qubit);

void applyZ_old(state_t* state, qubit_t qubit);

/*
 * =====================================================================================================================
 *                                                  Clifford gates
 * =====================================================================================================================
 */

void applyS_old(state_t* state, qubit_t qubit);

void applySdagger_old(state_t* state, qubit_t qubit);

void applyH_old(state_t* state, qubit_t qubit);

/*
 * =====================================================================================================================
 *                                                  Hadamard-Y gate
 * =====================================================================================================================
 */

void applyHy_old(state_t* state, qubit_t qubit);

/*
 * =====================================================================================================================
 *                                                      T gate
 * =====================================================================================================================
 */

void applyT_old(state_t* state, qubit_t qubit);

void applyTdagger_old(state_t* state, qubit_t qubit);

/*
 * =====================================================================================================================
 *                                                      P gate
 * =====================================================================================================================
 */

void applyP_old(state_t* state, qubit_t qubit, double angle);

void applyPdagger_old(state_t* state, qubit_t qubit, double angle);

/*
 * =====================================================================================================================
 *                                                  Rotation gates
 * =====================================================================================================================
 */

void applyRX_old(state_t* state, qubit_t qubit, double angle);

void applyRXdagger_old(state_t* state, qubit_t qubit, double angle);

void applyRY_old(state_t* state, qubit_t qubit, double angle);

void applyRYdagger_old(state_t* state, qubit_t qubit, double angle);

void applyRZ_old(state_t* state, qubit_t qubit, double angle);

void applyRZdagger_old(state_t* state, qubit_t qubit, double angle);

/*
 * =====================================================================================================================
 *                                              controlled Pauli gates
 * =====================================================================================================================
 */

void applyCX_old(state_t* state, qubit_t control, qubit_t target);

void applyCY_old(state_t* state, qubit_t control, qubit_t target);

void applyCZ_old(state_t* state, qubit_t control, qubit_t target);

/*
 * =====================================================================================================================
 *                                              Controlled Clifford gates
 * =====================================================================================================================
 */

void applyCS_old(state_t* state, qubit_t control, qubit_t target);

void applyCSdagger_old(state_t* state, qubit_t control, qubit_t target);

void applyCH_old(state_t* state, qubit_t control, qubit_t target);

/*
 * =====================================================================================================================
 *                                              Controlled Hadamard-Y gate
 * =====================================================================================================================
 */

void applyCHy_old(state_t* state, qubit_t control, qubit_t target);

/*
 * =====================================================================================================================
 *                                                  Controlled T gate
 * =====================================================================================================================
 */

void applyCT_old(state_t* state, qubit_t control, qubit_t target);

void applyCTdagger_old(state_t* state, qubit_t control, qubit_t target);

/*
 * =====================================================================================================================
 *                                                  Controlled P gate
 * =====================================================================================================================
 */

void applyCP_old(state_t* state, qubit_t control, qubit_t target, double angle);

void applyCPdagger_old(state_t* state, qubit_t control, qubit_t target, double angle);

/*
 * =====================================================================================================================
 *                                                      SWAP gate
 * =====================================================================================================================
 */

void applySWAP_old(state_t* state, qubit_t qubit1, qubit_t qubit2);

/*
 * =====================================================================================================================
 *                                              Toffoli gate
 * =====================================================================================================================
 */

void applyToffoli_old(state_t* state, qubit_t control1, qubit_t control2, qubit_t target);




#ifdef __cplusplus
}
#endif

#endif /* QHIPSTER_H */