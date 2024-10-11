//
// Created by Timo Ziegler on 09.10.24.
//

#ifndef QHIPBLAS_H
#define QHIPBLAS_H

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

#ifndef __VECLIB__
#include <vecLib/vecLib.h>
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

void applyX_blas(state_t* state, qubit_t qubit);

void applyY_blas(state_t* state, qubit_t qubit);

void applyZ_blas(state_t* state, qubit_t qubit);

#ifdef __cplusplus
extern "C" {
#endif

#endif //QHIPBLAS_H
