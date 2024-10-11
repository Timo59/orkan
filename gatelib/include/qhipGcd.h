//
// Created by Timo Ziegler on 09.10.24.
//

#ifndef QHIPGCD_H
#define QHIPGCD_H

/*
 * =====================================================================================================================
 *                                                      Includes
 * =====================================================================================================================
 */

#ifndef __DISPATCH_PUBLIC__
#include <dispatch/dispatch.h>
#endif

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

void applyX_gcd(state_t* state, qubit_t qubit);

void applyY_gcd(state_t* state, qubit_t qubit);

void applyZ_gcd(state_t* state, qubit_t qubit);

#ifdef __cplusplus
extern "C" {
#endif

#endif /* QHIPGCD_H */
