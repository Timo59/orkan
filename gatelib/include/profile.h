//
// Created by Timo Ziegler on 17.10.24.
//

#ifndef GATE_PROFILE_H
#define GATE_PROFILE_H

/*
 * =====================================================================================================================
 *                                                      Includes
 * =====================================================================================================================
 */
#ifndef __DISPATCH_PUBLIC__
#include <dispatch/dispatch.h>
#endif

#ifndef __OMP_H
#include <omp.h>
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
*                                                  Pauli X-gate
* =====================================================================================================================
*/
void applyX_blas(state_t* state, qubit_t qubit);

void applyX_blas_omp(state_t* state, qubit_t qubit);

void applyX_gcd(state_t* state, qubit_t qubit);

void applyX_omp(state_t* state, qubit_t qubit);

/*
* =====================================================================================================================
*                                                  Pauli Z-gate
* =====================================================================================================================
*/
void applyY_blas(state_t* state, qubit_t qubit);

void applyY_blas_omp(state_t* state, qubit_t qubit);

void applyY_gcd(state_t* state, qubit_t qubit);

void applyY_omp(state_t* state, qubit_t qubit);

/*
* =====================================================================================================================
*                                                  Pauli X-gate
* =====================================================================================================================
*/
void applyZ_blas(state_t* state, qubit_t qubit);

void applyZ_blas_omp(state_t* state, qubit_t qubit);

void applyZ_gcd(state_t* state, qubit_t qubit);

void applyZ_omp(state_t* state, qubit_t qubit);

#ifdef __cplusplus
extern "C" {
#endif

#endif /* GATE_PROFILE_H */
