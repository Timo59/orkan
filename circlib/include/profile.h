//
// Created by Timo Ziegler on 17.10.24.
//
#ifndef CIRC_PROFILE_H
#define CIRC_PROFILE_H

/*
 * =====================================================================================================================
 *                                                      includes
 * =====================================================================================================================
 */
#ifndef LINALG_H
#include "linalg.h"
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
*                                                  type definitions
* =====================================================================================================================
*/
#ifndef PAULI_STRUCTS
#define PAULI_STRUCTS
typedef struct pauliOp {
    pauli_t* comps;
    cplx_t* coeffs;
    complength_t length;
    qubit_t qubits;
} pauliOp_t;

typedef struct pauliObs {
    pauli_t* comps;
    double* coeffs;
    complength_t length;
    qubit_t qubits;
} pauliObs_t;
#endif

/*
 * =====================================================================================================================
 *                                                  Apply Pauli operator
 * =====================================================================================================================
 */
void applyOpPauli_blas(state_t* state, const pauliOp_t* op);

void applyOpPauli_blas_omp(state_t* state, const pauliOp_t* op);

void applyOpPauli_omp(state_t* state, const pauliOp_t* op);

/*
 * =====================================================================================================================
 *                                                  Apply Pauli observable
 * =====================================================================================================================
 */
void applyObsPauli_blas(state_t* state, const pauliObs_t* obs);

void applyObsPauli_blas_omp(state_t* state, const pauliObs_t* obs);

void applyObsPauli_omp(state_t* state, const pauliObs_t* obs);

/*
 * =====================================================================================================================
 *                                              Evolve with Pauli string
 * =====================================================================================================================
 */
void evolvePauliStr_blas(state_t* state, const pauli_t paulistr[], double angle);

void evolvePauliStr_blas_omp(state_t* state, const pauli_t paulistr[], double angle);

void evolvePauliStr_omp(state_t* state, const pauli_t paulistr[], double angle);

#ifdef __cplusplus
extern "C" {
#endif

#endif /* CIRC_PROFILE_H */
