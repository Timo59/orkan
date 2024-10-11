//
// Created by Timo Ziegler on 10.10.24.
//

#ifndef PAULIBLMP_H
#define PAULIBLMP_H

/*
 * =====================================================================================================================
 *                                                      includes
 * =====================================================================================================================
 */

#ifndef GATELIB_H
#include "gatelib.h"
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
#ifndef PAULI_STRUCT
#define PAULI_STRUCT
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
* ======================================================================================================================
*                                                   Apply Pauli structs
* ======================================================================================================================
*/
void applyPauliStr_blas_omp(state_t* state, const pauli_t paulistr[]);

void applyOpPauli_blas_omp(state_t* state, const pauliOp_t* op);

void applyObsPauli_blas_omp(state_t* state, const pauliObs_t* obs);

/*
 * =====================================================================================================================
 *                                              Evolve with Pauli string
 * =====================================================================================================================
 */

void evolvePauliStr_blas_omp(state_t* state, const pauli_t paulistr[], double angle);

#ifdef __cplusplus
extern "C" {
#endif

#endif /* PAULIBLMP_H */
