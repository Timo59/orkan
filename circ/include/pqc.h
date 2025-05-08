//
// Created by Timo Ziegler on 12.03.25.
//

#ifndef PQC_H
#define PQC_H
/*
 * =====================================================================================================================
 *                                                      Includes
 * =====================================================================================================================
 */
#ifndef PSTATE_H
#include "pstate.h"
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
 *                                                  Type definitions
 * =====================================================================================================================
 */

/*
 * Function: Apply Quantum Block (QB)
 * ------------------------------------
 * This function applies a fixed sequence of gates to a quantum state
 */
typedef void (*applyQB)(state_t* state);

/*
 * Function: Apply Parametrized Quantum Block (PQB)
 * -------------------------------------------------
 * This function applies a fixed sequence of (parametrized) gates to a quantum state
 */
typedef void (*applyPQB)(state_t* state, double par);

/*
 * Function: Apply Linear Combination of hermitian operators
 * ----------------------------------------------------------
 * This function applies a real valued linear combination of (parametrized) gates to a quantum state
 */
typedef void (*applyHerm)(state_t* state, double c[]);

/*
 * =====================================================================================================================
 *                                                  Direct application
 * =====================================================================================================================
 */

/*
 * This function applies the diagonal hermitian operator to the quantum state
 *
 * @param[in,out] state : Quantum state; On exit,
 */
void applyDiag(state_t* state, const double diag[]);
void lcQB(state_t* state, depth_t d, const applyQB qb[], const cplx_t c[]);

/*
 * =====================================================================================================================
 *                                              Imaginary time evolution
 * =====================================================================================================================
 */
void evoDiag(state_t* state, const double diag[], double par);

void evoQB(state_t* state, applyQB qb, double par);

void evoHerm(state_t* state, applyHerm herm, double par);

#define ite(X, Y, Z) _Generic((Y), \
    double*:    evoDiag,           \
    applyQB:    evoQB              \
    applyHerm:  evoHerm            \
    ) (X, Y, Z)

void applyPQC(state_t* state, depth_t d, const applyPQB pqbs[], const double par[]);

#ifdef __cplusplus
}
#endif
#endif //PQC_H
