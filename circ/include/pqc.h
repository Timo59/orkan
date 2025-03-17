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
 * Function: Apply Linear Combination of Quantum Blocks (LCQB)
 * -----------------------------------------------------------
 * This function applies a linear combination of (parametrized) gates to a quantum state
 */
typedef void (*applyLCQB)(state_t* state, double c[]);

/*
 * =====================================================================================================================
 *                                                  Function definitions
 * =====================================================================================================================
 */
void evoQB(state_t* state, applyQB qb, double par);
void lcQB(state_t* state, const depth_t d, applyQB qb[], const cplx_t c[]);
void applyPQC(state_t* state, depth_t d, const applyPQB pqbs[], const double par[]);

#ifdef __cplusplus
}
#endif
#endif //PQC_H
