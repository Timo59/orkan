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
 * Function: Apply Parametrized Quantum Block (PQB)
 * -------------------------------------------------
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
 * =====================================================================================================================
 *                                                  Function definitions
 * =====================================================================================================================
 */
void applyPQC(state_t* state, const depth_t d, applyPQB pqbs[], double par[]);

#ifdef __cplusplus
}
#endif
#endif //PQC_H
