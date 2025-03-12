//
// Created by Timo Ziegler on 12.03.25.
//
/*
 * =====================================================================================================================
 *                                                      Includes
 * =====================================================================================================================
 */
#ifndef PQC_H
#include "pqc.h"
#endif

/*
 * =====================================================================================================================
 *                                                  Function definitions
 * =====================================================================================================================
 */
/*
 * This function applies a Parametrized Quantum Circuit (PQC), which is a sequence of parametrized blocks, to a quantum
 * state
 *
 * @param[in,out]   state   Input quatum state to the PQC. On exit, output quantum state of the PQC
 * @param[in]       d       Number of parametrized quantum blocks
 * @param[in]       pqbs[]  Array of functions; resembling the parametrized quantum blocks
 * @param[in]       par[]   Parameter setting; par[i] is parameter for pqbs[i]
 */
void applyPQC(state_t* state, const depth_t d, const applyPQB pqbs[], const double par[]) {
    for (depth_t i = 0; i < d; i++) {
        pqbs[i](state, par[i]);
    }
}