//
// Created by Timo Ziegler on 12.03.25.
//
/*
 * =====================================================================================================================
 *                                                      Includes
 * =====================================================================================================================
 */
#ifndef __MATH__
#include <math.h>
#endif

#ifndef PQC_H
#include "pqc.h"
#endif

#ifndef _STDIO_H_
#include <stdio.h>
#endif

#ifndef _STDLIB_H_
#include <stdlib.h>
#endif

#ifdef MACOS
#include <vecLib/blas_new.h>
#include <vecLib/lapack.h>
#else
#include <openblas-pthread/f77blas.h>
#endif

/*
 * =====================================================================================================================
 *                                                  Function definitions
 * =====================================================================================================================
 */
/*
 * This function applies an imaginary time evolution of the quantum state by the specified quantum block (sequence of
 * gates) for the time specified by par
 *
 * @param[in,out]   state   Input quantum state to the PQB. On exit, imaginary time evolution of the quantum state
 * @param[in]       qb      Quantum block; i.e, a function applying quantum gates to the quantum state
 * @param[in]       par     Parameter; i.e., time
 *
 * WARNING: This function only gives appropriate results if the QB is involutory; i.e., QB*QB = Id
 */
void evolveQB(state_t* state, applyQB qb, double par) {
    const dim_t incr = 1;
    const cplx_t c = cos(par);
    const cplx_t s = - I * sin(par);
    cplx_t* tmp = malloc(state->dim * sizeof (cplx_t));
    if(!tmp) {
        fprintf(stderr, "applyPQB: tmp allocation failed\n");
        return;
    }

    zcopy_(&state->dim, state->vec, &incr, tmp, &incr);
    qb(state);
    zaxpby_(&state->dim, &c, tmp, &incr, &s, state->vec, &incr);    // cos(par)*state - i*sin(par)*qb(state)
}

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