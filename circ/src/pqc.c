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

#ifndef __OMP_H
#include <omp.h>
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
 *                                                  Direct application
 * =====================================================================================================================
 */

void applyDiag(state_t* state, const double diag[]) {
#pragma omp parallel for default(none) shared(state, diag)
    for (dim_t i = 0; i < state->dim; ++i) {
        state->vec[i] *= diag[i];
    }
}

void applyHerm(state_t* )

void lcQB(state_t* state, const depth_t d, const applyQB qb[], const cplx_t c[]) {
    const dim_t incr = 1;
    cplx_t* out = calloc(state->dim, sizeof (cplx_t));
    if(!out) {
        fprintf(stderr, "lcQB: out allocation failed\n");
        return;
    }
    for (depth_t i = 0; i < d; ++i) {
        cplx_t* tmp = malloc(state->dim * sizeof (cplx_t));
        if(!tmp) {
            fprintf(stderr, "lcQB: tmp allocation failed\n");
            free(out);
            return;
        }
        zcopy_(&state->dim, state->vec, &incr, tmp, &incr);
        qb[i](state);
        zaxpy_(&state->dim, c + i, state->vec, &incr, out, &incr);
        stateFreeVector(state);
        state->vec = tmp;
    }
    stateFreeVector(state);
    state->vec = out;
}

/*
 * =====================================================================================================================
 *                                              Unitary time evolution
 * =====================================================================================================================
 */
void evoDiag(state_t* state, const double diag[], const double par) {
#pragma omp parallel for default(none) shared(state, diag, par)
    for (dim_t i = 0; i < state->dim; ++i) {
        state->vec[i] *= cexp(-I * par * diag[i]);
    }
}

void evoQB(state_t* state, const applyQB qb, const double par) {
    const dim_t incr = 1;
    const cplx_t c = cos(par);
    const cplx_t s = - I * sin(par);
    cplx_t* tmp = malloc(state->dim * sizeof (cplx_t));
    if(!tmp) {
        fprintf(stderr, "evoPQB: tmp allocation failed\n");
        return;
    }

    zcopy_(&state->dim, state->vec, &incr, tmp, &incr);
    qb(state);
    zaxpby_(&state->dim, &c, tmp, &incr, &s, state->vec, &incr);    // cos(par)*state - i*sin(par)*qb(state)
    free(tmp);
}

void evoPQB(state_t* state, const applyPQB pqb, const double par) {
    pqb(state, par);
}

void applyPQC(state_t* state, const pqc_t pqc, const double par[]) {
    for (depth_t i = 0; i < pqc.len; ++i) {
        pqc.pqb[i](state, *(par + pqc.parIdx[i]));
    }
}