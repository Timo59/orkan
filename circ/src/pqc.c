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
#include <vecLib/cblas_new.h>
#include <vecLib/lapack.h>
#else
#include <cblas.h>
#endif

/*
 * =====================================================================================================================
 *                                                  Direct application
 * =====================================================================================================================
 */

void applyDiag(const state_t* state, const double diag[]) {
#pragma omp parallel for default(none) shared(state, diag)
    for (dim_t i = 0; i < state->dim; ++i) {
        state->vec[i] *= diag[i];
    }
}

void applyCGwrapper(state_t* state, const applyCG cg) {
    cg(state);
}

void applyHerm(state_t* state, const herm_t* herm) {
    cplx_t* out = calloc(state->dim, sizeof(cplx_t));
    if (out == NULL) {
        fprintf(stderr, "applyHerm(): out allocation failed\n");
        exit(EXIT_FAILURE);
    }

    // Creates a team of threads each applying one of the hermitian's composite gates to a copy of the input state
#pragma omp parallel for default(none) shared(state, herm, stderr, out)
    for (unsigned int i = 0; i < herm->len; ++i) {
        state_t tmp;
        tmp.qubits = state->qubits;
        tmp.dim = state->dim;
        if ((tmp.vec = malloc(tmp.dim * sizeof(cplx_t))) == NULL) {
            fprintf(stderr, "applyHerm(): tmp allocation failed\n");
            free(out);
            exit(EXIT_FAILURE);
        }
        stateInitVector(&tmp, state->vec);
        herm->comp[i](&tmp);
        cplx_t ALPHA = herm->weight[i];            // zaxpy_ takes in double complex variable

        // To avoid a race condition the resulting state vector is added to the output one by one
#pragma omp critical(sum)
        {
            cblas_zaxpy(state->dim, &ALPHA, tmp.vec, 1, out, 1);
        }
        stateFreeVector(&tmp);
    }
    stateFreeVector(state);
    state->vec = out;
}

void applyLCCG(state_t* state, const lccg_t* lccg) {

    cplx_t* out = calloc(state->dim, sizeof(cplx_t));
    if (out == NULL) {
        fprintf(stderr, "applyLQCB(): out allocation failed\n");
        exit(EXIT_FAILURE);
    }

    // Creates a team of threads each applying one of the linear combination's composite gates to a copy of the input
    // state
#pragma omp parallel for default(none) shared(state, lccg, stderr, out)
    for (depth_t i = 0; i < lccg->len; ++i) {
        state_t tmp;
        tmp.qubits = state->qubits;
        tmp.dim = state->dim;
        if ((tmp.vec = malloc(tmp.dim * sizeof(cplx_t))) == NULL) {
            fprintf(stderr, "applyLCCG(): tmp allocation failed\n");
            free(out);
            exit(EXIT_FAILURE);
        }
        stateInitVector(&tmp, state->vec);
        lccg->comp[i](&tmp);

        // To avoid a race condition the resulting state vector is added to the output one by one
#pragma omp critical(sum)
        {
            cblas_zaxpy(state->dim, lccg->weight + i, tmp.vec, 1, out, 1);
        }
        stateFreeVector(&tmp);
    }
    stateFreeVector(state);
    state->vec = out;
}

/*
 * =====================================================================================================================
 *                                              Unitary time evolution
 * =====================================================================================================================
 */
void evoDiag(const state_t* state, const double diag[], double par) {
    par /= 2.;
#pragma omp parallel for default(none) shared(state, diag, par)
    for (dim_t i = 0; i < state->dim; ++i) {
        state->vec[i] *= cexp(-I * par * diag[i]);
    }
}

void evoCG(state_t* state, const applyCG cg, double par) {
    par /= 2.;
    const cplx_t c = cos(par);
    const cplx_t s = - I * sin(par);
    // Initialize temporary state the composite quantum gate is applied to
    state_t tmp;
    stateInitEmpty(&tmp, state->qubits);
    stateInitVector(&tmp, state->vec);

    // Scale the input state with the cosine, apply the composite quantum gate to the temporary state and add the latter
    // times the negative sine to the former
    cblas_zscal(state->dim, &c, state->vec, 1);
    cg(&tmp);
    cblas_zaxpy(state->dim, &s, tmp.vec, 1, state->vec, 1);
    stateFreeVector(&tmp);
}

void evoHerm(state_t* state, const herm_t* herm, double par) {
    for (unsigned int i = 0; i < herm->len; ++i) {
        evoCG(state, herm->comp[i], herm->weight[i] * par);
    }
}

void evoPCGwrapper(state_t* state, const applyPCG pcg, double par) {
    par /= 2.;
    pcg(state, par);
}

void applyPQC(state_t* state, const pqc_t pqc, const double par[]) {
    for (depth_t i = 0; i < pqc.len; ++i) {
        pqc.pqb[i](state, *(par + pqc.parIdx[i]));
    }
}