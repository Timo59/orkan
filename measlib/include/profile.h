//
// Created by Timo Ziegler on 17.10.24.
//
#ifndef EXACT_PROFILE_H
#define EXACT_PROFILE_H

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
*                                                  Mean value operator
* =====================================================================================================================
*/
cplx_t expValPauli_blas(const state_t* state, const pauliOp_t* op);

cplx_t expValPauli_blas_omp(const state_t* state, const pauliOp_t* op);

cplx_t expValPauli_omp(const state_t* state, const pauliOp_t* op);

/*
 * =====================================================================================================================
 *                                                  Mean value observable
 * =====================================================================================================================
 */
double expValObsPauli_blas(const state_t* state, const pauliObs_t* obs);

double expValObsPauli_blas_omp(const state_t* state, const pauliObs_t* obs);

double expValObsPauli_omp(const state_t* state, const pauliObs_t* obs);

/*
 * =====================================================================================================================
 *                                                      gradient PQC
 * =====================================================================================================================
 */
double* gradPQC_blas(const state_t* state, const double par[], const obs_t* obs, const obs_t *evoOps[],
                     depth_t circdepth);

double* gradPQC_blas_omp(const state_t* state, const double par[], const obs_t* obs, const obs_t *evoOps[],
                         depth_t circdepth);

double* gradPQC_omp(const state_t* state, const double par[], const obs_t* obs, const obs_t *evoOps[],
                    depth_t circdepth);

#ifdef __cplusplus
extern "C" {
#endif

#endif /* eXACT_PROFILE_H */
