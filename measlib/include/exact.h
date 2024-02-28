#ifndef EXACT_H
#define EXACT_H

/*
 * =================================================================================================
 *                                              includes
 * =================================================================================================
 */

#include "cplxutil.h"
#include "circlib.h"

/*
 * =================================================================================================
 *                                              C++ check
 * =================================================================================================
 */

#ifdef __cplusplus
extern "C" {
#endif

/*
 * =================================================================================================
 *                                              expectation value
 * =================================================================================================
 */

cplx_t expValPauli(const state_t* state, const pauliOp_t* op);

double expValObsPauli(const state_t* state, const pauliObs_t* observable);

cplx_t expValDiag(const state_t* state, const cplx_t op[]);

double expValObsDiag(const state_t* state, const double observable[]);

double expValObs(const state_t* state, const obs_t* observable);

/*
 * =================================================================================================
 *                                              apply PQC
 * =================================================================================================
 */

void applyPQC(state_t* state, const double params[], const obs_t* evoOps[], depth_t circdepth);

/*
 * =================================================================================================
 *                                              expectation value PQC
 * =================================================================================================
 */

double expValPQC(const state_t* state, const double params[], const obs_t* observable, \
                 const obs_t* evoOps[], depth_t circdepth);

/*
 * =================================================================================================
 *                                              gradient PQC
 * =================================================================================================
 */

double* gradientPQC(const state_t* state, const double params[], const obs_t* observable, \
                   const obs_t* evoOps[], depth_t circdepth);

/*
 * =================================================================================================
 *                                              hessian PQC
 * =================================================================================================
 */

double* hessianPQC(const state_t* state, const double params[], const obs_t* observable, \
                   const obs_t* evoOps[], depth_t circdepth);

/*
 * =================================================================================================
 *                                    Moment matrices
 * =================================================================================================
 */

cplx_t* momentMatrices(const state_t* state, const obs_t* searchHerms, double* angles, \
                       depth_t searchDim, obs_t* observable);

#ifdef __cplusplus
}
#endif

#endif
