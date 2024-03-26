#ifndef EXACT_H
#define EXACT_H

/*
 * =====================================================================================================================
 *                                                      includes
 * =====================================================================================================================
 */

#include "cplxutil.h"
#include "circlib.h"

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

typedef double* (*gradPQC)(const state_t* state, const double params[], const obs_t* observable, \
                    const obs_t* evoOps[], depth_t circdepth);

typedef void (*srchU)(const state_t* state, const obs_t* srchOp);

/*
 * =====================================================================================================================
 *                                                  Mean value operator
 * =====================================================================================================================
 */
cplx_t expValPauli(const state_t* state, const pauliOp_t* op);

cplx_t expValDiag(const state_t* state, const cplx_t op[]);

cplx_t expVal(const state_t* state, const op_t* operator);

/*
 * =====================================================================================================================
 *                                                  Mean value observable
 * =====================================================================================================================
 */
double expValObsPauli(const state_t* state, const pauliObs_t* observable);

double expValObsDiag(const state_t* state, const double observable[]);

double expValObs(const state_t* state, const obs_t* observable);

/*
 * =====================================================================================================================
 *                                                  Mean value PQC
 * =====================================================================================================================
 */

double expValObsPauliPQC(const state_t* state, const double params[], const pauliObs_t* observable, \
                         const obs_t* evoOps[], depth_t circdepth);

double expValObsDiagPQC(const state_t* state, const double params[], const double observable[], \
                        const obs_t* evoOps[], depth_t circdepth);

double expValObsPQC(const state_t* state, const double params[], const obs_t* observable, \
                    const obs_t* evoOps[], depth_t circdepth);
/*
 * =====================================================================================================================
 *                                                      gradient PQC
 * =====================================================================================================================
 */

double* gradientPQC(const state_t* state, const double params[], const obs_t* observable, \
                   const obs_t* evoOps[], depth_t circdepth);

double* approxGradientPQC(const state_t* state, const double params[], const obs_t* observable, \
                               const obs_t* evoOps[], depth_t circdepth, double epsilon);

/*
 * =====================================================================================================================
 *                                                      hessian PQC
 * =====================================================================================================================
 */

double* hessianPQC(const state_t* state, const double params[], const obs_t* observable, \
                   const obs_t* evoOps[], depth_t circdepth);

double* approxHessianPQC(const state_t* state, const double params[], const obs_t* observable, \
                         const obs_t* evoOps[], depth_t circdepth, gradPQC grad, double epsilon);

/*
 * =====================================================================================================================
 *                                                  Moment matrices
 * =====================================================================================================================
 */

cplx_t** momMat(const state_t* state, const obs_t* obs[], depth_t obsc, const obs_t* srchOps[], depth_t matDim, \
                srchU applyU);


#ifdef __cplusplus
}
#endif

#endif
