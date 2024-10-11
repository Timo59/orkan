#ifndef EXACT_H
#define EXACT_H

/*
 * =====================================================================================================================
 *                                                      includes
 * =====================================================================================================================
 */

#if defined(EXACT_BLAS)
#include "exactBlas.h"
#elif defined(EXACT_OMP)
#include "exactOmp.h"
#elif defined(EXACT_BLAS_OMP)
#include "exactBlMp.h"
#elif defined(TEST)
#include "exactBlas.h"
#include "exactOmp.h"
#include "exactBlMp.h"
#include "exactOld.h"
#else
#include "exactOld.h"
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
*                                                      Macros
* =====================================================================================================================
*/

/*#define expValObs(X, Y) _Generic((Y), \
                        double*:        expValObsDiag, \
                        pauliObs_t*:    expValObsPauli \
                        ) (X, Y)*/
/*
 * =====================================================================================================================
 *                                                  type definitions
 * =====================================================================================================================
 */

typedef double* (*gradientPQC)(const state_t* state, const double par[], const obs_t* obs, const obs_t* evoOps[],
                               depth_t circdepth);

typedef void (*srchU)(state_t* state, const obs_t* srchOp);

/*
 * =====================================================================================================================
 *                                                  Mean value operator
 * =====================================================================================================================
 */
cplx_t expValPauli(const state_t* state, const pauliOp_t* op);

cplx_t expValDiag(const state_t* state, const cplx_t op[]);

cplx_t expVal(const state_t* state, const op_t* op);

/*
 * =====================================================================================================================
 *                                                  Mean value observable
 * =====================================================================================================================
 */
double expValObsPauli(const state_t* state, const pauliObs_t* obs);

double expValObsDiag(const state_t* state, const double observable[]);

double expValObs(const state_t* state, const obs_t* obs);

/*
 * =====================================================================================================================
 *                                                  Mean value PQC
 * =====================================================================================================================
 */

double expValObsPqcPauli(const state_t* state, const double par[], const pauliObs_t* obs, \
                         const obs_t *evoOps[], depth_t circdepth);

double expValObsPqcDiag(const state_t* state, const double par[], const double obs[], \
                        const obs_t *evoOps[], depth_t circdepth);

double expValObsPQC(const state_t* state, const double params[], const obs_t* obs, \
                    const obs_t* evoOps[], depth_t circdepth);
/*
 * =====================================================================================================================
 *                                                      gradient PQC
 * =====================================================================================================================
 */

double* gradPQC(const state_t* state, const double par[], const obs_t* obs, const obs_t *evoOps[],
                depth_t circdepth);

double* approxGradientPQC(const state_t* state, const double par[], const obs_t* observable,
                               const obs_t* evoOps[], depth_t circdepth, double epsilon);

/*
 * =====================================================================================================================
 *                                                      hessian PQC
 * =====================================================================================================================
 */

double* hessianPQC(const state_t* state, const double params[], const obs_t* observable, \
                   const obs_t* evoOps[], depth_t circdepth);

double* approxHessianPQC(const state_t* state, const double params[], const obs_t* observable, \
                         const obs_t* evoOps[], depth_t circdepth, gradientPQC grad, double epsilon);

/*
 * =====================================================================================================================
 *                                                  Moment matrices
 * =====================================================================================================================
 */

cplx_t** momMat(const state_t* state, const obs_t* obs[], depth_t obsc, const obs_t* srchOps[], depth_t matDim, \
                srchU applyU);

cplx_t** momMatPQC(const state_t* state, const obs_t* obs[], depth_t obsc, const obs_t* evoOps[], depth_t circdepth, \
                   double angles[]);


#ifdef __cplusplus
}
#endif

#endif /* EXACT_H */
