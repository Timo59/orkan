//
// Created by Timo Ziegler on 12.03.25.
//

#ifndef EXACT_H
#define EXACT_H
/*
 * =====================================================================================================================
 *                                                      Includes
 * =====================================================================================================================
 */
#ifndef PSTATE_H
#include "pstate.h"
#endif

#ifndef PQC_H
#include "pqc.h"
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
 *                                                  Function definitions
 * =====================================================================================================================
 */
double meanDiagObs(const state_t* state, const double obs[]);
double meanObsQb(const state_t* state, applyQB obs);
#define meanObs(X, Y) _Generic((Y), \
    double*: meanDiagObs,           \
    applyQB: meanObsQb              \
    ) (X, Y)
void gradPQC(state_t* state, depth_t d, const applyPQB pqc[], const double par[], const applyQB qb[],
    const double obs[], double* grad);
void mmseq(state_t* state, depth_t obsc, const applyQB obs[], depth_t circdepth, const cplx_t c[], depth_t uc,
           const applyQB u[], depth_t link, cplx_t* momMat[]);

#ifdef __cplusplus
}
#endif

#endif //EXACT_H
