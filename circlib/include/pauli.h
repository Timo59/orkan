#ifndef PAULI_H
#define PAULI_H

/*
 * =====================================================================================================================
 *                                                      includes
 * =====================================================================================================================
 */

#if defined(PAULI_BLAS)
#include "pauliBlas.h"
#elif defined(PAULI_OMP)
#include "pauliOmp.h"
#elif defined(PAULI_BLAS_OMP)
#include "pauliBlMp.h"
#elif defined(TEST)
#include "pauliBlas.h"
#include "pauliOmp.h"
#include "pauliBlMp.h"
#include "pauliOld.h"
#else
#include "pauliOld.h"
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
 *                                                  component applications
 * =====================================================================================================================
 */

void applyPauliStr(state_t* state, const pauli_t paulistr[]);

void applyOpPauli(state_t* state, const pauliOp_t* op);

void applyObsPauli(state_t* state, const pauliObs_t* obs);

void evolvePauliStr(state_t* state, const pauli_t paulistr[], double angle);

void evolveObsPauliTrotter(state_t* state, const pauliObs_t* obs, double angle);

void evolveObsPauliTrotterDagger(state_t* state, const pauliObs_t* obs, double angle);

#ifdef __cplusplus
}
#endif

#endif /* PAULI_H */