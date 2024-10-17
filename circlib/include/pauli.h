#ifndef PAULI_H
#define PAULI_H

/*
 * =====================================================================================================================
 *                                                      includes
 * =====================================================================================================================
 */
#ifndef GATELIB_H
#include "gatelib.h"
#endif

#if defined(PAULI_BLAS)
#ifndef __VECLIB__
#include <vecLib/vecLib.h>
#endif

#elif defined(PAULI_BLAS_OMP)
#ifndef __OMP_H
#include <omp.h>
#endif

#ifndef __VECLIB__
#include <vecLib/vecLib.h>
#endif

#elif defined(PAULI_OMP)
#ifndef LINALG_H
#include "linalg.h"
#endif

#ifndef __OMP_H
#include <omp.h>
#endif

#elif defined(PAULI_PROFILE)
#ifndef CIRC_PROFILE_H
#include "profile.h"
#endif

#else
#ifndef LINALG_H
#include "linalg.h"
#endif

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
*                                                  type definitions
* =====================================================================================================================
*/
#ifndef PAULI_STRUCTS
#define PAULI_STRUCTS
typedef struct pauliOp {
    pauli_t* comps;
    cplx_t* coeffs;
    complength_t length;
    qubit_t qubits;
} pauliOp_t;

typedef struct pauliObs {
    pauli_t* comps;
    double* coeffs;
    complength_t length;
    qubit_t qubits;
} pauliObs_t;
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