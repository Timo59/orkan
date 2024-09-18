#ifndef PAULI_H
#define PAULI_H

/*
 * =====================================================================================================================
 *                                                      includes
 * =====================================================================================================================
 */

#ifndef GATELIB_H
#include "../../gatelib/gatelib.h"
#endif

#ifndef LINALG_H
#include "linalg.h"
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

/*
 * =====================================================================================================================
 *                                                  component applications
 * =====================================================================================================================
 */

void applyPauliStr(state_t* state, const pauli_t paulistr[]);

void applyOpPauli(state_t* state, const pauliOp_t* op);

void applyObsPauli(state_t* state, const pauliObs_t* obs);

void applyObsPauliBlas(state_t* state, const pauliObs_t* observable);

void applyObsPauliOmp(state_t* state, const pauliObs_t* observable);

void applyObsPauliBlas_omp(state_t* state, const pauliObs_t* observable);

void applyObsPauliBlas_omp_new(state_t* state, const pauliObs_t* observable);

void evolvePauliStr(state_t* state, const pauli_t paulistr[], double angle);

void evolvePauliStrBlas(state_t* state, const pauli_t paulistr[], double angle);

void evolvePauliStrOmp(state_t* state, const pauli_t paulistr[], double angle);

void evolveObsPauliTrotter(state_t* state, const pauliObs_t* obs, double angle);

void evolveObsPauliTrotterDagger(state_t* state, const pauliObs_t* obs, double angle);

#ifdef __cplusplus
}
#endif

#endif /* PAULI_H */