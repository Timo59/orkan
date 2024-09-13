#ifndef QLIB_H
#define QLIB_H

/*
 * =================================================================================================
 *                                              includes
 * =================================================================================================
 */

#include <complex.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

/*
 * =================================================================================================
 *                                         macros
 * =================================================================================================
 */

#define INVSQRT2                0.7071067811865475
	
#define MIN(a,b)                ((a) < (b) ? (a) : (b))
#define MAX(a,b)                ((a) > (b) ? (a) : (b))
#define SWAP(a, b, T)           do { register T q; q = *(a); *(a) = *(b); *(b) = q; } while(0)
#define POW2(a, T)              (1 << (T) (a))

#define SETREAL(a, b)           (*((double*) &(a)) = (b))
#define SETIMAG(a, b)           (*(((double*) &(a)) + 1) = (b))

#define numberPrint(X) _Generic((X), \
    double: rnumberPrint, \
    cplx_t: cnumberPrint \
)

#define vectorPrint(X, Y) _Generic((X), \
    double*: rvectorPrint, \
    cplx_t*: cvectorPrint \
) ((X), (Y))

#define matrixPrint(X, Y) _Generic((X), \
    double*: rmatrixPrint, \
    cplx_t*: cmatrixPrint \
) (X, Y)

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
 *                                              type definitions
 * =================================================================================================
 */

typedef unsigned        complength_t;
typedef unsigned        qubit_t;
typedef unsigned        dim_t;
typedef unsigned        mat_t;
typedef unsigned        depth_t;
typedef unsigned        iter_t;
typedef double complex  cplx_t;

typedef enum pauli {
    ID,
    X,
    Y,
    Z
} pauli_t;

/*
 * Struct:      state_t
 * --------------------
 * Description: This struct represents a multi-qubit pure state.
 * Contents:    
 *      vector:             Vector representation of the state.
 *      qubits:             Number of qubits of the underlying system.
 *      length:             Length of the state's vector representation.
 */
typedef struct state {
    cplx_t* vector;
    qubit_t qubits;
    dim_t dimension;
} state_t;

typedef struct pauliOp {
    pauli_t* components;
    cplx_t* coefficients;
    complength_t length;
    qubit_t qubits;
} pauliOp_t;

typedef struct pauliObs {
    pauli_t* components;
    double* coefficients;
    complength_t length;
    qubit_t qubits;
} pauliObs_t;

typedef enum opType {
    DIAG,
    PAULI
} opType_t;

typedef struct op {
    opType_t type;
    union {
        cplx_t* diagOp;
        pauliOp_t* pauliOp;
    };
} op_t;

typedef struct obs {
    opType_t type;
    union {
        double* diagObs;
        pauliObs_t* pauliObs;
    };
} obs_t;

/*
 * =================================================================================================
 *                                              imaginary power multiplier
 * =================================================================================================
 */

void multiplyByIPower(cplx_t* number, int8_t power);

/*
 * =================================================================================================
 *                                              complex almost equalities
 * =================================================================================================
 */

bool calmostEqual(cplx_t a, cplx_t b, double precision);

bool cvectorAlmostEqual(cplx_t v[], cplx_t w[], dim_t dim, double precision);

/*
 * =================================================================================================
 *                                              double prints
 * =================================================================================================
 */

void rnumberPrint(double number);

void rvectorPrint(double vector[], dim_t dim);

void rmatrixPrint(double matrix[], dim_t dim);

/*
 * =================================================================================================
 *                                              complex prints
 * =================================================================================================
 */

void cnumberPrint(cplx_t number);

void cvectorPrint(cplx_t vector[], dim_t dim);

void cmatrixPrint(cplx_t matrix[], dim_t dim);

/*
 * =================================================================================================
 *                                              purely real matrix/vector operations
 * =================================================================================================
 */

double* matAdd(const double a[], const double b[], dim_t dim);

void matAddInPlace(double a[], const double b[], dim_t dim);

double* vecAdd(const double a[], const double b[], dim_t dim);

void vecAddInPlace(double a[], const double b[], dim_t dim);

double vecProduct(const double a[], const double b[], dim_t dim);

double* scalarMatMul(double a, const double b[], dim_t dim);

void scalarMatMulInPlace(double a, double b[], dim_t dim);

double* scalarVecMul(double a, const double b[], dim_t dim);

void scalarVecMulInPlace(double a, double b[], dim_t dim);

double* matMul(const double a[], const double b[], dim_t dim);

void matMulInPlace(double a[], const double b[], dim_t dim);

double* matVecMul(const double a[], const double b[], dim_t dim);

double* kronecker(const double a[], const double b[], dim_t dimA, dim_t dimB);

/* 
 * =================================================================================================
 *                                              purely complex matrix/vector operations
 * =================================================================================================
 */

cplx_t* cmatAdd(const cplx_t a[], const cplx_t b[], dim_t dim);

void cmatAddInPlace(cplx_t a[], const cplx_t b[], dim_t dim);

cplx_t* cvecAdd(const cplx_t a[], const cplx_t b[], dim_t dim);

void cvecAddInPlace(cplx_t a[], const cplx_t b[], dim_t dim);

cplx_t cvecProduct(const cplx_t a[], const cplx_t b[], dim_t dim);

cplx_t cinnerProduct(const cplx_t a[], const cplx_t b[], dim_t dim);

cplx_t* cscalarMatMul(cplx_t a, const cplx_t b[], dim_t dim);

void cscalarMatMulInPlace(cplx_t a, cplx_t b[], dim_t dim);

cplx_t* cscalarVecMul(cplx_t a, const cplx_t b[], dim_t dim);

void cscalarVecMulInPlace(cplx_t a, cplx_t b[], dim_t dim);

cplx_t* cmatMul(const cplx_t a[], const cplx_t b[], dim_t dim);

void cmatMulInPlace(cplx_t a[], const cplx_t b[], dim_t dim);

cplx_t* cmatVecMul(const cplx_t a[], const cplx_t b[], dim_t dim);

void cmatVecMulInPlace(const cplx_t a[], cplx_t v[], dim_t dim);

cplx_t* ckronecker(const cplx_t a[], const cplx_t b[], dim_t dimA, dim_t dimB);

/*
 * =================================================================================================
 *                                              state construction
 * =================================================================================================
 */

void stateInitEmpty(state_t* state, qubit_t qubits);

void stateInitZero(state_t* state, qubit_t qubits);

void stateInitPlus(state_t* state, qubit_t qubits);

void stateInitVector(state_t* state, cplx_t vector[], qubit_t qubits);

void stateCopyVector(state_t* state, const cplx_t vector[]);

/*
 * =================================================================================================
 *                                  		free state
 * =================================================================================================
 */

void stateFreeVector(state_t* state);

void stateFree(state_t* state);

/*
 * =================================================================================================
 *                                          binary operations
 * =================================================================================================
 */

cplx_t stateOverlap(const state_t* state1, const state_t* state2);

/*
 * =================================================================================================
 *                                              Pauli gates
 * =================================================================================================
 */

void applyX(state_t* state, qubit_t qubit);

void applyY(state_t* state, qubit_t qubit);

void applyZ(state_t* state, qubit_t qubit);

/*
 * =================================================================================================
 *                                              Clifford gates
 * =================================================================================================
 */

void applyS(state_t* state, qubit_t qubit);

void applySdagger(state_t* state, qubit_t qubit);

void applyH(state_t* state, qubit_t qubit);

/*
 * =================================================================================================
 *                                              Hadamard-Y gate
 * =================================================================================================
 */

void applyHy(state_t* state, qubit_t qubit);

/*
 * =================================================================================================
 *                                              T gate
 * =================================================================================================
 */

void applyT(state_t* state, qubit_t qubit);

void applyTdagger(state_t* state, qubit_t qubit);

/*
 * =================================================================================================
 *                                              P gate
 * =================================================================================================
 */

void applyP(state_t* state, qubit_t qubit, double angle);

void applyPdagger(state_t* state, qubit_t qubit, double angle);

/*
 * =================================================================================================
 *                                              rotation gates
 * =================================================================================================
 */

void applyRX(state_t* state, qubit_t qubit, double angle);

void applyRXdagger(state_t* state, qubit_t qubit, double angle);

void applyRY(state_t* state, qubit_t qubit, double angle);

void applyRYdagger(state_t* state, qubit_t qubit, double angle);

void applyRZ(state_t* state, qubit_t qubit, double angle);

void applyRZdagger(state_t* state, qubit_t qubit, double angle);

/*
 * =================================================================================================
 *                                              controlled Pauli gates
 * =================================================================================================
 */

void applyCX(state_t* state, qubit_t control, qubit_t target);

void applyCY(state_t* state, qubit_t control, qubit_t target);

void applyCZ(state_t* state, qubit_t control, qubit_t target);

/*
 * =================================================================================================
 *                                              controlled Clifford gates
 * =================================================================================================
 */

void applyCS(state_t* state, qubit_t control, qubit_t target);

void applyCSdagger(state_t* state, qubit_t control, qubit_t target);

void applyCH(state_t* state, qubit_t control, qubit_t target);

/*
 * =================================================================================================
 *                                              controlled Hadamard-Y gate
 * =================================================================================================
 */

void applyCHy(state_t* state, qubit_t control, qubit_t target);

/*
 * =================================================================================================
 *                                              controlled T gate
 * =================================================================================================
 */

void applyCT(state_t* state, qubit_t control, qubit_t target);

void applyCTdagger(state_t* state, qubit_t control, qubit_t target);

/*
 * =================================================================================================
 *                                              controlled P gate
 * =================================================================================================
 */

void applyCP(state_t* state, qubit_t control, qubit_t target, double angle);

void applyCPdagger(state_t* state, qubit_t control, qubit_t target, double angle);

/*
 * =================================================================================================
 *                                              SWAP gate
 * =================================================================================================
 */

void applySWAP(state_t* state, qubit_t qubit1, qubit_t qubit2);

/*
 * =================================================================================================
 *                                              Toffoli gate
 * =================================================================================================
 */

void applyToffoli(state_t* state, qubit_t control1, qubit_t control2, qubit_t target);

/*
 * =================================================================================================
 *                                              component applications
 * =================================================================================================
 */

void applyPauliStr(state_t* state, const pauli_t paulistr[]);

void evolvePauliStr(state_t* state, const pauli_t paulistr[], double angle);

/*
 * =================================================================================================
 *                                              apply operator
 * =================================================================================================
 */

void applyOperator(state_t* state, const op_t* op);

/*
 * =================================================================================================
 *                                              apply observable
 * =================================================================================================
 */

void applyObservable(state_t* state, const obs_t* observable);

/*
 * =================================================================================================
 *                                              evolve with trotterized observable
 * =================================================================================================
 */

void evolveWithTrotterizedObservable(state_t* state, const obs_t* observable, double angle);

/*
 * =================================================================================================
 *                                              expectation value
 * =================================================================================================
 */

cplx_t expValPauli(const state_t* state, const pauliOp_t* op);

double expValObsPauli(const state_t* state, const pauliObs_t* observable);

cplx_t expValDiag(const state_t* state, const cplx_t opr[]);

double expValObsDiag(const state_t* state, const double observable[]);

double expValObs(const state_t* state, const obs_t* observable);

/*
 * =================================================================================================
 *                                              apply PQC
 * =================================================================================================
 */

void applyPQC(state_t* state, const double params[], const obs_t evoOps[], depth_t circdepth);

/*
 * =================================================================================================
 *                                              expectation value PQC
 * =================================================================================================
 */

double expValPQC(const state_t* state, const double params[], const obs_t* observable, \
                 const obs_t evoOps[], depth_t circdepth);

/*
 * =================================================================================================
 *                                              gradient PQC
 * =================================================================================================
 */

double* gradientPQC(const state_t* state, const double params[], const obs_t* observable, \
                   const obs_t evoOps[], depth_t circdepth);

/*
 * =================================================================================================
 *                                              hessian PQC
 * =================================================================================================
 */

double* hessianPQC(const state_t* state, const double params[], const obs_t* observable, \
                   const obs_t evoOps[], depth_t circdepth);

#ifdef __cplusplus
}
#endif

#endif
