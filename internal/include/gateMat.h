//
// Created by Timo Ziegler on 18.09.24.
//

#ifndef GATEMAT_H
#define GATEMAT_H

#ifndef LINALG_H
#include "linalg.h"
#endif

#ifndef __MATH__
#include <math.h>
#endif

#ifndef QTYPES_H
#include "qTypes.h"
#endif

#ifndef _STDIO_H_
#include <stdio.h>
#endif

/*
 * =====================================================================================================================
 *                                                  Static matrices
 * =====================================================================================================================
 */

static const cplx_t IDMAT[4] = {1.0 + 0.0 * I, \
                                0.0 + 0.0 * I, \
                                0.0 + 0.0 * I, \
                                1.0 + 0.0 * I};

static const cplx_t XMAT[4] = {0.0 + 0.0 * I, \
                               1.0 + 0.0 * I, \
                               1.0 + 0.0 * I, \
                               0.0 + 0.0 * I};

static const cplx_t YMAT[4] = {0.0 + 0.0 * I, \
                               0.0 - 1.0 * I, \
                               0.0 + 1.0 * I, \
                               0.0 + 0.0 * I};

static const cplx_t ZMAT[4] = {1.0 + 0.0 * I, \
                               0.0 + 0.0 * I, \
                               0.0 + 0.0 * I, \
                               -1.0 + 0.0 * I};


static const cplx_t SMAT[4] = {1.0 + 0.0 * I, \
                               0.0 + 0.0 * I, \
                               0.0 + 0.0 * I, \
                               0.0 + 1.0 * I};

static const cplx_t HMAT[4] = {INVSQRT2 + 0.0 * I, \
                               INVSQRT2 + 0.0 * I, \
                               INVSQRT2 + 0.0 * I, \
                               -INVSQRT2 + 0.0 * I};

static const cplx_t TMAT[4] = {1.0 + 0.0 * I, \
                               0.0 + 0.0 * I, \
                               0.0 + 0.0 * I, \
                               INVSQRT2 + INVSQRT2 * I};

static const cplx_t P0MAT[4] = {1.0 + 0.0 * I, \
                                0.0 + 0.0 * I, \
                                0.0 + 0.0 * I, \
                                0.0 + 0.0 * I};

static const cplx_t P1MAT[4] = {0.0 + 0.0 * I, \
                                0.0 + 0.0 * I, \
                                0.0 + 0.0 * I, \
                                1.0 + 0.0 * I};

static const cplx_t SWAPMAT[16] = {1.0 + 0.0 * I, 0.0 + 0.0 * I, 0.0 + 0.0 * I, 0.0 + 0.0 * I, \
                                   0.0 + 0.0 * I, 0.0 + 0.0 * I, 1.0 + 0.0 * I, 0.0 + 0.0 * I, \
                                   0.0 + 0.0 * I, 1.0 + 0.0 * I, 0.0 + 0.0 * I, 0.0 + 0.0 * I, \
                                   0.0 + 0.0 * I, 0.0 + 0.0 * I, 0.0 + 0.0 * I, 1.0 + 0.0 * I};

static cplx_t ONE[1] = {1.0 + 0.0 * I};

/*
* =====================================================================================================================
*                                                  Gate constructions
* =====================================================================================================================
*/

cplx_t* identityMat(qubit_t qubits);

cplx_t* singleQubitGateMat(const cplx_t mat[], qubit_t qubits, qubit_t pos);

cplx_t* RXGateMat(qubit_t qubits, qubit_t pos, double angle);

cplx_t* RYGateMat(qubit_t qubits, qubit_t pos, double angle);

cplx_t* RZGateMat(qubit_t qubits, qubit_t pos, double angle);

cplx_t* controlledGateMat(const cplx_t mat[], qubit_t qubits, qubit_t control, qubit_t target);

cplx_t* neighboringSwapGateMat(qubit_t qubits, qubit_t target);

cplx_t* swapGateMat(qubit_t qubits, qubit_t qubit1, qubit_t qubit2);

cplx_t* RswapGateMat(qubit_t qubits, qubit_t qubit1, qubit_t qubit2, double angle);

cplx_t* doublyControlledGateMat(const cplx_t mat[],qubit_t qubits, qubit_t control1, qubit_t control2, qubit_t target);

#endif /* GATEMAT_H */
