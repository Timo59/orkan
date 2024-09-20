//
// Created by Timo Ziegler on 19.09.24.
//

#ifndef FINDIFF_H
#define FINDIFF_H

/*
 * =====================================================================================================================
 *                                                      Includes
 * =====================================================================================================================
 */

#ifndef PAULIMAT_H
#include "pauliMat.h"
#endif

#ifndef QTYPES_H
#include "qTypes.h"
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
 *                                              Finite difference methods
 * =====================================================================================================================
 */

double finiteExpValPQC(cplx_t* statevector,
                       qubit_t qubits,
                       dim_t dim,
                       depth_t circdepth,
                       pauli_t* compObs,
                       double* coeffObs,
                       complength_t lengthObs,
                       pauli_t* compEvoOps,
                       double* coeffEvoOps,
                       complength_t lengthEvoOps,
                       double* par);

double* finiteGradientPQC(cplx_t* statevector,
                          qubit_t qubits,
                          dim_t dim,
                          depth_t circdepth,
                          pauli_t* compObs,
                          double* coeffObs,
                          complength_t lengthObs,
                          pauli_t* compEvoOps,
                          double* coeffEvoOps,
                          complength_t lengthEvoOps,
                          double* par,
                          double epsilon);

double* finiteHessianPQC(cplx_t* statevector,
                         qubit_t qubits,
                         dim_t dim,
                         depth_t circdepth,
                         pauli_t* compObs,
                         double* coeffObs,
                         complength_t lengthObs,
                         pauli_t* compEvoOps,
                         double* coeffEvoOps,
                         complength_t lengthEvoOps,
                         double* par,
                         double epsilon);

#ifdef __cplusplus
}
#endif

#endif /* FINDIFF_H */
