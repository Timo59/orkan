//
// Created by Timo Ziegler on 19.09.24.
//

#include "finDiff.h"



/*
 * =====================================================================================================================
 *                                              Finite difference methods
 * =====================================================================================================================
 */

/*
 * This function computes an observable's mean value on a state vector after applying a PQC
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
                       double* par)
{
    /*
     * 1. Calculate the observable's matrix
     */
    cplx_t* observableMat = pauliObsMat(compObs, coeffObs, lengthObs, qubits);

    /*
     * 2. Initialize the matrices for the evolution operators
     */
    cplx_t** evoOpsMat = (cplx_t**) malloc(circdepth * sizeof(cplx_t*));
    for (depth_t i = 0; i < circdepth; ++i) {
        evoOpsMat[i] = expPauliObsMatTrotter(compEvoOps,
                                             coeffEvoOps + (i * lengthEvoOps),
                                             lengthEvoOps,
                                             par[i],
                                             qubits);
    }

    /*
     * 3. Evolve copies of the statevector with all evolution operators
     */
    cplx_t* bra = cmatVecMul(evoOpsMat[0], statevector, dim);
    for (depth_t i = 1; i < circdepth; ++i) {
        cmatVecMulInPlace(evoOpsMat[i], bra, dim);
    }

    cplx_t* ket = cmatVecMul(evoOpsMat[0], statevector, dim);
    for (depth_t i = 1; i < circdepth; ++i) {
        cmatVecMulInPlace(evoOpsMat[i], ket, dim);
    }

    /*
     * 4. Apply the observable to one of the copies and compute their inner product
     */
    cmatVecMulInPlace(observableMat, ket, dim);

    double result = creal(cInner(bra, ket, dim));

    free(bra);
    free(ket);
    free(observableMat);
    for (depth_t i = 0; i < circdepth; ++i) {
        free(evoOpsMat[i]);
    }
    free(evoOpsMat);

    return result;
}

/*
 * This function computes the gradient of an observable's mean value on a state vector after applying a PQC based on a
 * finite difference method at the current parameter setting
 */
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
                          double epsilon)
{
    double* result = malloc(circdepth * sizeof(double));    // The resulting array holding the gradients components

    double expVal = finiteExpValPQC(statevector,                // Calculate the obserable's expectation value at the
                                    qubits,                     // current parameter setting
                                    dim,
                                    circdepth,
                                    compObs,
                                    coeffObs,
                                    lengthObs,
                                    compEvoOps,
                                    coeffEvoOps,
                                    lengthEvoOps,
                                    par);

    /*
     * Iterate the entries of the parameter vector, shift it by epsilon, calculate the observable's expectation value
     * and enter the difference quotient to result's component
     */
    for (depth_t i = 0; i < circdepth; ++i) {
        par[i] += epsilon;

        double value = finiteExpValPQC(statevector,
                                       qubits,
                                       dim,
                                       circdepth,
                                       compObs,
                                       coeffObs,
                                       lengthObs,
                                       compEvoOps,
                                       coeffEvoOps,
                                       lengthEvoOps,
                                       par);

        result[i] = (1. / epsilon) * (value - expVal);
        par[i] -= epsilon;
    }

    return result;
}

/*
 * This function computes the hessian of an observable's mean value on a state vector after applying a PQC based on a
 * finite difference method at the current parameter setting
 */

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
                         double epsilon)
{
    double* result = malloc(circdepth*circdepth * sizeof(double));

    /*
     * 1. Calculate the gradient at the current parameter setting
     */
    double* refGradient = finiteGradientPQC(statevector,
                                            qubits,
                                            dim,
                                            circdepth,
                                            compObs,
                                            coeffObs,
                                            lengthObs,
                                            compEvoOps,
                                            coeffEvoOps,
                                            lengthEvoOps,
                                            par,
                                            epsilon);

    /*
     * 2. Shift the parameter setting by epsilon in all directions; the entrywise difference to the prvious gradient
     * gives the corresponding column of the hessian
     */
    for (depth_t i = 0; i < circdepth; ++i) {
        par[i] += epsilon;
        double* gradient = finiteGradientPQC(statevector,
                                             qubits,
                                             dim,
                                             circdepth,
                                             compObs,
                                             coeffObs,
                                             lengthObs,
                                             compEvoOps,
                                             coeffEvoOps,
                                             lengthEvoOps,
                                             par,
                                             epsilon);

        for (depth_t j = 0; j < circdepth; ++j) {
            result[j * circdepth + i] = (1. / epsilon) * (gradient[j] - refGradient[j]);
        }
        par[i] -= epsilon;
        free(gradient);
    }
    free(refGradient);
    return result;
}
