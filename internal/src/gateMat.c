//
// Created by Timo Ziegler on 18.09.24.
//

#include "../../test/gateMat.h"

/*
 * =====================================================================================================================
 *                                                  Gate constructions
 * =====================================================================================================================
 */

cplx_t* identityMat(qubit_t qubits) {
    cplx_t* result = ONE;
    cplx_t* tmp;
    for (qubit_t i = 0; i < qubits; ++i) {
        tmp = ckronecker(result, IDMAT, POW2(i, dim_t), 2);
        if (i != 0) {
            free(result);
        }
        result = tmp;
    }
    return result;
}

/*
 * This function returns the matrix representation of single qubit gate applied at pos in a system of qubits.
 *
 * Input:
 *      cplx_t mat[]:       Array representing the single qubit gate
 *      qubit_t qubits:     Number of qubits in the system
 *      qubit_t pos:        Position of qubit the gate is applied to
 *
 * Output:
 *      Kronecker product of identity matrices and the single qubit matrix at the specified position.
 */
cplx_t* singleQubitGateMat(const cplx_t mat[], qubit_t qubits, qubit_t pos) {
    /*
     * If either no qubits are involved or the position of the affected qubit lies outside the range of available
     * qubits, the function terminates with an error message.
    */
    if (qubits == 0 || pos >= qubits) {
        printf("Couldn't create tensorized single qubit gate!\n");
        return NULL;
    }

    cplx_t* result = ONE;       // Array to hold the resulting matrix initialized to 1.0+0.0i
    cplx_t* tmp;                // Temporary array holding all intermediate Kronecker products

    for (qubit_t i = 0; i < qubits; ++i) {                                              // Iterate qubits from the left:
        if (i == qubits - pos - 1) {                                                    // If current is the targeted
                                                                                        // position seen from the right:
            tmp = ckronecker(result, mat, POW2(i, dim_t), 2);           // Kronecker product of
        }                                                                               // current matrix with single
                                                                                        // qubit matrix

        else {                                                                          // Otherwise: Kronecker product
            tmp = ckronecker(result, IDMAT, POW2(i, dim_t), 2);         // of the current matrix with
        }                                                                               // with another identity
        if (i != 0) {           // Free current result pointer except for in the first iteration(not dynamically alloc.)
            free(result);
        }
        result = tmp;           // Identify result with tmp's pointer address
    }
    return result;
}

cplx_t* RXGateMat(qubit_t qubits, qubit_t pos, double angle) {
    angle /= 2.;
    cplx_t mat[4] = {cos(angle) + 0.0 * I, \
                     0.0 - sin(angle) * I, \
                     0.0 - sin(angle) * I, \
                     cos(angle) + 0.0 * I};
    return singleQubitGateMat(mat, qubits, pos);
}

cplx_t* RYGateMat(qubit_t qubits, qubit_t pos, double angle) {
    angle /= 2.;
    cplx_t mat[4] = {cos(angle) + 0.0 * I, \
                     -sin(angle) + 0.0 * I, \
                     sin(angle) + 0.0 * I, \
                     cos(angle) + 0.0 * I};
    return singleQubitGateMat(mat, qubits, pos);
}

cplx_t* RZGateMat(qubit_t qubits, qubit_t pos, double angle) {
    angle /= 2.;
    cplx_t mat[4] = {cos(angle) - sin(angle) * I, \
                     0.0 + 0.0 * I, \
                     0.0 + 0.0 * I, \
                     cos(angle) + sin(angle) * I};
    return singleQubitGateMat(mat, qubits, pos);
}

/*
 * This function computes the matrix representation of a single qubit gate with conditional execution depending on the
 * state of control.
 *
 * Input:
 *      const cplx_t mat[]:     Array representing the single qubit gate
 *      qubit_t qubits:         Number of qubits involved in a state
 *      qubit_t control:        Position of the control qubit
 *      qubit_t target:         Position of the target qubit
 * Output:
 *      The matrix representation of a controlled single qubit gate built up by a sum of Kronecker products
*/
cplx_t* controlledGateMat(const cplx_t mat[], qubit_t qubits, qubit_t control, qubit_t target) {
    /*
     * If either no qubits are involved or the position of the affected qubit lies outside the range of available
     * qubits, the function terminates with an error message.
     */
    if (qubits == 0 || control >= qubits || target >= qubits || control == target) {
        printf("Couldn't create controlled single qubit gate!\n");
        return NULL;
    }

    cplx_t* active = ONE;        // Array holding the matrix for control qubit equals 1
    cplx_t* inactive = ONE;      // Array holding the matrix for control qubit equals 0
    cplx_t* tmpActive;           // Temporary holding each intermediate Kronecker product of active control
    cplx_t* tmpInactive;         // Temporary holding each intermediate Kronecker product of inactive control

    for (qubit_t i = 0; i < qubits; ++i) {      // Iterate qubits from the left:
        if (i == qubits - target - 1) {         // If current is the target (counted from the right)
            tmpActive = ckronecker(active, mat, POW2(i, dim_t), 2);       // Kronecker product of
                                                                                            // the active control matrix
                                                                                            // with the single qubit
                                                                                            // gate
            tmpInactive = ckronecker(inactive, IDMAT, POW2(i, dim_t), 2);   // Kronecker product of
                                                                                                // the inactive control
                                                                                                // matrix with identity
        }
        else if (i == qubits - control - 1) {   // If else the current is the control (counted from the right)
            tmpActive = ckronecker(active, P1MAT, POW2(i, dim_t), 2);       // Kronecker product of
                                                                                            // the active control matrix
                                                                                            // with the projector onto
                                                                                            // the one state
            tmpInactive = ckronecker(inactive, P0MAT, POW2(i, dim_t), 2);   // Kronecker product of
                                                                                            // the inactive control
                                                                                            // matrix with the projector
                                                                                            // onto the zero state
        }
        else {                                                                                  // If it corresponds to
            tmpActive = ckronecker(active, IDMAT, POW2(i, dim_t), 2);       // neither, Kronecker
            tmpInactive = ckronecker(inactive, IDMAT, POW2(i, dim_t), 2);   // product with identity
        }

        if (i != 0) {       // Free current result pointers except for in the first iteration(not dynamically alloc.)
            free(active);
            free(inactive);
        }
        active = tmpActive;          // active is now pointing to the outcome of the Kronecker product in tmpActive
        inactive = tmpInactive;      // inactive is now pointing to the outcome of the Kronecker product in tmpInactive
    }
    cmatAddInPlace(inactive, active, POW2(qubits, dim_t));  // Add the two matrices (they are projectors)
    free(active);
    return inactive;
}

/*
 * This function computes the matrix corresponding to the swap of the specified qubit with its neighbor to the left
 *
 * Input:
 *      qubit_t qubits:     Number of qubits
 *      qubit_t target:     Index of the qubit that is to be swapped with its neighbor to the right
 *
 * Output:
 *      The Kronecker product of identity matrices with a swap matrix at the position of target and qubit to its left
 *      counted from the right.
 */

cplx_t* neighboringSwapGateMat(qubit_t qubits, qubit_t target) {
    /*
    If either there are no qubits at all or one of the qubits to be swapped is out of range, return an error.
    */
    if (qubits < 2 || target >= qubits - 1) {
        fprintf(stderr, "Couldn't create a neighboring swap gate for qubit %d!\n", target);
        return NULL;
    }
    cplx_t* result = ONE;               // Array to hold the resulting matrix initialized to 1.0+0.0i
    cplx_t* tmp;                        // double complex pointer used for temporarily holding the adress
    // for the intermediate Kronecker products
    for (qubit_t i = 0; i < qubits; ++i){       // Iterate all qubits:
        if (i == target) {                          // If the current qubit is targeted
            tmp = ckronecker(SWAPMAT, result, 4, POW2(i, dim_t));   // Kronecker product with the SWAP
            // matrix from the left (qubits are
            // counted from the right)
            if (i != 0) {
                free(result);
            }
            ++i;
        }
        else {                                                                      // Otherwise:
            tmp = ckronecker(IDMAT, result, 2, POW2(i, dim_t));     // Kronecker product with ID
                                                                                    // matrix from the left
            if (i != 0) {
                free(result);
            }
        }

        result = tmp;
    }
    return result;
}

/*
 * This function computes the matrix for the swap of two arbitrary qubits by concatenating swaps of neighboring qubits.
 */
cplx_t* swapGateMat(qubit_t qubits, qubit_t qubit1, qubit_t qubit2) {
    dim_t dim = POW2(qubits, dim_t);                                    // Hilbert space dimension
    qubit_t left_qubit = MAX(qubit1, qubit2);                           // index of the qubit to the left to be swapped
    qubit_t right_qubit = MIN(qubit1, qubit2);                          // index of the qubit to the right to be swapped
    cplx_t* result = neighboringSwapGateMat(qubits, right_qubit);   // Array to later hold the swap gate
                                                                        // initialized to swap the right qubit with its
                                                                        // left neighbor
                                                                        // To swap right qubit to the left position:
    for (qubit_t i = right_qubit + 1; i < left_qubit; ++i) {            // Iterate qubits between left and right qubit
        cplx_t* neighborSwap = neighboringSwapGateMat(qubits, i);       // Calculate the matrix for swapping with
                                                                                // the next left neighbor
        cplx_t* tmp = cmatMul(neighborSwap, result, dim);               // Multiply the neighbor swap matrix from
                                                                                // the left
        free(result);
        free(neighborSwap);
        result = tmp;
    }
                                                                        // To swap the left qubit to the right position
    for (qubit_t i = left_qubit - 1; i > right_qubit; --i) {            // Iterate qubits from left to right
        cplx_t* neighborSwap = neighboringSwapGateMat(qubits, i - 1);   // Calculate the matrix for swapping with
                                                                                // the next left neighbor
        cplx_t* tmp = cmatMul(neighborSwap, result, dim);         // Multiply the neighbor swap matrix from
                                                                        // the left
        free(result);
        free(neighborSwap);
        result = tmp;
    }

    return result;
}

/*
 * This function computes the matrix representation of any single qubit gate with conditional execution depending on the
 * state of two other qubits.
 * Input:
 *      cplx_t mat[]:       Matrix representation of the gate to be applied to target
 *      qubit_t qubits:     The number of qubits
 *      qubit_t control1:   Index of the first control from the right
 *      qubit_t control2:   Index of the second control from the right
 *      qubit_t target:     Index of the target from the right
 * Output:
 *      The matrix representation of a single qubit gate controlled by two qubits from a Kronecker product
*/
cplx_t* doublyControlledGateMat(const cplx_t mat[],qubit_t qubits, qubit_t control1, qubit_t control2, qubit_t target) {
    /*
    If either no qubits are involved or the position of the affected qubit lies outside the range of
    available qubits, the function terminates with an error message.
    */
    if (qubits == 0 || control1 >= qubits || control2 >= qubits || target >= qubits || \
        control1 == target || control1 == control2 || control2 == target) {
        printf("Couldn't create controlled single qubit gate!\n");
        return NULL;
    }

    cplx_t* inactiveInactiveControl = ONE;
    cplx_t* activeInactiveControl = ONE;
    cplx_t* inactiveActiveControl = ONE;
    cplx_t* activeActiveControl = ONE;
    cplx_t* tmpInactiveInactive;
    cplx_t* tmpInactiveActive;
    cplx_t* tmpActiveInactive;
    cplx_t* tmpActiveActive;
    for (qubit_t i = 0; i < qubits; ++i) {

        if (i == qubits - target - 1) {
            tmpInactiveInactive = ckronecker(inactiveInactiveControl, IDMAT, POW2(i, dim_t), 2);
            tmpActiveInactive = ckronecker(activeInactiveControl, IDMAT, POW2(i, dim_t), 2);
            tmpInactiveActive = ckronecker(inactiveActiveControl, IDMAT, POW2(i, dim_t), 2);
            tmpActiveActive = ckronecker(activeActiveControl, mat, POW2(i, dim_t), 2);
        }
        else if (i == qubits - control1 - 1) {
            tmpInactiveInactive = ckronecker(inactiveInactiveControl, P0MAT, POW2(i, dim_t), 2);
            tmpActiveInactive = ckronecker(activeInactiveControl, P1MAT, POW2(i, dim_t), 2);
            tmpInactiveActive = ckronecker(inactiveActiveControl, P0MAT, POW2(i, dim_t), 2);
            tmpActiveActive = ckronecker(activeActiveControl, P1MAT, POW2(i, dim_t), 2);
        }
        else if (i == qubits - control2 - 1) {
            tmpInactiveInactive = ckronecker(inactiveInactiveControl, P0MAT, POW2(i, dim_t), 2);
            tmpActiveInactive = ckronecker(activeInactiveControl, P0MAT, POW2(i, dim_t), 2);
            tmpInactiveActive = ckronecker(inactiveActiveControl, P1MAT, POW2(i, dim_t), 2);
            tmpActiveActive = ckronecker(activeActiveControl, P1MAT, POW2(i, dim_t), 2);
        }
        else {
            tmpInactiveInactive = ckronecker(inactiveInactiveControl, IDMAT, POW2(i, dim_t), 2);
            tmpActiveInactive = ckronecker(activeInactiveControl, IDMAT, POW2(i, dim_t), 2);
            tmpInactiveActive = ckronecker(inactiveActiveControl, IDMAT, POW2(i, dim_t), 2);
            tmpActiveActive = ckronecker(activeActiveControl, IDMAT, POW2(i, dim_t), 2);
        }
        if (i != 0) {
            free(inactiveInactiveControl);
            free(activeInactiveControl);
            free(inactiveActiveControl);
            free(activeActiveControl);
        }
        inactiveInactiveControl = tmpInactiveInactive;
        activeInactiveControl = tmpActiveInactive;
        inactiveActiveControl = tmpInactiveActive;
        activeActiveControl = tmpActiveActive;
    }
    cmatAddInPlace(inactiveInactiveControl, activeInactiveControl, POW2(qubits, dim_t));
    cmatAddInPlace(inactiveInactiveControl, inactiveActiveControl, POW2(qubits, dim_t));
    cmatAddInPlace(inactiveInactiveControl, activeActiveControl, POW2(qubits, dim_t));
    free(activeInactiveControl);
    free(inactiveActiveControl);
    free(activeActiveControl);
    return inactiveInactiveControl;
}