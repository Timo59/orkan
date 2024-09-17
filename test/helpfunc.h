#ifndef HELPFUNC_H
#define HELPFUNC_H

/*
 * =====================================================================================================================
 *                                                      includes
 * =====================================================================================================================
 */
#ifndef CPLXUTIL_H
#include "cplxutil.h"
#endif

#ifndef LINALG_H
#include "linalg.h"
#endif

/*
 * =====================================================================================================================
 *                                                  static variables
 * =====================================================================================================================
 */

static const cplx_t IDMAT[4] = {1.0 + 0.0 * I, \
                                0.0 + 0.0 * I, \
                                0.0 + 0.0 * I, \
                                1.0 + 0.0 * I};         // Identity matrix as a 1d array stored as a
                                                        // read-only variable in the data segment

static const cplx_t XMAT[4] = {0.0 + 0.0 * I, \
                               1.0 + 0.0 * I, \
                               1.0 + 0.0 * I, \
                               0.0 + 0.0 * I};          // Pauli-X matrix as a 1d array stored as a
                                                        // read-only variable in the data segment

static const cplx_t YMAT[4] = {0.0 + 0.0 * I, \
                               0.0 - 1.0 * I, \
                               0.0 + 1.0 * I, \
                               0.0 + 0.0 * I};          // Pauli-Y matrix as a 1d arry stored as a
                                                        // read-only variable in the data segment

static const cplx_t ZMAT[4] = {1.0 + 0.0 * I, \
                               0.0 + 0.0 * I, \
                               0.0 + 0.0 * I, \
                               -1.0 + 0.0 * I};         // Pauli-Z matrix as a 1d array stored as a
                                                        // read-only variable in the data segment


static const cplx_t SMAT[4] = {1.0 + 0.0 * I, \
                               0.0 + 0.0 * I, \
                               0.0 + 0.0 * I, \
                               0.0 + 1.0 * I};          // Phase gate matrix as a 1d array stored as
                                                        // a read-only variable in the data segment

static const cplx_t HMAT[4] = {INVSQRT2 + 0.0 * I, \
                               INVSQRT2 + 0.0 * I, \
                               INVSQRT2 + 0.0 * I, \
                               -INVSQRT2 + 0.0 * I};    // Hadamard gate matrix as a 1d array stored
                                                        // as a read-only variable in the data
                                                        //segment

static const cplx_t TMAT[4] = {1.0 + 0.0 * I, \
                               0.0 + 0.0 * I, \
                               0.0 + 0.0 * I, \
                               INVSQRT2 + INVSQRT2 * I};// T-gate matrix as a 1d array stored as a
                                                        // read-only variable in the data segment

static const cplx_t P0MAT[4] = {1.0 + 0.0 * I, \
                                0.0 + 0.0 * I, \
                                0.0 + 0.0 * I, \
                                0.0 + 0.0 * I};         // Projector onto the first computational
                                                        // eigenstate as a 1d array stored as a read
                                                        // only variable in the data segment

static const cplx_t P1MAT[4] = {0.0 + 0.0 * I, \
                                0.0 + 0.0 * I, \
                                0.0 + 0.0 * I, \
                                1.0 + 0.0 * I};         // Projector onto the second computational
                                                        // eigenstate as a 1d array stored as a read
                                                        // only variable in the data segment

static const cplx_t SWAPMAT[16] = {1.0 + 0.0 * I, 0.0 + 0.0 * I, 0.0 + 0.0 * I, 0.0 + 0.0 * I, \
                                   0.0 + 0.0 * I, 0.0 + 0.0 * I, 1.0 + 0.0 * I, 0.0 + 0.0 * I, \
                                   0.0 + 0.0 * I, 1.0 + 0.0 * I, 0.0 + 0.0 * I, 0.0 + 0.0 * I, \
                                   0.0 + 0.0 * I, 0.0 + 0.0 * I, 0.0 + 0.0 * I, 1.0 + 0.0 * I};
                                						// SWAP matrix of two neighboring qubits as
                                                        // a 1d array stored as a read-only variable
                                                        // in the data segment

static cplx_t ONE[1] = {1.0 + 0.0 * I};                 // 1d array holding the value 1 as a double
                                                        // complex

/*
 * =====================================================================================================================
 *                                                  Helper functions
 * =====================================================================================================================
 */

/*
This function generates the identity matrix on the number of specified qubits

Input:
    qubit_t qubits: unsigned integer holding the number of qubits involved in a state

Output:
    Double complex pointer to the array holding the identity matrix on qubits qubits

*/
cplx_t* identityMat(qubit_t qubits) {
    cplx_t* result = ONE;                                               // double complex pointer
                                                                        // used to store the final
                                                                        // Kronecker product of
                                                                        // of Paulis on the qubits

    cplx_t* tmp;                                                        // double complex pointer
                                                                        // hodling each intermediate
                                                                        // Kronecker product of
                                                                        // identities                      
    /*
    Starting at the first qubit add an identity matrix to the Kronecker product until all qubits
    contribute.
    */
    for (qubit_t i = 0; i < qubits; ++i) {
        tmp = ckronecker(result, IDMAT, POW2(i, dim_t), 2);             // the Kronecker product of
                                                                        // the current intermediate
                                                                        // matrix with the next 
                                                                        // identity matrix is stored
                                                                        // to tmp
        /*
        Except for the first iteration step, result is freed.
        */
        if (i != 0) {
            free(result);
        }
        /*
        result is now pointing to the same adress tmp is pointing to.
        */
        result = tmp;                                                   
    }
    return result;
}

/*
This function generates a set of vectors used to unittest the elementary operations on pure states.

Input:
    qubit_t qubits: the number of qubits involved in a state

Output:
    The output is a pointer to an array of double complex pointers each representing a state vector.
    In total the comprises 2^{qubits} + 1 different states, namely all computational basis states
    plus one vector with alternating 1 and im in the entries.
*/
cplx_t** generateTestVectors(qubit_t qubits) {
    dim_t dim = POW2(qubits, dim_t);                                    // unsigned integer holding
                                                                        // the dimension of the 
                                                                        // Hilbert space according
                                                                        // to the number of qubits

    cplx_t** vectors = calloc(dim + 1, sizeof(cplx_t*));     // pointer to an array of
                                                                        // dim + 1 double complex
                                                                        // pointers holding the 
                                                                        // test vectors

    vectors[dim] = calloc(dim, sizeof(*(vectors[dim])));     // initialise the last entry
                                                                        // of the vector array as a
                                                                        // pointer to a 1d array
                                                                        // with dim entries
                                                                        // representing the
                                                                        // additional test vector
    
    /*
    Starting at zero, the integer i runs through all integers up to dim - 1, initialising the
    corresponding basis vector in the z basis.  
    */                                                                    
    for (dim_t i = 0; i < dim; ++i) {
        
        vectors[i] = calloc(dim, sizeof(*(vectors[i])));     // initialise the i-th entry
                                                                        // of the vectors array as a
                                                                        // pointer to a 1d array of
                                                                        // size dim with all entries
                                                                        // set to zero

        vectors[i][i] = 1.0 + 0.0 * I;                                  // set the i-th entry of the
                                                                        // i-th array pointed to
                                                                        // from the vectors array to
                                                                        // one representing the i-th
                                                                        // basis vector

        vectors[dim][i] = (i % 2) ? 0.5 + 0.0 * I : 0.0 + 0.5 * I;      // set the i-th entry of the
                                                                        // additional test vector to
                                                                        // 1 or im, depending on
                                                                        // whether i is even or odd
    }
    return vectors;                                                     // return the pointer to the
                                                                        // array of pointers to the
                                                                        // test vectors
}

/*
This function frees the allocated memory for the set of test vectors used to unittest the elementary
operations

Input:
    cplx_t** vectors: pointer to an array of pointers holding the memory adresses of the test
    vectors
    qubit_t qubits: the number of qubits involved in a state
Output:
    There is no output value, instead the allocated memory for the pointers to the test vectors as
    well as the allocated memory for the test vectors themeselves are freed.

*/
void freeTestVectors(cplx_t** vectors, qubit_t qubits) {
    dim_t dim = POW2(qubits, dim_t);        // unsigned integer holding the dimension of the
                                            // Hilbert space according to the number of qubits
    /*
    Starting at zero, delete the test vector stored at position of the integer i and increase i by
    one. After dim + 1 iterations all basis vectors and the additional vector are deleted, aka their
    memory is freed.
    */   
    for (dim_t i = 0; i < dim + 1; ++i) {
        free(vectors[i]);                   // free the allocated memory for the i-th test vector
    }
    free(vectors);                          // free the allocated memory for the pointers to the
                                            // test vectors
}

/*
 * This function returns an array of all possible Pauli strings on the specified number of qubits
 */
pauli_t* allPauliStrings(qubit_t qubits) {
    dim_t stringc =  (1 << (2 * qubits));
    pauli_t* result = calloc(qubits * stringc, sizeof(pauli_t));

    for (dim_t i = 1; i < stringc; ++i) {
        for (qubit_t j = 0; j < qubits; ++j) {
            result[i * qubits + j] = (i / (1 << (2 * j))) % 4;
        }
    }
    return result;
}

pauli_t* allPauliStringsDiag(qubit_t qubits) {
    dim_t stringc =  (1 << qubits);
    pauli_t* result = calloc(qubits * stringc, sizeof(pauli_t));

    for (dim_t i = 1; i < stringc; ++i) {
        for (qubit_t j = 0; j < qubits; ++j) {
            result[i * qubits + j] = 3 * ((i / (1 << j)) % 2);
        }
    }
    return result;
}

pauli_t* allPauliStringsX(qubit_t qubits) {
    dim_t stringc =  (1 << qubits);
    pauli_t* result = calloc(qubits * stringc, sizeof(pauli_t));

    for (dim_t i = 1; i < stringc; ++i) {
        for (qubit_t j = 0; j < qubits; ++j) {
            result[i * qubits + j] = (i / (1 << j)) % 2;
        }
    }
    return result;
}

pauli_t* allPauliStringsY(qubit_t qubits) {
    dim_t stringc =  (1 << qubits);
    pauli_t* result = calloc(qubits * stringc, sizeof(pauli_t));

    for (dim_t i = 1; i < stringc; ++i) {
        for (qubit_t j = 0; j < qubits; ++j) {
            result[i * qubits + j] = 2 * ((i / (1 << j)) % 2);
        }
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
            tmpActive = ckronecker(active, mat, POW2(i, dim_t), 2);         // Kronecker product of
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
        printf("Couldn't create a neighboring swap gate for qubit %d!\n", target);
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
        else {                                      // Otherwise
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
        cplx_t* neighborSwap = neighboringSwapGateMat(qubits, i - 1);       // Calculate the matrix for swapping with
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
This function computes the matrix representation of any single qubit gate with conditional execution
depending on the state of another qubit.

Input:
    const cplx_t math[]: read only complex double of the gate applied conditioned on control in
    matrix representation
    qubit_t qubits: the number of qubits involved in a state
    qubit_t control1: unsigned integer holding the position of the first control qubit from the
    right
    qubit_t control2: unsigned integer holding the position of the second control qubit from the
    right
    qubit_t target: unsigned integer holding the position of the target qubit from the right

Output:
    Complex double pointer to the array holding the sum of two kronecker products corresponding to
    the cases where the control qubit is active and inactive, respectively. 
*/
cplx_t* doublyControlledGateMat(const cplx_t mat[], qubit_t qubits, qubit_t control1, \
                                qubit_t control2, qubit_t target) {
    /*
    If either no qubits are involved or the position of the affected qubit lies outside the range of
    available qubits, the function terminates with an error message.
    */
    if (qubits == 0 || control1 >= qubits || control2 >= qubits || target >= qubits || \
        control1 == target || control1 == control2 || control2 == target) {
        printf("Couldn't create controlled single qubit gate!\n");
        return NULL;
    }

    cplx_t* inactiveInactiveControl = ONE;      // double complex pointer to the array holding the
                                                // matrix corresponding to the case where both
                                                // control qubits are in the zero state

    cplx_t* activeInactiveControl = ONE;        // double complex pointer to the array holding the
                                                // matrix corresponding to the case where the first
                                                // qubit is in the one state but the second is in
                                                // the zero state

    cplx_t* inactiveActiveControl = ONE;        // double complex pointer to the array holding the
                                                // matrix corresponding to the case where the first
                                                // qubit is in the zero state but the second is in
                                                // the one state     

    cplx_t* activeActiveControl = ONE;          // double complex pointer to the array holding the
                                                // matrix corresponding to the case where both 
                                                // control qubits are in the one state

    cplx_t* tmpInactiveInactive;                // double complex pointer holding each intermediate
                                                // Kronecker products of matrices where both control
                                                // qubits are inactive

    cplx_t* tmpActiveInactive;                  // double complex pointer holding each intermediate
                                                // Kronecker products of matrices where the first
                                                // and the second qubit is active and inactive,
                                                // respectively

    cplx_t* tmpInactiveActive;                  // double complex pointer holding each intermediate
                                                // Kronecker products of matrices where the first
                                                // and the second qubit is inactive and active,
                                                // respectively

    cplx_t* tmpActiveActive;                    // double complex pointer holding each intermediate
                                                // Kronecker products of matrices where both control
                                                // qubits are active

    /*
    Starting at the leftmost position, identity matrices are appended to Kronecker product of the
    all intermediate matrices until the index reaches either one of the control qubits or the target
    qubit. In the former case, projectors onto the second and first computational basis states are
    appended according to the state of the respective control qubit. In the latter case, the
    specified gate is appended only to the all active case and the identity matrix otherwise. 
    */
    for (qubit_t i = 0; i < qubits; ++i) {
        /*
        If the current integer corresponds to the position of the target qubit from the left, the
        kronecker product with the single qubit gate is computed and its output pointer is
        identified with tmpActiveActive. Moreover, the kronecker product with the identity matrix is
        computed for all other cases and its output pointer is identified with the respective
        temporary
        */
        if (i == qubits - target - 1) {
            tmpInactiveInactive = ckronecker(inactiveInactiveControl, IDMAT, POW2(i, dim_t), 2);
                                                // Kronecker product of the current matrix of both
                                                // qubits inactive with another identity matrix;
                                                // pointer address is stored to the corresponding
                                                // temporary                             

            tmpActiveInactive = ckronecker(activeInactiveControl, IDMAT, POW2(i, dim_t), 2);
                                                // Kronecker product of the current matrix of the
                                                // first qubit active and the second inactive with
                                                // another identity matrix; pointer address is
                                                // stored to the corresponding temporary

            tmpInactiveActive = ckronecker(inactiveActiveControl, IDMAT, POW2(i, dim_t), 2);
                                                // Kronecker product of the current matrix of the
                                                // first qubit inactive and the second active with
                                                // another identity matrix; pointer address is
                                                // stored to the corresponding temporary
            
            tmpActiveActive = ckronecker(activeActiveControl, mat, POW2(i, dim_t), 2);      
                                                // Kronecker product of the current matrix of both
                                                // qubits active with the specified gate matrix;
                                                // pointer address is stored to the corresponding
                                                // temporary
        } else if (i == qubits - control1 - 1) {
        /*
        If otherwise, the current integer corresponds to the position of the first control qubit
        from the left, the projector onto the first and second computational basis state is appended
        to those temporaries corresponding to the first control qubit inactive and active,
        respectively. 
        */
            tmpInactiveInactive = ckronecker(inactiveInactiveControl, P0MAT, POW2(i, dim_t), 2);
                                                // Kronecker product of the current matrix of both
                                                // qubits inactive with the projector onto the zero
                                                // state since the first control qubit is inactive

            tmpActiveInactive = ckronecker(activeInactiveControl, P1MAT, POW2(i, dim_t), 2);
                                                // Kronecker product of the current matrix of the
                                                // first control qubit active and the second control
                                                // qubit inactive with the projector onto the one
                                                // state since the first control qubit is active

            tmpInactiveActive = ckronecker(inactiveActiveControl, P0MAT, POW2(i, dim_t), 2);
                                                // Kronecker product of the current matrix of the
                                                // first control qubit inactive and the second control
                                                // qubit active with the projector onto the zero
                                                // state since the first control qubit is inactive
            
            tmpActiveActive = ckronecker(activeActiveControl, P1MAT, POW2(i, dim_t), 2);
                                                // Kronecker product of the current matrix of both
                                                // qubits active with the projector onto the one 
                                                // state since the first control qubit is active
        } else if (i == qubits - control2 - 1) {
        /*
        If otherwise, the current integer corresponds to the position of the second control qubit
        from the left, the projector onto the first and second computational basis state is appended
        to those temporaries corresponding to the first control qubit inactive and active,
        respectively. 
        */
            tmpInactiveInactive = ckronecker(inactiveInactiveControl, P0MAT, POW2(i, dim_t), 2);   
                                                // Kronecker product of the current matrix of both
                                                // qubits inactive with the projector onto the zero
                                                // state since the second control qubit is inactive
            
            tmpActiveInactive = ckronecker(activeInactiveControl, P0MAT, POW2(i, dim_t), 2);
                                                // Kronecker product of the current matrix of the
                                                // first control qubit active and the second control
                                                // qubit inactive with the projector onto the zero
                                                // state since the second control qubit isinactive

            tmpInactiveActive = ckronecker(inactiveActiveControl, P1MAT, POW2(i, dim_t), 2);
                                                // Kronecker product of the current matrix of the
                                                // first control qubit active and the second control
                                                // qubit active with the projector onto the
                                                // one state since the second control qubit is
                                                // active

            tmpActiveActive = ckronecker(activeActiveControl, P1MAT, POW2(i, dim_t), 2);       
                                                // Kronecker product of the current matrix of both
                                                // qubits inactive with the projector onto the one
                                                // state since the second control qubit is active
        } else {
        /*
        If the current integer corresponds to neither position of the control qubits nor the target
        qubit's position from the left, an identity matrix is appended to all intermediate Kronecker
        products.
        */
            tmpInactiveInactive = ckronecker(inactiveInactiveControl, IDMAT, POW2(i, dim_t), 2);   
                                                // Kronecker product of the current matrix with both
                                                // control qubits inactive with another identity
                                                // matrix

            tmpActiveInactive = ckronecker(activeInactiveControl, IDMAT, POW2(i, dim_t), 2);
                                                // Kronecker product of the current matrix of the
                                                // first control qubit active and the second control
                                                // qubit inactive with another identity matrix

            tmpInactiveActive = ckronecker(inactiveActiveControl, IDMAT, POW2(i, dim_t), 2);
                                                // Kronecker product of the current matrix of the
                                                // first control qubit active and the second control
                                                // qubit active with another identity matrix
            
            tmpActiveActive = ckronecker(activeActiveControl, IDMAT, POW2(i, dim_t), 2);       
                                                // Kronecker product of the current matrix with both
                                                // control qubits active with another identity
                                                // matrix
        }
        /*
        Except for the first iterational step, the pointers to both current result matrices, active
        and inactive control qubit are freed.
        */
        if (i != 0) {
            free(inactiveInactiveControl);
            free(activeInactiveControl);
            free(inactiveActiveControl);
            free(activeActiveControl);
        }
        inactiveInactiveControl = tmpInactiveInactive; 
                                            // inactiveControl is now pointing to the outcome of the
                                            // Kronecker product in tmpInactive

        activeInactiveControl = tmpActiveInactive;

        inactiveActiveControl = tmpInactiveActive;

        activeActiveControl = tmpActiveActive;
                                            // activeControl is now pointing to the outcome of the
                                            // Kronecker product in tmpActive
    }
    /*
    When all position are executed the matrices in inactiveControl and activeControl are merged via
    an in place addition to inactiveControl.
    */
    cmatAddInPlace(inactiveInactiveControl, activeInactiveControl, POW2(qubits, dim_t));
    cmatAddInPlace(inactiveInactiveControl, inactiveActiveControl, POW2(qubits, dim_t));
    cmatAddInPlace(inactiveInactiveControl, activeActiveControl, POW2(qubits, dim_t));
    free(activeInactiveControl);
    free(inactiveActiveControl);
    free(activeActiveControl);
    return inactiveInactiveControl;
}

/*
 * =====================================================================================================================
 *                                                  Pauli string matrices
 * =====================================================================================================================
 */

/*
 * This function returns the matrix representation of a single Pauli string.
 *
 * Input:
 *      pauli_t paulistr[]      Array holding the Pauli string's components
 *      qubit_t qubits          Number of qubits; i.e., number of components
 *
 * Output:
 *      The Kronecker product of qubit Pauli matrices stored in a dynamical memory block in row major form
 */
cplx_t* pauliStrMat(const pauli_t paulistr[], qubit_t qubits) {
    cplx_t* result = ONE;       // Output initialized to the array {1.0+0.0I}
    cplx_t* tmp;                // Temporary array holding each intermediate Kronecker product

    /*
     * Iterate through the Pauli string from the right and depending on the Pauli matrix store the current Kronecker
     * product in tmp.
     */
    for (qubit_t i = 0; i < qubits; ++i) {
        switch (paulistr[qubits - i - 1]) {
            case ID: {
                tmp = ckronecker(result, IDMAT, POW2(i, dim_t), 2);
                break;
            }
            case X: {
                tmp = ckronecker(result, XMAT, POW2(i, dim_t), 2);
                break;
            }
            case Y: {
                tmp = ckronecker(result, YMAT, POW2(i, dim_t), 2);
                break;
            }
            case Z: {
                tmp = ckronecker(result, ZMAT, POW2(i, dim_t), 2);
                break;
            }
            default: {
                printf("Error in constructing Pauli string matrix.\n");
                exit(1);
            }
        }

        if (i != 0) {           // Except for the first iteration(no memory allocated), free the memory allocated for
            free(result);       // result
        }
        result = tmp;           // Change result's address to that of tmp
    }
    return result;
}

/*
 * This function returns the matrix representation of an observable built from Pauli strings; i.e., a real linear
 * combination of Pauli strings.
 *
 * Input:
 *      pauli_t comps[]         Concatenation of the observable's Pauli strings
 *      double coeffs[]         Array of the coefficients for each observable's Pauli string
 *      complength_t length     Number of constituents; i.e, number of Pauli strings
 *      qubit_t qubits          Number of qubits; i.e., number of components for each Pauli string
 *
 * Output:
 *      The sum of Kronecker products of qubit Pauli matrices weighted with their respective coefficients stored in a
 *      dynamical memory block in row major form
 */
cplx_t* pauliObsMat(const pauli_t comps[], const double coeffs[], complength_t length, qubit_t qubits) {
    dim_t dim = POW2(qubits, dim_t);                                            // Hilbert space dimension
    cplx_t* result = (cplx_t*) calloc(dim*dim, sizeof(cplx_t));     // The resulting square matrix
    cplx_t* tmp;                                                                // Temporary holding each intermediate

    /*
     * Iterate the Pauli strings, create their matrices and add them times the coefficient to the resulting matrix
     */
    for(complength_t i = 0; i < length; ++i) {
        tmp = pauliStrMat(comps + (i * qubits), qubits);
        cscalarMatMulInPlace(coeffs[i], tmp, dim);
        cmatAddInPlace(result, tmp, dim);
        free(tmp);
    }

    return result;
}

/*
 * This function returns the diagonal of an observable built only from diagonal Pauli strings; i.e., a real linear
 * combination of diagonal Pauli strings.
 *
 * Input:
 *      pauli_t comps[]         Concatenation of the observable's Pauli strings
 *      double coeffs[]         Array of the coefficients for each observable's Pauli string
 *      complength_t length     Number of constituents; i.e, number of Pauli strings
 *      qubit_t qubits          Number of qubits; i.e., number of components for each Pauli string
 *
 * Output:
 *      The diagonal entries of the observable matrix generated by pauliObservableMat stored in a
 *      dynamical memory block
 */
double* pauliObsMatDiag(const pauli_t comps[], const double coeffs[], complength_t length, qubit_t qubits) {
    dim_t dim = POW2(qubits, dim_t);                                    // Hilbert space dimension
    double* result = (double*) malloc(dim * sizeof(double));        // The resulting array of diagonal entries
    cplx_t* tmp = pauliObsMat(comps, coeffs, length, qubits);    // Temporary holding the observable matrix

    for (dim_t i = 0; i <dim; ++i) {                                    // Iterate through diagonal entries and store
        result[i] = creal(tmp[i*dim + i]);                              // them to result
    }
    free(tmp);
    return result;
}

/*
 * This function returns the matrix representation of the evolution with a Pauli string by an angle.
 *
 * Input:
 *      pauli_t paulistr[]      Array holding the Pauli string's components
 *      double angle            Angle of the evolution; i.e., factor in the exponential expression
 *      qubit_t qubits          Number of qubits; i.e., number of components
 *
 * Output:
 *      The matrix exponential of the Pauli string matrix generated by pauliStringMat, calculated with cosine and sine
 *      since any Pauli string is involutive, stored in a dynamical memory block in row major form
 */

cplx_t* expPauliStrMat(const pauli_t paulistr[], double angle, qubit_t qubits) {
    dim_t dim = POW2(qubits, dim_t);                            // Hilbert space dimension
    cplx_t* result = identityMat(qubits);                       // The resulting matrix initialized to the all identity
                                                                // matrix

    cscalarMatMulInPlace(cos(angle), result, dim);      // Multiply the all identity matrix with cos(angle)
    cplx_t* tmp = pauliStrMat(paulistr, qubits);             // tmp holds the Pauli string matrix
    cscalarMatMulInPlace(-I * sin(angle), tmp, dim);    // Multiply the Pauli string matrix by -I * sin(angle)
    cmatAddInPlace(result, tmp, dim);                   // Add tmp entry-wise to result

    free(tmp);
    return result;
}

/*
 * This function returns the matrix representation of the evolution with an observable by an angle in Lie-Trotter
 * approximation to first order.
 *
 * Input:
 *      pauli_t comps[]         Concatenation of the observable's Pauli strings
 *      double coeffs[]         Array of the coefficients for each observable's Pauli string
 *      complength_t length     Number of constituents; i.e, number of Pauli strings
 *      double angle            Angle of the evolution; i.e., factor in the exponential expression
 *      qubit_t qubits          Number of qubits; i.e., number of components for each Pauli string
 *
 * Output:
 *      The matrix exponential of the observable times (-I * angle), calculated as the product of matrix exponentials
 *      from the observable's Pauli strings(assuming they commute) generated by expPauliStringMat, stored in a dynamical
 *      memory block in row major form
 */

cplx_t* expPauliObsMatTrotter(const pauli_t comps[],
                              const double coeffs[],
                              complength_t length,
                              double angle,
                              qubit_t qubits)
{
    dim_t dim = POW2(qubits, dim_t);            // Hilbert space dimension
    cplx_t* result = identityMat(qubits);       // The resulting matrix
    cplx_t* tmp;                                // Temporary holding each Pauli string's matrix exponential

    /*
     * Iterate through the observable's Pauli strings from the right, store their matrix exponential to tmp and multiply
     * it to the resulting matrix
     */
    for (complength_t i = 0; i < length; ++i) {
        tmp = expPauliStrMat(comps + ((length - i - 1) * qubits),
                             coeffs[length - i - 1] * angle,
                             qubits);
        cmatMulInPlace(result, tmp, dim);
        free(tmp);
    }

    return result;
}

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
                       double* par
                       ) {
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

    double result = creal(cinnerProduct(bra, ket, dim));

    free(bra);
    free(ket);
    free(observableMat);
    for (depth_t i = 0; i < circdepth; ++i) {
        free(evoOpsMat[i]);
    }
    free(evoOpsMat);

    return result;
}

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
                         double epsilon
                         )
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

#endif /* HELPFUNC_H */