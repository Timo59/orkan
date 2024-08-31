/*
 * =================================================================================================
 *                                              includes
 * =================================================================================================
 */

//#include "cplxutil.h"
#include "linalg.h"
#include <stdio.h>
//#include "unity.h"

/*
 * =================================================================================================
 *                                              type definitions
 * =================================================================================================
 */

typedef long                 squbit_t;

/*
 * =================================================================================================
 *                                              static variables
 * =================================================================================================
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
 * =================================================================================================
 *                                              helper functions
 * =================================================================================================
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

/*
This function computes the matrix of applying one single qubit gate to the qubit at position in a
system of qubits qubits.

Input:
    const cplex_t mat[]: matrix representation of the gate applied as complex double read-only
    1d array
    qubit_t qubits: the number of qubits involved in a state
    qubit_t position: the position of the qubit the gate acts on from the right in a qubit string

Output:
    Complex double pointer to the array holding the kronecker product of qubits - 1 identity matrix
    with one single qubit gate at the position specified by position.

*/
cplx_t* singleQubitGateMat(const cplx_t mat[], qubit_t qubits, qubit_t position) {
    /*
    If either no qubits are involved or the position of the affected qubit lies outside the range of
    available qubits, the function terminates with an error message.
    */
    if (qubits == 0 || position >= qubits) {
        printf("Couldn't create tensorized single qubit gate!\n");
        return NULL;
    }

    cplx_t* result = ONE;                                       // double complex pointer holding
                                                                // the resulting state vector
                                                                // initialized as a double complex
                                                                // array holding only a 1

    cplx_t* tmp;                                                // double complex pointer holding
                                                                // a temporary state vector

    /*
    Starting a the leftmost position, the kronecker with the next matrix depending on the position
    of the single qubit gate is computed until all qubits have contributed.
    */
    for (qubit_t i = 0; i < qubits; ++i) {
        /*
        If the current integer corresponds to the position of the single qubit gate from the left,
        the kronecker product of the matrix stored in result with the single qubit gate is computed.
        The output, a pointer to the resulting array, is stored in tmp. 
        */
        if (i == qubits - position - 1) {
            tmp = ckronecker(result, mat, POW2(i, dim_t), 2);   // kronecker product of the current
                                                                // result with 2^{i} rows and
                                                                // columns with the single qubit
                                                                // gate in matrix representation
                                                                // Output is a pointer to the array
                                                                // which is stored in tmp 
        }
        /*
        Otherwise, the kronecker product of the matrix stored in result with the identity is computed.
        The output, a pointer to the resulting array, is stored in tmp. 
        */
        else {
            tmp = ckronecker(result, IDMAT, POW2(i, dim_t), 2); // kronecker product of the current
                                                                // result with 2^{i} rows and
                                                                // columns with the identity matrix
                                                                // Output is a pointer to the array
                                                                // which is stored in tmp
        }
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

cplx_t* RXGateMat(qubit_t qubits, qubit_t position, double angle) {
    angle /= 2.;
    cplx_t mat[4] = {cos(angle) + 0.0 * I, \
                     0.0 - sin(angle) * I, \
                     0.0 - sin(angle) * I, \
                     cos(angle) + 0.0 * I};
    return singleQubitGateMat(mat, qubits, position);
}

cplx_t* RYGateMat(qubit_t qubits, qubit_t position, double angle) {
    angle /= 2.;
    cplx_t mat[4] = {cos(angle) + 0.0 * I, \
                     -sin(angle) + 0.0 * I, \
                     sin(angle) + 0.0 * I, \
                     cos(angle) + 0.0 * I};
    return singleQubitGateMat(mat, qubits, position);
}

cplx_t* RZGateMat(qubit_t qubits, qubit_t position, double angle) {
    angle /= 2.;
    cplx_t mat[4] = {cos(angle) - sin(angle) * I, \
                     0.0 + 0.0 * I, \
                     0.0 + 0.0 * I, \
                     cos(angle) + sin(angle) * I};
    return singleQubitGateMat(mat, qubits, position);
}

/*
This function computes the matrix representation of any single qubit gate with conditional execution
depending on the state of another qubit.

Input:
    const cplx_t math[]: read only complex double of the gate applied conditioned on control in
    matrix representation
    qubit_t qubits: the number of qubits involved in a state
    qubit_t control: unsigned integer holding the position of the control qubit from the right
    qubit_t target: unsigned integer holding the position of the target qubit from the right

Output:
    Complex double pointer to the array holding the sum of two kronecker products corresponding to
    the cases where the control qubit is active and inactive, respectively. 
*/
cplx_t* controlledGateMat(const cplx_t mat[], qubit_t qubits, qubit_t control, qubit_t target) {
    /*
    If either no qubits are involved or the position of the affected qubit lies outside the range of
    available qubits, the function terminates with an error message.
    */
    if (qubits == 0 || control >= qubits || target >= qubits || control == target) {
        printf("Couldn't create controlled single qubit gate!\n");
        return NULL;
    }

    cplx_t* inactiveControl = ONE;      // double complex pointer to the array holding the matrix 
                                        // corresponding to the case where the control qubit is
                                        // in the zero state and therefore inactive
    cplx_t* activeControl = ONE;        // double complex pointer to the array holding the matrix 
                                        // corresponding to the case where the control qubit is
                                        // in the one state and therefore active
    cplx_t* tmpInactive;                // double complex pointer holding the temporary matrix of
                                        // inactive control qubit
    cplx_t* tmpActive;                  // double complex pointer holding the temporary matrix of
                                        // sactive control qubit

    /*
    Starting at the leftmost position, identity matrices are appended to Kronecker product of the
    active as well as the inactive intermediate matrices until the index reaches either the control
    or the target qubit. In that case the projector onto the second and first computational basis
    states, respectively are appended or the specified matrix and the identity, respectively are
    appended. After that identity matrices keep being added until the number of constituents equals 
    the number qubits. 
    */
    for (qubit_t i = 0; i < qubits; ++i) {
        /*
        If the current integer corresponds to the position of the targeted qubit from the left, the
        kronecker product with the single qubit gate is computed and its output pointer is
        identified with tmpActive. Moreover, the kronecker product with the identity matrix is
        computed and its output pointer is identified with tmpInactive.
        */
        if (i == qubits - target - 1) {
            tmpInactive = ckronecker(inactiveControl, IDMAT, POW2(i, dim_t), 2);// kronecker product
                                                                                // of the current
                                                                                // inactive control
                                                                                // with 2^{i} rows
                                                                                // and columns with
                                                                                // the identity
                                                                                // matrix.
                                                                                // Output is pointer
                                                                                // to the array
                                                                                // which is stored
                                                                                // in tmpInactive
            
            tmpActive = ckronecker(activeControl, mat, POW2(i, dim_t), 2);      // kronecker product
                                                                                // of the current
                                                                                // result with 2^{i}
                                                                                // rows and columns
                                                                                // with the single
                                                                                // qubit gate.
                                                                                // Output is pointer
                                                                                // to the array
                                                                                // which is stored
                                                                                // in tmpActive
        }
        /*
        If otherwise, the current integer corresponds to the position of the control qubit from the
        left, the kronecker product with either the projector on the first and second computational
        eigenstate is computed. 
        */
        else if (i == qubits - control - 1) {
            tmpInactive = ckronecker(inactiveControl, P0MAT, POW2(i, dim_t), 2);// kronecker product
                                                                                // of the current
                                                                                // inactive control
                                                                                // with 2^{i} rows
                                                                                // and columns with
                                                                                // the projector
                                                                                // onto the first
                                                                                // computational
                                                                                // basis state.
                                                                                // Output is pointer
                                                                                // to the array
                                                                                // which is stored
                                                                                // in tmpInactive
            
            tmpActive = ckronecker(activeControl, P1MAT, POW2(i, dim_t), 2);    // kronecker product
                                                                                // of the current
                                                                                // active control
                                                                                // with 2^{i} rows
                                                                                // and columns with
                                                                                // the projector
                                                                                // onto the second
                                                                                // computational
                                                                                // basis state.
                                                                                // Output is pointer
                                                                                // to the array
                                                                                // which is stored
                                                                                // in tmpActive
        }
        /*
        If the current integer neither corresponds the position of the control nor the target qubit
        from the left both, the inactive and the active control matrices, are kronecker multiplied
        with the identity matrix.
        */
        else {
            tmpInactive = ckronecker(inactiveControl, IDMAT, POW2(i, dim_t), 2);// kronecker product
                                                                                // of the current
                                                                                // inactive control
                                                                                // with 2^{i} rows
                                                                                // and columns with
                                                                                // the identity
                                                                                // matrix.
                                                                                // Output is pointer
                                                                                // to the array
                                                                                // which is stored
                                                                                // in tmpInactive
            
            tmpActive = ckronecker(activeControl, IDMAT, POW2(i, dim_t), 2);    // kronecker product
                                                                                // of the current
                                                                                // active control
                                                                                // with 2^{i} rows
                                                                                // and columns with
                                                                                // the identity
                                                                                // matrix.
                                                                                // Output is pointer
                                                                                // to the array
                                                                                // which is stored
                                                                                // in tmpActive
        }
        /*
        Except for the first iterational step, the pointers to both current result matrices, active
        and inactive control qubit are freed.
        */
        if (i != 0) {
            free(inactiveControl);
            free(activeControl);
        }
        inactiveControl = tmpInactive;      // inactiveControl is now pointing to the outcome of the
                                            // Kronecker product in tmpInactive
        activeControl = tmpActive;          // activeControl is now pointing to the outcome of the
                                            // Kronecker product in tmpActive
    }
    /*
    When all position are executed the matrices in inactiveControl and activeControl are merged via
    an in place addition to inactiveControl.
    */
    cmatAddInPlace(inactiveControl, activeControl, POW2(qubits, dim_t));
    free(activeControl);
    return inactiveControl;
}

/*
This function computes the matrix corresponding to the swap of the specified qubit with its
neighbor to the left.

Input:
    qubit_t qubits: the number of qubits involved in a state
    qubit_t qubit1: unsigned integer corresponding to the one qubit to be swapped with its left
    neighbor

Output:
    Double complex pointer to the matrix representing the swap of the specified qubit with its
    neighbor to the left.
*/
cplx_t* neighboringSwapGateMat(squbit_t qubits, squbit_t qubit) {
    /*
    If either there are no qubits at all, the qubit to be swapped is the outmost qubit to the left
    (the last in qhipster notation) or even outside the rage of qubits, return an error.
    */
    if (qubits == 0 || qubit >= qubits - 1) {
        printf("Couldn't create a neighboring swap gate for qubit %ld!\n", qubit);
        return NULL;
    }
    cplx_t* result = ONE;       // double complex pointer to the 1d array representing the 
                                // all intermediate steps of the Kronecker product matrix;
                                // initialised as complex 1
    cplx_t* tmp;                // double complex pointer used for temporarily holding the adress
                                // for the intermediate Kronecker products
    /*
    Starting at the leftmost position, identity matrices are appended to the Kronecker product until
    the index reaches the position ahead of the the specified. Then, the 4x4 SWAP matrix is appended
    while the index is shifted by to steps instead of one. After that identity matrices keep being
    added until the number of constituents equals the number qubits.
    */
    for (squbit_t i = 0; i < qubits; ++i){
        /*
        If the index matches the position of the qubit that is to be swapped with the specified,
        i.e. the one right to it (the one of smaller index in qhipster notation) one appends the 4x4
        SWAP gate to the Kronecker product. In this case the index needs to be incremented by one
        manually to skip the qubit that took part in the SWAP.
        */
        if (i + 1 == qubits - qubit - 1) {
            tmp = ckronecker(result, SWAPMAT, POW2(i, dim_t), 4);   // the output pointer to the
                                                                    // intermediate Kronecker
                                                                    // product is stored to tmp

            ++i;                                                    // index is raised by one since
                                                                    // the SWAP operates also on the
                                                                    // next qubit

            /*
            Except for the first iterational step, i.e., the index is one since it was incremented
            after appending the SWAP gate, free the pointer to result.
            */
            if (i != 1) {
                free(result);
            }
        }
        /*
        Otherwise, another identity matrix is appended to the Kronecker product.
        */
        else {
            tmp = ckronecker(result, IDMAT, POW2(i, dim_t), 2);     // the output pointer to the
                                                                    // intermediate Kronecker
                                                                    // product is stored to tmp

            /*
            Except for the first iterational step, free the pointer to result.
            */
            if (i != 0) {
                free(result);
            }
        }

        result = tmp;       // result is now pointing to the outcome of the intermediate Kronecker
                            // product matrix

    }
    return result;
}

/*
This function computes the matrix for the swap of two arbitrary qubits by concatenating swaps of 
neighboring qubits.

Input:
    qubit_t qubits: the number of qubits involved in a state
    qubit_t qubit1: unsigned integer corresponding to the one qubit to be swapped
    qubit_t qubit2: unsigned integer corresponding to the other qubit to be swapped

Output:
    Double complex pointer to the matrix product of neighboring swaps representing the matrix of
    swapping two arbitrary qubits.
*/
cplx_t* swapGateMat(squbit_t qubits, squbit_t qubit1, squbit_t qubit2) {
    dim_t dimension = POW2(qubits, dim_t);                                  // unsigned integer
                                                                            // holding the dimension
    squbit_t left_qubit = MAX(qubit1, qubit2);                              // unsigned integer 
                                                                            // representing the
                                                                            // leftmost (lower index
                                                                            // in qhipster) qubit in
                                                                            // the SWAP

    squbit_t right_qubit = MIN(qubit1, qubit2);                             // unsigned integer
                                                                            // representing the
                                                                            // rightmost (larger
                                                                            // index in qhipster)
                                                                            // qubit in the SWAP

    cplx_t* result = neighboringSwapGateMat(qubits, left_qubit - 1);        // double complex
                                                                            // pointer to the
                                                                            // intermediate product
                                                                            // matrix; initialised
                                                                            // as the leftmost (last
                                                                            // in qhipster) swap

    cplx_t* tmp;                                                            // double complex
                                                                            // pointer to the
                                                                            // temporary 
                                                                            // intermediate product
                                                                            // matrices

    cplx_t* neighborSwap;                                                  // double complex
                                                                            // pointer to the
                                                                            // neighboring swap
                                                                            // matrices
    /*
    Starting at the position right to the left qubit (larger index in qhipster), multiply the swap
    between the next neighbors to the right from the left and the right to the intermediate matrix
    until arriving at the righmost qubit.
    */
    for (squbit_t i = left_qubit - 2; i >= right_qubit; --i) {

        neighborSwap = neighboringSwapGateMat(qubits, i);       // compute the matrix for the next
                                                                // neighboring swap

        tmp = cmatMul(neighborSwap, result, dimension);         // multiply the next neighboring
                                                                // swap to the current permutation
                                                                // matrix from the left and store 
                                                                // the result in the double complex
                                                                // pointer tmp

        tmp = cmatMul(tmp, neighborSwap, dimension);           // multiply the next neighboring
                                                                // swap from the right to previously
                                                                // stored double complex pointer tmp

        free(result);                                           // free the allocated memory for
                                                                // result

        result = tmp;                                           // attach the pointer to the
                                                                // inermediate permutation matrix
                                                                // to result 
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
This function generates the matrix corresponding to a Pauli string

Input:
    const paulistr_t paulistr: user defined struct holding an array of Paulis and the number of
    qubits they act on

Output:
    Double complex pointer to the array holding the Kronecker product of single qubit Pauli matrices
    according to the array stored in paulistr
*/
cplx_t* pauliStringMat(const pauli_t paulistr[], qubit_t qubits) {
    cplx_t* result = ONE;                                               // double complex pointer
                                                                        // used to store the final
                                                                        // Kronecker product of
                                                                        // of Paulis on the qubits

    cplx_t* tmp;                                                        // double complex pointer
                                                                        // hodling each intermediate
                                                                        // Kronecker product of
                                                                        // Paulis

    /*
    Starting at the first position determines which of the Pauli matrices has to contribute to the
    Kronecker product with the current intermediate matrix stored in result until all qubits are
    evaluated.
    */
    for (qubit_t i = 0; i < qubits; ++i) {
        /*
        Case distinction: Which static matrix defined in the beginning has to contribute to the
        Kronecker product according to paulistr.
        */
        switch (paulistr[qubits - i - 1]) {
            /*
            For each case of a Pauli the Kronecker product of the intermediate matrix with the
            corresponding Pauli is stored in tmp.
            */
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
                printf("Error in constructing Pauli string matrix.\n");     // If the current Pauli
                                                                            // mentioned in paulistr
                                                                            // is not one out of the
                                                                            // above return an error
                exit(1);
            }
        }
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
This function generates the matrix corresponding to an observable that is a sum of weighted Pauli
strings

Input:
    const paulistr_t paulistr: user defined struct holding an array of Paulis and the number of
    qubits they act on

Output:
    Double complex pointer to the array holding the Kronecker product of single qubit Pauli matrices
    according to the array stored in paulistr
*/
cplx_t* pauliObservableMat(const pauli_t components[], const double coefficients[], complength_t length, qubit_t qubits)
{
    dim_t dim = POW2(qubits, dim_t);
    cplx_t* result = (cplx_t*) calloc(dim*dim, sizeof(cplx_t));
    cplx_t* tmp;

    for(complength_t i = 0; i < length; ++i) {
        tmp = pauliStringMat(components + (i*qubits), qubits);
        cscalarMatMulInPlace(coefficients[i], tmp, dim);
        cmatAddInPlace(result, tmp, dim);
    }
    free(tmp);
    return result;
}

/*
 * This function generates an array of the diagonal entries of an observable's matrix representation.
 *
 * Input:
 *      pauli_t components[]:       Pauli strings corresponding to the observable
 *      double coefficients[]:      Coefficients of the Pauli strings
 *
 */
double* pauliObservableDiag(const pauli_t components[], const double coefficients[], complength_t length, \
                            qubit_t qubits) {
    dim_t dim = POW2(qubits, dim_t);
    double* result = (double*) malloc(dim * sizeof(double));
    cplx_t* tmp = pauliObservableMat(components, coefficients, length, qubits);
    for (dim_t i = 0; i <dim; ++i) {
        result[i] = creal(tmp[i*dim + i]);
    }
    return result;
}

/*
This function generates the matrix for an evolution of a Pauli string by an angle.

Input:
    const paulistr_t paulistr: user defined struct holding an array of Paulis and the number of
    qubits they act on
Output:
    Double complex pointer to an array holding the exponential of a Pauli string times an angle and
    the imaginary unit. Since any Pauli string is involutive, one computes this as
    cos(angle) * ID - I * sin(angle) * PauliMat
*/
cplx_t* expPauliStringMat(const pauli_t paulistr[], double angle, qubit_t qubits) {
    dim_t dim = POW2(qubits, dim_t);                    // unsigned integer holding the Hilbert
                                                        // space's dimension according to the number
                                                        // of qubits

    cplx_t* result = identityMat(qubits);               // complex double pointer holding the final
                                                        // evolution matrix, initialised to be the 
                                                        // identity matrix on the specified number  
                                                        // of qubits

    cscalarMatMulInPlace(cos(angle), result, dim);      // multiply all entries of the identity 
                                                        // matrix stored in result by cos(angle)

    cplx_t* tmp = pauliStringMat(paulistr, qubits);     // complex double pointer holding the matrix
                                                        // corresponding to the Pauli string given 
                                                        // by paulistr

    cscalarMatMulInPlace(-I * sin(angle), tmp, dim);    // multiply all entries of the combined
                                                        // Pauli matrix stored in tmp by
                                                        // -I * sin(angle)

    cmatAddInPlace(result, tmp, dim);                   // add all entries in the result array by 
                                                        // the entries in tmp

    free(tmp);                                          // free the memory allocated to tmp
    return result;
}

cplx_t* expTrotterizedPauliObservableMat(const pauli_t components[], const double coefficients[], \
                                         complength_t length, double angle, qubit_t qubits) {
    dim_t dim = POW2(qubits, dim_t);

    cplx_t* result = identityMat(qubits);

    cplx_t* tmp;

    for (complength_t i = 0; i < length; ++i)
    {
        tmp = expPauliStringMat(components + (i * qubits), \
                                coefficients[i] * angle, qubits);
        cmatMulInPlace(result, tmp, dim);
    }

    return result;
}

/*
 * Finite difference methods
 */

double finiteExpValPQC(cplx_t* statevector, qubit_t qubits, dim_t dim, depth_t circdepth, \
                          pauli_t* compObs, double* coeffObs, complength_t lengthObs, \
                          pauli_t* compEvoOps, double* coeffEvoOps, complength_t lengthEvoOps, double* parameters
                       ) {
    /*
     * 1. Calculate the observable's matrix
     */
    cplx_t* observableMat = pauliObservableMat(compObs, coeffObs, lengthObs, qubits);

    /*
     * 2. Initialize the matrices for the evolution operators
     */
    cplx_t** evoOpsMat = (cplx_t**) malloc(circdepth * sizeof(cplx_t*));
    for (depth_t i = 0; i < circdepth; ++i) {
        evoOpsMat[i] = expTrotterizedPauliObservableMat(compEvoOps, \
                                                            coeffEvoOps + (i * lengthEvoOps), \
                                                            lengthEvoOps, \
                                                            parameters[i], qubits);
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
    for (depth_t i = 0; i < circdepth; ++i) {
        free(evoOpsMat[i]);
    }
    free(evoOpsMat);

    return result;
}

double* finiteGradientPQC(cplx_t* statevector, qubit_t qubits, dim_t dim, depth_t circdepth, \
                          pauli_t* compObs, double* coeffObs, complength_t lengthObs, \
                          pauli_t* compEvoOps, double* coeffEvoOps, complength_t lengthEvoOps, double* parameters, \
                          double epsilon) {
    double* result = malloc(circdepth * sizeof(double));
    /*
     * 1. Calculate the observable's matrix
     */
    cplx_t* observableMat = pauliObservableMat(compObs, coeffObs, lengthObs, qubits);

    /*
     * 2. Calculate the mean value at the current parameter setting
     */
    double expVal = finiteExpValPQC(statevector, qubits, dim, circdepth, compObs, coeffObs, lengthObs, \
                                    compEvoOps, coeffEvoOps, lengthEvoOps, parameters);

    /*
     * 3. Shift the parameter setting by epsilon in all directions; the difference of the mean and the previous mean
     * give the finite difference's entry
     */
    for (depth_t i = 0; i < circdepth; ++i) {
        parameters[i] += epsilon;

        double value = finiteExpValPQC(statevector, qubits, dim, circdepth, compObs, coeffObs, lengthObs, \
                                       compEvoOps, coeffEvoOps, lengthEvoOps, parameters);

        result[i] = (1. / epsilon) * (value - expVal);
        parameters[i] -= epsilon;
    }

    return result;
}

double* finiteHessianPQC(cplx_t* statevector, qubit_t qubits, dim_t dim, depth_t circdepth, \
                         pauli_t* compObs, double* coeffObs, complength_t lengthObs, \
                         pauli_t* compEvoOps, double* coeffEvoOps, complength_t lengthEvoOps, double* parameters, \
                         double epsilon
                         ) {
    double* result = malloc(circdepth*circdepth * sizeof(double));

    /*
     * 1. Calculate the gradient at the current parameter setting
     */
    double* refGradient = finiteGradientPQC(statevector, qubits, dim, circdepth, compObs, coeffObs, lengthObs, \
                                            compEvoOps, coeffEvoOps, lengthEvoOps, parameters, epsilon);

    /*
     * 2. Shift the parameter setting by epsilon in all directions; the entrywise difference to the prvious gradient
     * gives the corresponding column of the hessian
     */
    for (depth_t i = 0; i < circdepth; ++i) {
        parameters[i] += epsilon;
        double* gradient = finiteGradientPQC(statevector, qubits, dim, circdepth, compObs, coeffObs, lengthObs, \
                                            compEvoOps, coeffEvoOps, lengthEvoOps, parameters, epsilon);
        for (depth_t j = 0; j < circdepth; ++j) {
            result[j*circdepth + i] = (1./epsilon) * (gradient[j] - refGradient[j]);
        }
        parameters[i] -= epsilon;
    }
    return result;
}
