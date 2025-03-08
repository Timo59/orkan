//
// Created by Timo Ziegler on 19.09.24.
//

#include "../../test/pauliMat.h"

/*
 * =====================================================================================================================
 *                                              Pauli string constructions
 * =====================================================================================================================
 */

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
 * This function returns an array of all possible diagonal Pauli strings; i.e. combinations of id and Z on the specified
 * number of qubits
 */
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

/*
 * This function returns an array of all possible Pauli strings consisting of id and Pauli X on the specified number of
 * qubits
 */
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

/*
 * This function returns an array of all possible Pauli strings consisting of id and Pauli Y on the specified number of
 * qubits
 */
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
 * =====================================================================================================================
 *                                              Pauli matrix constructions
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