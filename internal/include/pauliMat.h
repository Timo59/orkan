//
// Created by Timo Ziegler on 19.09.24.
//

#ifndef PAULIMAT_H
#define PAULIMAT_H

/*
 * =====================================================================================================================
 *                                                      Includes
 * =====================================================================================================================
 */

#ifndef GATEMAT_H
#include "gateMat.h"
#endif

#ifndef QTYPES_H
#include "qTypes.h"
#endif

#ifndef _STDLIB_H_
#include <stdlib.h>
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
*                                              Pauli string constructions
* =====================================================================================================================
*/

pauli_t* allPauliStrings(qubit_t qubits);

pauli_t* allPauliStringsDiag(qubit_t qubits);

pauli_t* allPauliStringsX(qubit_t qubits);

pauli_t* allPauliStringsY(qubit_t qubits);

/*
 * =====================================================================================================================
 *                                              Pauli matrix constructions
 * =====================================================================================================================
 */

cplx_t* pauliStrMat(const pauli_t paulistr[], qubit_t qubits);

cplx_t* pauliObsMat(const pauli_t comps[], const double coeffs[], complength_t length, qubit_t qubits);

double* pauliObsMatDiag(const pauli_t comps[], const double coeffs[], complength_t length, qubit_t qubits);

cplx_t* expPauliStrMat(const pauli_t paulistr[], double angle, qubit_t qubits);

cplx_t* expPauliObsMatTrotter(const pauli_t comps[],
                              const double coeffs[],
                              complength_t length,
                              double angle,
                              qubit_t qubits);

#ifdef __cplusplus
}
#endif

#endif /* PAULIMAT_H */
