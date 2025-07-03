//
// Created by Timo Ziegler on 03.07.25.
//

#ifndef MSTATE_H
#define MSTATE_H
/*
 * =====================================================================================================================
 *                                                      includes
 * =====================================================================================================================
 */
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
#endif //MSTATE_H

    /*
 * =====================================================================================================================
 *                                                  type definitions
 * =====================================================================================================================
 */

/*
 * Struct:      state_t
 * --------------------
 * This struct represents a multi-qubit mixed state.
 * Contents:
 *      vec:        Vector representation of the state.
 *      qubits:     Number of qubits of the underlying system.
 *      dim:        Length of the state's vector representation.
 */
typedef struct state {
    cplx_t* vec;
    qubit_t qubits;
    dim_t dim;
} state_t;