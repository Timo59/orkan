/*
 * =====================================================================================================================
 *                                                      includes
 * =====================================================================================================================
 */
#ifndef GATELIB_H
#include "gatelib.h"
#endif

/*
 * =====================================================================================================================
 *                                                  Pauli gates
 * =====================================================================================================================
 */

/*
 * This function executes the Pauli-X gate on a specified qubit of a quantum state in a qubit system
 *
 * Input:
 *      state_t* state:     State of a qubits system
 *      qubit_t qubit:      Qubit the gate is executed on
 *
 * Output:
 *      There is no return value; the state vector is changed in place.
 */

void applyX(state_t* state, qubit_t qubit) {
    dim_t blockDistance = POW2(qubit + 1, dim_t);       // Range of indices whose state vector entries after the
                                                        // execution of the single qubit gate will be linear
                                                        // combinations of entries with indices only within the range
    dim_t flipDistance = POW2(qubit, dim_t);            // Distance of indices whose state vector entries after the
                                                        // execution of the single qubit gate will be linear
                                                        // combinations of their entries prior to the execution

    for (dim_t i = 0; i < state->dim; i += blockDistance) {     // Iterate the computational basis states with all zeros
                                                                // to the right and including the target qubit
        for (dim_t j = i; j < i + flipDistance; ++j) {              // Iterate the computational basis states whose
                                                                    // qubits left to and including the target qubit
                                                                    // are constant
            SWAP(state->vec + j, state->vec + (j + flipDistance), cplx_t);      // Swap entry with entry at flipDistance
        }
    }
}

void applyY(state_t* state, qubit_t qubit) {
    dim_t blockDistance = POW2(qubit + 1, dim_t);
    dim_t flipDistance = POW2(qubit, dim_t);
    for (dim_t i = 0; i < state->dim; i += blockDistance) {
        for (dim_t j = i; j < i + flipDistance; ++j) {
            SWAP(state->vec + j, state->vec + (j + flipDistance), cplx_t);
            state->vec[j] *= -I;
            state->vec[j + flipDistance] *= I;
        }
    }
}

void applyZ(state_t* state, qubit_t qubit) {
    dim_t blockDistance = POW2(qubit + 1, dim_t);
    dim_t flipDistance = POW2(qubit, dim_t);
    for (dim_t i = 0; i < state->dim; i += blockDistance) {
        for (dim_t j = i; j < i + flipDistance; ++j) {
            state->vec[j + flipDistance] *= -1;
        }
    }
}