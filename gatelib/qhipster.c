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

void applyXblas(state_t* state, qubit_t qubit) {
    __LAPACK_int flipDistance = POW2(qubit, __LAPACK_int);
    __LAPACK_int blockDistance = POW2(qubit + 1, __LAPACK_int);

    for (__LAPACK_int i = 0; i < state->dim; i += blockDistance) {
        cblas_zswap(flipDistance, state->vec + i, 1, state->vec + (i + flipDistance), 1);
    }
}

void applyXomp(state_t* state, qubit_t qubit) {
    dim_t flipDistance = POW2(qubit, dim_t);
    dim_t indices = POW2(state->qubits - 1, dim_t);

#pragma omp parallel for default(none) shared(flipDistance, indices, qubit, state)
    for (dim_t i = 0; i < indices; ++i) {
        dim_t left = (i & (~0U << qubit)) << 1;
        dim_t right = i & ~(~0U << qubit);
        dim_t j = (left | right);
        SWAP(state->vec + j, state->vec + (j + flipDistance), cplx_t);
    }
}

void applyXomp_blas(state_t* state, qubit_t qubit) {
    __LAPACK_int flipDistance = POW2(qubit, __LAPACK_int);
    __LAPACK_int blockDistance = POW2(qubit + 1, __LAPACK_int);

#pragma omp parallel for default(none) shared(blockDistance, flipDistance, qubit, state)
    for (__LAPACK_int i = 0; i < state->dim; i += blockDistance) {
        cblas_zswap(flipDistance, state->vec + i, 1, state->vec + (i + flipDistance), 1);
    }
}

void applyXgcd(state_t* state, qubit_t qubit) {
    dim_t flipDistance = POW2(qubit, dim_t);
    dim_t indices = POW2(state->qubits - 1, dim_t);
    dim_t stride = (state->qubits > 6) ? POW2(state->qubits - 5, dim_t) : indices;

    /* Get a concurrent dispatch queue */
    dispatch_queue_t queue = dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0);

    dispatch_apply(indices / stride, queue, ^(size_t i) {
        size_t j = i * stride;
        size_t j_stop = j + stride;

        do {
            dim_t left = (j & (~0U << qubit)) << 1;
            dim_t right = j & ~(~0U << qubit);
            dim_t k = (left | right);
            SWAP(state->vec + k, state->vec + (k + flipDistance), cplx_t);
            ++j;
        }while (j < j_stop);
    });
}

/**********************************************************************************************************************/

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

void applyYblas(state_t* state, qubit_t qubit) {
    __LAPACK_int flipDistance = POW2(qubit, __LAPACK_int);
    __LAPACK_int blockDistance = POW2(qubit + 1, __LAPACK_int);
    __LAPACK_double_complex plus = I;
    __LAPACK_double_complex minus = -I;

    for (__LAPACK_int i = 0; i < state->dim; i += blockDistance) {
        cblas_zswap(flipDistance, state->vec + i, 1, state->vec + (i + flipDistance), 1);
        cblas_zscal(flipDistance, &minus, state->vec + i, 1);
        cblas_zscal(flipDistance, &plus, state->vec + (i + flipDistance), 1);
    }
}

void applyYomp(state_t* state, qubit_t qubit) {
    dim_t flipDistance = POW2(qubit, dim_t);
    dim_t indices = POW2(state->qubits - 1, dim_t);

#pragma omp parallel for default(none) shared(flipDistance, indices, qubit, state)
    for (dim_t i = 0; i < indices; ++i) {
        dim_t left = (i & (~0U << qubit)) << 1;
        dim_t right = i & ~(~0U << qubit);
        dim_t j = (left | right);
        SWAP(state->vec + j, state->vec + (j + flipDistance), cplx_t);
        state->vec[j] *= -I;
        state->vec[j + flipDistance] *= I;
    }
}

void applyYomp_blas(state_t* state, qubit_t qubit) {
    __LAPACK_int flipDistance = POW2(qubit, __LAPACK_int);
    __LAPACK_int blockDistance = POW2(qubit + 1, __LAPACK_int);
    __LAPACK_double_complex plus = I;
    __LAPACK_double_complex minus = -I;

#pragma omp parallel for default(none) shared(blockDistance, flipDistance, minus, plus, qubit, state)
    for (__LAPACK_int i = 0; i < state->dim; i += blockDistance) {
        cblas_zswap(flipDistance, state->vec + i, 1, state->vec + (i + flipDistance), 1);
        cblas_zscal(flipDistance, &minus, state->vec + i, 1);
        cblas_zscal(flipDistance, &plus, state->vec + (i + flipDistance), 1);
    }
}

void applyYgcd(state_t* state, qubit_t qubit) {
    dim_t flipDistance = POW2(qubit, dim_t);
    dim_t indices = POW2(state->qubits - 1, dim_t);
    dim_t stride = (state->qubits > 6) ? POW2(state->qubits - 5, dim_t) : indices;

    /* Get a concurrent dispatch queue */
    dispatch_queue_t queue = dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0);

    dispatch_apply(indices / stride, queue, ^(size_t i) {
        size_t j = i * stride;
        size_t j_stop = j + stride;

        do {
            dim_t left = (j & (~0U << qubit)) << 1;
            dim_t right = j & ~(~0U << qubit);
            dim_t k = (left | right);
            SWAP(state->vec + k, state->vec + (k + flipDistance), cplx_t);
            state->vec[k] *= -I;
            state->vec[k + flipDistance] *= I;
            ++j;
        }while (j < j_stop);
    });
}

/**********************************************************************************************************************/

void applyZ(state_t* state, qubit_t qubit) {
    dim_t blockDistance = POW2(qubit + 1, dim_t);
    dim_t flipDistance = POW2(qubit, dim_t);
    for (dim_t i = 0; i < state->dim; i += blockDistance) {
        for (dim_t j = i; j < i + flipDistance; ++j) {
            state->vec[j + flipDistance] *= -1;
        }
    }
}

void applyZblas(state_t* state, qubit_t qubit) {
    __LAPACK_int flipDistance = POW2(qubit, __LAPACK_int);
    __LAPACK_int blockDistance = POW2(qubit + 1, __LAPACK_int);
    __LAPACK_double_complex minus = -1;

    for (__LAPACK_int i = 0; i < state->dim; i += blockDistance) {
        cblas_zscal(flipDistance, &minus, state->vec + (i + flipDistance), 1);
    }
}

void applyZomp(state_t* state, qubit_t qubit) {
    dim_t flipDistance = POW2(qubit, dim_t);
    dim_t indices = POW2(state->qubits - 1, dim_t);

#pragma omp parallel for default(none) shared(flipDistance, indices, qubit, state)
    for (dim_t i = 0; i < indices; ++i) {
        dim_t left = (i & (~0U << qubit)) << 1;
        left |= (1 << qubit);
        dim_t right = i & ~(~0U << qubit);
        dim_t j = (left | right);
        state->vec[j] *= -1;
    }
}

void applyZomp_blas(state_t* state, qubit_t qubit) {
    __LAPACK_int flipDistance = POW2(qubit, __LAPACK_int);
    __LAPACK_int blockDistance = POW2(qubit + 1, __LAPACK_int);
    __LAPACK_double_complex minus = -1;

#pragma omp parallel for default(none) shared(blockDistance, flipDistance, minus, qubit, state)
    for (__LAPACK_int i = 0; i < state->dim; i += blockDistance) {
        cblas_zscal(flipDistance, &minus, state->vec + (i + flipDistance), 1);
    }
}

void applyZgcd(state_t* state, qubit_t qubit) {
    dim_t flipDistance = POW2(qubit, dim_t);
    dim_t indices = POW2(state->qubits - 1, dim_t);
    dim_t stride = (state->qubits > 6) ? POW2(state->qubits - 5, dim_t) : indices;

    /* Get a concurrent dispatch queue */
    dispatch_queue_t queue = dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0);

    dispatch_apply(indices / stride, queue, ^(size_t i) {
        size_t j = i * stride;
        size_t j_stop = j + stride;

        do {
            dim_t left = (j & (~0U << qubit)) << 1;
            left |= (1 << qubit);
            dim_t right = j & ~(~0U << qubit);
            dim_t k = (left | right);
            state->vec[k] *= -1;
            ++j;
        }while (j < j_stop);
    });
}

/*
 * =====================================================================================================================
 *                                                  Clifford gates
 * =====================================================================================================================
 */

void applyS(state_t* state, qubit_t qubit) {
    dim_t blockDistance = POW2(qubit + 1, dim_t);
    dim_t flipDistance = POW2(qubit, dim_t);
    for (dim_t i = 0; i < state->dim; i += blockDistance) {
        for (dim_t j = i; j < i + flipDistance; ++j) {
            state->vec[j + flipDistance] *= I;
        }
    }
}

void applySdagger(state_t* state, qubit_t qubit) {
    dim_t blockDistance = POW2(qubit + 1, dim_t);
    dim_t flipDistance = POW2(qubit, dim_t);
    for (dim_t i = 0; i < state->dim; i += blockDistance) {
        for (dim_t j = i; j < i + flipDistance; ++j) {
            state->vec[j + flipDistance] *= -I;
        }
    }
}

void applyH(state_t* state, qubit_t qubit) {
    dim_t blockDistance = POW2(qubit + 1, dim_t);
    dim_t flipDistance = POW2(qubit, dim_t);
    cplx_t tmp;
    for (dim_t i = 0; i < state->dim; i += blockDistance) {
        for (dim_t j = i; j < i + flipDistance; ++j) {
            tmp = state->vec[j];
            state->vec[j] = INVSQRT2 * (tmp + state->vec[j + flipDistance]);
            state->vec[j + flipDistance] = INVSQRT2 * (tmp - state->vec[j + flipDistance]);
        }
    }
}

/*
 * =====================================================================================================================
 *                                                  Hadamard-Y gate
 * =====================================================================================================================
 */

void applyHy(state_t* state, qubit_t qubit) {
    dim_t blockDistance = POW2(qubit + 1, dim_t);
    dim_t flipDistance = POW2(qubit, dim_t);
    cplx_t tmp;
    for (dim_t i = 0; i < state->dim; i += blockDistance) {
        for (dim_t j = i; j < i + flipDistance; ++j) {
            tmp = state->vec[j];
            state->vec[j] = INVSQRT2 * (tmp - I * state->vec[j + flipDistance]);
            state->vec[j + flipDistance] = INVSQRT2 * (I * tmp - state->vec[j + flipDistance]);
        }
    }
}

/*
 * =====================================================================================================================
 *                                                      T gate
 * =====================================================================================================================
 */

void applyT(state_t* state, qubit_t qubit) {
    dim_t blockDistance = POW2(qubit + 1, dim_t);
    dim_t flipDistance = POW2(qubit, dim_t);
    for (dim_t i = 0; i < state->dim; i += blockDistance) {
        for (dim_t j = i; j < i + flipDistance; ++j) {
            state->vec[j + flipDistance] *= (INVSQRT2 + INVSQRT2 * I);
        }
    }
}

void applyTdagger(state_t* state, qubit_t qubit) {
    dim_t blockDistance = POW2(qubit + 1, dim_t);
    dim_t flipDistance = POW2(qubit, dim_t);
    for (dim_t i = 0; i < state->dim; i += blockDistance) {
        for (dim_t j = i; j < i + flipDistance; ++j) {
            state->vec[j + flipDistance] *= (INVSQRT2 - INVSQRT2 * I);
        }
    }
}

/*
 * =====================================================================================================================
 *                                                      P gate
 * =====================================================================================================================
 */

void applyP(state_t* state, qubit_t qubit, double angle) {
    dim_t blockDistance = POW2(qubit + 1, dim_t);
    dim_t flipDistance = POW2(qubit, dim_t);
    for (dim_t i = 0; i < state->dim; i += blockDistance) {
        for (dim_t j = i; j < i + flipDistance; ++j) {
            state->vec[j + flipDistance] *= cos(angle) + sin(angle) * I;
        }
    }
}

void applyPdagger(state_t* state, qubit_t qubit, double angle) {
    applyP(state, qubit, -angle);
}

/*
 * =====================================================================================================================
 *                                                  Rotation gates
 * =====================================================================================================================
 */

void applyRX(state_t* state, qubit_t qubit, double angle) {
    dim_t blockDistance = POW2(qubit + 1, dim_t);
    dim_t flipDistance = POW2(qubit, dim_t);
    angle /= 2.;
    cplx_t tmp;
    for (dim_t i = 0; i < state->dim; i += blockDistance) {
        for (dim_t j = i; j < i + flipDistance; ++j) {
            tmp = state->vec[j];
            state->vec[j] = cos(angle) * tmp - I * sin(angle) * state->vec[j + flipDistance];
            state->vec[j + flipDistance] = -I * sin(angle) * tmp + cos(angle) \
 * state->vec[j + flipDistance];
        }
    }
}

void applyRXdagger(state_t* state, qubit_t qubit, double angle) {
    applyRX(state, qubit, -angle);
}

void applyRY(state_t* state, qubit_t qubit, double angle) {
    dim_t blockDistance = POW2(qubit + 1, dim_t);
    dim_t flipDistance = POW2(qubit, dim_t);
    angle /= 2.;
    cplx_t tmp;
    for (dim_t i = 0; i < state->dim; i += blockDistance) {
        for (dim_t j = i; j < i + flipDistance; ++j) {
            tmp = state->vec[j];
            state->vec[j] = cos(angle) * tmp - sin(angle) * state->vec[j + flipDistance];
            state->vec[j + flipDistance] = sin(angle) * tmp + cos(angle) \
 * state->vec[j + flipDistance];
        }
    }
}

void applyRYdagger(state_t* state, qubit_t qubit, double angle) {
    applyRY(state, qubit, -angle);
}

void applyRZ(state_t* state, qubit_t qubit, double angle) {
    dim_t blockDistance = POW2(qubit + 1, dim_t);
    dim_t flipDistance = POW2(qubit, dim_t);
    angle /= 2.;
    for (dim_t i = 0; i < state->dim; i += blockDistance) {
        for (dim_t j = i; j < i + flipDistance; ++j) {
            state->vec[j] *= cos(angle) - I * sin(angle);
            state->vec[j + flipDistance] *= cos(angle) + I * sin(angle);
        }
    }
}

void applyRZdagger(state_t* state, qubit_t qubit, double angle) {
    applyRZ(state, qubit, -angle);
}

/*
 * =====================================================================================================================
 *                                              controlled Pauli gates
 * =====================================================================================================================
 */

/*
 * This function applies a Pauli-X gate on the target qubit conditioned on the control qubit
 *
 * Input:
 *      state_t* state:     State of qubit system
 *      qubit_t control:    Qubit the execution is conditioned on
 *      qubit_t target:     Qubit the single qubit gate is applied to conditioned on control
 *
 * Output:
 *      The function has no output, but alters the state vector of state according to the execution of the gate.
 */
void applyCX(state_t* state, qubit_t control, qubit_t target) {
    dim_t blockDistance = POW2(target + 1, dim_t);      // Range of indices whose state vector entries after the
                                                        // execution of the single qubit gate will be linear
                                                        // combinations of entries with indices only within the range
    dim_t flipDistance = POW2(target, dim_t);           // Distance of indices whose state vector entries after the
                                                        // execution of the single qubit gate will be linear
                                                        // combinations of their entries prior to the execution


    for (dim_t i = 0; i < state->dim; i += blockDistance) {   // Iterate the computational basis states with all zeros
                                                                // to the right and including the target qubit

        for (dim_t j = i; j < i + flipDistance; ++j) {          // Iterate the computational basis states whose
                                                                // qubits left to and including the qubit specified
                                                                    // by target are constant

            if (j & POW2(control ,dim_t)) {                                             // If the current iterate has a
                                                                                        // one at the control position:
                SWAP(state->vec + j, state->vec + (j + flipDistance), cplx_t);    // Execute Pauli-X on target
            }
        }
    }
}

void applyCY(state_t* state, qubit_t control, qubit_t target) {
    dim_t blockDistance = POW2(target + 1, dim_t);
    dim_t flipDistance = POW2(target, dim_t);
    for (dim_t i = 0; i < state->dim; i += blockDistance) {
        for (dim_t j = i; j < i + flipDistance; ++j) {
            if (j & POW2(control ,dim_t)) {
                SWAP(state->vec + j, state->vec + (j + flipDistance), cplx_t);
                state->vec[j] *= -I;
                state->vec[j + flipDistance] *= I;
            }
        }
    }
}

void applyCZ(state_t* state, qubit_t control, qubit_t target) {
    dim_t blockDistance = POW2(target + 1, dim_t);
    dim_t flipDistance = POW2(target, dim_t);
    for (dim_t i = 0; i < state->dim; i += blockDistance) {
        for (dim_t j = i; j < i + flipDistance; ++j) {
            if (j & POW2(control ,dim_t)) {
                state->vec[j + flipDistance] *= -1;
            }
        }
    }
}

/*
 * =====================================================================================================================
 *                                              Controlled Clifford gates
 * =====================================================================================================================
 */

void applyCS(state_t* state, qubit_t control, qubit_t target) {
    dim_t blockDistance = POW2(target + 1, dim_t);
    dim_t flipDistance = POW2(target, dim_t);
    for (dim_t i = 0; i < state->dim; i += blockDistance) {
        for (dim_t j = i; j < i + flipDistance; ++j) {
            if (j & POW2(control ,dim_t)) {
                state->vec[j + flipDistance] *= I;
            }
        }
    }
}

void applyCSdagger(state_t* state, qubit_t control, qubit_t target) {
    dim_t blockDistance = POW2(target + 1, dim_t);
    dim_t flipDistance = POW2(target, dim_t);
    for (dim_t i = 0; i < state->dim; i += blockDistance) {
        for (dim_t j = i; j < i + flipDistance; ++j) {
            if (j & POW2(control ,dim_t)) {
                state->vec[j + flipDistance] *= -I;
            }
        }
    }
}

void applyCH(state_t* state, qubit_t control, qubit_t target) {
    dim_t blockDistance = POW2(target + 1, dim_t);
    dim_t flipDistance = POW2(target, dim_t);
    cplx_t tmp;
    for (dim_t i = 0; i < state->dim; i += blockDistance) {
        for (dim_t j = i; j < i + flipDistance; ++j) {
            if (j & POW2(control ,dim_t)) {
                tmp = state->vec[j];
                state->vec[j] = INVSQRT2 * (tmp + state->vec[j + flipDistance]);
                state->vec[j + flipDistance] = INVSQRT2 * (tmp - state->vec[j + flipDistance]);
            }
        }
    }
}

/*
 * =====================================================================================================================
 *                                              Controlled Hadamard-Y gate
 * =====================================================================================================================
 */

void applyCHy(state_t* state, qubit_t control, qubit_t target) {
    dim_t blockDistance = POW2(target + 1, dim_t);
    dim_t flipDistance = POW2(target, dim_t);
    cplx_t tmp;
    for (dim_t i = 0; i < state->dim; i += blockDistance) {
        for (dim_t j = i; j < i + flipDistance; ++j) {
            if (j & POW2(control ,dim_t)) {
                tmp = state->vec[j];
                state->vec[j] = INVSQRT2 * (tmp - I * state->vec[j + flipDistance]);
                state->vec[j + flipDistance] = INVSQRT2 * (I * tmp - state->vec[j + flipDistance]);
            }
        }
    }
}

/*
 * =====================================================================================================================
 *                                                  Controlled T gate
 * =====================================================================================================================
 */

void applyCT(state_t* state, qubit_t control, qubit_t target) {
    dim_t blockDistance = POW2(target + 1, dim_t);
    dim_t flipDistance = POW2(target, dim_t);
    for (dim_t i = 0; i < state->dim; i += blockDistance) {
        for (dim_t j = i; j < i + flipDistance; ++j) {
            if (j & POW2(control ,dim_t)) {
                state->vec[j + flipDistance] *= INVSQRT2 + INVSQRT2 * I;
            }
        }
    }
}

void applyCTdagger(state_t* state, qubit_t control, qubit_t target) {
    dim_t blockDistance = POW2(target + 1, dim_t);
    dim_t flipDistance = POW2(target, dim_t);
    for (dim_t i = 0; i < state->dim; i += blockDistance) {
        for (dim_t j = i; j < i + flipDistance; ++j) {
            if (j & POW2(control ,dim_t)) {
                state->vec[j + flipDistance] *= INVSQRT2 - INVSQRT2 * I;
            }
        }
    }
}

/*
 * =====================================================================================================================
 *                                                  Controlled P gate
 * =====================================================================================================================
 */

void applyCP(state_t* state, qubit_t control, qubit_t target, double angle) {
    dim_t blockDistance = POW2(target + 1, dim_t);
    dim_t flipDistance = POW2(target, dim_t);
    for (dim_t i = 0; i < state->dim; i += blockDistance) {
        for (dim_t j = i; j < i + flipDistance; ++j) {
            if (j & POW2(control ,dim_t)) {
                state->vec[j + flipDistance] *= cos(angle) + sin(angle) * I;
            }
        }
    }
}

void applyCPdagger(state_t* state, qubit_t control, qubit_t target, double angle) {
    applyCP(state, control, target, -angle);
}

/*
 * =====================================================================================================================
 *                                                      SWAP gate
 * =====================================================================================================================
 */

/*
 * This function executes a SWAP of qubit1 and qubit2
 *
 * Input:
 *      state_t* state:     State of qubit system
 *      qubit_t qubit1:     Position of one qubit to swap
 *      qubit_t qubit2:     Position of the other qubit to swap
 *
 */
void applySWAP(state_t* state, qubit_t qubit1, qubit_t qubit2) {

    qubit_t left_qubit = MAX(qubit1, qubit2);                       // index of the qubit to the left to be swapped
    qubit_t right_qubit = MIN(qubit1, qubit2);                      // index of the qubit to the right to be swapped
    dim_t majorBlockDistance = POW2(left_qubit + 1, dim_t);         // Range of indices whose state vector entries after
                                                                    // execution of the two qubit gate will be linear
                                                                    // combinations of entries with indices only within
                                                                    // the range
    dim_t minorBlockDistance = POW2(right_qubit + 1, dim_t);        //
    dim_t leftOffset = POW2(left_qubit, dim_t);                     // Range of computational basis states that leave
                                                                    // the left qubit at zero when started from all
                                                                    // qubits left to it set to zero
    dim_t rightOffset = POW2(right_qubit, dim_t);                // Distance of computational basis states with zero
                                                                    // and one at the right qubit, respectively
    dim_t flipDistance = POW2(left_qubit, dim_t) - POW2(right_qubit, dim_t);    // Distance of two computational basis
                                                                                // with zero at the left and one at the
                                                                                // right qubit and with one at the left
                                                                                // and zero at the right qubit
    for (dim_t i = rightOffset; i < state->dim; i += majorBlockDistance) {        // Iterate the computational basis states with
                                                                        // all zeros to the right of and including the
                                                                        // left qubit and a one at the right qubit
        for (dim_t j = i ; j < i + leftOffset; j += minorBlockDistance) {   // Iterate the computational basis
                                                                                // states with constant qubits to the
                                                                                // left and including the left qubit
                                                                                // starting at the state where the right
                                                                                // qubit is one
            for (dim_t k = j ; k < j + rightOffset; ++k) {                                     // If the right qubit is one:
                SWAP(state->vec + k, state->vec + (k + flipDistance), cplx_t);      // Swap entry with the one
                                                                                    // corresponding to zero at the
                                                                                    // right and one at the left qubit
            }
        }
    }
}

/*
This function applies the Toffoli gate to the state.

Input
    state_t state: a user defined structure holding the pointer to a double complex valued array
    corresponding to the state vector, the number of qubits as an unisgned integer of type qubits_t
    and the Hilbert space dimension aka the number of entries in the state vector as an unsigned
    integer of type dim_t
    qubit_t control1: an unsigned integer holding the first of the two qubits' positions from the
    right in a bit string that control the execution of the Pauli-X gate on the target qubit
    qubit_t control2: an unsigned integer holding the second of the two qubits' positions from the
    right in a bit string that control the execution of the Pauli-X gate on the target qubit
    qubit_t target: an unsigned integer holding the qubit's position from the right in a bit string
    the Pauli-X gate acts on if the control qubits both are in the |1> state
    
*/

/*
 * =====================================================================================================================
 *                                              Toffoli gate
 * =====================================================================================================================
 */

void applyToffoli(state_t* state, qubit_t control1, qubit_t control2, qubit_t target)
{
    dim_t i;                                            // index for the outer loop which runs
                                                        // through the blocks, i.e., increseases
                                                        // the bits left to target by one

    dim_t j;                                            // index of the inner loop running through
                                                        // all indices within the block, i.e.,
                                                        // increases bits right to target by one

    dim_t blockDistance = POW2(target + 1, dim_t);      // distance of two integers whose binary
                                                        // representation differs only by one seen
                                                        // from all bits left to target

    dim_t flipDistance = POW2(target, dim_t);           // distance of two integers whose binary
                                                        // representation only differs in the bit 
                                                        // at position target

    /*
    Starting with all bits set to zero the integer i is the starting index for the inner loop and
    is increased by blockDistance after each iteration, thereby incrementing the bits left to target
    by one.
    */
    for (i = 0; i < state->dim; i += blockDistance)
    {
        /*
        The integer j starts at the same configuration of bits in the index's binary representation
        gradually changing all bits to the right of target until all are set to one
        */
        for (j = i; j < i + flipDistance; ++j)
        {
            /*
            The Pauli-X gate is applied to target only if both control bits in the index's binary 
            representation are set to one.
            */
            if (j & POW2(control1, dim_t) && j & POW2(control2, dim_t)) 
            {
                /*
                Meeting the condition, the value of the j-th entry originating from the pointer
                stored in state is stored to a temporary double complex variable. Subesequently, one
                sets the value the j-th entry's pointer holds to be the value of the entry whose
                index's binary representation only differs in qubit. Finally, the latter pointer is
                set to the value stored under the temporary and j.
                */
                SWAP(state->vec + j, state->vec + (j + flipDistance), cplx_t);
            }
        }
    }
}
