//
// Created by Timo Ziegler on 09.10.24.
//
/*
 * =====================================================================================================================
 *                                                      includes
 * =====================================================================================================================
 */
#ifndef __MATH__
#include <math.h>
#endif

#ifndef PSTATE_H
#include "pstate.h"
#endif

#ifndef _STDLIB_H_
#include <stdlib.h>
#endif

#ifdef MACOS
#include <vecLib/blas_new.h>
#include <vecLib/lapack.h>
#else
#include <f77blas.h>
#endif
#include <stdio.h>

/*
 * =====================================================================================================================
 *                                                  Internal functions
 * =====================================================================================================================
 */
/*
 * This function inserts a one at the specified position in the bit string corresponding to the given long integer.
 *
 * @param[in] a     Integer the one is inserted into
 * @param[in] pos   Position in the bit string representation to insert the one
 *
 * @return          Integer where one is inserted at pos
 */
static inline dim_t insone(dim_t a, dim_t pos) {
    dim_t left = a & (~0U << pos);                                  // Create a mask of all qubits left to control
    left = (left << 1) | 1;                                         // Append another bit and set it to one
    const dim_t right = a & ~(~0U << pos);                          // Create a mask of all qubits right to control
    const dim_t out = (left | right);                               // Concatenate the two bit strings
    return out;
}

/*
 * =====================================================================================================================
 *                                                  Pauli gates
 * =====================================================================================================================
 */
/*
 * This function applies the Pauli-X gate to the specified qubit in the specified pure quantum state.
 *
 * @param[in,out] state     Address of the state struct the gate is applied to
 * @param[in] qubit         Qubit the gate acts on
 */
void applyX(state_t* state, const qubit_t qubit) {
    const dim_t flipStride = POW2(qubit, dim_t);                    // Distance between vector elements differing only
                                                                    // at position qubit
    const dim_t blockStride = POW2(qubit + 1, dim_t);               // Number of contiguous elements the tensor product
                                                                    // of gate and identities acts on
    const dim_t incr = 1;                                           // Increment (needed for BLAS function)

    for (dim_t i = 0; i < state->dim; i += blockStride) {           // Iterate blocks with all qubits left to the
                                                                    // addressed qubit fixed
        zswap_(&flipStride, state->vec + i, &incr, state->vec + i + flipStride, &incr);
    }
}

/*
 * This function applies the Pauli-Y gate to the specified qubit in the specified pure quantum state.
 *
 * @param[in,out] state     Address of the state struct the gate is applied to
 * @param[in] qubit         Qubit the gate acts on
 */
void applyY(state_t* state, const qubit_t qubit) {
    const dim_t flipStride = POW2(qubit, dim_t);
    const dim_t blockStride = POW2(qubit + 1, dim_t);
    const dim_t incr = 1;
    const cplx_t plus = I;
    const cplx_t minus = -I;

    for (dim_t i = 0; i < state->dim; i += blockStride) {
        zswap_(&flipStride, state->vec + i, &incr, state->vec + i + flipStride, &incr);
        zscal_(&flipStride, &minus, state->vec + i, &incr);
        zscal_(&flipStride, &plus, state->vec + i + flipStride, &incr);
    }
}

/*
 * This function applies the Pauli-Z gate to the specified qubit in the specified pure quantum state.
 *
 * @param[in,out] state     Address of the state struct the gate is applied to
 * @param[in] qubit         Qubit the gate acts on
 */
void applyZ(state_t* state, const qubit_t qubit) {
    const dim_t flipStride = POW2(qubit, dim_t);
    const dim_t blockStride = POW2(qubit + 1, dim_t);
    const dim_t incr = 1;
    const cplx_t minus = -1;

    for (dim_t i = 0; i < state->dim; i += blockStride) {
        zscal_(&flipStride, &minus, state->vec + i + flipStride, &incr);
    }
}

/*
 * =====================================================================================================================
 *                                                  Clifford gates
 * =====================================================================================================================
 */
/*
 * This function applies the phase gate to the specified qubit in the specified pure quantum state.
 *
 * @param[in,out] state    Address of the state struct the gate is applied to
 * @param[in,out] qubit    Qubit the gate acts on
 */
void applyS(state_t* state, const qubit_t qubit) {
    const dim_t offset = POW2(qubit, dim_t);
    const dim_t blockStride = POW2(qubit + 1, dim_t);
    const dim_t incr = 1;
    const cplx_t phase = I;
    for (dim_t i = 0; i < state->dim; i += blockStride) {
        zscal_(&offset, &phase, state->vec + i + offset, &incr);
    }
}

/*
 * This function applies the inverse phase gate to the specified qubit in the specified pure quantum state.
 *
 * @param[in,out] state    Address of the state struct the gate is applied to
 * @param[in,out] qubit    Qubit the gate acts on
 */
void applySdagger(state_t* state, const qubit_t qubit) {
    const dim_t offset = POW2(qubit, dim_t);
    const dim_t blockStride = POW2(qubit + 1, dim_t);
    const dim_t incr = 1;
    const cplx_t phase = -I;
    for (dim_t i = 0; i < state->dim; i += blockStride) {
        zscal_(&offset, &phase, state->vec + i + offset, &incr);
    }
}

/*
 * This function applies the Hadamard gate to the specified qubit in the specified pure quantum state.
 *
 * @param[in,out] state    Address of the state struct the gate is applied to
 * @param[in,out] qubit    Qubit the gate acts on
 */
void applyH(state_t* state, const qubit_t qubit) {
    const dim_t flipStride = POW2(qubit, dim_t);
    const dim_t blockStride = POW2(qubit + 1, dim_t);
    const dim_t incr = 1;
    const double minus = -1;
    const double c = INVSQRT2;
    const double s = -INVSQRT2;
    for (dim_t i = 0; i < state->dim; i += blockStride) {
        zscal_(&flipStride, &minus, state->vec + i + flipStride, &incr);
        zdrot_(&flipStride, state->vec + i, &incr, state->vec + i + flipStride, &incr, &c, &s);
    }                   // Perform a pi/4 rotation in the plane spanned by the vectors with 0 and 1 at position qubit
}

/*
 * =====================================================================================================================
 *                                                  Hadamard-Y gate
 * =====================================================================================================================
 */
/*
 * This function applies the Hadamard gate to the specified qubit in the specified pure quantum state.
 *
 * @param[in,out] state    Address of the state struct the gate is applied to
 * @param[in,out] qubit    Qubit the gate acts on
 */
void applyHy(state_t* state, const qubit_t qubit) {
    const dim_t blockStride = POW2(qubit + 1, dim_t);
    const dim_t flipStride = POW2(qubit, dim_t);
    const dim_t incr = 1;
    const cplx_t minus = -1;
    const double c = INVSQRT2;
    const cplx_t s = INVSQRT2 * I;

    for (dim_t i = 0; i < state->dim; i += blockStride) {
        zscal_(&flipStride, &minus, state->vec + i + flipStride, &incr);
        zrot_(&flipStride, state->vec + i, &incr, state->vec + i + flipStride, &incr, &c, &s);
    }
}

/*
 * =====================================================================================================================
 *                                                      T gate
 * =====================================================================================================================
 */
/*
 * This function applies the T-gate to the specified qubit in the specified pure quantum state.
 *
 * @param[in,out] state    Address of the state struct the gate is applied to
 * @param[in,out] qubit    Qubit the gate acts on
 */
void applyT(state_t* state, const qubit_t qubit) {
    const dim_t blockStride = POW2(qubit + 1, dim_t);
    const dim_t flipStride = POW2(qubit, dim_t);
    const dim_t incr = 1;
    const cplx_t alpha = INVSQRT2 + INVSQRT2 * I;

    for (dim_t i = 0; i < state->dim; i += blockStride) {
        zscal_(&flipStride, &alpha, state->vec + i + flipStride, &incr);
    }
}

/*
 * This function applies the inverse T-gate to the specified qubit in the specified pure quantum state.
 *
 * @param[in,out] state    Address of the state struct the gate is applied to
 * @param[in,out] qubit    Qubit the gate acts on
 */
void applyTdagger(state_t* state, const qubit_t qubit) {
    const dim_t blockStride = POW2(qubit + 1, dim_t);
    const dim_t flipStride = POW2(qubit, dim_t);
    const dim_t incr = 1;
    const cplx_t alpha = INVSQRT2 - INVSQRT2 * I;

    for (dim_t i = 0; i < state->dim; i += blockStride) {
        zscal_(&flipStride, &alpha, state->vec + i + flipStride, &incr);
    }
}
}

/*
 * =====================================================================================================================
 *                                                      P gate
 * =====================================================================================================================
 */

void applyP(state_t* state, const qubit_t qubit, const double angle) {
    const dim_t blockStride = POW2(qubit + 1, dim_t);
    const dim_t flipStride = POW2(qubit, dim_t);
    const dim_t incr = 1;
    const cplx_t alpha = cos(angle) + sin(angle) * I;

    for (dim_t i = 0; i < state->dim; i += blockStride) {
        zscal_(&flipStride, &alpha, state->vec + i + flipStride, &incr);
    }
}

void applyPdagger(state_t* state, const qubit_t qubit, const double angle) {
    applyP(state, qubit, -angle);
}

/*
 * =====================================================================================================================
 *                                                  Pauli rotation gates
 * =====================================================================================================================
 */
/*
 * This function applies the Pauli-X rotation on the specified qubit in the specified quantum state by the specified
 * angle
 *
 * @param[in,out] state     Address of the state struct the gate is applied to
 * @param[in] qubit         Qubit the gate acts on
 * @param[in] angle         Rotation angle on the Bloch sphere
 */
void applyRX(state_t* state, const qubit_t qubit, const double angle) {
    const dim_t flipStride = POW2(qubit, dim_t);
    const dim_t blockStride = POW2(qubit + 1, dim_t);
    const dim_t incr = 1;
    const double c = cos(angle / 2.);
    const cplx_t s = - I * sin(angle / 2.);

    for (dim_t i = 0; i < state->dim; i += blockStride) {
        zrot_(&flipStride, state->vec + i, &incr, state->vec + i + flipStride, &incr, &c, &s);
    }
}

void applyRXdagger(state_t* state, const qubit_t qubit, const double angle) {
    applyRX(state, qubit, -angle);
}

/*
 * This function applies the Pauli-Y rotation on the specified qubit in the specified quantum state by the specified
 * angle
 *
 * @param[in,out] state     Address of the state struct the gate is applied to
 * @param[in] qubit         Qubit the gate acts on
 * @param[in] angle         Rotation angle on the Bloch sphere
 */
void applyRY(state_t* state, const qubit_t qubit, const double angle) {
    const dim_t flipStride = POW2(qubit, dim_t);
    const dim_t blockStride = POW2(qubit + 1, dim_t);
    const dim_t incr = 1;
    const double c = cos(angle / 2.);
    const double s = sin(angle / 2.);

    for (dim_t i = 0; i < state->dim; i += blockStride) {
        zdrot_(&flipStride, state->vec + i + flipStride, &incr, state->vec + i, &incr, &c, &s);
    }
}

void applyRYdagger(state_t* state, const qubit_t qubit, double const angle) {
    applyRY(state, qubit, -angle);
}

/*
 * This function applies the Pauli-Z rotation on the specified qubit in the specified quantum state by the specified
 * angle
 *
 * @param[in,out] state     Address of the state struct the gate is applied to
 * @param[in] qubit         Qubit the gate acts on
 * @param[in] angle         Rotation angle on the Bloch sphere
 */
void applyRZ(state_t* state, const qubit_t qubit, const double angle) {
    const dim_t flipStride = POW2(qubit, dim_t);
    const dim_t blockStride = POW2(qubit + 1, dim_t);
    const dim_t incr = 1;
    const cplx_t minus = cos(angle / 2.) - I * sin(angle / 2.);
    const cplx_t plus = cos(angle / 2.) + I * sin(angle / 2.);
    for (dim_t i = 0; i < state->dim; i += blockStride) {
        zscal_(&flipStride, &minus, state->vec + i, &incr);
        zscal_(&flipStride, &plus, state->vec + i + flipStride, &incr);
    }
}

void applyRZdagger(state_t* state, const qubit_t qubit, const double angle) {
    applyRZ(state, qubit, -angle);
}

/*
 * =====================================================================================================================
 *                                              controlled Pauli gates
 * =====================================================================================================================
 */
/*
 * This function applies a Pauli-X gate to the target qubit conditioned on the control qubit.
 *
 * @param[in,out] state     Address of the state struct the gate is applied to
 * @param[in] target        Qubit the gate is applied on if control is active
 * @param[in] control       Qubit controlling the application of the gate
 *
 */
void applyCX(state_t* state, const qubit_t target, const qubit_t control) {
    const dim_t flipStride = POW2(target, dim_t);                   // Distance between vector elements differing only
                                                                    // at position qubit
    const dim_t blockStride = POW2(target + 1, dim_t);              // Number of contiguous elements the tensor product
                                                                    // of gate and identities acts on
    const dim_t indices = POW2(state->qubits - 1, dim_t);           // Number of qubits not controlling the gate (free)
    const dim_t incr = 1;

    for (dim_t i = 0; i < indices; i += blockStride) {              // Iterate basis states of free qubits and insert an
                                                                    // active bit at the position control:
        const dim_t j = insone(i, control);

        zswap_(&flipStride, state->vec + j, &incr, state->vec + j + flipStride, &incr);
    }
}

/*
 * This function applies a Pauli-Y gate to the target qubit conditioned on the control qubit.
 *
 * @param[in,out] state     Address of the state struct the gate is applied to
 * @param[in] target        Qubit the gate is applied on if control is active
 * @param[in] control       Qubit controlling the application of the gate
 *
 */
void applyCY(state_t* state, const qubit_t target, const qubit_t control) {
    const dim_t flipStride = POW2(target, dim_t);
    const dim_t blockStride = POW2(target + 1, dim_t);
    const dim_t indices = POW2(state->qubits - 1, dim_t);
    const dim_t incr = 1;
    const cplx_t plus = I;
    const cplx_t minus = -I;

    for (dim_t i = 0; i < indices; i += blockStride) {
        const dim_t j = insone(i, control);

        zswap_(&flipStride, state->vec + j, &incr, state->vec + j + flipStride, &incr);
        zscal_(&flipStride, &minus, state->vec + j, &incr);
        zscal_(&flipStride, &plus, state->vec + j + flipStride, &incr);
    }
}

/*
 * This function applies a Pauli-Z gate to the target qubit conditioned on the control qubit.
 *
 * @param[in,out] state     Address of the state struct the gate is applied to
 * @param[in] target        Qubit the gate is applied on if control is active
 * @param[in] control       Qubit controlling the application of the gate
 *
 */
void applyCZ(state_t* state, const qubit_t target, const qubit_t control) {
    const dim_t flipStride = POW2(target, dim_t);
    const dim_t blockStride = POW2(target + 1, dim_t);
    const dim_t indices = POW2(state->qubits - 1, dim_t);
    const dim_t incr = 1;
    const cplx_t minus = -1;

    for (dim_t i = 0; i < indices; i += blockStride) {
        const dim_t j = insone(i, control);

        zscal_(&flipStride, &minus, state->vec + j + flipStride, &incr);
    }
}

/*
 * =====================================================================================================================
 *                                              Controlled Clifford gates
 * =====================================================================================================================
 */
/*
 * This function applies an S gate to the target qubit conditioned on the control qubit.
 *
 * @param[in,out] state     Address of the state struct the gate is applied to
 * @param[in] target        Qubit the gate is applied on if control is active
 * @param[in] control       Qubit controlling the application of the gate
 *
 */
void applyCS(state_t* state, const qubit_t target, const qubit_t control) {
    const dim_t flipStride = POW2(target, dim_t);
    const dim_t blockStride = POW2(target + 1, dim_t);
    const dim_t indices = POW2(state->qubits - 1, dim_t);
    const dim_t incr = 1;
    const cplx_t alpha = I;

    for (dim_t i = 0; i < indices; i += blockStride) {
        const dim_t j = insone(i, control);

        zscal_(&flipStride, &alpha, state->vec + j + flipStride, &incr);
    }
}

/*
 * This function applies in inverse S gate to the target qubit conditioned on the control qubit.
 *
 * @param[in,out] state     Address of the state struct the gate is applied to
 * @param[in] target        Qubit the gate is applied on if control is active
 * @param[in] control       Qubit controlling the application of the gate
 *
 */
void applyCSdagger(state_t* state, const qubit_t target, const qubit_t control) {
    const dim_t flipStride = POW2(target, dim_t);
    const dim_t blockStride = POW2(target + 1, dim_t);
    const dim_t indices = POW2(state->qubits - 1, dim_t);
    const dim_t incr = 1;
    const cplx_t alpha = -I;

    for (dim_t i = 0; i < indices; i += blockStride) {
        const dim_t j = insone(i, control);

        zscal_(&flipStride, &alpha, state->vec + j + flipStride, &incr);
    }
}

/*
 * This function applies a Hadamard gate to the target qubit conditioned on the control qubit.
 *
 * @param[in,out] state     Address of the state struct the gate is applied to
 * @param[in] target        Qubit the gate is applied on if control is active
 * @param[in] control       Qubit controlling the application of the gate
 *
 */
void applyCH(state_t* state, const qubit_t target, const qubit_t control) {
    const dim_t flipStride = POW2(target, dim_t);
    const dim_t blockStride = POW2(target + 1, dim_t);
    const dim_t indices = POW2(state->qubits - 1, dim_t);
    const dim_t incr = 1;
    const cplx_t minus = -1;
    const double c = INVSQRT2;
    const double s = -INVSQRT2;

    for (dim_t i = 0; i < indices; i += blockStride) {
        const dim_t j = insone(i, control);

        zscal_(&flipStride, &minus, state->vec + j + flipStride, &incr);
        zdrot_(&flipStride, state->vec + j, &incr, state->vec + j + flipStride, &incr, &c, &s);
    }
}

/*
 * =====================================================================================================================
 *                                              Controlled Hadamard-Y gate
 * =====================================================================================================================
 */
/*
 * This function applies a Hadamard-Y gate to the target qubit conditioned on the control qubit.
 *
 * @param[in,out] state     Address of the state struct the gate is applied to
 * @param[in] target        Qubit the gate is applied on if control is active
 * @param[in] control       Qubit controlling the application of the gate
 *
 */
void applyCHy(state_t* state, const qubit_t target, const qubit_t control) {
    const dim_t flipStride = POW2(target, dim_t);
    const dim_t blockStride = POW2(target + 1, dim_t);
    const dim_t indices = POW2(state->qubits - 1, dim_t);
    const dim_t incr = 1;
    const cplx_t minus = -1;
    const double c = INVSQRT2;
    const cplx_t s = INVSQRT2 * I;

    for (dim_t i = 0; i < indices; i += blockStride) {
        const dim_t j = insone(i, control);

        zscal_(&flipStride, &minus, state->vec + j + flipStride, &incr);
        zrot_(&flipStride, state->vec + j, &incr, state->vec + j + flipStride, &incr, &c, &s);
    }
}

/*
 * =====================================================================================================================
 *                                                  Controlled T gate
 * =====================================================================================================================
 */
/*
 * This function applies a T gate gate to the target qubit conditioned on the control qubit.
 *
 * @param[in,out] state     Address of the state struct the gate is applied to
 * @param[in] target        Qubit the gate is applied on if control is active
 * @param[in] control       Qubit controlling the application of the gate
 *
 */
void applyCT(state_t* state, const qubit_t target, const qubit_t control) {
    const dim_t flipStride = POW2(target, dim_t);
    const dim_t blockStride = POW2(target + 1, dim_t);
    const dim_t indices = POW2(state->qubits - 1, dim_t);
    const dim_t incr = 1;
    const cplx_t alpha = INVSQRT2 + INVSQRT2 * I;

    for (dim_t i = 0; i < indices; i += blockStride) {
        const dim_t j = insone(i, control);

        zscal_(&flipStride, &alpha, state->vec + j + flipStride, &incr);
    }
}

/*
 * This function applies an inverse T gate to the target qubit conditioned on the control qubit.
 *
 * @param[in,out] state     Address of the state struct the gate is applied to
 * @param[in] target        Qubit the gate is applied on if control is active
 * @param[in] control       Qubit controlling the application of the gate
 *
 */
void applyCTdagger(state_t* state, const qubit_t target, const qubit_t control) {
    const dim_t flipStride = POW2(target, dim_t);
    const dim_t blockStride = POW2(target + 1, dim_t);
    const dim_t indices = POW2(state->qubits - 1, dim_t);
    const dim_t incr = 1;
    const cplx_t alpha = INVSQRT2 - INVSQRT2 * I;

    for (dim_t i = 0; i < indices; i += blockStride) {
        const dim_t j = insone(i, control);

        zscal_(&flipStride, &alpha, state->vec + j + flipStride, &incr);
    }
}

/*
 * =====================================================================================================================
 *                                                  Controlled P gate
 * =====================================================================================================================
 */
/*
 * This function applies a P gate to the target qubit conditioned on the control qubit.
 *
 * @param[in,out] state     Address of the state struct the gate is applied to
 * @param[in] target        Qubit the gate is applied on if control is active
 * @param[in] control       Qubit controlling the application of the gate
 *
 */
void applyCP(state_t* state, const qubit_t target, const qubit_t control, const double angle) {
    const dim_t flipStride = POW2(target, dim_t);
    const dim_t blockStride = POW2(target + 1, dim_t);
    const dim_t indices = POW2(state->qubits - 1, dim_t);
    const dim_t incr = 1;
    const cplx_t alpha = cos(angle) + sin(angle) * I;

    for (dim_t i = 0; i < indices; i += blockStride) {
        const dim_t j = insone(i, control);

        zscal_(&flipStride, &alpha, state->vec + j + flipStride, &incr);
    }
}

void applyCPdagger(state_t* state, const qubit_t target, const qubit_t control, const double angle) {
    applyCP(state, control, target, -angle);
}

/*
 * =====================================================================================================================
 *                                                      SWAP
 * =====================================================================================================================
 */
/*
 * This function swaps the specified qubits.
 *
 * @param[in,out] state     Address of the state struct the gate is applied to
 * @param[in] qubit1        Qubit to be swapped
 * @param[in] qubit2        Qubits to be swapped
 */
void applySWAP(state_t* state, const qubit_t qubit1, const qubit_t qubit2) {
    const qubit_t left = MAX(qubit1, qubit2);                       // index of the qubit to the left
    const qubit_t right = MIN(qubit1, qubit2);                      // index of the qubit to the right
    const dim_t majorStride = POW2(left + 1, dim_t);                // Number of contiguous elements the tensor product
                                                                    // of gate and identities acts on
    const dim_t ubound = POW2(left, dim_t);                         // Number of contiguous elements actually affected
                                                                    // by the swap (left=0, right=1)
    const dim_t minorStride = POW2(right + 1, dim_t);               // Range of computational basis states that leave
    const dim_t lbound = POW2(right, dim_t);                        // Number of contiguous elements swapped
    const dim_t flipStride = POW2(left, dim_t) - POW2(right, dim_t);    // Stride of basis states with left=0 and
                                                                    // right=1
    const dim_t incr = 1;
    for (dim_t i = lbound; i < state->dim; i += majorStride) {      // Iterate basis states with all zeros to the right
                                                                    // of and including left and right=1
        for (dim_t j = i ; j < i + ubound; j += minorStride) {      // Iterate basis states with fixed qubits to the
                                                                    // left of and including left (left=0)
            zswap_(&lbound, state->vec + j, &incr, state->vec + j + flipStride, &incr);
        }
    }
}

/*
 * =====================================================================================================================
 *                                              Toffoli gate
 * =====================================================================================================================
 */
/*
 * This function applies the Toffoli gate; i.e., a Pauli-X gate to the target controlled by the specified qubits.
 *
 * @param[in,out] state     Address of the state struct the gate is applied to
 * @param[in] target        Qubit the gate is applied to if both control qubits are active
 * @param[in] control1      Control qubit one
 * @param[in] control2      Control qubit two
 */

void applyToffoli(state_t* state, const qubit_t target, const qubit_t control1, const qubit_t control2) {
    const dim_t flipStride = POW2(target, dim_t);                   // Distance between vector elements differing only
                                                                    // at position qubit
    const dim_t blockStride = POW2(target + 1, dim_t);              // Number of contiguous elements the tensor product
                                                                    // of gate and identities acts on
    for (dim_t i = 0; i < state->dim; i += blockStride) {

        for (dim_t j = i; j < i + flipStride; ++j)
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
                SWAP(state->vec + j, state->vec + (j + flipStride), cplx_t);
            }
        }
    }
}