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
#include "../../include/state.h"
#endif

#ifdef MACOS
#include <vecLib/cblas_new.h>
#include <vecLib/lapack.h>
#else
#include <cblas.h>
#endif

#ifndef _STDIO_H_
#include <stdio.h>
#endif

#ifndef QHIPSTER_H
#include "qhipster.h"
#endif

/*
 * =====================================================================================================================
 *                                                  Pauli gates
 * =====================================================================================================================
 */
/*
 * This function applies the Pauli-X gate to the specified qubit in the specified pure quantum state.
 *                                          | 0   1 |
 *                                      X = |       |
 *                                          | 1   0 |
 * @param[in,out] state     Address of the state struct the gate is applied to
 * @param[in] qubit         Qubit the gate acts on
 */
void applyX(state_t* state, const qubit_t qubit) {
    const dim_t stride = POW2(qubit, dim_t);                        // Stride between elements the tensor product builds
                                                                    // linear combinations from
    const dim_t stepSize = POW2(qubit + 1, dim_t);                  // Number of contiguous elements the linear
                                                                    // combinations induced by the tensor product are
                                                                    // restricted to
    for (dim_t i = 0; i < state->dim; i += stepSize) {              // Iterate blocks with all qubits left to the
                                                                    // addressed qubit fixed
        cblas_zswap(stride, state->vec + i, 1, state->vec + i + stride, 1);
    }
}

/*
 * This function applies the Pauli-Y gate to the specified qubit in the specified pure quantum state.
 *                                          | 0  -i |
 *                                      Y = |       |
 *                                          | i   0 |
 * @param[in,out] state     Address of the state struct the gate is applied to
 * @param[in] qubit         Qubit the gate acts on
 */
void applyY(state_t* state, const qubit_t qubit) {
    const dim_t stride = POW2(qubit, dim_t);
    const dim_t stepSize = POW2(qubit + 1, dim_t);
    const dim_t incr = 1;
    const cplx_t minus = -1;
    const double c = 0;
    const cplx_t s = I;

    for (dim_t i = 0; i < state->dim; i += stepSize) {
        zrot_(&stride, state->vec + i, &incr, state->vec + i + stride, &incr, &c, &s);
        cblas_zscal(stride, &minus, state->vec + i, 1);
    }
}

/*
 * This function applies the Pauli-Z gate to the specified qubit in the specified pure quantum state.
 *                                          | 1   0 |
 *                                      Z = |       |
 *                                          | 0  -1 |
 * @param[in,out] state     Address of the state struct the gate is applied to
 * @param[in] qubit         Qubit the gate acts on
 */
void applyZ(state_t* state, const qubit_t qubit) {
    const dim_t stride = POW2(qubit, dim_t);
    const dim_t stepSize = POW2(qubit + 1, dim_t);
    const cplx_t minus = -1;

    for (dim_t i = 0; i < state->dim; i += stepSize) {
        cblas_zscal(stride, &minus, state->vec + i + stride, 1);
    }
}

/*
 * =====================================================================================================================
 *                                                  Clifford gates
 * =====================================================================================================================
 */
/*
 * This function applies the phase gate to the specified qubit in the specified pure quantum state.
 *                                          | 1   0 |
 *                                      S = |       |
 *                                          | 0   i |
 *
 * @param[in,out] state    Address of the state struct the gate is applied to
 * @param[in,out] qubit    Qubit the gate acts on
 */
void applyS(state_t* state, const qubit_t qubit) {
    const dim_t vecSize = POW2(qubit, dim_t);                       // Number of elements multiplied by I
    const dim_t stepSize = POW2(qubit + 1, dim_t);
    const cplx_t phase = I;
    for (dim_t i = vecSize; i < state->dim; i += stepSize) {        // Only entries with t=1 obtain a phase
        cblas_zscal(vecSize, &phase, state->vec + i, 1);
    }
}

/*
 * This function applies the inverse phase gate to the specified qubit in the specified pure quantum state.
 *
 *                                          | 1   0 |
 *                                   S^-1 = |       |
 *                                          | 0  -i |
 * @param[in,out] state    Address of the state struct the gate is applied to
 * @param[in,out] qubit    Qubit the gate acts on
 */
void applySdagger(state_t* state, const qubit_t qubit) {
    const dim_t vecSize = POW2(qubit, dim_t);
    const dim_t stride = POW2(qubit + 1, dim_t);
    const cplx_t phase = -I;
    for (dim_t i = vecSize; i < state->dim; i += stride) {
        cblas_zscal(vecSize, &phase, state->vec + i, 1);
    }
}

/*
 * This function applies the Hadamard gate to the specified qubit in the specified pure quantum state.
 *                                             1    | 1   1 |
 *                                      H = ------- |       |
 *                                          sqrt(2) | 1  -1 |
 * @param[in,out] state    Address of the state struct the gate is applied to
 * @param[in,out] qubit    Qubit the gate acts on
 */
void applyH(state_t* state, const qubit_t qubit) {
    const dim_t stride = POW2(qubit, dim_t);
    const dim_t stepSize = POW2(qubit + 1, dim_t);
    const cplx_t minus = -1;
    const double c = INVSQRT2;
    const double s = -INVSQRT2;
    for (dim_t i = 0; i < state->dim; i += stepSize) {
        cblas_zscal(stride, &minus, state->vec + i + stride, 1);
        cblas_zdrot(stride, state->vec + i, 1, state->vec + i + stride, 1, c, s);
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
    const dim_t stride = POW2(qubit, dim_t);
    const dim_t stepSize = POW2(qubit + 1, dim_t);
    const dim_t incr = 1;
    const cplx_t minus = -1;
    const double c = INVSQRT2;
    const cplx_t s = INVSQRT2 * I;

    for (dim_t i = 0; i < state->dim; i += stepSize) {
        cblas_zscal(stride, &minus, state->vec + i + stride, 1);
        zrot_(&stride, state->vec + i, &incr, state->vec + i + stride, &incr, &c, &s);
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
    const dim_t stride = POW2(qubit, dim_t);
    const dim_t stepSize = POW2(qubit + 1, dim_t);
    const cplx_t alpha = INVSQRT2 + INVSQRT2 * I;

    for (dim_t i = 0; i < state->dim; i += stepSize) {
        cblas_zscal(stride, &alpha, state->vec + i + stride, 1);
    }
}

/*
 * This function applies the inverse T-gate to the specified qubit in the specified pure quantum state.
 *
 * @param[in,out] state    Address of the state struct the gate is applied to
 * @param[in,out] qubit    Qubit the gate acts on
 */
void applyTdagger(state_t* state, const qubit_t qubit) {
    const dim_t stride = POW2(qubit, dim_t);
    const dim_t stepSize = POW2(qubit + 1, dim_t);
    const cplx_t alpha = INVSQRT2 - INVSQRT2 * I;

    for (dim_t i = 0; i < state->dim; i += stepSize) {
        cblas_zscal(stride, &alpha, state->vec + i + stride, 1);
    }
}

/*
 * =====================================================================================================================
 *                                                      P gate
 * =====================================================================================================================
 */
/*
 * This function applies the P-gate to the specified qubit in the specified pure quantum state.
 *
 * @param[in,out] state     Address of the state struct the gate is applied to
 * @param[in] qubit         Qubit the gate acts on
 * @param[in] angle
 */
void applyP(state_t* state, const qubit_t qubit, const double angle) {
    const dim_t stride = POW2(qubit, dim_t);
    const dim_t stepSize = POW2(qubit + 1, dim_t);
    const cplx_t alpha = cos(angle) + sin(angle) * I;

    for (dim_t i = 0; i < state->dim; i += stepSize) {
        cblas_zscal(stride, &alpha, state->vec + i + stride, 1);
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
    const dim_t stride = POW2(qubit, dim_t);
    const dim_t stepSize = POW2(qubit + 1, dim_t);
    const dim_t incr = 1;
    const double c = cos(angle / 2.);
    const cplx_t s = - I * sin(angle / 2.);

    for (dim_t i = 0; i < state->dim; i += stepSize) {
        zrot_(&stride, state->vec + i, &incr, state->vec + i + stride, &incr, &c, &s);
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
    const dim_t stride = POW2(qubit, dim_t);
    const dim_t stepSize = POW2(qubit + 1, dim_t);
    const dim_t incr = 1;
    const double c = cos(angle / 2.);
    const double s = sin(angle / 2.);

    for (dim_t i = 0; i < state->dim; i += stepSize) {
        cblas_zdrot(stride, state->vec + i + stride, 1, state->vec + i, 1, c, s);
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
    const dim_t stride = POW2(qubit, dim_t);
    const dim_t stepSize = POW2(qubit + 1, dim_t);
    const dim_t incr = 1;
    const cplx_t minus = cos(angle / 2.) - I * sin(angle / 2.);
    const cplx_t plus = cos(angle / 2.) + I * sin(angle / 2.);

    for (dim_t i = 0; i < state->dim; i += stepSize) {
        cblas_zscal(stride, &minus, state->vec + i, 1);
        cblas_zscal(stride, &plus, state->vec + i + stride, 1);
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
    dim_t outerStep, innerStep, innerBound, vecSize;
    if (target > control) {                                         // If t is left of c, outer loop iterates qubit con-
        outerStep = POW2(target + 1, dim_t);                        // figurations left of t, inner loop iterates qubit
        innerStep = POW2(control + 1, dim_t);                       // configurations right of t and left of c and
        innerBound = POW2(target, dim_t);                           // vector over elements right of c are swapped
        vecSize = POW2(control, dim_t);
    }
    else if (target < control) {                                    // If t is right of c, outer loop iterates qubit
        outerStep = POW2(control + 1, dim_t);                       // configurations left of c, inner loop interates
        innerStep = POW2(target + 1, dim_t);                        // qubit configurations right of c and left of t and
        innerBound = POW2(control, dim_t);                          // vector over elements right of t are swapped
        vecSize = POW2(target, dim_t);
    }
    else {
        fprintf(stderr, "applyCX: Target and control qubit coincide\n");
        return;
    }
    const dim_t stride = POW2(target, dim_t);                       // Stride between elements to be swapped

    for (dim_t i = 1L << control; i < state->dim; i += outerStep) { // Outer loop ensures that control qubit is active
        for (dim_t j = i; j < i + innerBound; j += innerStep) {
            cblas_zswap(vecSize, state->vec + j, 1, state->vec + j + stride, 1);
        }
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
    dim_t outerStep, innerStep, innerBound, vecSize;
    if (target > control) {
        outerStep = POW2(target + 1, dim_t);
        innerStep = POW2(control + 1, dim_t);
        innerBound = POW2(target, dim_t);
        vecSize = POW2(control, dim_t);
    }
    else if (target < control) {
        outerStep = POW2(control + 1, dim_t);
        innerStep = POW2(target + 1, dim_t);
        innerBound = POW2(control, dim_t);
        vecSize = POW2(target, dim_t);
    }
    else {
        fprintf(stderr, "applyCY: Target and control qubit coincide\n");
        return;
    }
    const dim_t stride = POW2(target, dim_t);
    const dim_t incr = 1;
    const cplx_t minus = -1;
    const double c = 0;
    const cplx_t s = I;

    for (dim_t i = 1L << control; i < state->dim; i += outerStep) {
        for (dim_t j = i; j < i + innerBound; j += innerStep) {
            zrot_(&vecSize, state->vec + j, &incr, state->vec + j + stride, &incr, &c, &s);
            cblas_zscal(vecSize, &minus, state->vec + j, 1);
        }
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
    dim_t outerStep, innerStep, innerBound, vecSize;
    if (target > control) {
        outerStep = POW2(target + 1, dim_t);
        innerStep = POW2(control + 1, dim_t);
        innerBound = POW2(target, dim_t);
        vecSize = POW2(control, dim_t);
    }
    else if (target < control) {
        outerStep = POW2(control + 1, dim_t);
        innerStep = POW2(target + 1, dim_t);
        innerBound = POW2(control, dim_t);
        vecSize = POW2(target, dim_t);
    }
    else {
        fprintf(stderr, "applyCZ: Target and control qubit coincide\n");
        return;
    }
    const dim_t stride = POW2(target, dim_t);
    const cplx_t minus = -1;

    for (dim_t i = 1L << control; i < state->dim; i += outerStep) {
        for (dim_t j = i; j < i + innerBound; j += innerStep) {
            cblas_zscal(vecSize, &minus, state->vec + j + stride, 1);
        }
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
    dim_t outerStep, innerStep, innerBound, vecSize;
    if (target > control) {
        outerStep = POW2(target + 1, dim_t);
        innerStep = POW2(control + 1, dim_t);
        innerBound = POW2(target, dim_t);
        vecSize = POW2(control, dim_t);
    }
    else if (target < control) {
        outerStep = POW2(control + 1, dim_t);
        innerStep = POW2(target + 1, dim_t);
        innerBound = POW2(control, dim_t);
        vecSize = POW2(target, dim_t);
    }
    else {
        fprintf(stderr, "applyCS: Target and control qubit coincide\n");
        return;
    }
    const dim_t stride = POW2(target, dim_t);
    const cplx_t alpha = I;

    for (dim_t i = 1L << control; i < state->dim; i += outerStep) {
        for (dim_t j = i; j < i + innerBound; j += innerStep) {
            cblas_zscal(vecSize, &alpha, state->vec + j + stride, 1);
        }
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
    dim_t outerStep, innerStep, innerBound, vecSize;
    if (target > control) {
        outerStep = POW2(target + 1, dim_t);
        innerStep = POW2(control + 1, dim_t);
        innerBound = POW2(target, dim_t);
        vecSize = POW2(control, dim_t);
    }
    else if (target < control) {
        outerStep = POW2(control + 1, dim_t);
        innerStep = POW2(target + 1, dim_t);
        innerBound = POW2(control, dim_t);
        vecSize = POW2(target, dim_t);
    }
    else {
        fprintf(stderr, "applySdagger: Target and control qubit coincide\n");
        return;
    }
    const dim_t stride = POW2(target, dim_t);
    const cplx_t alpha = -I;

    for (dim_t i = 1L << control; i < state->dim; i += outerStep) {
        for (dim_t j = i; j < i + innerBound; j += innerStep) {
            cblas_zscal(vecSize, &alpha, state->vec + j + stride, 1);
        }
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
    dim_t outerStep, innerStep, innerBound, vecSize;
    if (target > control) {
        outerStep = POW2(target + 1, dim_t);
        innerStep = POW2(control + 1, dim_t);
        innerBound = POW2(target, dim_t);
        vecSize = POW2(control, dim_t);
    }
    else if (target < control) {
        outerStep = POW2(control + 1, dim_t);
        innerStep = POW2(target + 1, dim_t);
        innerBound = POW2(control, dim_t);
        vecSize = POW2(target, dim_t);
    }
    else {
        fprintf(stderr, "applyCH: Target and control qubit coincide\n");
        return;
    }
    const dim_t stride = POW2(target, dim_t);
    const cplx_t minus = -1;

    for (dim_t i = 1L << control; i < state->dim; i += outerStep) {
        for (dim_t j = i; j < i + innerBound; j += innerStep) {
            cblas_zscal(vecSize, &minus, state->vec + j + stride, 1);
            cblas_zdrot(vecSize, state->vec + j, 1, state->vec + j + stride, 1, INVSQRT2, -INVSQRT2);
        }
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
    dim_t outerStep, innerStep, innerBound, vecSize;
    if (target > control) {
        outerStep = POW2(target + 1, dim_t);
        innerStep = POW2(control + 1, dim_t);
        innerBound = POW2(target, dim_t);
        vecSize = POW2(control, dim_t);
    }
    else if (target < control) {
        outerStep = POW2(control + 1, dim_t);
        innerStep = POW2(target + 1, dim_t);
        innerBound = POW2(control, dim_t);
        vecSize = POW2(target, dim_t);
    }
    else {
        fprintf(stderr, "applyCHy: Target and control qubit coincide\n");
        return;
    }
    const dim_t stride = POW2(target, dim_t);
    const dim_t incr = 1;
    const cplx_t minus = -1;
    const double c = INVSQRT2;
    const cplx_t s = INVSQRT2 * I;

    for (dim_t i = 1L << control; i < state->dim; i += outerStep) {
        for (dim_t j = i; j < i + innerBound; j += innerStep) {
            cblas_zscal(vecSize, &minus, state->vec + j + stride, 1);
            zrot_(&vecSize, state->vec + j, &incr, state->vec + j + stride, &incr, &c, &s);
        }
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
    dim_t outerStep, innerStep, innerBound, vecSize;
    if (target > control) {
        outerStep = POW2(target + 1, dim_t);
        innerStep = POW2(control + 1, dim_t);
        innerBound = POW2(target, dim_t);
        vecSize = POW2(control, dim_t);
    }
    else if (target < control) {
        outerStep = POW2(control + 1, dim_t);
        innerStep = POW2(target + 1, dim_t);
        innerBound = POW2(control, dim_t);
        vecSize = POW2(target, dim_t);
    }
    else {
        fprintf(stderr, "applyCT: Target and control qubit coincide\n");
        return;
    }
    const dim_t stride = POW2(target, dim_t);
    const cplx_t alpha = INVSQRT2 + INVSQRT2 * I;

    for (dim_t i = 1L << control; i < state->dim; i += outerStep) {
        for (dim_t j = i; j < i + innerBound; j += innerStep) {
            cblas_zscal(vecSize, &alpha, state->vec + j + stride, 1);
        }
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
    dim_t outerStep, innerStep, innerBound, vecSize;
    if (target > control) {
        outerStep = POW2(target + 1, dim_t);
        innerStep = POW2(control + 1, dim_t);
        innerBound = POW2(target, dim_t);
        vecSize = POW2(control, dim_t);
    }
    else if (target < control) {
        outerStep = POW2(control + 1, dim_t);
        innerStep = POW2(target + 1, dim_t);
        innerBound = POW2(control, dim_t);
        vecSize = POW2(target, dim_t);
    }
    else {
        fprintf(stderr, "applyCTdagger: Target and control qubit coincide\n");
        return;
    }
    const dim_t stride = POW2(target, dim_t);
    const cplx_t alpha = INVSQRT2 - INVSQRT2 * I;

    for (dim_t i = 1L << control; i < state->dim; i += outerStep) {
        for (dim_t j = i; j < i + innerBound; j += innerStep) {
            cblas_zscal(vecSize, &alpha, state->vec + j + stride, 1);
        }
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
    dim_t outerStep, innerStep, innerBound, vecSize;
    if (target > control) {
        outerStep = POW2(target + 1, dim_t);
        innerStep = POW2(control + 1, dim_t);
        innerBound = POW2(target, dim_t);
        vecSize = POW2(control, dim_t);
    } else if (target < control) {
        outerStep = POW2(control + 1, dim_t);
        innerStep = POW2(target + 1, dim_t);
        innerBound = POW2(control, dim_t);
        vecSize = POW2(target, dim_t);
    } else {
        fprintf(stderr, "applyCP: Target and control qubit coincide\n");
        return;
    }
    const dim_t stride = POW2(target, dim_t);
    const cplx_t alpha = cos(angle) + sin(angle) * I;

    for (dim_t i = 1L << control; i < state->dim; i += outerStep) {
        for (dim_t j = i; j < i + innerBound; j += innerStep) {
            cblas_zscal(vecSize, &alpha, state->vec + j + stride, 1);
        }
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
    const dim_t outerStep = POW2(left + 1, dim_t);                  // Number of contiguous elements the tensor product
                                                                    // of gate and identities acts on
    const dim_t innerBound = POW2(left, dim_t);                     // Number of contiguous elements actually affected
                                                                    // by the swap (left=0, right=1)
    const dim_t innerStep = POW2(right + 1, dim_t);                 // Range of computational basis states that leave
    const dim_t stride = POW2(right, dim_t);                        // Number of contiguous elements swapped

    for (dim_t i = 0; i < state->dim; i += outerStep) {             // Iterate basis states with all zeros to the right
                                                                    // of and including left and right=1
        for (dim_t j = i ; j < i + innerBound; j += innerStep) {        // Iterate basis states with fixed qubits to the
                                                                        // left of and including left (left=0)
            cblas_zswap(stride, state->vec + j + stride, 1, state->vec + j + innerBound, 1);
        }
    }
}
/*
 * This function performs an imaginary time evolution of the specified state by a swap gate the specified qubits.
 *
 * @param[in,out] state     Address of the state struct the gate is applied to
 * @param[in] qubit1        Qubit to be swapped
 * @param[in] qubit2        Qubits to be swapped
 */
void applyRSWAP(state_t* state, qubit_t qubit1, qubit_t qubit2, double angle) {
    const qubit_t left = MAX(qubit1, qubit2);                       // index of the qubit to the left
    const qubit_t right = MIN(qubit1, qubit2);                      // index of the qubit to the right
    const dim_t outerStep = POW2(left + 1, dim_t);                  // Number of contiguous elements the tensor product
                                                                    // of gate and identities acts on
    const dim_t innerBound = POW2(left, dim_t);                     // Number of contiguous elements actually affected
                                                                    // by the swap (left=0, right=1)
    const dim_t innerStep = POW2(right + 1, dim_t);                 // Range of computational basis states that leave
    const dim_t stride = POW2(right, dim_t);                        // Number of contiguous elements swapped
    const dim_t incr = 1;
    const double c = cos(angle / 2.);
    const cplx_t s = -I * sin(angle / 2.);
    const cplx_t e = (cplx_t) c + s;                                // cos(angle/2) - i sin(angle/2)

    for (dim_t i = 0; i < state->dim; i += outerStep) {             // Iterate basis states with all zeros to the right
                                                                    // of and including left and right=1
        for (dim_t j = i ; j < i + innerBound; j += innerStep) {    // Iterate basis states with fixed qubits to the
                                                                    // left of and including left (left=0)
           cblas_zscal(stride, &e, state->vec + j, 1);                          // Multiply bases with q1=0 and q2=0 by
                                                                                // e^-i*angle
           cblas_zscal(stride, &e, state->vec + j + stride + innerBound, 1);    // Multiply bases with q1=1 and q2=1 by
                                                                                // e^-i*angle
           zrot_(&stride, state->vec + j + stride, &incr, state->vec + j + innerBound, &incr, &c, &s);
                                                                    // Perform a rotation in the complex plain of bases
                                                                    // with q1=0, q2=1 and q1=1, q2=0, respectively
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
    const dim_t stride = POW2(target, dim_t);                   // Distance between vector elements differing only
                                                                    // at position qubit
    const dim_t stepSize = POW2(target + 1, dim_t);              // Number of contiguous elements the tensor product
                                                                    // of gate and identities acts on
    for (dim_t i = 0; i < state->dim; i += stepSize) {
        for (dim_t j = i; j < i + stride; ++j) {
            if (j & POW2(control1, dim_t) && j & POW2(control2, dim_t)) {
                SWAP(state->vec + j, state->vec + j + stride, cplx_t);
            }
        }
    }
}