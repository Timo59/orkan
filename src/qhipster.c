// qhipster.c - Functions representing native (locally acting) quantum gates

#ifndef QHIPSTER_H
#include "qhipster.h"
#endif

#include <math.h>

#include <stdio.h>

#if defined(__APPLE__)
    #include <vecLib/cblas_new.h>
    #include <vecLib/lapack.h>
#elif defined(__linux__)
    #include <cblas.h>
#endif

/*
 * =====================================================================================================================
 * Pauli gates
 * =====================================================================================================================
 */

void applyX(state_t* state, const qubit_t target) {
    // Check that target is in the system's scope
    if (target >= state->qubits) {
        fprintf(stderr, "applyX(): Target is out of scope; Expected %u, Was %u", state->qubits, target);
        state_free(state);
    }

    // Swap entries with target = 0 and 1 but all other qubits fixed
    const dim_t dim = POW2(state->qubits, dim_t);   // Number of entries in the state vector
    const dim_t stride = POW2(target, dim_t);    // Stride between elements only differing in the targeted qubit
    const dim_t step = POW2(target + 1, dim_t);  // Size of the mutually independent blocks

    // Iterate blocks with all qubits left to the addressed qubit fixed
    for (dim_t i = 0; i < dim; i += step) {
        cblas_zswap(stride, state->data + i, 1, state->data + i + stride, 1);
    }
}


void applyY(state_t* state, const qubit_t target) {
    // Check that target is in the system's scope
    if (target >= state->qubits) {
        fprintf(stderr, "applyY(): Target is out of scope; Expected %u, Was %u", state->qubits, target);
        state_free(state);
    }

    const dim_t dim = POW2(state->qubits, dim_t);
    const dim_t stride = POW2(target, dim_t);
    const dim_t step = POW2(target + 1, dim_t);

    const dim_t incr = 1;   // Increment for qubits right to the target
    const cplx_t minus = -1;
    const double c = 0;
    const cplx_t s = I;

    // |...0...> --> -i |...1...> and |...1...> --> i |...0...>
    for (dim_t i = 0; i < dim; i += step) {
        zrot_(&stride, state->data + i, &incr, state->data + i + stride, &incr, &c, &s);
        cblas_zscal(stride, &minus, state->data + i, 1);
    }
}


void applyZ(state_t* state, const qubit_t target) {
    // Check that target is in the system's scope
    if (target >= state->qubits) {
        fprintf(stderr, "applyX(): Target is out of scope; Expected %u, Was %u", state->qubits, target);
        state_free(state);
    }

    const dim_t dim = POW2(state->qubits, dim_t);
    const dim_t stride = POW2(target, dim_t);
    const dim_t step = POW2(target + 1, dim_t);
    const cplx_t minus = -1;

    // |...0...> --> |...0...> and |...1...> --> - |...1...>
    for (dim_t i = 0; i < dim; i += step) {
        cblas_zscal(stride, &minus, state->data + i + stride, 1);
    }
}


/*
 * =====================================================================================================================
 * Clifford gates
 * =====================================================================================================================
 */

/*
void applyS(state_t* state, const qubit_t qubit) {
    const dim_t vecSize = POW2(qubit, dim_t);                       // Number of elements multiplied by I
    const dim_t stepSize = POW2(qubit + 1, dim_t);
    const cplx_t phase = I;
    for (dim_t i = vecSize; i < state->dim; i += stepSize) {        // Only entries with t=1 obtain a phase
        cblas_zscal(vecSize, &phase, state->vec + i, 1);
    }
}

void applySdagger(state_t* state, const qubit_t qubit) {
    const dim_t vecSize = POW2(qubit, dim_t);
    const dim_t stride = POW2(qubit + 1, dim_t);
    const cplx_t phase = -I;
    for (dim_t i = vecSize; i < state->dim; i += stride) {
        cblas_zscal(vecSize, &phase, state->vec + i, 1);
    }
}


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
 */
/*
 * =====================================================================================================================
 *                                                  Hadamard-Y gate
 * =====================================================================================================================
 */

/*
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
*/

/*
 * =====================================================================================================================
 *                                                      T gate
 * =====================================================================================================================
 */

/*
void applyT(state_t* state, const qubit_t qubit) {
    const dim_t stride = POW2(qubit, dim_t);
    const dim_t stepSize = POW2(qubit + 1, dim_t);
    const cplx_t alpha = INVSQRT2 + INVSQRT2 * I;

    for (dim_t i = 0; i < state->dim; i += stepSize) {
        cblas_zscal(stride, &alpha, state->vec + i + stride, 1);
    }
}

void applyTdagger(state_t* state, const qubit_t qubit) {
    const dim_t stride = POW2(qubit, dim_t);
    const dim_t stepSize = POW2(qubit + 1, dim_t);
    const cplx_t alpha = INVSQRT2 - INVSQRT2 * I;

    for (dim_t i = 0; i < state->dim; i += stepSize) {
        cblas_zscal(stride, &alpha, state->vec + i + stride, 1);
    }
}
 */

/*
 * =====================================================================================================================
 *                                                      P gate
 * =====================================================================================================================
 */

/*
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
 */

/*
 * =====================================================================================================================
 *                                                  Pauli rotation gates
 * =====================================================================================================================
 */

/*
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

 */

/*
 * =====================================================================================================================
 *                                              controlled Pauli gates
 * =====================================================================================================================
 */

/*
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
 */

/*
 * =====================================================================================================================
 *                                              Controlled Clifford gates
 * =====================================================================================================================
 */

/*
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
 */

/*
 * =====================================================================================================================
 *                                              Controlled Hadamard-Y gate
 * =====================================================================================================================
 */

/*
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
 */

/*
 * =====================================================================================================================
 *                                                  Controlled T gate
 * =====================================================================================================================
 */

/*
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
 */

/*
 * =====================================================================================================================
 *                                                  Controlled P gate
 * =====================================================================================================================
 */

/*
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
 */

/*
 * =====================================================================================================================
 *                                                      SWAP
 * =====================================================================================================================
 */

/*
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
 */

/*
 * =====================================================================================================================
 *                                              Toffoli gate
 * =====================================================================================================================
 */

/*
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
 */