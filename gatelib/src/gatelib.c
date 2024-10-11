//
// Created by Timo Ziegler on 09.10.24.
//

#ifndef GATELIB_H
#include "gatelib.h"
#endif

/*
 * =====================================================================================================================
 *                                                  Pauli gates
 * =====================================================================================================================
 */

void applyX(state_t* state, qubit_t qubit) {
#if defined(GATE_BLAS)
    return applyX_blas(state, qubit);
#elif defined(GATE_OMP)
    return applyX_omp(state, qubit);
#elif defined(GATE_GCD)
    return applyX_gcd(state, qubit);
#elif defined(GATE_BLAS_OMP)
    return applyX_blas_omp(state, qubit);
#else
    return applyX_old(state, qubit);
#endif
}

void applyY(state_t* state, qubit_t qubit) {
#if defined(GATE_BLAS)
    return applyY_blas(state, qubit);
#elif defined(GATE_OMP)
    return applyY_omp(state, qubit);
#elif defined(GATE_GCD)
    return applyY_gcd(state, qubit);
#elif defined(GATE_BLAS_OMP)
    return applyY_blas_omp(state, qubit);
#else
    return applyY_old(state, qubit);
#endif
}

void applyZ(state_t* state, qubit_t qubit) {
#if defined(GATE_BLAS)
    return applyZ_blas(state, qubit);
#elif defined(GATE_OMP)
    return applyZ_omp(state, qubit);
#elif defined(GATE_GCD)
    return applyZ_gcd(state, qubit);
#elif defined(GATE_BLAS_OMP)
    return applyZ_blas_omp(state, qubit);
#else
    return applyZ_old(state, qubit);
#endif
}

/*
 * =====================================================================================================================
 *                                                  Clifford gates
 * =====================================================================================================================
 */

void applyS(state_t* state, qubit_t qubit) {
    return applyS_old(state, qubit);
}

void applySdagger(state_t* state, qubit_t qubit) {
    return applySdagger_old(state, qubit);
}

void applyH(state_t* state, qubit_t qubit) {
    return applyH_old(state, qubit);
}

/*
 * =====================================================================================================================
 *                                                  Hadamard-Y gate
 * =====================================================================================================================
 */

void applyHy(state_t* state, qubit_t qubit) {
    return applyHy_old(state, qubit);
}

/*
 * =====================================================================================================================
 *                                                      T gate
 * =====================================================================================================================
 */

void applyT(state_t* state, qubit_t qubit) {
    return applyT_old(state, qubit);
}

void applyTdagger(state_t* state, qubit_t qubit) {
    return applyTdagger_old(state, qubit);
}

/*
 * =====================================================================================================================
 *                                                      P gate
 * =====================================================================================================================
 */

void applyP(state_t* state, qubit_t qubit, double angle) {
    return applyP_old(state, qubit, angle);
}

void applyPdagger(state_t* state, qubit_t qubit, double angle) {
    return applyPdagger_old(state, qubit, angle);
}

/*
 * =====================================================================================================================
 *                                                  Rotation gates
 * =====================================================================================================================
 */

void applyRX(state_t* state, qubit_t qubit, double angle) {
    return applyRX_old(state, qubit, angle);
}

void applyRXdagger(state_t* state, qubit_t qubit, double angle) {
    return applyRXdagger_old(state, qubit, angle);
}

void applyRY(state_t* state, qubit_t qubit, double angle) {
    return applyRY_old(state, qubit, angle);
}

void applyRYdagger(state_t* state, qubit_t qubit, double angle) {
    return applyRYdagger_old(state, qubit, angle);
}

void applyRZ(state_t* state, qubit_t qubit, double angle) {
    return applyRZ_old(state, qubit, angle);
}

void applyRZdagger(state_t* state, qubit_t qubit, double angle) {
    return applyRZdagger_old(state, qubit, angle);
}

/*
 * =====================================================================================================================
 *                                              Controlled Pauli gates
 * =====================================================================================================================
 */

void applyCX(state_t* state, qubit_t control, qubit_t target) {
    return applyCX_old(state, control, target);
}

void applyCY(state_t* state, qubit_t control, qubit_t target) {
    return applyCY_old(state, control, target);
}

void applyCZ(state_t* state, qubit_t control, qubit_t target) {
    return applyCZ_old(state, control, target);
}

/*
 * =====================================================================================================================
 *                                              Controlled Clifford gates
 * =====================================================================================================================
 */

void applyCS(state_t* state, qubit_t control, qubit_t target) {
    return applyCS_old(state, control, target);
}

void applyCSdagger(state_t* state, qubit_t control, qubit_t target) {
    return applyCSdagger_old(state, control, target);
}

void applyCH(state_t* state, qubit_t control, qubit_t target) {
    return applyCH_old(state, control, target);
}

/*
 * =====================================================================================================================
 *                                              Controlled Hadamard-Y gate
 * =====================================================================================================================
 */

void applyCHy(state_t* state, qubit_t control, qubit_t target) {
    return applyCHy_old(state, control, target);
}

/*
 * =====================================================================================================================
 *                                                  Controlled T gate
 * =====================================================================================================================
 */

void applyCT(state_t* state, qubit_t control, qubit_t target) {
    return applyCT_old(state, control, target);
}

void applyCTdagger(state_t* state, qubit_t control, qubit_t target) {
    return applyCTdagger_old(state, control, target);
}

/*
 * =====================================================================================================================
 *                                                  Controlled P gate
 * =====================================================================================================================
 */

void applyCP(state_t* state, qubit_t control, qubit_t target, double angle) {
    return applyCP_old(state, control, target, angle);
}

void applyCPdagger(state_t* state, qubit_t control, qubit_t target, double angle) {
    return applyCPdagger_old(state, control, target, angle);
}

/*
 * =====================================================================================================================
 *                                                      SWAP gate
 * =====================================================================================================================
 */

void applySWAP(state_t* state, qubit_t qubit1, qubit_t qubit2) {
    return applySWAP_old(state, qubit1, qubit2);
}

/*
 * =====================================================================================================================
 *                                              Toffoli gate
 * =====================================================================================================================
 */

void applyToffoli(state_t* state, qubit_t control1, qubit_t control2, qubit_t target) {
    return applyToffoli_old(state, control1, control2, target);
}
