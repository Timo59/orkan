// qhipster.c - Functions representing native (locally acting) quantum gates

#ifndef GATE_H
#define GATE_H

/*
 * =====================================================================================================================
 * Includes
 * =====================================================================================================================
 */
#ifndef STATE_H
#include "state.h"
#endif

/*
 * =====================================================================================================================
 * C++ check
 * =====================================================================================================================
 */
#ifdef __cplusplus
extern "C" {
#endif

/*
 * =====================================================================================================================
 * Fortran BLAS function zrot_
 * =====================================================================================================================
 */
#ifdef LINUX
void zrot_(
    const blasint *n,
    openblas_complex_double *x, const blasint *incx,
    openblas_complex_double *y, const blasint *incy,
    const double *c, const openblas_complex_double *s);
#endif

/*
 * =====================================================================================================================
 * Pauli gates
 * =====================================================================================================================
 */

/*
 * @brief   Applies the Pauli-X gate to the qubit <qubit> in the quantum state <state>
 *                                          | 0   1 |
 *                                      X = |       |
 *                                          | 1   0 |
 * @param[in,out]   state     Quantum state
 * @param[in]       qubit     Index of the qubit (counted from the tight)
 */
void x_pure(state_t *state, qubit_t target);
void x_mixed(state_t *state, qubit_t target);
inline void x(state_t* state, const qubit_t target) {
    if (state->type == PURE) x_pure(state, target);
    else x_mixed(state, target);
}


/*
 * @brief   Applies the Pauli-Y gate to the qubit <qubit> in the quantum state <state>
 *                                          | 0  -i |
 *                                      Y = |       |
 *                                          | i   0 |
 * @param[in,out]   state     Quantum state
 * @param[in]       qubit     Index of the qubit (counted from the tight)
 */
void applyY(state_t* state, qubit_t target);


/*
 * @brief   Applies the Pauli-Z gate to the qubit <qubit> in the quantum state <state>
 *                                          | 1   0 |
 *                                      Z = |       |
 *                                          | 0  -1 |
 * @param[in,out]   state     Quantum state
 * @param[in]       qubit     Index of the qubit (counted from the tight)
 */
void applyZ(state_t* state, qubit_t target);

/*
 * =====================================================================================================================
 *                                                  Clifford gates
 * =====================================================================================================================
 */

/*
 * @brief   Applies the phase gate to the qubit <qubit> in the pure quantum state <state>
 *                                          | 1   0 |
 *                                      S = |       |
 *                                          | 0   i |
 *
 * @param[in,out]   state    Quantum state
 * @param[in]       qubit    Index of the qubit (counted from the tight)
 */
void applyS(state_t* state, qubit_t qubit);

/*
 * @brief   Applies the inverse phase gate to the qubit <qubit> in the pure quantum state <state>
 *                                          | 1   0 |
 *                                   S^-1 = |       |
 *                                          | 0  -i |
 *
 * @param[in,out]   state    Quantum state
 * @param[in]       qubit    Index of the qubit (counted from the tight)
 */
void applySdagger(state_t* state, qubit_t qubit);


/*
 * @brief   Applies the inverse phase gate to the qubit <qubit> in the pure quantum state <state>
 *                                          | 1   1 |
 *                            H = 1/sqrt(2) |       |
 *                                          | 1  -1 |
 *
 * @param[in,out]   state    Quantum state
 * @param[in]       qubit    Index of the qubit (counted from the tight)
 */
void applyH(state_t* state, qubit_t qubit);

/*
 * =====================================================================================================================
 *                                                  Hadamard-Y gate
 * =====================================================================================================================
 */
void applyHy(state_t* state, qubit_t qubit);

/*
 * =====================================================================================================================
 *                                                      T gate
 * =====================================================================================================================
 */
void applyT(state_t* state, qubit_t qubit);
void applyTdagger(state_t* state, qubit_t qubit);

/*
 * =====================================================================================================================
 *                                                      P gate
 * =====================================================================================================================
 */
void applyP(state_t* state, qubit_t qubit, double angle);
void applyPdagger(state_t* state, qubit_t qubit, double angle);

/*
 * =====================================================================================================================
 *                                                  Rotation gates
 * =====================================================================================================================
 */
void applyRX(state_t* state, qubit_t qubit, double angle);
void applyRXdagger(state_t* state, qubit_t qubit, double angle);
void applyRY(state_t* state, qubit_t qubit, double angle);
void applyRYdagger(state_t* state, qubit_t qubit, double angle);
void applyRZ(state_t* state, qubit_t qubit, double angle);
void applyRZdagger(state_t* state, qubit_t qubit, double angle);

/*
 * =====================================================================================================================
 *                                              Controlled Pauli gates
 * =====================================================================================================================
 */
void applyCX(state_t* state, qubit_t target, qubit_t control);
void applyCY(state_t* state, qubit_t target, qubit_t control);
void applyCZ(state_t* state, qubit_t target, qubit_t control);

/*
 * =====================================================================================================================
 *                                              Controlled Clifford gates
 * =====================================================================================================================
 */
void applyCS(state_t* state, qubit_t target, qubit_t control);
void applyCSdagger(state_t* state, qubit_t target, qubit_t control);
void applyCH(state_t* state, qubit_t target, qubit_t control);

/*
 * =====================================================================================================================
 *                                              Controlled Hadamard-Y gate
 * =====================================================================================================================
 */
void applyCHy(state_t* state, qubit_t target, qubit_t control);

/*
 * =====================================================================================================================
 *                                                  Controlled T gate
 * =====================================================================================================================
 */
void applyCT(state_t* state, qubit_t target, qubit_t control);
void applyCTdagger(state_t* state, qubit_t target, qubit_t control);

/*
 * =====================================================================================================================
 *                                                  Controlled P gate
 * =====================================================================================================================
 */
void applyCP(state_t* state, qubit_t target, qubit_t control, double angle);
void applyCPdagger(state_t* state, qubit_t target, qubit_t control, double angle);

/*
 * =====================================================================================================================
 *                                                      SWAP gate
 * =====================================================================================================================
 */
void applySWAP(state_t* state, qubit_t qubit1, qubit_t qubit2);
void applyRSWAP(state_t* state, qubit_t qubit1, qubit_t qubit2, double angle);

/*
 * =====================================================================================================================
 *                                              Toffoli gate
 * =====================================================================================================================
 */
void applyToffoli(state_t* state, qubit_t target, qubit_t control1, qubit_t control2);

#ifdef __cplusplus
}
#endif

#endif // GATE_H