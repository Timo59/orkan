#ifndef QHIPSTER_H
#define QHIPSTER_H

/*
 * =====================================================================================================================
 *                                                      Includes
 * =====================================================================================================================
 */
#ifndef PSTATE_H
#include "pstate.h"
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
 *                                                  Pauli gates
 * =====================================================================================================================
 */
void applyX(state_t* state, qubit_t qubit);
void applyY(state_t* state, qubit_t qubit);
void applyZ(state_t* state, qubit_t qubit);

/*
 * =====================================================================================================================
 *                                                  Clifford gates
 * =====================================================================================================================
 */
void applyS(state_t* state, qubit_t qubit);
void applySdagger(state_t* state, qubit_t qubit);
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
void applyRswap(state_t* state, qubit_t qubit1, qubit_t qubit2, double angle);

/*
 * =====================================================================================================================
 *                                              Toffoli gate
 * =====================================================================================================================
 */
void applyToffoli(state_t* state, qubit_t target, qubit_t control1, qubit_t control2);

#ifdef __cplusplus
}
#endif

#endif /* QHIPSTER_H */