// gate.h - Functions representing the action of local quantum gates to a general quantum state

#ifndef GATE_H
#define GATE_H

/*
 * =====================================================================================================================
 * Includes
 * =====================================================================================================================
 */
#include "state.h"

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
 * Pauli gates
 * =====================================================================================================================
 */

void x(state_t* state, const qubit_t target);
void y(state_t* state, const qubit_t target);
void z(state_t* state, const qubit_t target);

/*
 * =====================================================================================================================
 * Clifford gates
 * =====================================================================================================================
 */

void h(state_t* state, const qubit_t target);
void s(state_t* state, const qubit_t target);
void sdg(state_t* state, const qubit_t target);
void t(state_t* state, const qubit_t target);
void tdg(state_t* state, const qubit_t target);
void hy(state_t* state, const qubit_t target);

/*
 * =====================================================================================================================
 * Rotation gates
 * =====================================================================================================================
 */

void rx(state_t* state, const qubit_t target, const double theta);
void ry(state_t* state, const qubit_t target, const double theta);
void rz(state_t* state, const qubit_t target, const double theta);
void p(state_t* state, const qubit_t target, const double theta);

/*
 * =====================================================================================================================
 * Two-qubit gates
 * =====================================================================================================================
 */

void cx(state_t* state, const qubit_t control, const qubit_t target);
void cy(state_t* state, const qubit_t control, const qubit_t target);
void cz(state_t* state, const qubit_t control, const qubit_t target);
void swap_gate(state_t* state, const qubit_t q1, const qubit_t q2);

/*
 * =====================================================================================================================
 * Three-qubit gates
 * =====================================================================================================================
 */

void ccx(state_t* state, const qubit_t ctrl1, const qubit_t ctrl2, const qubit_t target);

#ifdef __cplusplus
}
#endif

#endif // GATE_H
