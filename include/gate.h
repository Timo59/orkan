// gate.h   - Functions representing the action of local quantum gates to a general quantum state

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

#include <stdio.h>
#include <stdlib.h>

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
 * Error handling
 * =====================================================================================================================
 */

#define GATE_VALIDATE(cond, msg) do {                           \
    if (!(cond)) {                                              \
        fprintf(stderr, "qlib: gate: %s\n", (msg));            \
        exit(EXIT_FAILURE);                                     \
    }                                                           \
} while(0)

/*
 * =====================================================================================================================
 * Pauli gates
 * =====================================================================================================================
 */

void x(state_t* state, qubit_t target);
void y(state_t* state, qubit_t target);
void z(state_t* state, qubit_t target);

/*
 * =====================================================================================================================
 * Clifford gates
 * =====================================================================================================================
 */

void h(state_t* state, qubit_t target);
void s(state_t* state, qubit_t target);
void sdg(state_t* state, qubit_t target);
void t(state_t* state, qubit_t target);
void tdg(state_t* state, qubit_t target);
void hy(state_t* state, qubit_t target);

/*
 * =====================================================================================================================
 * Rotation gates
 * =====================================================================================================================
 */

void rx(state_t* state, qubit_t target, double theta);
void ry(state_t* state, qubit_t target, double theta);
void rz(state_t* state, qubit_t target, double theta);
void p(state_t* state, qubit_t target, double theta);

/*
 * =====================================================================================================================
 * Two-qubit gates
 * =====================================================================================================================
 */

void cx(state_t* state, qubit_t control, qubit_t target);
void cy(state_t* state, qubit_t control, qubit_t target);
void cz(state_t* state, qubit_t control, qubit_t target);
void cs(state_t* state, qubit_t control, qubit_t target);
void csdg(state_t* state, qubit_t control, qubit_t target);
void ch(state_t* state, qubit_t control, qubit_t target);
void chy(state_t* state, qubit_t control, qubit_t target);
void ct(state_t* state, qubit_t control, qubit_t target);
void ctdg(state_t* state, qubit_t control, qubit_t target);
void cp(state_t* state, qubit_t control, qubit_t target, double theta);
void cpdg(state_t* state, qubit_t control, qubit_t target, double theta);
void swap_gate(state_t* state, qubit_t q1, qubit_t q2);

/*
 * =====================================================================================================================
 * Three-qubit gates
 * =====================================================================================================================
 */

void ccx(state_t* state, qubit_t ctrl1, qubit_t ctrl2, qubit_t target);

#ifdef __cplusplus
}
#endif

#endif // GATE_H
