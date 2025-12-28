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


#ifdef __cplusplus
}
#endif

#endif // GATE_H