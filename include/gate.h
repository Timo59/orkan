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
 * @brief   Applies the Pauli-X gate to the qubit <target> in the quantum state <state>
 *                                          | 0   1 |
 *                                      X = |       |
 *                                          | 1   0 |
 * @param[in,out]   state     Quantum state
 * @param[in]       target    Index of the qubit (counted from the right)
 * @return  QS_OK on success, or:
 *          QS_ERR_NULL if state or state->data is NULL
 *          QS_ERR_QUBIT if target >= state->qubits
 */
qs_error_t x(state_t* state, qubit_t target);

/*
 * @brief   Applies the Pauli-Y gate to the qubit <target> in the quantum state <state>
 *                                          | 0  -i |
 *                                      Y = |       |
 *                                          | i   0 |
 * @param[in,out]   state     Quantum state
 * @param[in]       target    Index of the qubit (counted from the right)
 * @return  QS_OK on success, or:
 *          QS_ERR_NULL if state or state->data is NULL
 *          QS_ERR_QUBIT if target >= state->qubits
 */
qs_error_t y(state_t* state, qubit_t target);

/*
 * @brief   Applies the Pauli-Z gate to the qubit <target> in the quantum state <state>
 *                                          | 1   0 |
 *                                      Z = |       |
 *                                          | 0  -1 |
 * @param[in,out]   state     Quantum state
 * @param[in]       target    Index of the qubit (counted from the right)
 * @return  QS_OK on success, or:
 *          QS_ERR_NULL if state or state->data is NULL
 *          QS_ERR_QUBIT if target >= state->qubits
 */
qs_error_t z(state_t* state, qubit_t target);

/*
 * =====================================================================================================================
 * Clifford gates
 * =====================================================================================================================
 */

/*
 * @brief   Applies the Hadamard gate to the qubit <target> in the quantum state <state>
 *                                               1  | 1   1 |
 *                                      H = -------|       |
 *                                              √2  | 1  -1 |
 * @param[in,out]   state     Quantum state
 * @param[in]       target    Index of the qubit (counted from the right)
 * @return  QS_OK on success, or:
 *          QS_ERR_NULL if state or state->data is NULL
 *          QS_ERR_QUBIT if target >= state->qubits
 */
qs_error_t h(state_t* state, qubit_t target);

/*
 * @brief   Applies the S (phase) gate to the qubit <target> in the quantum state <state>
 *                                          | 1   0 |
 *                                      S = |       |
 *                                          | 0   i |
 * @param[in,out]   state     Quantum state
 * @param[in]       target    Index of the qubit (counted from the right)
 * @return  QS_OK on success, or:
 *          QS_ERR_NULL if state or state->data is NULL
 *          QS_ERR_QUBIT if target >= state->qubits
 */
qs_error_t s(state_t* state, qubit_t target);

/*
 * @brief   Applies the S-dagger gate to the qubit <target> in the quantum state <state>
 *                                           | 1   0 |
 *                                      S† = |       |
 *                                           | 0  -i |
 * @param[in,out]   state     Quantum state
 * @param[in]       target    Index of the qubit (counted from the right)
 * @return  QS_OK on success, or:
 *          QS_ERR_NULL if state or state->data is NULL
 *          QS_ERR_QUBIT if target >= state->qubits
 */
qs_error_t sdg(state_t* state, qubit_t target);

/*
 * @brief   Applies the T (π/8) gate to the qubit <target> in the quantum state <state>
 *                                          | 1       0      |
 *                                      T = |                |
 *                                          | 0   e^(iπ/4)   |
 * @param[in,out]   state     Quantum state
 * @param[in]       target    Index of the qubit (counted from the right)
 * @return  QS_OK on success, or:
 *          QS_ERR_NULL if state or state->data is NULL
 *          QS_ERR_QUBIT if target >= state->qubits
 */
qs_error_t t(state_t* state, qubit_t target);

/*
 * @brief   Applies the T-dagger gate to the qubit <target> in the quantum state <state>
 *                                           | 1        0      |
 *                                      T† = |                 |
 *                                           | 0   e^(-iπ/4)   |
 * @param[in,out]   state     Quantum state
 * @param[in]       target    Index of the qubit (counted from the right)
 * @return  QS_OK on success, or:
 *          QS_ERR_NULL if state or state->data is NULL
 *          QS_ERR_QUBIT if target >= state->qubits
 */
qs_error_t tdg(state_t* state, qubit_t target);


#ifdef __cplusplus
}
#endif

#endif // GATE_H