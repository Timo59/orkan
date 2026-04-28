#ifndef STATE_H
#define STATE_H

/*
 * =====================================================================================================================
 * includes
 * =====================================================================================================================
 */

#include "q_types.h"

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
 * Type definitions
 * =====================================================================================================================
 */

/*
 * @brief   Constant indicating the storage type for the quantum state
 */
typedef enum state_type {
    PURE,           /* Pure quantum state: statevector of dimension 2^n */
    MIXED_PACKED,   /* Mixed state: packed lower-triangular column-major storage */
    MIXED_TILED     /* Mixed state: tiled lower-triangular storage for cache efficiency */
} state_type_t;

/*
 * @brief   Structure for multi-qubit quantum state
 */
typedef struct state {
    state_type_t type;  /* Storage type indicator */
    cplx_t *data;       /* State representation array */
    qubit_t qubits;     /* Number of qubits in the quantum system */
} state_t;

/*
 * @brief   Function pointer type for quantum operations
 *
 * A quantum operation modifies a state in-place.  Operations with additional
 * parameters (target qubit, rotation angle, ...) must be wrapped in a
 * function matching this signature.
 */
typedef void (*qop_t)(state_t *state);

/*
 * =====================================================================================================================
 * Public API (dispatches based on state->type)
 * =====================================================================================================================
 */

/*
 * @brief   Returns the size of the array representing the quantum state
 *
 * @param[in]   state   Quantum state (type and qubits must be initialized)
 * @return      Number of complex elements in the state representation
 */
idx_t state_len(const state_t *state);

/*
 * @brief   Initialize a quantum state with the specified number of qubits
 *
 * If data is non-NULL, ownership is transferred (source pointer nullified).
 * If data is NULL, allocates zero-initialized storage.
 *
 * @param[in,out]   state   Quantum state (type must be set)
 * @param[in]       qubits  Number of qubits
 * @param[in,out]   data    Pointer to data array (ownership transferred), or NULL
 */
void state_init(state_t *state, qubit_t qubits, cplx_t **data);

/*
 * @brief   Free allocated memory and reset the quantum state
 *
 * Safe to call multiple times (handles NULL gracefully).
 *
 * @param[in,out]   state   Quantum state
 */
void state_free(state_t *state);

/*
 * @brief   Initialize a quantum state to the uniform superposition (plus state)
 *
 * For PURE: |+>^n with amplitude 1/sqrt(2^n) for each basis state
 * For MIXED_PACKED/MIXED_TILED: maximally mixed state rho = (1/2^n) * |+><+|
 *                                where all elements equal 1/2^n
 *
 * @param[in,out]   state   Quantum state (type must be set)
 * @param[in]       qubits  Number of qubits
 */
void state_plus(state_t *state, qubit_t qubits);

/*
 * @brief   Create a deep copy of the quantum state
 *
 * @param[in]   state   Source quantum state
 * @return      Deep copy; on allocation failure, returned state has data == NULL
 */
state_t state_cp(const state_t *state);

/*
 * @brief   Print the quantum state information and values
 *
 * @param[in]   state   Quantum state
 */
void state_print(const state_t *state);

/*
 * @brief   Get an element from the state representation
 *
 * For PURE: returns state->data[row] (col is ignored)
 * For MIXED: returns element (row, col) from the density matrix,
 *            handling Hermitian symmetry (upper triangle returns conjugate)
 *
 * @param[in]   state   Quantum state
 * @param[in]   row     Row index (or amplitude index for PURE)
 * @param[in]   col     Column index (ignored for PURE)
 * @return      The requested element
 */
cplx_t state_get(const state_t *state, idx_t row, idx_t col);

/*
 * @brief   Set an element in the state representation
 *
 * For PURE: sets state->data[row] = val (col is ignored)
 * For MIXED: sets element (row, col) in the density matrix,
 *            handling Hermitian symmetry (upper triangle stores conjugate)
 *
 * @param[in,out]   state   Quantum state
 * @param[in]       row     Row index (or amplitude index for PURE)
 * @param[in]       col     Column index (ignored for PURE)
 * @param[in]       val     Value to set
 */
void state_set(state_t *state, idx_t row, idx_t col, cplx_t val);

#ifdef __cplusplus
}
#endif

#endif  /* STATE_H */
