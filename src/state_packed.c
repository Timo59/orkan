/**
 * @file state_packed.c
 * @brief Packed lower-triangular storage for mixed quantum states (density matrices)
 *
 * Mixed states are represented as density matrices rho of dimension 2^n x 2^n.
 * Since density matrices are Hermitian (rho = rho^dagger), we only store the
 * lower triangle in column-major (LAPACK-compatible) packed format.
 *
 * Storage layout for an n x n Hermitian matrix:
 * - Column 0: elements (0,0), (1,0), (2,0), ..., (n-1,0)
 * - Column 1: elements (1,1), (2,1), ..., (n-1,1)
 * - Column j: elements (j,j), (j+1,j), ..., (n-1,j)
 *
 * Index formula for element (row, col) where row >= col:
 *   idx = col * n - col * (col + 1) / 2 + row
 *
 * This is equivalent to the LAPACK 'L' (lower) packed format for zhpmv, etc.
 *
 * Reference: thesis sec4.tex (Eq 6-10, Alg 5)
 */

#include "state.h"
#include "utils.h"

#include <assert.h>
#include <complex.h>
#include <stdio.h>

/*
 * =====================================================================================================================
 * Helper functions
 * =====================================================================================================================
 */

/**
 * @brief Compute index for element (row, col) in packed lower-triangular storage
 *
 * Assumes row >= col (lower triangle).
 *
 * @param dim   Matrix dimension (2^qubits)
 * @param row   Row index
 * @param col   Column index (col <= row)
 * @return      Index into packed storage array
 */
static inline dim_t packed_idx(dim_t dim, dim_t row, dim_t col) {
    return col * dim - col * (col + 1) / 2 + row;
}

/*
 * =====================================================================================================================
 * Packed state implementations
 * =====================================================================================================================
 */

dim_t state_packed_len(qubit_t qubits) {
    /*
     * Storage length for packed lower-triangular format: dim*(dim+1)/2
     *
     * To avoid overflow in dim*(dim+1), we compute (dim/2)*(dim+1) since
     * dim = 2^qubits is always even for qubits >= 1.
     *
     * Maximum supported qubits depends on sizeof(dim_t):
     * - 32-bit dim_t: max ~30 qubits (storage ~2^60 elements won't fit anyway)
     * - 64-bit dim_t: max ~62 qubits (memory limited long before overflow)
     */
    assert(qubits < sizeof(dim_t) * 8 - 1);  /* Prevent dim overflow */

    const dim_t dim = POW2(qubits, dim_t);
    if (qubits == 0) {
        return 1;  /* dim=1 is odd, so (dim>>1)*(dim+1) = 0*2 = 0, wrong */
    }
    /* dim is even for qubits >= 1, so (dim>>1)*(dim+1) = dim*(dim+1)/2 */
    return (dim >> 1) * (dim + 1);
}


void state_packed_plus(state_t *state, qubit_t qubits) {
    /* Initialize empty state (allocates zero-initialized storage) */
    state_init(state, qubits, NULL);
    if (!state->data) return;

    /*
     * The plus state |+>^n = (1/sqrt(2^n)) * sum_i |i> has density matrix:
     *   rho = |+><+| = (1/2^n) * sum_{i,j} |i><j|
     *
     * All matrix elements are equal to 1/2^n.
     */
    const dim_t dim = POW2(qubits, dim_t);
    const dim_t len = state_packed_len(qubits);
    const cplx_t prefactor = 1.0 / (double)dim;

    /* Set all stored elements to the uniform value */
    for (dim_t i = 0; i < len; ++i) {
        state->data[i] = prefactor;
    }
}


cplx_t state_packed_get(const state_t *state, dim_t row, dim_t col) {
    const dim_t dim = POW2(state->qubits, dim_t);
    assert(row < dim && col < dim);

    if (row >= col) {
        /* Lower triangle: direct access */
        return state->data[packed_idx(dim, row, col)];
    } else {
        /* Upper triangle: return conjugate of (col, row) */
        return conj(state->data[packed_idx(dim, col, row)]);
    }
}


void state_packed_set(state_t *state, dim_t row, dim_t col, cplx_t val) {
    const dim_t dim = POW2(state->qubits, dim_t);
    assert(row < dim && col < dim);

    if (row >= col) {
        /* Lower triangle: direct storage */
        state->data[packed_idx(dim, row, col)] = val;
    } else {
        /* Upper triangle: store conjugate at (col, row) */
        state->data[packed_idx(dim, col, row)] = conj(val);
    }
}


void state_packed_print(const state_t *state) {
    printf("Type: MIXED_PACKED\nQubits: %u\n", state->qubits);
    mprint_packed(state->data, POW2(state->qubits, dim_t));
}
