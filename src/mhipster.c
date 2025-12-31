// mhipster.c - Functions representing fundamental quantum gates applied to mixed states, i.e., density matrices on
// Hilbert spaces

#ifndef GATE_H
#include "gate.h"
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

void x_mixed(state_t *state, const qubit_t target) {
    // Check that target is in the system's scope
    if (target >= state->qubits) {
        fprintf(stderr, "x(): Target is out of scope; Expected %u, Was %u", state->qubits, target);
        state_free(state);
        return;
    }

    const dim_t dim = POW2(state->qubits, dim_t);   // Hilbert space dimension; Order of density matrix
    const dim_t incr = POW2(target, dim_t); // Index distance between elements only differing in the targeted qubit
    const dim_t subdim = POW2(target + 1, dim_t);   // Invariant subspace dimension

    // Iterate the invariant column blocks; col_block is the first column in the block
    for (dim_t col_block = 0; col_block < dim; col_block += subdim) {
        // Number of elements before current block
        dim_t offset = col_block * (2 * dim - col_block + 1) / 2;

        // Distance between elements (i,j) and (i+incr,j+incr)
        dim_t stride = incr * (2 * (dim - col_block) - incr + 1) / 2;

        // Iterate the single columns with j_a=0
        for (dim_t col = col_block; col < col_block + incr; ++col) {

            // First row block of each column needs special treatment
            cplx_t *data = state->data + offset;    // First element in the current column

            // Due to row>=col, not all entries with row_a=col_a=0 are part of the array
            dim_t n_rows = incr - (col - col_block);    // Number of rows in col with row<col_block+2**a

            // Swap entries with row_a=col_a=0 and row_a=col_a=1 and row<col+2**a, starting from row=col
            cblas_zswap(n_rows, data, 1, data + stride, 1);

            // For the first block, i<j+2**a, therefore take the conjugate of the transpose value
            dim_t step = 0; // Distance between entries (col+2**a,row) and (row+2**a,col)
            data += n_rows; // Data points to the first element in col with row_a=1
            // Iterate rows with row-2**a >= col
            for (dim_t k = col - col_block; k < incr; ++k) {
                cplx_t tmp = data[k];
                data[k] = conj(data[k + step]);
                data[k + step] = conj(tmp);
                step += dim - (k + col_block + 2);
            }
            data += incr;

            // Iterate the invariant row blocks; row_block is the first row in the block
            dim_t row_block = col_block + subdim;
            while (row_block < dim) {
                // Swap entries in the row block with i_a=j_a=0 and i_a=j_a=1
                cblas_zswap(incr, data, 1, data + stride, 1);

                // Complex conjugation of entries with i_a=1 and j_a=0
                data += incr;
                cblas_zswap(incr, data, 1, data + stride - subdim, 1);

                data += incr;
                row_block += subdim;
            }

            // Add to the offset/Subtract from the stride the number of entries in current column
            offset += dim - col;
            stride -= incr;
        }
    }
}