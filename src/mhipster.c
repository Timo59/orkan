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

    printf("\nInput state:\n");
    state_print(state);
    printf("Hilbert space dimension: %ld\n", dim);
    printf("Invariant subspace dimension: %ld\n", subdim);
    printf("Increment: %ld\n", incr);

    // Iterate the invariant column blocks; col_block is the first column in the block
    for (dim_t col_block = 0; col_block < dim; col_block += subdim) {
        // Number of elements before current block
        dim_t offset = col_block * (2 * dim - col_block + 1) / 2;

        // Distance between elements (i,j) and (i+incr,j+incr)
        dim_t stride = incr * (2 * (dim - col_block) - incr + 1) / 2;

        // Iterate the single columns with j_a=0
        for (dim_t col = col_block; col < col_block + incr; ++col) {
            printf("\nColumn: %ld\n", col);
            printf("Offset: %ld\n", offset);
            printf("Stride: %ld\n", stride);

            // Iterate the invariant row blocks; row_block is the first row in the block
            for (dim_t row_block = 0; row_block < dim - col_block; row_block += subdim) {
                cplx_t *data = state->data + row_block + offset;    // First element in column block and row block

                // Due to row >= col (i_a=j_a=0), only a portion of the first row block contributes
                dim_t n_rows;
                if (row_block == 0) n_rows = incr - (col - col_block);
                else n_rows = incr;

                // Swap entries in the row block with i_a=j_a=0 and i_a=j_a=1
                cblas_zswap(n_rows, data, 1, data + stride, 1);

                // Complex conjugation of entries with i_a=1 and j_a=0
                for (dim_t k = incr; k < subdim; ++k) {
                    *(data + k) = conj(*(data + k));
                }
            }

            // Add to the offset/Subtract from the stride the number of entries in current column
            offset += dim - col;
            stride -= incr;
        }
    }

    printf("\nOutput state:\n");
    state_print(state);
}