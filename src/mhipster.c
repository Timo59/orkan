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
    const dim_t dim = POW2(state->qubits, dim_t);   // Hilbert space dimension; Order of density matrix
    const dim_t incr = POW2(target, dim_t); // Distance between indices only differing in the targeted qubit
    const dim_t subdim = POW2(target + 1, dim_t);   // Invariant subspace dimension

    // Iterate invariant column blocks, i.e., start with all bits set to zero and increment only bits left to the
    // targeted bit in the column's index binary representation. col_block is the first column in the block.
    for (dim_t col_block = 0; col_block < dim; col_block += subdim) {
        // Number of elements in data array prior to current block; sum over numbers of entries in the first `col_block`
        // columns in the lower triangle of a square matrix of order `dim` in column major format
        dim_t offset = col_block * (2 * dim - col_block + 1) / 2;

        // Distance between the element i in column `col_block` and the element i+`incr` in column `col_block`+`incr` in
        // contiguous array representing the lower triangle of a square matrix of order `dim`in column major format. The
        // sum over numbers of entries in columns [`col_block`,`col_block`+`incr`] minus the number of entries up to i
        // in column `col_block` plus the number of entries up to i+`incr` in column `col_block`+`incr`.
        dim_t stride = incr * (2 * (dim - col_block) - incr + 1) / 2;

        // Iterate the single columns whose index binary representation is zero at the targted bit, i.e., increment only
        // bits right to the targeted bit. Updated entries are superpositions of entries in columns whose index only
        // differs in the targeted bit.
        for (dim_t col = col_block; col < col_block + incr; ++col) {
            cplx_t *data = state->data + offset;    // First element in the current column

            // First row block, i.e., row in [`col`,`col_block`+`incr`], of each column needs special treatment:
            // Only those entries whose row index binary representation is zero in the targeted bit and >=`col` is in
            // the range of the data array
            dim_t n_rows = incr - (col - col_block);    // Number of rows in [`col`,`col_block`+`incr`]

            // Swap entries whose row AND `col` index binary representation is zero in the targeted bit and entries
            // whose row and `col` representation is one in the targeted bit.
            cblas_zswap(n_rows, data, 1, data + stride, 1);

            // Entries whose row and `col` index binary representation is one and zero at the targeted position,
            // respectively, are swapped with entries whose row and `col` index binary representation is zero and one at
            // the targeted position, respectively. For row<`col`+`incr`, the latter is not in the range of the data
            // array: take the conjugate of the entry (`col`+`incr`,row) instead.
            data += n_rows; // First element in `col` whose index binary representation is one at the targeted bit
            dim_t step = 0; // First element updated to complex conjugate of itself (row=`col`)

            // Iterate rows whose index binary representation is one at the targeted bit and coincides with `col`'s
            // binary representation elsewhere (rows up to `col`+`incr` were updated with previous columns).
            for (dim_t k = col - col_block; k < incr; ++k) {
                // Swap entries (row,`col`+`incr`) and (`col`+`incr`,row) and complex conjugate them
                cplx_t tmp = data[k];
                data[k] = conj(data[k + step]);
                data[k + step] = conj(tmp);

                // The distance between elements (j+d,i) and (i+d,j) in the lower triangle of an N-by-N matrix stored in
                // column major format is computed via simple index calculation to be: (i-j)*(2N-i-j-3)/2, where i>=j.
                // The update is difference of distances ((j+d,i+1),(i+1+d,j)) and ((j+d,i),(i+d,j))
                step += dim - (k + col_block + 2);
            }

            // First element of the second invariant row block, i.e., row `col_block`+`subdim` in current column
            data += incr;

            // Iterate the invariant row blocks, i.e., start with `col_block`+`subdim` and increment only bits left to
            // the targeted bit. `row_block` is the first row in the block
            dim_t row_block = col_block + subdim;
            while (row_block < dim) {
                // Swap entries whose row AND `col` index binary representation is zero in the targeted bit and entries
                // whose row and `col` index binary representation is one in the targeted bit.
                cblas_zswap(incr, data, 1, data + stride, 1);

                // Swap entries whose row AND `col` index binary representation is one and zero, respectively, in the
                // targeted bit and entries whose row and `col` index binary representation is zero and one,
                // respectively, in the targeted bit
                data += incr;   // First element in `row_block` whose index is one at the targeted bit

                // The distance between elements (i+d,j) and (i,j+d) is the sum of numbers of entries in columns [j,j+d]
                // minus entries up to the (i+d)-th in column j and plus entries up to the i-th entry in column j+d in
                // the lower triangle of the N-by-N square matrix.
                cblas_zswap(incr, data, 1, data + stride - subdim, 1);

                data += incr;   // First element in subsequent row_block
                row_block += subdim;    // Index of the first element in subsequent row_block
            }

            // Add the number of entries in the current column to the offset for the next column
            offset += dim - col;

            // Subtract the difference between the number of entries in the current column and the number of entries in
            // the (current + `incr`) column to update the stride
            stride -= incr;
        }
    }
}