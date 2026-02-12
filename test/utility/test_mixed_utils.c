// test_mixed_utils.c - Utility functions for converting between packed and full density matrix formats

/*
 * =====================================================================================================================
 * includes
 * =====================================================================================================================
 */

#include "test_mixed_utils.h"

#include <stdlib.h>
#include <stdio.h>

/*
 * =====================================================================================================================
 * Packed storage utilities
 * =====================================================================================================================
 */

cplx_t* density_unpack(unsigned n, const cplx_t *packed) {
    cplx_t *full = NULL;

    // Allocate memory for full n×n matrix
    if (!((full = malloc(n * n * sizeof(*full))))) {
        fprintf(stderr, "density_unpack(): full allocation failed\n");
        return full;
    }

    // Unpack lower triangle from packed storage (column-major)
    // Packed format: column 0 (n elements), column 1 (n-1 elements), ..., column n-1 (1 element)
    unsigned packed_idx = 0;
    for (unsigned col = 0; col < n; ++col) {
        for (unsigned row = col; row < n; ++row) {
            full[row + col * n] = packed[packed_idx];  // Lower triangle
            if (row != col) {
                // Hermitian property: H[col,row] = conj(H[row,col])
                full[col + row * n] = conj(packed[packed_idx]);  // Upper triangle
            }
            ++packed_idx;
        }
    }

    return full;
}


cplx_t* density_pack(unsigned n, const cplx_t *full) {
    cplx_t *packed = NULL;

    // Allocate memory for packed storage
    const unsigned packed_len = n * (n + 1) / 2;
    if (!((packed = malloc(packed_len * sizeof(*packed))))) {
        fprintf(stderr, "density_pack(): packed allocation failed\n");
        return packed;
    }

    // Pack lower triangle into packed storage (column-major)
    unsigned packed_idx = 0;
    for (unsigned col = 0; col < n; ++col) {
        for (unsigned row = col; row < n; ++row) {
            packed[packed_idx] = full[row + col * n];
            ++packed_idx;
        }
    }

    return packed;
}
