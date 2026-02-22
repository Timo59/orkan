// test_mixed_utils.c - Utility functions for mixed state test harnesses (packed and tiled)

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
    for (unsigned col = 0; col < n; ++col) {
        for (unsigned row = col; row < n; ++row) {
            full[row + col * n] = packed[packed_index(n, row, col)];  // Lower triangle
            if (row != col) {
                // Hermitian property: H[col,row] = conj(H[row,col])
                full[col + row * n] = conj(packed[packed_index(n, row, col)]);  // Upper triangle
            }
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
    for (unsigned col = 0; col < n; ++col) {
        for (unsigned row = col; row < n; ++row) {
            packed[packed_index(n, row, col)] = full[row + col * n];
        }
    }

    return packed;
}

/*
 * =====================================================================================================================
 * Tiled storage utilities
 * =====================================================================================================================
 */

void tiled_state_from_full(state_t *state, unsigned nqubits, const cplx_t *full) {
    const unsigned dim = POW2(nqubits, dim_t);

    state->type = MIXED_TILED;
    state_init(state, nqubits, NULL);  // Allocate zero-initialized tiled storage

    if (!state->data) return;

    // Populate via state_set -- only need lower triangle (Hermitian)
    for (unsigned col = 0; col < dim; ++col) {
        for (unsigned row = col; row < dim; ++row) {
            state_set(state, row, col, full[row + col * dim]);
        }
    }
}

void assert_tiled_equals_full(const state_t *state, const cplx_t *ref_full, unsigned dim) {
    for (unsigned col = 0; col < dim; ++col) {
        for (unsigned row = 0; row < dim; ++row) {
            cplx_t got = state_get(state, row, col);
            cplx_t exp = ref_full[row + col * dim];
            double re_diff = creal(got) - creal(exp);
            double im_diff = cimag(got) - cimag(exp);
            if (re_diff > PRECISION || re_diff < -PRECISION ||
                im_diff > PRECISION || im_diff < -PRECISION) {
                char msg[256];
                snprintf(msg, sizeof(msg),
                    "Tiled mismatch at (%u,%u): got (%.15e, %.15e) expected (%.15e, %.15e)",
                    row, col, creal(got), cimag(got), creal(exp), cimag(exp));
                TEST_FAIL_MESSAGE(msg);
            }
        }
    }
}
