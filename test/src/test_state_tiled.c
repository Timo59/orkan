// test_state_tiled.c - Unit tests for MIXED_TILED quantum states (tiled/blocked density matrices)

#include "test.h"
#include "state.h"
#include "utils.h"

#include <math.h>
#include <stdlib.h>
#include <complex.h>

/*
 * =====================================================================================================================
 * Test: state_tiled_len() - small systems (1-4 qubits, fit in single tile)
 * =====================================================================================================================
 */

void test_tiled_len_small(void) {
    // For small systems that fit in a single tile (TILE_DIM=32 by default)
    // 1 qubit: dim=2, n_tiles=1, len = 1*TILE_SIZE
    // 2 qubits: dim=4, n_tiles=1, len = 1*TILE_SIZE
    // 3 qubits: dim=8, n_tiles=1, len = 1*TILE_SIZE
    // 4 qubits: dim=16, n_tiles=1, len = 1*TILE_SIZE

    // With TILE_DIM=32, TILE_SIZE=1024
    // Single tile systems store full TILE_SIZE elements
    dim_t tile_size = TILE_SIZE;

    TEST_ASSERT_EQUAL_INT64(tile_size, state_tiled_len(1));
    TEST_ASSERT_EQUAL_INT64(tile_size, state_tiled_len(2));
    TEST_ASSERT_EQUAL_INT64(tile_size, state_tiled_len(3));
    TEST_ASSERT_EQUAL_INT64(tile_size, state_tiled_len(4));
}

/*
 * =====================================================================================================================
 * Test: state_tiled_len() - 8 qubits (crosses tile boundaries)
 * =====================================================================================================================
 */

void test_tiled_len_8_qubits(void) {
    // 8 qubits: dim = 256
    // With TILE_DIM=32: n_tiles = ceil(256/32) = 8
    // Lower triangle: 8*(8+1)/2 = 36 tiles
    // Total elements: 36 * TILE_SIZE

    dim_t dim = 256;
    dim_t n_tiles = (dim + TILE_DIM - 1) / TILE_DIM;  // ceil(256/32) = 8
    dim_t expected_tiles = n_tiles * (n_tiles + 1) / 2;  // 8*9/2 = 36
    dim_t expected_len = expected_tiles * TILE_SIZE;

    TEST_ASSERT_EQUAL_INT64(expected_len, state_tiled_len(8));
}

/*
 * =====================================================================================================================
 * Test: state_init() with NULL data - zero initialization
 * =====================================================================================================================
 */

void test_tiled_init_null_data(void) {
    for (qubit_t n = 1; n <= 4; n++) {
        state_t state = {.type = MIXED_TILED, .data = NULL, .qubits = 0};

        state_init(&state, n, NULL);

        // Verify initialization
        TEST_ASSERT_EQUAL_UINT8(n, state.qubits);
        TEST_ASSERT_NOT_NULL(state.data);

        // Verify all elements are zero
        dim_t len = state_tiled_len(n);
        for (dim_t i = 0; i < len; i++) {
            TEST_ASSERT_COMPLEX_WITHIN(0.0 + 0.0*I, state.data[i], PRECISION);
        }

        state_free(&state);
    }
}

/*
 * =====================================================================================================================
 * Test: state_init() with data - ownership transfer
 * =====================================================================================================================
 */

void test_tiled_init_ownership(void) {
    state_t state = {.type = MIXED_TILED, .data = NULL, .qubits = 0};

    // Allocate data externally
    dim_t len = state_tiled_len(2);
    cplx_t *data = malloc(len * sizeof(cplx_t));
    TEST_ASSERT_NOT_NULL(data);

    // Set some values
    for (dim_t i = 0; i < len; i++) {
        data[i] = (double)(i % 100) + (double)(i % 50) * 0.1 * I;
    }

    cplx_t *original_ptr = data;

    // Transfer ownership
    state_init(&state, 2, &data);

    // Verify ownership transferred
    TEST_ASSERT_NULL(data);  // Source pointer nullified
    TEST_ASSERT_EQUAL_PTR(original_ptr, state.data);  // State owns the data

    // Verify values preserved
    for (dim_t i = 0; i < len; i++) {
        cplx_t expected = (double)(i % 100) + (double)(i % 50) * 0.1 * I;
        TEST_ASSERT_COMPLEX_WITHIN(expected, state.data[i], PRECISION);
    }

    state_free(&state);
}

/*
 * =====================================================================================================================
 * Test: state_free() - safe to call twice (double-free safety)
 * =====================================================================================================================
 */

void test_tiled_free_double_free(void) {
    state_t state = {.type = MIXED_TILED, .data = NULL, .qubits = 0};

    state_init(&state, 2, NULL);
    TEST_ASSERT_NOT_NULL(state.data);
    TEST_ASSERT_EQUAL_UINT8(2, state.qubits);

    // First free
    state_free(&state);
    TEST_ASSERT_NULL(state.data);
    TEST_ASSERT_EQUAL_UINT8(0, state.qubits);

    // Second free - should be safe
    state_free(&state);
    TEST_ASSERT_NULL(state.data);
    TEST_ASSERT_EQUAL_UINT8(0, state.qubits);
}

/*
 * =====================================================================================================================
 * Test: state_tiled_plus() - ALL elements = 1/dim (|+⟩⟨+|^⊗n density matrix)
 * =====================================================================================================================
 */

void test_tiled_plus_values(void) {
    // Test for 1-4 qubits
    // The plus state |+>^n has density matrix rho = |+><+|^⊗n where ALL elements equal 1/dim
    // This matches the behavior of state_packed_plus()
    for (qubit_t n = 1; n <= 4; n++) {
        state_t state = {.type = MIXED_TILED, .data = NULL, .qubits = 0};

        state_tiled_plus(&state, n);

        TEST_ASSERT_NOT_NULL(state.data);
        TEST_ASSERT_EQUAL_UINT8(n, state.qubits);

        dim_t dim = POW2(n, dim_t);
        double expected_value = 1.0 / (double)dim;

        // Check ALL elements = 1/dim (not just diagonal!)
        // The |+⟩⟨+| density matrix is a completely uniform matrix
        for (dim_t row = 0; row < dim; row++) {
            for (dim_t col = 0; col < dim; col++) {
                cplx_t val = state_tiled_get(&state, row, col);
                TEST_ASSERT_DOUBLE_WITHIN(PRECISION, expected_value, creal(val));
                TEST_ASSERT_DOUBLE_WITHIN(PRECISION, 0.0, cimag(val));
            }
        }

        state_free(&state);
    }
}

/*
 * =====================================================================================================================
 * Test: state_tiled_get() - within single tile (small system)
 * =====================================================================================================================
 */

void test_tiled_get_within_tile(void) {
    state_t state = {.type = MIXED_TILED, .data = NULL, .qubits = 0};

    state_init(&state, 2, NULL);  // 4x4 matrix, fits in single tile

    // Set some known values in lower triangle using tiled_set
    state_tiled_set(&state, 0, 0, 1.0 + 0.0*I);
    state_tiled_set(&state, 1, 0, 2.0 + 1.0*I);
    state_tiled_set(&state, 1, 1, 3.0 + 0.0*I);
    state_tiled_set(&state, 2, 0, 4.0 + 2.0*I);
    state_tiled_set(&state, 2, 1, 5.0 + 3.0*I);
    state_tiled_set(&state, 2, 2, 6.0 + 0.0*I);
    state_tiled_set(&state, 3, 0, 7.0 + 4.0*I);
    state_tiled_set(&state, 3, 1, 8.0 + 5.0*I);
    state_tiled_set(&state, 3, 2, 9.0 + 6.0*I);
    state_tiled_set(&state, 3, 3, 10.0 + 0.0*I);

    // Verify get retrieves correct values
    TEST_ASSERT_COMPLEX_WITHIN(1.0 + 0.0*I, state_tiled_get(&state, 0, 0), PRECISION);
    TEST_ASSERT_COMPLEX_WITHIN(2.0 + 1.0*I, state_tiled_get(&state, 1, 0), PRECISION);
    TEST_ASSERT_COMPLEX_WITHIN(3.0 + 0.0*I, state_tiled_get(&state, 1, 1), PRECISION);
    TEST_ASSERT_COMPLEX_WITHIN(4.0 + 2.0*I, state_tiled_get(&state, 2, 0), PRECISION);
    TEST_ASSERT_COMPLEX_WITHIN(5.0 + 3.0*I, state_tiled_get(&state, 2, 1), PRECISION);
    TEST_ASSERT_COMPLEX_WITHIN(6.0 + 0.0*I, state_tiled_get(&state, 2, 2), PRECISION);
    TEST_ASSERT_COMPLEX_WITHIN(7.0 + 4.0*I, state_tiled_get(&state, 3, 0), PRECISION);
    TEST_ASSERT_COMPLEX_WITHIN(8.0 + 5.0*I, state_tiled_get(&state, 3, 1), PRECISION);
    TEST_ASSERT_COMPLEX_WITHIN(9.0 + 6.0*I, state_tiled_get(&state, 3, 2), PRECISION);
    TEST_ASSERT_COMPLEX_WITHIN(10.0 + 0.0*I, state_tiled_get(&state, 3, 3), PRECISION);

    state_free(&state);
}

/*
 * =====================================================================================================================
 * Test: state_tiled_get() - across tile boundaries (8 qubits)
 * =====================================================================================================================
 */

void test_tiled_get_cross_tile(void) {
    state_t state = {.type = MIXED_TILED, .data = NULL, .qubits = 0};

    state_init(&state, 8, NULL);  // 256x256 matrix, crosses multiple 32x32 tiles

    // Set some elements in different tiles
    // Diagonal tiles
    state_tiled_set(&state, 0, 0, 1.0 + 0.0*I);      // Tile (0,0), element (0,0)
    state_tiled_set(&state, 31, 31, 2.0 + 0.0*I);    // Tile (0,0), element (31,31)
    state_tiled_set(&state, 32, 32, 3.0 + 0.0*I);    // Tile (1,1), element (0,0)
    state_tiled_set(&state, 63, 63, 4.0 + 0.0*I);    // Tile (1,1), element (31,31)
    state_tiled_set(&state, 255, 255, 5.0 + 0.0*I);  // Tile (7,7), element (31,31)

    // Off-diagonal tiles
    state_tiled_set(&state, 32, 0, 6.0 + 1.0*I);     // Tile (1,0), element (0,0)
    state_tiled_set(&state, 63, 31, 7.0 + 2.0*I);    // Tile (1,0), element (31,31)
    state_tiled_set(&state, 64, 0, 8.0 + 3.0*I);     // Tile (2,0), element (0,0)
    state_tiled_set(&state, 128, 64, 9.0 + 4.0*I);   // Tile (4,2), element (0,0)

    // Verify retrieval
    TEST_ASSERT_COMPLEX_WITHIN(1.0 + 0.0*I, state_tiled_get(&state, 0, 0), PRECISION);
    TEST_ASSERT_COMPLEX_WITHIN(2.0 + 0.0*I, state_tiled_get(&state, 31, 31), PRECISION);
    TEST_ASSERT_COMPLEX_WITHIN(3.0 + 0.0*I, state_tiled_get(&state, 32, 32), PRECISION);
    TEST_ASSERT_COMPLEX_WITHIN(4.0 + 0.0*I, state_tiled_get(&state, 63, 63), PRECISION);
    TEST_ASSERT_COMPLEX_WITHIN(5.0 + 0.0*I, state_tiled_get(&state, 255, 255), PRECISION);
    TEST_ASSERT_COMPLEX_WITHIN(6.0 + 1.0*I, state_tiled_get(&state, 32, 0), PRECISION);
    TEST_ASSERT_COMPLEX_WITHIN(7.0 + 2.0*I, state_tiled_get(&state, 63, 31), PRECISION);
    TEST_ASSERT_COMPLEX_WITHIN(8.0 + 3.0*I, state_tiled_get(&state, 64, 0), PRECISION);
    TEST_ASSERT_COMPLEX_WITHIN(9.0 + 4.0*I, state_tiled_get(&state, 128, 64), PRECISION);

    state_free(&state);
}

/*
 * =====================================================================================================================
 * Test: state_tiled_get() - upper triangle returns conjugate
 * =====================================================================================================================
 */

void test_tiled_get_upper(void) {
    state_t state = {.type = MIXED_TILED, .data = NULL, .qubits = 0};

    state_init(&state, 2, NULL);  // 4x4 matrix

    // Set some lower triangle values with non-zero imaginary parts
    state_tiled_set(&state, 1, 0, 2.0 + 1.0*I);
    state_tiled_set(&state, 2, 0, 4.0 + 2.0*I);
    state_tiled_set(&state, 2, 1, 5.0 + 3.0*I);
    state_tiled_set(&state, 3, 0, 7.0 + 4.0*I);

    // Test upper triangle returns conjugate (row < col)
    TEST_ASSERT_COMPLEX_WITHIN(2.0 - 1.0*I, state_tiled_get(&state, 0, 1), PRECISION);
    TEST_ASSERT_COMPLEX_WITHIN(4.0 - 2.0*I, state_tiled_get(&state, 0, 2), PRECISION);
    TEST_ASSERT_COMPLEX_WITHIN(5.0 - 3.0*I, state_tiled_get(&state, 1, 2), PRECISION);
    TEST_ASSERT_COMPLEX_WITHIN(7.0 - 4.0*I, state_tiled_get(&state, 0, 3), PRECISION);

    state_free(&state);
}

/*
 * =====================================================================================================================
 * Test: state_tiled_set() and state_tiled_get() - roundtrip
 * =====================================================================================================================
 */

void test_tiled_set_roundtrip(void) {
    state_t state = {.type = MIXED_TILED, .data = NULL, .qubits = 0};

    state_init(&state, 3, NULL);  // 8x8 matrix

    // Test various complex values at different positions
    cplx_t test_values[] = {
        1.0 + 0.0*I,
        0.0 + 1.0*I,
        -0.5 + 0.5*I,
        0.7 - 0.3*I,
        1.234 + 5.678*I,
        -9.876 - 5.432*I,
        0.0 + 0.0*I,
        100.0 + 100.0*I
    };

    // Set values at various positions (lower triangle only)
    dim_t positions[][2] = {
        {0, 0}, {1, 0}, {1, 1}, {2, 1}, {3, 0}, {4, 3}, {5, 2}, {7, 7}
    };

    for (int i = 0; i < 8; i++) {
        dim_t row = positions[i][0];
        dim_t col = positions[i][1];
        state_tiled_set(&state, row, col, test_values[i]);
    }

    // Verify roundtrip
    for (int i = 0; i < 8; i++) {
        dim_t row = positions[i][0];
        dim_t col = positions[i][1];
        cplx_t retrieved = state_tiled_get(&state, row, col);
        TEST_ASSERT_COMPLEX_WITHIN(test_values[i], retrieved, PRECISION);
    }

    state_free(&state);
}

/*
 * =====================================================================================================================
 * Test: state_cp() - independence of copy
 * =====================================================================================================================
 */

void test_tiled_cp_independence(void) {
    state_t original = {.type = MIXED_TILED, .data = NULL, .qubits = 0};

    state_tiled_plus(&original, 2);  // Maximally mixed state

    // Create copy
    state_t copy = state_cp(&original);

    // Verify metadata copied
    TEST_ASSERT_EQUAL_UINT8(original.type, copy.type);
    TEST_ASSERT_EQUAL_UINT8(original.qubits, copy.qubits);

    // Verify different memory allocations
    TEST_ASSERT_NOT_NULL(copy.data);
    TEST_ASSERT_NOT_EQUAL_PTR(original.data, copy.data);

    // Verify values initially identical
    dim_t len = state_tiled_len(2);
    for (dim_t i = 0; i < len; i++) {
        TEST_ASSERT_COMPLEX_WITHIN(original.data[i], copy.data[i], PRECISION);
    }

    // Modify original by setting a value
    state_tiled_set(&original, 0, 0, 99.0 + 99.0*I);
    state_tiled_set(&original, 1, 0, -50.0 - 50.0*I);

    // Verify original modified
    TEST_ASSERT_COMPLEX_WITHIN(99.0 + 99.0*I, state_tiled_get(&original, 0, 0), PRECISION);
    TEST_ASSERT_COMPLEX_WITHIN(-50.0 - 50.0*I, state_tiled_get(&original, 1, 0), PRECISION);

    // Verify copy unchanged (should still be |+><+|: ALL elements = 0.25)
    double expected_val = 0.25;
    TEST_ASSERT_COMPLEX_WITHIN(expected_val + 0.0*I, state_tiled_get(&copy, 0, 0), PRECISION);
    TEST_ASSERT_COMPLEX_WITHIN(expected_val + 0.0*I, state_tiled_get(&copy, 1, 0), PRECISION);

    state_free(&original);
    state_free(&copy);
}

/*
 * =====================================================================================================================
 * Test: Hermitian symmetry - ρ[i,j] = conj(ρ[j,i])
 * =====================================================================================================================
 */

void test_tiled_hermitian_symmetry(void) {
    state_t state = {.type = MIXED_TILED, .data = NULL, .qubits = 0};

    state_init(&state, 2, NULL);  // 4x4 matrix

    // Set some complex values in lower triangle
    state_tiled_set(&state, 1, 0, 2.0 + 1.5*I);
    state_tiled_set(&state, 2, 0, 3.0 + 2.5*I);
    state_tiled_set(&state, 2, 1, 4.0 + 3.5*I);
    state_tiled_set(&state, 3, 1, 5.0 + 4.5*I);

    // Verify Hermitian symmetry: ρ[i,j] = conj(ρ[j,i])
    cplx_t val_10 = state_tiled_get(&state, 1, 0);
    cplx_t val_01 = state_tiled_get(&state, 0, 1);
    TEST_ASSERT_COMPLEX_WITHIN(conj(val_10), val_01, PRECISION);

    cplx_t val_20 = state_tiled_get(&state, 2, 0);
    cplx_t val_02 = state_tiled_get(&state, 0, 2);
    TEST_ASSERT_COMPLEX_WITHIN(conj(val_20), val_02, PRECISION);

    cplx_t val_21 = state_tiled_get(&state, 2, 1);
    cplx_t val_12 = state_tiled_get(&state, 1, 2);
    TEST_ASSERT_COMPLEX_WITHIN(conj(val_21), val_12, PRECISION);

    cplx_t val_31 = state_tiled_get(&state, 3, 1);
    cplx_t val_13 = state_tiled_get(&state, 1, 3);
    TEST_ASSERT_COMPLEX_WITHIN(conj(val_31), val_13, PRECISION);

    state_free(&state);
}

/*
 * =====================================================================================================================
 * Test: Edge case - single qubit
 * =====================================================================================================================
 */

void test_tiled_single_qubit(void) {
    state_t state = {.type = MIXED_TILED, .data = NULL, .qubits = 0};

    state_tiled_plus(&state, 1);  // 1 qubit, 2x2 matrix

    TEST_ASSERT_EQUAL_UINT8(1, state.qubits);
    TEST_ASSERT_NOT_NULL(state.data);

    // |+><+| density matrix for single qubit: ALL elements = 1/2
    double expected = 0.5;
    TEST_ASSERT_COMPLEX_WITHIN(expected + 0.0*I, state_tiled_get(&state, 0, 0), PRECISION);
    TEST_ASSERT_COMPLEX_WITHIN(expected + 0.0*I, state_tiled_get(&state, 1, 1), PRECISION);
    TEST_ASSERT_COMPLEX_WITHIN(expected + 0.0*I, state_tiled_get(&state, 1, 0), PRECISION);
    TEST_ASSERT_COMPLEX_WITHIN(expected + 0.0*I, state_tiled_get(&state, 0, 1), PRECISION);

    state_free(&state);
}

/*
 * =====================================================================================================================
 * Test: Tile boundary conditions (8 qubits) - verify tile transitions
 * =====================================================================================================================
 */

void test_tiled_tile_boundaries(void) {
    state_t state = {.type = MIXED_TILED, .data = NULL, .qubits = 0};

    state_init(&state, 8, NULL);  // 256x256 matrix

    // Test elements at tile boundaries (with TILE_DIM=32)
    // Boundary positions: 31, 32, 63, 64, etc.

    // Set values at boundaries
    state_tiled_set(&state, 31, 0, 1.0 + 0.1*I);   // Last row of tile (0,0)
    state_tiled_set(&state, 32, 0, 2.0 + 0.2*I);   // First row of tile (1,0)
    state_tiled_set(&state, 31, 31, 3.0 + 0.3*I);  // Corner of tile (0,0)
    state_tiled_set(&state, 32, 31, 4.0 + 0.4*I);  // Boundary between tiles
    state_tiled_set(&state, 32, 32, 5.0 + 0.5*I);  // Corner of tile (1,1)
    state_tiled_set(&state, 63, 32, 6.0 + 0.6*I);  // Last row of tile (1,1)
    state_tiled_set(&state, 64, 32, 7.0 + 0.7*I);  // First row of tile (2,1)

    // Verify values are correctly stored and retrieved
    TEST_ASSERT_COMPLEX_WITHIN(1.0 + 0.1*I, state_tiled_get(&state, 31, 0), PRECISION);
    TEST_ASSERT_COMPLEX_WITHIN(2.0 + 0.2*I, state_tiled_get(&state, 32, 0), PRECISION);
    TEST_ASSERT_COMPLEX_WITHIN(3.0 + 0.3*I, state_tiled_get(&state, 31, 31), PRECISION);
    TEST_ASSERT_COMPLEX_WITHIN(4.0 + 0.4*I, state_tiled_get(&state, 32, 31), PRECISION);
    TEST_ASSERT_COMPLEX_WITHIN(5.0 + 0.5*I, state_tiled_get(&state, 32, 32), PRECISION);
    TEST_ASSERT_COMPLEX_WITHIN(6.0 + 0.6*I, state_tiled_get(&state, 63, 32), PRECISION);
    TEST_ASSERT_COMPLEX_WITHIN(7.0 + 0.7*I, state_tiled_get(&state, 64, 32), PRECISION);

    state_free(&state);
}
