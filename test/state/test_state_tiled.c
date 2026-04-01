// test_state_tiled.c - Unit tests for MIXED_TILED quantum states (tiled/blocked density matrices)

#include "test.h"
#include "state.h"
#include "utils.h"

#include <math.h>
#include <stdlib.h>
#include <complex.h>

/*
 * =====================================================================================================================
 * Test: state_tiled_len() - sub-tile and exact-tile systems (1 to LOG_TILE_DIM qubits)
 * =====================================================================================================================
 */

void test_tiled_len_small(void) {
    // Sub-tile (dim < TILE_DIM) and exact-tile (dim == TILE_DIM) cases
    for (qubit_t n = 1; n <= LOG_TILE_DIM; n++) {
        idx_t dim = POW2(n, idx_t);
        idx_t n_tiles = (dim + TILE_DIM - 1) / TILE_DIM;
        idx_t expected = n_tiles * (n_tiles + 1) / 2 * TILE_SIZE;
        TEST_ASSERT_EQUAL_INT64(expected, state_tiled_len(n));
    }
}

/*
 * =====================================================================================================================
 * Test: state_tiled_len() - multi-tile system (LOG_TILE_DIM+1 qubits)
 * =====================================================================================================================
 */

void test_tiled_len_multi_tile(void) {
    // LOG_TILE_DIM+1 qubits: dim = 2*TILE_DIM, n_tiles = 2, 3 lower-triangular tiles
    const qubit_t n = LOG_TILE_DIM + 1;
    idx_t dim = POW2(n, idx_t);
    idx_t n_tiles = (dim + TILE_DIM - 1) / TILE_DIM;
    idx_t expected_len = n_tiles * (n_tiles + 1) / 2 * TILE_SIZE;

    TEST_ASSERT_EQUAL_INT64(expected_len, state_tiled_len(n));
}

/*
 * =====================================================================================================================
 * Test: state_init() with NULL data - zero initialization
 * =====================================================================================================================
 */

void test_tiled_init_null_data(void) {
    for (qubit_t n = 1; n <= LOG_TILE_DIM + 1; n++) {
        state_t state = {.type = MIXED_TILED, .data = NULL, .qubits = 0};

        state_init(&state, n, NULL);

        // Verify initialization
        TEST_ASSERT_EQUAL_UINT8(n, state.qubits);
        TEST_ASSERT_NOT_NULL(state.data);

        // Verify all elements are zero
        idx_t len = state_tiled_len(n);
        for (idx_t i = 0; i < len; i++) {
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

    // Allocate data externally (64-byte aligned)
    idx_t len = state_tiled_len(2);
    size_t size = (len * sizeof(cplx_t) + 63) & ~(size_t)63;
    cplx_t *data = aligned_alloc(64, size);
    TEST_ASSERT_NOT_NULL(data);

    // Set some values
    for (idx_t i = 0; i < len; i++) {
        data[i] = (double)(i % 100) + (double)(i % 50) * 0.1 * I;
    }

    cplx_t *original_ptr = data;

    // Transfer ownership
    state_init(&state, 2, &data);

    // Verify ownership transferred
    TEST_ASSERT_NULL(data);  // Source pointer nullified
    TEST_ASSERT_EQUAL_PTR(original_ptr, state.data);  // State owns the data

    // Verify values preserved
    for (idx_t i = 0; i < len; i++) {
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
    // The plus state |+>^n has density matrix rho = |+><+|^⊗n where ALL elements equal 1/dim
    // Covers sub-tile, exact-tile, and multi-tile cases
    for (qubit_t n = 1; n <= LOG_TILE_DIM + 1; n++) {
        state_t state = {.type = MIXED_TILED, .data = NULL, .qubits = 0};

        state_tiled_plus(&state, n);

        TEST_ASSERT_NOT_NULL(state.data);
        TEST_ASSERT_EQUAL_UINT8(n, state.qubits);

        idx_t dim = POW2(n, idx_t);
        double expected_value = 1.0 / (double)dim;

        // Check ALL elements = 1/dim (not just diagonal!)
        // The |+⟩⟨+| density matrix is a completely uniform matrix
        for (idx_t row = 0; row < dim; row++) {
            for (idx_t col = 0; col < dim; col++) {
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
 * Test: state_tiled_get() - across tile boundaries (LOG_TILE_DIM+1 qubits)
 * =====================================================================================================================
 */

void test_tiled_get_cross_tile(void) {
    const qubit_t n = LOG_TILE_DIM + 1;
    const idx_t T = TILE_DIM;
    state_t state = {.type = MIXED_TILED, .data = NULL, .qubits = 0};

    state_init(&state, n, NULL);  // dim = 2*T, 3 lower-triangular tiles

    // Diagonal tile T(0,0)
    state_tiled_set(&state, 0, 0, 1.0 + 0.0*I);
    state_tiled_set(&state, T-1, T-1, 2.0 + 0.0*I);

    // Diagonal tile T(1,1)
    state_tiled_set(&state, T, T, 3.0 + 0.0*I);
    state_tiled_set(&state, 2*T-1, 2*T-1, 4.0 + 0.0*I);

    // Off-diagonal tile T(1,0)
    state_tiled_set(&state, T, 0, 5.0 + 1.0*I);
    state_tiled_set(&state, 2*T-1, T-1, 6.0 + 2.0*I);
    state_tiled_set(&state, T+1, 1, 7.0 + 3.0*I);

    // Verify retrieval
    TEST_ASSERT_COMPLEX_WITHIN(1.0 + 0.0*I, state_tiled_get(&state, 0, 0), PRECISION);
    TEST_ASSERT_COMPLEX_WITHIN(2.0 + 0.0*I, state_tiled_get(&state, T-1, T-1), PRECISION);
    TEST_ASSERT_COMPLEX_WITHIN(3.0 + 0.0*I, state_tiled_get(&state, T, T), PRECISION);
    TEST_ASSERT_COMPLEX_WITHIN(4.0 + 0.0*I, state_tiled_get(&state, 2*T-1, 2*T-1), PRECISION);
    TEST_ASSERT_COMPLEX_WITHIN(5.0 + 1.0*I, state_tiled_get(&state, T, 0), PRECISION);
    TEST_ASSERT_COMPLEX_WITHIN(6.0 + 2.0*I, state_tiled_get(&state, 2*T-1, T-1), PRECISION);
    TEST_ASSERT_COMPLEX_WITHIN(7.0 + 3.0*I, state_tiled_get(&state, T+1, 1), PRECISION);

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
    idx_t positions[][2] = {
        {0, 0}, {1, 0}, {1, 1}, {2, 1}, {3, 0}, {4, 3}, {5, 2}, {7, 7}
    };

    for (int i = 0; i < 8; i++) {
        idx_t row = positions[i][0];
        idx_t col = positions[i][1];
        state_tiled_set(&state, row, col, test_values[i]);
    }

    // Verify roundtrip
    for (int i = 0; i < 8; i++) {
        idx_t row = positions[i][0];
        idx_t col = positions[i][1];
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

    state_tiled_plus(&original, 2);  // |+><+| state

    // Create copy
    state_t copy = state_cp(&original);

    // Verify metadata copied
    TEST_ASSERT_EQUAL_UINT8(original.type, copy.type);
    TEST_ASSERT_EQUAL_UINT8(original.qubits, copy.qubits);

    // Verify different memory allocations
    TEST_ASSERT_NOT_NULL(copy.data);
    TEST_ASSERT_NOT_EQUAL_PTR(original.data, copy.data);

    // Verify values initially identical
    idx_t len = state_tiled_len(2);
    for (idx_t i = 0; i < len; i++) {
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
 * Test: Tile boundary conditions (LOG_TILE_DIM+1 qubits) - verify tile transitions
 * =====================================================================================================================
 */

void test_tiled_tile_boundaries(void) {
    const qubit_t n = LOG_TILE_DIM + 1;
    const idx_t T = TILE_DIM;
    state_t state = {.type = MIXED_TILED, .data = NULL, .qubits = 0};

    state_init(&state, n, NULL);  // dim = 2*T

    // Set values at tile boundaries: T-1 is last index of tile 0, T is first of tile 1
    state_tiled_set(&state, T-1, 0, 1.0 + 0.1*I);       // Last row of T(0,0)
    state_tiled_set(&state, T, 0, 2.0 + 0.2*I);          // First row of T(1,0)
    state_tiled_set(&state, T-1, T-1, 3.0 + 0.3*I);      // Corner of T(0,0)
    state_tiled_set(&state, T, T-1, 4.0 + 0.4*I);        // Boundary: T(1,0) adjacent to T(0,0)
    state_tiled_set(&state, T, T, 5.0 + 0.5*I);          // Corner of T(1,1)
    state_tiled_set(&state, 2*T-1, T, 6.0 + 0.6*I);      // Last row of T(1,1)

    // Verify values are correctly stored and retrieved
    TEST_ASSERT_COMPLEX_WITHIN(1.0 + 0.1*I, state_tiled_get(&state, T-1, 0), PRECISION);
    TEST_ASSERT_COMPLEX_WITHIN(2.0 + 0.2*I, state_tiled_get(&state, T, 0), PRECISION);
    TEST_ASSERT_COMPLEX_WITHIN(3.0 + 0.3*I, state_tiled_get(&state, T-1, T-1), PRECISION);
    TEST_ASSERT_COMPLEX_WITHIN(4.0 + 0.4*I, state_tiled_get(&state, T, T-1), PRECISION);
    TEST_ASSERT_COMPLEX_WITHIN(5.0 + 0.5*I, state_tiled_get(&state, T, T), PRECISION);
    TEST_ASSERT_COMPLEX_WITHIN(6.0 + 0.6*I, state_tiled_get(&state, 2*T-1, T), PRECISION);

    state_free(&state);
}

/*
 * =====================================================================================================================
 * Test: state_tiled_set() - upper triangle stores conjugate at transposed position
 * =====================================================================================================================
 */

void test_tiled_set_upper(void) {
    state_t state = {.type = MIXED_TILED, .data = NULL, .qubits = 0};

    state_init(&state, 2, NULL);  // 4x4 matrix

    // Set upper triangle elements (row < col)
    // These should be stored as conjugates at the transposed position
    state_tiled_set(&state, 0, 1, 2.0 + 1.0*I);
    state_tiled_set(&state, 0, 2, 4.0 + 2.0*I);
    state_tiled_set(&state, 1, 2, 5.0 + 3.0*I);
    state_tiled_set(&state, 0, 3, 7.0 + 4.0*I);

    // Verify that state_tiled_get on the lower triangle position returns the conjugate
    TEST_ASSERT_COMPLEX_WITHIN(2.0 - 1.0*I, state_tiled_get(&state, 1, 0), PRECISION);
    TEST_ASSERT_COMPLEX_WITHIN(4.0 - 2.0*I, state_tiled_get(&state, 2, 0), PRECISION);
    TEST_ASSERT_COMPLEX_WITHIN(5.0 - 3.0*I, state_tiled_get(&state, 2, 1), PRECISION);
    TEST_ASSERT_COMPLEX_WITHIN(7.0 - 4.0*I, state_tiled_get(&state, 3, 0), PRECISION);

    // Verify get on the original upper position returns the original value (double conjugation roundtrip)
    TEST_ASSERT_COMPLEX_WITHIN(2.0 + 1.0*I, state_tiled_get(&state, 0, 1), PRECISION);
    TEST_ASSERT_COMPLEX_WITHIN(4.0 + 2.0*I, state_tiled_get(&state, 0, 2), PRECISION);
    TEST_ASSERT_COMPLEX_WITHIN(5.0 + 3.0*I, state_tiled_get(&state, 1, 2), PRECISION);
    TEST_ASSERT_COMPLEX_WITHIN(7.0 + 4.0*I, state_tiled_get(&state, 0, 3), PRECISION);

    state_free(&state);
}

/*
 * =====================================================================================================================
 * Test: Multi-tile end-to-end (LOG_TILE_DIM+1 qubits) - set/get in all 3 tiles
 * =====================================================================================================================
 */

void test_tiled_multi_tile(void) {
    const qubit_t n = LOG_TILE_DIM + 1;
    const idx_t T = TILE_DIM;
    state_t state = {.type = MIXED_TILED, .data = NULL, .qubits = 0};

    state_init(&state, n, NULL);  // dim = 2*T, 3 lower-triangular tiles

    // Verify storage length
    idx_t dim = POW2(n, idx_t);
    idx_t nt = (dim + TILE_DIM - 1) / TILE_DIM;
    idx_t expected_len = nt * (nt + 1) / 2 * TILE_SIZE;
    TEST_ASSERT_EQUAL_INT64(expected_len, state_tiled_len(n));
    TEST_ASSERT_EQUAL_UINT8(n, state.qubits);
    TEST_ASSERT_NOT_NULL(state.data);

    // Tile T(0,0): rows 0..T-1, cols 0..T-1
    state_tiled_set(&state, 0, 0, 1.0 + 0.0*I);
    state_tiled_set(&state, T/2, T/4, 2.0 + 1.0*I);
    state_tiled_set(&state, T-1, T-1, 3.0 + 2.0*I);

    // Tile T(1,0): rows T..2T-1, cols 0..T-1
    state_tiled_set(&state, T, 0, 4.0 + 3.0*I);
    state_tiled_set(&state, T + T/2, T/4, 5.0 + 4.0*I);
    state_tiled_set(&state, 2*T-1, T-1, 6.0 + 5.0*I);

    // Tile T(1,1): rows T..2T-1, cols T..2T-1
    state_tiled_set(&state, T, T, 7.0 + 6.0*I);
    state_tiled_set(&state, T + T/4, T + T/8, 8.0 + 7.0*I);
    state_tiled_set(&state, 2*T-1, 2*T-1, 9.0 + 8.0*I);

    // Verify get returns correct values
    TEST_ASSERT_COMPLEX_WITHIN(1.0 + 0.0*I, state_tiled_get(&state, 0, 0), PRECISION);
    TEST_ASSERT_COMPLEX_WITHIN(2.0 + 1.0*I, state_tiled_get(&state, T/2, T/4), PRECISION);
    TEST_ASSERT_COMPLEX_WITHIN(3.0 + 2.0*I, state_tiled_get(&state, T-1, T-1), PRECISION);
    TEST_ASSERT_COMPLEX_WITHIN(4.0 + 3.0*I, state_tiled_get(&state, T, 0), PRECISION);
    TEST_ASSERT_COMPLEX_WITHIN(5.0 + 4.0*I, state_tiled_get(&state, T + T/2, T/4), PRECISION);
    TEST_ASSERT_COMPLEX_WITHIN(6.0 + 5.0*I, state_tiled_get(&state, 2*T-1, T-1), PRECISION);
    TEST_ASSERT_COMPLEX_WITHIN(7.0 + 6.0*I, state_tiled_get(&state, T, T), PRECISION);
    TEST_ASSERT_COMPLEX_WITHIN(8.0 + 7.0*I, state_tiled_get(&state, T + T/4, T + T/8), PRECISION);
    TEST_ASSERT_COMPLEX_WITHIN(9.0 + 8.0*I, state_tiled_get(&state, 2*T-1, 2*T-1), PRECISION);

    state_free(&state);
}
