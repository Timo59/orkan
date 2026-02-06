// test_state_packed.c - Unit tests for MIXED_PACKED quantum states (packed density matrices)

#include "test.h"
#include "state.h"
#include "utils.h"

#include <math.h>
#include <stdlib.h>
#include <complex.h>

/*
 * =====================================================================================================================
 * Test: state_packed_len()
 * =====================================================================================================================
 */

void test_packed_len(void) {
    // Test returns dim*(dim+1)/2 for 1-4 qubits
    // 1 qubit: dim=2, len = 2*3/2 = 3
    TEST_ASSERT_EQUAL_INT64(3, state_packed_len(1));

    // 2 qubits: dim=4, len = 4*5/2 = 10
    TEST_ASSERT_EQUAL_INT64(10, state_packed_len(2));

    // 3 qubits: dim=8, len = 8*9/2 = 36
    TEST_ASSERT_EQUAL_INT64(36, state_packed_len(3));

    // 4 qubits: dim=16, len = 16*17/2 = 136
    TEST_ASSERT_EQUAL_INT64(136, state_packed_len(4));
}

/*
 * =====================================================================================================================
 * Test: state_init() with NULL data - zero initialization
 * =====================================================================================================================
 */

void test_packed_init_null_data(void) {
    for (qubit_t n = 1; n <= 4; n++) {
        state_t state = {.type = MIXED_PACKED, .data = NULL, .qubits = 0};

        state_init(&state, n, NULL);

        // Verify initialization
        TEST_ASSERT_EQUAL_UINT8(n, state.qubits);
        TEST_ASSERT_NOT_NULL(state.data);

        // Verify all elements are zero
        dim_t len = state_packed_len(n);
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

void test_packed_init_ownership(void) {
    state_t state = {.type = MIXED_PACKED, .data = NULL, .qubits = 0};

    // Allocate data externally (2 qubits = 10 elements, 64-byte aligned)
    dim_t len = state_packed_len(2);
    size_t size = (len * sizeof(cplx_t) + 63) & ~(size_t)63;
    cplx_t *data = aligned_alloc(64, size);
    TEST_ASSERT_NOT_NULL(data);

    // Set values
    for (dim_t i = 0; i < len; i++) {
        data[i] = (double)i + (double)i * 0.1 * I;
    }

    cplx_t *original_ptr = data;

    // Transfer ownership
    state_init(&state, 2, &data);

    // Verify ownership transferred
    TEST_ASSERT_NULL(data);  // Source pointer nullified
    TEST_ASSERT_EQUAL_PTR(original_ptr, state.data);  // State owns the data

    // Verify values preserved
    for (dim_t i = 0; i < len; i++) {
        cplx_t expected = (double)i + (double)i * 0.1 * I;
        TEST_ASSERT_COMPLEX_WITHIN(expected, state.data[i], PRECISION);
    }

    state_free(&state);
}

/*
 * =====================================================================================================================
 * Test: state_free() - safe to call twice (double-free safety)
 * =====================================================================================================================
 */

void test_packed_free_double_free(void) {
    state_t state = {.type = MIXED_PACKED, .data = NULL, .qubits = 0};

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
 * Test: state_packed_plus() - all elements = 1/dim (|+><+| state)
 * =====================================================================================================================
 */

void test_packed_plus_values(void) {
    // Test for 1-4 qubits
    for (qubit_t n = 1; n <= 4; n++) {
        state_t state = {.type = MIXED_PACKED, .data = NULL, .qubits = 0};

        state_packed_plus(&state, n);

        TEST_ASSERT_NOT_NULL(state.data);
        TEST_ASSERT_EQUAL_UINT8(n, state.qubits);

        dim_t dim = POW2(n, dim_t);
        dim_t len = state_packed_len(n);
        double expected_value = 1.0 / (double)dim;

        // For |+><+| state, all elements should be 1/dim
        for (dim_t i = 0; i < len; i++) {
            double real_part = creal(state.data[i]);
            double imag_part = cimag(state.data[i]);

            TEST_ASSERT_DOUBLE_WITHIN(PRECISION, expected_value, real_part);
            TEST_ASSERT_DOUBLE_WITHIN(PRECISION, 0.0, imag_part);
        }

        state_free(&state);
    }
}

/*
 * =====================================================================================================================
 * Test: state_packed_get() - lower triangle (row >= col) direct access
 * =====================================================================================================================
 */

void test_packed_get_lower(void) {
    state_t state = {.type = MIXED_PACKED, .data = NULL, .qubits = 0};

    state_init(&state, 2, NULL);  // 4x4 matrix, 10 elements in packed storage

    // Manually set lower triangle elements to known values
    // Column-major packed storage (LAPACK 'L' format):
    // Column 0: (0,0), (1,0), (2,0), (3,0) -> indices 0, 1, 2, 3
    // Column 1: (1,1), (2,1), (3,1)       -> indices 4, 5, 6
    // Column 2: (2,2), (3,2)              -> indices 7, 8
    // Column 3: (3,3)                     -> index 9

    state.data[0] = 1.0 + 0.0*I;   // (0,0)
    state.data[1] = 2.0 + 1.0*I;   // (1,0)
    state.data[2] = 3.0 + 2.0*I;   // (2,0)
    state.data[3] = 4.0 + 3.0*I;   // (3,0)
    state.data[4] = 5.0 + 0.0*I;   // (1,1)
    state.data[5] = 6.0 + 4.0*I;   // (2,1)
    state.data[6] = 7.0 + 5.0*I;   // (3,1)
    state.data[7] = 8.0 + 0.0*I;   // (2,2)
    state.data[8] = 9.0 + 6.0*I;   // (3,2)
    state.data[9] = 10.0 + 0.0*I;  // (3,3)

    // Test lower triangle access (row >= col)
    TEST_ASSERT_COMPLEX_WITHIN(1.0 + 0.0*I, state_packed_get(&state, 0, 0), PRECISION);
    TEST_ASSERT_COMPLEX_WITHIN(2.0 + 1.0*I, state_packed_get(&state, 1, 0), PRECISION);
    TEST_ASSERT_COMPLEX_WITHIN(5.0 + 0.0*I, state_packed_get(&state, 1, 1), PRECISION);
    TEST_ASSERT_COMPLEX_WITHIN(3.0 + 2.0*I, state_packed_get(&state, 2, 0), PRECISION);
    TEST_ASSERT_COMPLEX_WITHIN(6.0 + 4.0*I, state_packed_get(&state, 2, 1), PRECISION);
    TEST_ASSERT_COMPLEX_WITHIN(8.0 + 0.0*I, state_packed_get(&state, 2, 2), PRECISION);
    TEST_ASSERT_COMPLEX_WITHIN(4.0 + 3.0*I, state_packed_get(&state, 3, 0), PRECISION);
    TEST_ASSERT_COMPLEX_WITHIN(7.0 + 5.0*I, state_packed_get(&state, 3, 1), PRECISION);
    TEST_ASSERT_COMPLEX_WITHIN(9.0 + 6.0*I, state_packed_get(&state, 3, 2), PRECISION);
    TEST_ASSERT_COMPLEX_WITHIN(10.0 + 0.0*I, state_packed_get(&state, 3, 3), PRECISION);

    state_free(&state);
}

/*
 * =====================================================================================================================
 * Test: state_packed_get() - upper triangle (row < col) returns conjugate
 * =====================================================================================================================
 */

void test_packed_get_upper(void) {
    state_t state = {.type = MIXED_PACKED, .data = NULL, .qubits = 0};

    state_init(&state, 2, NULL);  // 4x4 matrix

    // Set some lower triangle values with non-zero imaginary parts
    // Column-major packed storage:
    // Column 0: indices 0-3 -> (0,0), (1,0), (2,0), (3,0)
    // Column 1: indices 4-6 -> (1,1), (2,1), (3,1)
    state.data[1] = 2.0 + 1.0*I;   // (1,0) at index 1
    state.data[2] = 4.0 + 2.0*I;   // (2,0) at index 2
    state.data[5] = 5.0 + 3.0*I;   // (2,1) at index 5
    state.data[3] = 7.0 + 4.0*I;   // (3,0) at index 3

    // Test upper triangle returns conjugate (row < col)
    // (0,1) should return conj((1,0))
    TEST_ASSERT_COMPLEX_WITHIN(2.0 - 1.0*I, state_packed_get(&state, 0, 1), PRECISION);

    // (0,2) should return conj((2,0))
    TEST_ASSERT_COMPLEX_WITHIN(4.0 - 2.0*I, state_packed_get(&state, 0, 2), PRECISION);

    // (1,2) should return conj((2,1))
    TEST_ASSERT_COMPLEX_WITHIN(5.0 - 3.0*I, state_packed_get(&state, 1, 2), PRECISION);

    // (0,3) should return conj((3,0))
    TEST_ASSERT_COMPLEX_WITHIN(7.0 - 4.0*I, state_packed_get(&state, 0, 3), PRECISION);

    state_free(&state);
}

/*
 * =====================================================================================================================
 * Test: state_packed_set() - lower triangle (row >= col) direct storage
 * =====================================================================================================================
 */

void test_packed_set_lower(void) {
    state_t state = {.type = MIXED_PACKED, .data = NULL, .qubits = 0};

    state_init(&state, 2, NULL);  // 4x4 matrix

    // Set lower triangle elements
    state_packed_set(&state, 0, 0, 1.0 + 0.0*I);
    state_packed_set(&state, 1, 0, 2.0 + 1.0*I);
    state_packed_set(&state, 2, 0, 3.0 + 2.0*I);
    state_packed_set(&state, 3, 0, 4.0 + 3.0*I);
    state_packed_set(&state, 1, 1, 5.0 + 0.0*I);
    state_packed_set(&state, 2, 1, 6.0 + 4.0*I);
    state_packed_set(&state, 3, 1, 7.0 + 5.0*I);
    state_packed_set(&state, 2, 2, 8.0 + 0.0*I);
    state_packed_set(&state, 3, 2, 9.0 + 6.0*I);
    state_packed_set(&state, 3, 3, 10.0 + 0.0*I);

    // Verify stored values in packed array (column-major order)
    // Column 0: (0,0), (1,0), (2,0), (3,0) -> indices 0, 1, 2, 3
    // Column 1: (1,1), (2,1), (3,1)       -> indices 4, 5, 6
    // Column 2: (2,2), (3,2)              -> indices 7, 8
    // Column 3: (3,3)                     -> index 9
    TEST_ASSERT_COMPLEX_WITHIN(1.0 + 0.0*I, state.data[0], PRECISION);   // (0,0)
    TEST_ASSERT_COMPLEX_WITHIN(2.0 + 1.0*I, state.data[1], PRECISION);   // (1,0)
    TEST_ASSERT_COMPLEX_WITHIN(3.0 + 2.0*I, state.data[2], PRECISION);   // (2,0)
    TEST_ASSERT_COMPLEX_WITHIN(4.0 + 3.0*I, state.data[3], PRECISION);   // (3,0)
    TEST_ASSERT_COMPLEX_WITHIN(5.0 + 0.0*I, state.data[4], PRECISION);   // (1,1)
    TEST_ASSERT_COMPLEX_WITHIN(6.0 + 4.0*I, state.data[5], PRECISION);   // (2,1)
    TEST_ASSERT_COMPLEX_WITHIN(7.0 + 5.0*I, state.data[6], PRECISION);   // (3,1)
    TEST_ASSERT_COMPLEX_WITHIN(8.0 + 0.0*I, state.data[7], PRECISION);   // (2,2)
    TEST_ASSERT_COMPLEX_WITHIN(9.0 + 6.0*I, state.data[8], PRECISION);   // (3,2)
    TEST_ASSERT_COMPLEX_WITHIN(10.0 + 0.0*I, state.data[9], PRECISION);  // (3,3)

    state_free(&state);
}

/*
 * =====================================================================================================================
 * Test: state_packed_set() - upper triangle (row < col) stores conj at transposed position
 * =====================================================================================================================
 */

void test_packed_set_upper(void) {
    state_t state = {.type = MIXED_PACKED, .data = NULL, .qubits = 0};

    state_init(&state, 2, NULL);  // 4x4 matrix

    // Set upper triangle elements (row < col)
    // These should be stored as conjugates at the transposed position
    state_packed_set(&state, 0, 1, 2.0 + 1.0*I);  // Should store conj(2+i) at (1,0) -> index 1
    state_packed_set(&state, 0, 2, 4.0 + 2.0*I);  // Should store conj(4+2i) at (2,0) -> index 2
    state_packed_set(&state, 1, 2, 5.0 + 3.0*I);  // Should store conj(5+3i) at (2,1) -> index 5
    state_packed_set(&state, 0, 3, 7.0 + 4.0*I);  // Should store conj(7+4i) at (3,0) -> index 3

    // Verify conjugates stored at lower triangle positions (column-major)
    // Column 0: indices 0-3 -> (0,0), (1,0), (2,0), (3,0)
    // Column 1: indices 4-6 -> (1,1), (2,1), (3,1)
    TEST_ASSERT_COMPLEX_WITHIN(2.0 - 1.0*I, state.data[1], PRECISION);  // (1,0) at index 1
    TEST_ASSERT_COMPLEX_WITHIN(4.0 - 2.0*I, state.data[2], PRECISION);  // (2,0) at index 2
    TEST_ASSERT_COMPLEX_WITHIN(5.0 - 3.0*I, state.data[5], PRECISION);  // (2,1) at index 5
    TEST_ASSERT_COMPLEX_WITHIN(7.0 - 4.0*I, state.data[3], PRECISION);  // (3,0) at index 3

    // Verify get retrieves original values (via double conjugation)
    TEST_ASSERT_COMPLEX_WITHIN(2.0 + 1.0*I, state_packed_get(&state, 0, 1), PRECISION);
    TEST_ASSERT_COMPLEX_WITHIN(4.0 + 2.0*I, state_packed_get(&state, 0, 2), PRECISION);
    TEST_ASSERT_COMPLEX_WITHIN(5.0 + 3.0*I, state_packed_get(&state, 1, 2), PRECISION);
    TEST_ASSERT_COMPLEX_WITHIN(7.0 + 4.0*I, state_packed_get(&state, 0, 3), PRECISION);

    state_free(&state);
}

/*
 * =====================================================================================================================
 * Test: state_cp() - independence of copy
 * =====================================================================================================================
 */

void test_packed_cp_independence(void) {
    state_t original = {.type = MIXED_PACKED, .data = NULL, .qubits = 0};

    state_packed_plus(&original, 2);  // |+><+| state

    // Create copy
    state_t copy = state_cp(&original);

    // Verify metadata copied
    TEST_ASSERT_EQUAL_UINT8(original.type, copy.type);
    TEST_ASSERT_EQUAL_UINT8(original.qubits, copy.qubits);

    // Verify different memory allocations
    TEST_ASSERT_NOT_NULL(copy.data);
    TEST_ASSERT_NOT_EQUAL_PTR(original.data, copy.data);

    // Verify values initially identical
    dim_t len = state_packed_len(2);
    for (dim_t i = 0; i < len; i++) {
        TEST_ASSERT_COMPLEX_WITHIN(original.data[i], copy.data[i], PRECISION);
    }

    // Modify original
    double expected_original = 0.25;  // 1/4
    original.data[0] = 99.0 + 99.0*I;
    original.data[1] = -50.0 - 50.0*I;

    // Verify original modified
    TEST_ASSERT_DOUBLE_WITHIN(PRECISION, 99.0, creal(original.data[0]));
    TEST_ASSERT_DOUBLE_WITHIN(PRECISION, 99.0, cimag(original.data[0]));
    TEST_ASSERT_DOUBLE_WITHIN(PRECISION, -50.0, creal(original.data[1]));
    TEST_ASSERT_DOUBLE_WITHIN(PRECISION, -50.0, cimag(original.data[1]));

    // Verify copy unchanged
    TEST_ASSERT_DOUBLE_WITHIN(PRECISION, expected_original, creal(copy.data[0]));
    TEST_ASSERT_DOUBLE_WITHIN(PRECISION, 0.0, cimag(copy.data[0]));
    TEST_ASSERT_DOUBLE_WITHIN(PRECISION, expected_original, creal(copy.data[1]));
    TEST_ASSERT_DOUBLE_WITHIN(PRECISION, 0.0, cimag(copy.data[1]));

    state_free(&original);
    state_free(&copy);
}

/*
 * =====================================================================================================================
 * Test: Hermitian symmetry - ρ[i,j] = conj(ρ[j,i])
 * =====================================================================================================================
 */

void test_packed_hermitian_symmetry(void) {
    state_t state = {.type = MIXED_PACKED, .data = NULL, .qubits = 0};

    state_init(&state, 2, NULL);  // 4x4 matrix

    // Set some complex values in lower triangle
    state_packed_set(&state, 1, 0, 2.0 + 1.5*I);
    state_packed_set(&state, 2, 0, 3.0 + 2.5*I);
    state_packed_set(&state, 2, 1, 4.0 + 3.5*I);
    state_packed_set(&state, 3, 1, 5.0 + 4.5*I);

    // Verify Hermitian symmetry
    cplx_t val_10 = state_packed_get(&state, 1, 0);
    cplx_t val_01 = state_packed_get(&state, 0, 1);
    TEST_ASSERT_COMPLEX_WITHIN(conj(val_10), val_01, PRECISION);

    cplx_t val_20 = state_packed_get(&state, 2, 0);
    cplx_t val_02 = state_packed_get(&state, 0, 2);
    TEST_ASSERT_COMPLEX_WITHIN(conj(val_20), val_02, PRECISION);

    cplx_t val_21 = state_packed_get(&state, 2, 1);
    cplx_t val_12 = state_packed_get(&state, 1, 2);
    TEST_ASSERT_COMPLEX_WITHIN(conj(val_21), val_12, PRECISION);

    cplx_t val_31 = state_packed_get(&state, 3, 1);
    cplx_t val_13 = state_packed_get(&state, 1, 3);
    TEST_ASSERT_COMPLEX_WITHIN(conj(val_31), val_13, PRECISION);

    state_free(&state);
}

/*
 * =====================================================================================================================
 * Test: Edge case - single qubit
 * =====================================================================================================================
 */

void test_packed_single_qubit(void) {
    state_t state = {.type = MIXED_PACKED, .data = NULL, .qubits = 0};

    state_packed_plus(&state, 1);  // 1 qubit, 2x2 matrix, 3 elements

    TEST_ASSERT_EQUAL_UINT8(1, state.qubits);
    TEST_ASSERT_EQUAL_INT64(3, state_packed_len(1));
    TEST_ASSERT_NOT_NULL(state.data);

    // |+><+| single qubit: all elements = 1/2
    double expected = 0.5;
    TEST_ASSERT_DOUBLE_WITHIN(PRECISION, expected, creal(state.data[0]));  // (0,0)
    TEST_ASSERT_DOUBLE_WITHIN(PRECISION, expected, creal(state.data[1]));  // (1,0)
    TEST_ASSERT_DOUBLE_WITHIN(PRECISION, expected, creal(state.data[2]));  // (1,1)

    state_free(&state);
}

/*
 * =====================================================================================================================
 * Test: Zero-qubit state - special case for 0 qubits
 * =====================================================================================================================
 */

void test_packed_zero_qubit(void) {
    // Test that state_packed_len(0) returns 1 (special case for 0-qubit = 1x1 scalar density matrix)
    TEST_ASSERT_EQUAL_INT64(1, state_packed_len(0));

    state_t state = {.type = MIXED_PACKED, .data = NULL, .qubits = 0};

    state_init(&state, 0, NULL);

    // Verify initialization
    TEST_ASSERT_EQUAL_UINT8(0, state.qubits);
    TEST_ASSERT_NOT_NULL(state.data);
    TEST_ASSERT_EQUAL_INT64(1, state_packed_len(0));

    // Set the single element to 1.0 (valid scalar density matrix)
    state_packed_set(&state, 0, 0, 1.0 + 0.0*I);

    // Verify get returns the value
    cplx_t val = state_packed_get(&state, 0, 0);
    TEST_ASSERT_COMPLEX_WITHIN(1.0 + 0.0*I, val, PRECISION);

    state_free(&state);
}
