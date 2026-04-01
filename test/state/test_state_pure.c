// test_state_pure.c - Unit tests for PURE quantum states (state vectors)

#include "test.h"
#include "state_internal.h"
#include "utils.h"

#include <math.h>
#include <stdlib.h>
#include <complex.h>

/*
 * =====================================================================================================================
 * Test: state_pure_len()
 * =====================================================================================================================
 */

void test_pure_len(void) {
    // Test returns 2^n for 1-4 qubits
    TEST_ASSERT_EQUAL_INT64(2, state_pure_len(1));      // 2^1 = 2
    TEST_ASSERT_EQUAL_INT64(4, state_pure_len(2));      // 2^2 = 4
    TEST_ASSERT_EQUAL_INT64(8, state_pure_len(3));      // 2^3 = 8
    TEST_ASSERT_EQUAL_INT64(16, state_pure_len(4));     // 2^4 = 16
}

/*
 * =====================================================================================================================
 * Test: state_init() with NULL data - zero initialization
 * =====================================================================================================================
 */

void test_pure_init_null_data(void) {
    for (qubit_t n = 1; n <= 4; n++) {
        state_t state = {.type = PURE, .data = NULL, .qubits = 0};

        state_init(&state, n, NULL);

        // Verify initialization
        TEST_ASSERT_EQUAL_UINT8(n, state.qubits);
        TEST_ASSERT_NOT_NULL(state.data);

        // Verify all elements are zero
        idx_t len = state_pure_len(n);
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

void test_pure_init_ownership(void) {
    state_t state = {.type = PURE, .data = NULL, .qubits = 0};

    // Allocate data externally (64-byte aligned for SIMD gate kernels)
    idx_t len = state_pure_len(3);  // 8 elements for 3 qubits
    size_t size = (len * sizeof(cplx_t) + 63) & ~(size_t)63;
    cplx_t *data = aligned_alloc(64, size);
    TEST_ASSERT_NOT_NULL(data);

    // Set values
    for (idx_t i = 0; i < len; i++) {
        data[i] = (double)i + (double)i * I;
    }

    cplx_t *original_ptr = data;

    // Transfer ownership
    state_init(&state, 3, &data);

    // Verify ownership transferred
    TEST_ASSERT_NULL(data);  // Source pointer nullified
    TEST_ASSERT_EQUAL_PTR(original_ptr, state.data);  // State owns the data

    // Verify values preserved
    for (idx_t i = 0; i < len; i++) {
        cplx_t expected = (double)i + (double)i * I;
        TEST_ASSERT_COMPLEX_WITHIN(expected, state.data[i], PRECISION);
    }

    state_free(&state);
}

/*
 * =====================================================================================================================
 * Test: state_free() - safe to call twice (double-free safety)
 * =====================================================================================================================
 */

void test_pure_free_double_free(void) {
    state_t state = {.type = PURE, .data = NULL, .qubits = 0};

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
 * Test: state_pure_plus() - all elements = 1/sqrt(dim)
 * =====================================================================================================================
 */

void test_pure_plus_values(void) {
    // Test for 1-4 qubits
    for (qubit_t n = 1; n <= 4; n++) {
        state_t state = {.type = PURE, .data = NULL, .qubits = 0};

        state_pure_plus(&state, n);

        TEST_ASSERT_NOT_NULL(state.data);
        TEST_ASSERT_EQUAL_UINT8(n, state.qubits);

        idx_t len = state_pure_len(n);
        double expected_amplitude = 1.0 / sqrt((double)len);

        // Check all amplitudes are equal to 1/sqrt(dim)
        for (idx_t i = 0; i < len; i++) {
            double real_part = creal(state.data[i]);
            double imag_part = cimag(state.data[i]);

            TEST_ASSERT_DOUBLE_WITHIN(PRECISION, expected_amplitude, real_part);
            TEST_ASSERT_DOUBLE_WITHIN(PRECISION, 0.0, imag_part);
        }

        // Verify normalization: Σ|cᵢ|² = 1
        double norm_sq = 0.0;
        for (idx_t i = 0; i < len; i++) {
            double mag = cabs(state.data[i]);
            norm_sq += mag * mag;
        }
        TEST_ASSERT_DOUBLE_WITHIN(PRECISION, 1.0, norm_sq);

        state_free(&state);
    }
}

/*
 * =====================================================================================================================
 * Test: state_pure_get() and state_pure_set() - roundtrip
 * =====================================================================================================================
 */

void test_pure_get_set_roundtrip(void) {
    state_t state = {.type = PURE, .data = NULL, .qubits = 0};

    state_init(&state, 3, NULL);  // 8 elements

    // Test various complex values
    cplx_t test_values[] = {
        0.5 + 0.0*I,
        0.0 + 0.5*I,
        -0.3 + 0.4*I,
        0.7 - 0.2*I,
        1.0 + 1.0*I,
        -1.0 - 1.0*I,
        0.0 + 0.0*I,
        0.123456789 + 0.987654321*I
    };

    idx_t len = state_pure_len(3);

    // Set all values
    for (idx_t i = 0; i < len; i++) {
        state_pure_set(&state, i, test_values[i]);
    }

    // Get and verify all values
    for (idx_t i = 0; i < len; i++) {
        cplx_t retrieved = state_pure_get(&state, i);
        TEST_ASSERT_COMPLEX_WITHIN(test_values[i], retrieved, PRECISION);
    }

    state_free(&state);
}

/*
 * =====================================================================================================================
 * Test: state_cp() - independence of copy
 * =====================================================================================================================
 */

void test_pure_cp_independence(void) {
    state_t original = {.type = PURE, .data = NULL, .qubits = 0};

    state_pure_plus(&original, 2);  // |++++⟩ state

    // Create copy
    state_t copy = state_cp(&original);

    // Verify metadata copied
    TEST_ASSERT_EQUAL_UINT8(original.type, copy.type);
    TEST_ASSERT_EQUAL_UINT8(original.qubits, copy.qubits);

    // Verify different memory allocations
    TEST_ASSERT_NOT_NULL(copy.data);
    TEST_ASSERT_NOT_EQUAL_PTR(original.data, copy.data);

    // Verify values initially identical
    idx_t len = state_pure_len(2);
    for (idx_t i = 0; i < len; i++) {
        TEST_ASSERT_COMPLEX_WITHIN(original.data[i], copy.data[i], PRECISION);
    }

    // Modify original
    double expected_original = 1.0 / sqrt(4.0);
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
 * Test: Edge case - single qubit
 * =====================================================================================================================
 */

void test_pure_single_qubit(void) {
    state_t state = {.type = PURE, .data = NULL, .qubits = 0};

    state_init(&state, 1, NULL);

    TEST_ASSERT_EQUAL_UINT8(1, state.qubits);
    TEST_ASSERT_EQUAL_INT64(2, state_pure_len(1));
    TEST_ASSERT_NOT_NULL(state.data);

    // Set to |0⟩ state
    state.data[0] = 1.0 + 0.0*I;
    state.data[1] = 0.0 + 0.0*I;

    // Verify values
    TEST_ASSERT_COMPLEX_WITHIN(1.0 + 0.0*I, state.data[0], PRECISION);
    TEST_ASSERT_COMPLEX_WITHIN(0.0 + 0.0*I, state.data[1], PRECISION);

    // Verify normalization
    double norm_sq = cabs(state.data[0]) * cabs(state.data[0]) +
                     cabs(state.data[1]) * cabs(state.data[1]);
    TEST_ASSERT_DOUBLE_WITHIN(PRECISION, 1.0, norm_sq);

    state_free(&state);
}

/*
 * =====================================================================================================================
 * Test: Boundary test - maximum test qubits
 * =====================================================================================================================
 */

void test_pure_four_qubits(void) {
    state_t state = {.type = PURE, .data = NULL, .qubits = 0};

    state_pure_plus(&state, 4);  // 16-element state

    TEST_ASSERT_EQUAL_UINT8(4, state.qubits);
    TEST_ASSERT_EQUAL_INT64(16, state_pure_len(4));
    TEST_ASSERT_NOT_NULL(state.data);

    double expected_amplitude = 1.0 / 4.0;  // 1/sqrt(16)

    // Verify all 16 elements
    for (idx_t i = 0; i < 16; i++) {
        TEST_ASSERT_DOUBLE_WITHIN(PRECISION, expected_amplitude, creal(state.data[i]));
        TEST_ASSERT_DOUBLE_WITHIN(PRECISION, 0.0, cimag(state.data[i]));
    }

    state_free(&state);
}

/*
 * =====================================================================================================================
 * Test: Reinitialize already-initialized state
 * =====================================================================================================================
 */

void test_pure_reinit(void) {
    state_t state = {.type = PURE, .data = NULL, .qubits = 0};

    // First initialization with 2 qubits
    state_init(&state, 2, NULL);
    TEST_ASSERT_EQUAL_UINT8(2, state.qubits);
    TEST_ASSERT_NOT_NULL(state.data);
    TEST_ASSERT_EQUAL_INT64(4, state_pure_len(2));

    // Set some values
    state.data[0] = 1.0 + 1.0*I;
    state.data[1] = 2.0 + 2.0*I;
    state.data[2] = 3.0 + 3.0*I;
    state.data[3] = 4.0 + 4.0*I;

    // Reinitialize the same state with 3 qubits
    state_init(&state, 3, NULL);

    // Verify the state now has 3 qubits
    TEST_ASSERT_EQUAL_UINT8(3, state.qubits);
    TEST_ASSERT_NOT_NULL(state.data);
    TEST_ASSERT_EQUAL_INT64(8, state_pure_len(3));

    // Verify all new elements are zero (old data was freed)
    idx_t len = state_pure_len(3);
    for (idx_t i = 0; i < len; i++) {
        TEST_ASSERT_COMPLEX_WITHIN(0.0 + 0.0*I, state.data[i], PRECISION);
    }

    state_free(&state);
}
