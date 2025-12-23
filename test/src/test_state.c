#ifndef Q_TEST_H
#include "test.h"
#endif

#ifndef STATE_H
#include "state.h"
#endif

#ifndef UTILS_H
#include "utils.h"
#endif

#include <math.h>
#include <stdlib.h>

/*
 * =====================================================================================================================
 * Unity: setUp and tearDown
 * =====================================================================================================================
 */

void setUp(void) {}

void tearDown(void) {}

/*
 * =====================================================================================================================
 * Test: state_len()
 * =====================================================================================================================
 */

void test_state_len_pure(void) {
    state_t state;
    state.type = PURE;

    state.qubits = 1;
    TEST_ASSERT_EQUAL_INT64(2, state_len(&state));

    state.qubits = 2;
    TEST_ASSERT_EQUAL_INT64(4, state_len(&state));

    state.qubits = 3;
    TEST_ASSERT_EQUAL_INT64(8, state_len(&state));

    state.qubits = 4;
    TEST_ASSERT_EQUAL_INT64(16, state_len(&state));
}

void test_state_len_mixed(void) {
    state_t state;
    state.type = MIXED;

    state.qubits = 1;
    TEST_ASSERT_EQUAL_INT64(3, state_len(&state));

    state.qubits = 2;
    TEST_ASSERT_EQUAL_INT64(10, state_len(&state));

    state.qubits = 3;
    TEST_ASSERT_EQUAL_INT64(36, state_len(&state));

    state.qubits = 4;
    TEST_ASSERT_EQUAL_INT64(136, state_len(&state));
}

/*
 * =====================================================================================================================
 * Test: state_init()
 * =====================================================================================================================
 */

void test_state_init_pure_null_data(void) {
    state_t state = {.type = PURE, .data = NULL, .qubits = 0};

    state_init(&state, 2, NULL);

    TEST_ASSERT_EQUAL_UINT8(2, state.qubits);
    TEST_ASSERT_NOT_NULL(state.data);

    // Verify zero initialization
    dim_t len = state_len(&state);
    for (dim_t i = 0; i < len; i++) {
        TEST_ASSERT_COMPLEX_WITHIN(0.0 + 0.0*I, state.data[i], PRECISION);
    }

    state_free(&state);
}

void test_state_init_pure_with_data(void) {
    state_t state = {.type = PURE, .data = NULL, .qubits = 0};

    // Allocate and initialize data
    cplx_t *data = malloc(4 * sizeof(cplx_t));
    TEST_ASSERT_NOT_NULL(data);
    data[0] = 0.25 + 0.0*I;
    data[1] = 0.0 + 0.25*I;
    data[2] = -0.25 + 0.0*I;
    data[3] = 0.0 - 0.25*I;

    cplx_t *original_ptr = data;

    state_init(&state, 2, &data);

    // Verify ownership transfer
    TEST_ASSERT_NULL(data);
    TEST_ASSERT_EQUAL_PTR(original_ptr, state.data);

    // Verify values preserved
    TEST_ASSERT_COMPLEX_WITHIN(0.25 + 0.0*I, state.data[0], PRECISION);
    TEST_ASSERT_COMPLEX_WITHIN(0.0 + 0.25*I, state.data[1], PRECISION);
    TEST_ASSERT_COMPLEX_WITHIN(-0.25 + 0.0*I, state.data[2], PRECISION);
    TEST_ASSERT_COMPLEX_WITHIN(0.0 - 0.25*I, state.data[3], PRECISION);

    state_free(&state);
}

void test_state_init_mixed_null_data(void) {
    state_t state = {.type = MIXED, .data = NULL, .qubits = 0};

    state_init(&state, 2, NULL);

    TEST_ASSERT_EQUAL_UINT8(2, state.qubits);
    TEST_ASSERT_NOT_NULL(state.data);

    // Verify zero initialization
    dim_t len = state_len(&state);
    for (dim_t i = 0; i < len; i++) {
        TEST_ASSERT_COMPLEX_WITHIN(0.0 + 0.0*I, state.data[i], PRECISION);
    }

    state_free(&state);
}

void test_state_init_mixed_with_data(void) {
    state_t state = {.type = MIXED, .data = NULL, .qubits = 0};

    // Allocate and initialize data (10 elements for 2-qubit mixed state)
    cplx_t *data = malloc(10 * sizeof(cplx_t));
    TEST_ASSERT_NOT_NULL(data);
    for (int i = 0; i < 10; i++) {
        data[i] = (double)i + 0.5*I;
    }

    cplx_t *original_ptr = data;

    state_init(&state, 2, &data);

    // Verify ownership transfer
    TEST_ASSERT_NULL(data);
    TEST_ASSERT_EQUAL_PTR(original_ptr, state.data);

    // Verify values preserved
    for (int i = 0; i < 10; i++) {
        TEST_ASSERT_COMPLEX_WITHIN((double)i + 0.5*I, state.data[i], PRECISION);
    }

    state_free(&state);
}

void test_state_init_reinitialization(void) {
    state_t state = {.type = PURE, .data = NULL, .qubits = 0};

    // First initialization
    state_init(&state, 2, NULL);
    TEST_ASSERT_EQUAL_UINT8(2, state.qubits);
    TEST_ASSERT_NOT_NULL(state.data);
    cplx_t *first_ptr = state.data;

    // Modify some data to verify it gets replaced
    state.data[0] = 99.0 + 99.0*I;

    // Re-initialize with different size
    state_init(&state, 3, NULL);
    TEST_ASSERT_EQUAL_UINT8(3, state.qubits);
    TEST_ASSERT_NOT_NULL(state.data);

    // Verify new allocation (different pointer or at least different size)
    dim_t len = state_len(&state);
    TEST_ASSERT_EQUAL_INT64(8, len);

    // Verify zero initialization of new array
    for (dim_t i = 0; i < len; i++) {
        TEST_ASSERT_COMPLEX_WITHIN(0.0 + 0.0*I, state.data[i], PRECISION);
    }

    state_free(&state);
}

/*
 * =====================================================================================================================
 * Test: state_free()
 * =====================================================================================================================
 */

void test_state_free_pure(void) {
    state_t state = {.type = PURE, .data = NULL, .qubits = 0};

    state_init(&state, 3, NULL);
    TEST_ASSERT_NOT_NULL(state.data);
    TEST_ASSERT_EQUAL_UINT8(3, state.qubits);

    state_free(&state);

    TEST_ASSERT_NULL(state.data);
    TEST_ASSERT_EQUAL_UINT8(0, state.qubits);

    // Safe to call again (double-free safety)
    state_free(&state);
    TEST_ASSERT_NULL(state.data);
}

void test_state_free_mixed(void) {
    state_t state = {.type = MIXED, .data = NULL, .qubits = 0};

    state_init(&state, 2, NULL);
    TEST_ASSERT_NOT_NULL(state.data);
    TEST_ASSERT_EQUAL_UINT8(2, state.qubits);

    state_free(&state);

    TEST_ASSERT_NULL(state.data);
    TEST_ASSERT_EQUAL_UINT8(0, state.qubits);

    // Safe to call again (double-free safety)
    state_free(&state);
    TEST_ASSERT_NULL(state.data);
}

/*
 * =====================================================================================================================
 * Test: state_plus()
 * =====================================================================================================================
 */

void test_state_plus_pure_normalization(void) {
    state_t state = {.type = PURE, .data = NULL, .qubits = 0};

    state_plus(&state, 3);

    dim_t len = state_len(&state);
    TEST_ASSERT_EQUAL_INT64(8, len);

    double expected_amplitude = 1.0 / sqrt(8.0);

    // Check all amplitudes
    for (dim_t i = 0; i < len; i++) {
        double real_part = creal(state.data[i]);
        double imag_part = cimag(state.data[i]);
        TEST_ASSERT_DOUBLE_WITHIN(PRECISION, expected_amplitude, real_part);
        TEST_ASSERT_DOUBLE_WITHIN(PRECISION, 0.0, imag_part);
    }

    // Verify normalization
    double norm_sq = 0.0;
    for (dim_t i = 0; i < len; i++) {
        double mag = cabs(state.data[i]);
        norm_sq += mag * mag;
    }
    TEST_ASSERT_DOUBLE_WITHIN(PRECISION, 1.0, norm_sq);

    state_free(&state);
}

void test_state_plus_mixed_trace(void) {
    state_t state = {.type = MIXED, .data = NULL, .qubits = 0};

    state_plus(&state, 2);

    dim_t dim = POW2(2, dim_t);  // Dimension = 4
    double expected_entry = 1.0 / (double)dim;  // 1/4

    // For mixed state, all entries should be 1/4
    dim_t len = state_len(&state);
    TEST_ASSERT_EQUAL_INT64(10, len);
    for (dim_t i = 0; i < len; i++) {
        double real_part = creal(state.data[i]);
        double imag_part = cimag(state.data[i]);
        TEST_ASSERT_DOUBLE_WITHIN(PRECISION, expected_entry, real_part);
        TEST_ASSERT_DOUBLE_WITHIN(PRECISION, 0.0, imag_part);
    }

    state_free(&state);
}

/*
 * =====================================================================================================================
 * Test: state_cp()
 * =====================================================================================================================
 */

void test_state_cp_pure_deep_copy(void) {
    state_t original = {.type = PURE, .data = NULL, .qubits = 0};

    state_plus(&original, 2);

    state_t copy = state_cp(&original);

    // Verify metadata copied
    TEST_ASSERT_EQUAL_UINT8(original.type, copy.type);
    TEST_ASSERT_EQUAL_UINT8(original.qubits, copy.qubits);

    // Verify different allocations
    TEST_ASSERT_NOT_EQUAL_PTR(original.data, copy.data);

    // Verify values identical
    dim_t len = state_len(&original);
    for (dim_t i = 0; i < len; i++) {
        TEST_ASSERT_COMPLEX_WITHIN(original.data[i], copy.data[i], PRECISION);
    }

    // Modify original
    original.data[0] = 99.0 + 99.0*I;

    // Verify copy unchanged
    double expected = 1.0 / sqrt(4.0);
    TEST_ASSERT_DOUBLE_WITHIN(PRECISION, expected, creal(copy.data[0]));

    state_free(&original);
    state_free(&copy);
}

void test_state_cp_mixed_deep_copy(void) {
    state_t original = {.type = MIXED, .data = NULL, .qubits = 0};

    state_plus(&original, 2);

    state_t copy = state_cp(&original);

    // Verify metadata copied
    TEST_ASSERT_EQUAL_UINT8(original.type, copy.type);
    TEST_ASSERT_EQUAL_UINT8(original.qubits, copy.qubits);

    // Verify different allocations
    TEST_ASSERT_NOT_EQUAL_PTR(original.data, copy.data);

    // Verify values identical
    dim_t len = state_len(&original);
    for (dim_t i = 0; i < len; i++) {
        TEST_ASSERT_COMPLEX_WITHIN(original.data[i], copy.data[i], PRECISION);
    }

    // Modify original
    original.data[0] = 99.0 + 99.0*I;

    // Verify copy unchanged (should be 0.25 from state_plus)
    TEST_ASSERT_DOUBLE_WITHIN(PRECISION, 0.25, creal(copy.data[0]));

    state_free(&original);
    state_free(&copy);
}

/*
 * =====================================================================================================================
 * Test: Edge Cases
 * =====================================================================================================================
 */

void test_state_single_qubit(void) {
    // Test pure state with 1 qubit
    state_t pure = {.type = PURE, .data = NULL, .qubits = 0};
    state_init(&pure, 1, NULL);

    TEST_ASSERT_EQUAL_UINT8(1, pure.qubits);
    TEST_ASSERT_EQUAL_INT64(2, state_len(&pure));
    TEST_ASSERT_NOT_NULL(pure.data);

    // Set to |0⟩ state
    pure.data[0] = 1.0 + 0.0*I;
    pure.data[1] = 0.0 + 0.0*I;

    // Verify normalization
    double norm_sq = cabs(pure.data[0]) * cabs(pure.data[0]) +
                      cabs(pure.data[1]) * cabs(pure.data[1]);
    TEST_ASSERT_DOUBLE_WITHIN(PRECISION, 1.0, norm_sq);

    state_free(&pure);

    // Test mixed state with 1 qubit
    state_t mixed = {.type = MIXED, .data = NULL, .qubits = 0};
    state_init(&mixed, 1, NULL);

    TEST_ASSERT_EQUAL_UINT8(1, mixed.qubits);
    TEST_ASSERT_EQUAL_INT64(3, state_len(&mixed));  // 2*3/2 = 3
    TEST_ASSERT_NOT_NULL(mixed.data);

    state_free(&mixed);
}

void test_state_many_qubits(void) {
    // Test with 10 qubits - mainly to verify no overflow
    state_t pure = {.type = PURE, .data = NULL, .qubits = 0};
    state_init(&pure, 10, NULL);

    TEST_ASSERT_EQUAL_UINT8(10, pure.qubits);

    dim_t expected_len = POW2(10, dim_t);  // 1024
    TEST_ASSERT_EQUAL_INT64(expected_len, state_len(&pure));

    // Verify allocation succeeded
    TEST_ASSERT_NOT_NULL(pure.data);

    state_free(&pure);

    // Test mixed state with smaller qubit count due to memory constraints
    state_t mixed = {.type = MIXED, .data = NULL, .qubits = 0};
    state_init(&mixed, 5, NULL);  // 5 qubits = 32x32 matrix = 528 elements

    TEST_ASSERT_EQUAL_UINT8(5, mixed.qubits);

    dim_t dim = POW2(5, dim_t);  // 32
    dim_t expected_mixed_len = dim * (dim + 1) / 2;  // 32*33/2 = 528
    TEST_ASSERT_EQUAL_INT64(expected_mixed_len, state_len(&mixed));

    // Verify allocation succeeded
    TEST_ASSERT_NOT_NULL(mixed.data);

    state_free(&mixed);
}

/*
 * =====================================================================================================================
 * Unity Test Runner
 * =====================================================================================================================
 */

int main(void) {
    UNITY_BEGIN();

    // state_len() tests
    RUN_TEST(test_state_len_pure);
    RUN_TEST(test_state_len_mixed);

    // state_init() tests
    RUN_TEST(test_state_init_pure_null_data);
    RUN_TEST(test_state_init_pure_with_data);
    RUN_TEST(test_state_init_mixed_null_data);
    RUN_TEST(test_state_init_mixed_with_data);
    RUN_TEST(test_state_init_reinitialization);

    // state_free() tests
    RUN_TEST(test_state_free_pure);
    RUN_TEST(test_state_free_mixed);

    // state_plus() tests
    RUN_TEST(test_state_plus_pure_normalization);
    RUN_TEST(test_state_plus_mixed_trace);

    // state_cp() tests
    RUN_TEST(test_state_cp_pure_deep_copy);
    RUN_TEST(test_state_cp_mixed_deep_copy);

    // Edge case tests
    RUN_TEST(test_state_single_qubit);
    RUN_TEST(test_state_many_qubits);

    return UNITY_END();
}
