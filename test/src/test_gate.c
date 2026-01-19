// test_gates.c  - Unit tests for quantum gate operations on pure and mixed states

/*
 * =====================================================================================================================
 * includes
 * =====================================================================================================================
 */

#ifndef TEST_GATE_H
#include "test_gate.h"
#endif

/*
 * =====================================================================================================================
 * Unity: setUp and tearDown
 * =====================================================================================================================
 */

void setUp(void) {}

void tearDown(void) {}


/*
 * =====================================================================================================================
 * Test functions: Pure states
 * =====================================================================================================================
 */

void test_X_pure(void) {testSingleQubitGate(x, XMAT);}
void test_Y_pure(void) {testSingleQubitGate(y, YMAT);}
void test_Z_pure(void) {testSingleQubitGate(z, ZMAT);}

/*
 * =====================================================================================================================
 * Test functions: Mixed states
 * =====================================================================================================================
 */

void test_X_mixed(void) {testSingleQubitGateMixed(x, XMAT);}
void test_Y_mixed(void) {testSingleQubitGateMixed(y, YMAT);}
void test_Z_mixed(void) {testSingleQubitGateMixed(z, ZMAT);}

/*
 * =====================================================================================================================
 * Test functions: Error handling
 * =====================================================================================================================
 */

void test_error_null_state(void) {
    qs_error_t err = x(NULL, 0);
    TEST_ASSERT_EQUAL_INT(QS_ERR_NULL, err);
}

void test_error_null_data(void) {
    state_t state = {0};
    state.type = PURE;
    state.qubits = 2;
    state.data = NULL;
    qs_error_t err = x(&state, 0);
    TEST_ASSERT_EQUAL_INT(QS_ERR_NULL, err);
}

void test_error_qubit_out_of_range(void) {
    state_t state = {0};
    state_plus(&state, 2);  // 2-qubit state, valid targets are 0 and 1

    qs_error_t err = x(&state, 2);  // Target 2 is out of range
    TEST_ASSERT_EQUAL_INT(QS_ERR_QUBIT, err);

    err = x(&state, 255);  // Way out of range
    TEST_ASSERT_EQUAL_INT(QS_ERR_QUBIT, err);

    state_free(&state);
}

/*
 * =====================================================================================================================
 * Unity body
 * =====================================================================================================================
 */

int main(void) {
    UNITY_BEGIN();

    // Error handling tests
    RUN_TEST(test_error_null_state);
    RUN_TEST(test_error_null_data);
    RUN_TEST(test_error_qubit_out_of_range);

    // Pure state tests
    RUN_TEST(test_X_pure);
    RUN_TEST(test_Y_pure);
    RUN_TEST(test_Z_pure);

    // Mixed state tests
    RUN_TEST(test_X_mixed);
    RUN_TEST(test_Y_mixed);
    RUN_TEST(test_Z_mixed);

    return UNITY_END();
}
