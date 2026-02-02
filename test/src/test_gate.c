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
void test_H_pure(void) {testSingleQubitGate(h, HMAT);}
void test_S_pure(void) {testSingleQubitGate(s, SMAT);}
void test_Sdg_pure(void) {testSingleQubitGate(sdg, SDGMAT);}
void test_T_pure(void) {testSingleQubitGate(t, TMAT);}
void test_Tdg_pure(void) {testSingleQubitGate(tdg, TDGMAT);}
void test_Rx_pure(void) {testRotationGate(rx, mat_rx);}
void test_Ry_pure(void) {testRotationGate(ry, mat_ry);}
void test_Rz_pure(void) {testRotationGate(rz, mat_rz);}

/*
 * =====================================================================================================================
 * Test functions: Mixed states
 * =====================================================================================================================
 */

void test_X_mixed(void) {testSingleQubitGateMixed(x, XMAT);}
void test_Y_mixed(void) {testSingleQubitGateMixed(y, YMAT);}
void test_Z_mixed(void) {testSingleQubitGateMixed(z, ZMAT);}
void test_H_mixed(void) {testSingleQubitGateMixed(h, HMAT);}
void test_S_mixed(void) {testSingleQubitGateMixed(s, SMAT);}
void test_Sdg_mixed(void) {testSingleQubitGateMixed(sdg, SDGMAT);}
void test_T_mixed(void) {testSingleQubitGateMixed(t, TMAT);}
void test_Tdg_mixed(void) {testSingleQubitGateMixed(tdg, TDGMAT);}

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
    RUN_TEST(test_H_pure);
    RUN_TEST(test_S_pure);
    RUN_TEST(test_Sdg_pure);
    RUN_TEST(test_T_pure);
    RUN_TEST(test_Tdg_pure);
    RUN_TEST(test_Rx_pure);
    RUN_TEST(test_Ry_pure);
    RUN_TEST(test_Rz_pure);

    // Mixed state tests
    RUN_TEST(test_X_mixed);
    RUN_TEST(test_Y_mixed);
    RUN_TEST(test_Z_mixed);
    RUN_TEST(test_H_mixed);
    RUN_TEST(test_S_mixed);
    RUN_TEST(test_Sdg_mixed);
    RUN_TEST(test_T_mixed);
    RUN_TEST(test_Tdg_mixed);

    return UNITY_END();
}
