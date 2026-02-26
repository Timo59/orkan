// test_circuit.c - Unity runner for circuit-module unit tests
//
// Forward declarations for all tests defined in test_pauli.c.
// setUp/tearDown are Unity required hooks (no per-test state needed here).

/*
 * =====================================================================================================================
 * Includes
 * =====================================================================================================================
 */

#ifndef Q_TEST_H
#include "test.h"
#endif

/*
 * =====================================================================================================================
 * Unity hooks
 * =====================================================================================================================
 */

void setUp(void)    {}
void tearDown(void) {}

/*
 * =====================================================================================================================
 * Forward declarations — construction and masks
 * =====================================================================================================================
 */

void test_pauli_masks_X(void);
void test_pauli_masks_Y(void);
void test_pauli_masks_Z(void);
void test_pauli_masks_str(void);
void test_pauli_from_str_negative(void);

/*
 * =====================================================================================================================
 * Forward declarations — pauli_phase
 * =====================================================================================================================
 */

void test_pauli_phase(void);

/*
 * =====================================================================================================================
 * Forward declarations — single-qubit application (n=1)
 * =====================================================================================================================
 */

void test_pauli_apply_I_q0(void);
void test_pauli_apply_X_q0(void);
void test_pauli_apply_Y_q0(void);
void test_pauli_apply_Z_q0(void);

/*
 * =====================================================================================================================
 * Forward declarations — two-qubit strings (n=2)
 * =====================================================================================================================
 */

void test_pauli_apply_XX(void);
void test_pauli_apply_ZZ(void);
void test_pauli_apply_XZ(void);
void test_pauli_apply_YY(void);
void test_pauli_apply_IZ(void);
void test_pauli_apply_IX(void);
void test_pauli_apply_IY(void);

/*
 * =====================================================================================================================
 * Forward declarations — three-qubit strings (n=3)
 * =====================================================================================================================
 */

void test_pauli_apply_XYZ(void);
void test_pauli_apply_ZZZ(void);
void test_pauli_apply_IXI(void);
void test_pauli_apply_XII(void);
void test_pauli_apply_IZI(void);
void test_pauli_apply_ZIZ(void);

/*
 * =====================================================================================================================
 * Forward declarations — four-qubit strings (n=4)
 * =====================================================================================================================
 */

void test_pauli_apply_4q(void);
void test_pauli_apply_4q_split1(void);
void test_pauli_apply_4q_split2(void);
void test_pauli_apply_4q_split3(void);

/*
 * =====================================================================================================================
 * Forward declarations — round-trip and involution
 * =====================================================================================================================
 */

void test_pauli_from_str_roundtrip(void);
void test_pauli_apply_involution(void);

/*
 * =====================================================================================================================
 * Forward declarations — identity / edge cases
 * =====================================================================================================================
 */

void test_pauli_apply_identity_noop(void);

/*
 * =====================================================================================================================
 * Unity body
 * =====================================================================================================================
 */

int main(void) {
    UNITY_BEGIN();

    printf("\n========== Pauli Construction Tests ==========\n");
    RUN_TEST(test_pauli_masks_X);
    RUN_TEST(test_pauli_masks_Y);
    RUN_TEST(test_pauli_masks_Z);
    RUN_TEST(test_pauli_masks_str);
    RUN_TEST(test_pauli_from_str_negative);

    printf("\n========== pauli_phase Tests ==========\n");
    RUN_TEST(test_pauli_phase);

    printf("\n========== Single-Qubit Application Tests (n=1) ==========\n");
    RUN_TEST(test_pauli_apply_I_q0);
    RUN_TEST(test_pauli_apply_X_q0);
    RUN_TEST(test_pauli_apply_Y_q0);
    RUN_TEST(test_pauli_apply_Z_q0);

    printf("\n========== Two-Qubit String Tests (n=2) ==========\n");
    RUN_TEST(test_pauli_apply_XX);
    RUN_TEST(test_pauli_apply_ZZ);
    RUN_TEST(test_pauli_apply_XZ);
    RUN_TEST(test_pauli_apply_YY);
    RUN_TEST(test_pauli_apply_IZ);
    RUN_TEST(test_pauli_apply_IX);
    RUN_TEST(test_pauli_apply_IY);

    printf("\n========== Three-Qubit String Tests (n=3) ==========\n");
    RUN_TEST(test_pauli_apply_XYZ);
    RUN_TEST(test_pauli_apply_ZZZ);
    RUN_TEST(test_pauli_apply_IXI);
    RUN_TEST(test_pauli_apply_XII);
    RUN_TEST(test_pauli_apply_IZI);
    RUN_TEST(test_pauli_apply_ZIZ);

    printf("\n========== Four-Qubit String Tests (n=4) ==========\n");
    RUN_TEST(test_pauli_apply_4q);
    RUN_TEST(test_pauli_apply_4q_split1);
    RUN_TEST(test_pauli_apply_4q_split2);
    RUN_TEST(test_pauli_apply_4q_split3);

    printf("\n========== Round-Trip and Involution Tests ==========\n");
    RUN_TEST(test_pauli_from_str_roundtrip);
    RUN_TEST(test_pauli_apply_involution);

    printf("\n========== Identity / Edge Case Tests ==========\n");
    RUN_TEST(test_pauli_apply_identity_noop);

    return UNITY_END();
}
