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
void test_pauli_new_n0(void);
void test_pauli_new_null_labels(void);
void test_pauli_new_too_many_qubits(void);
void test_pauli_new_n64(void);
void test_pauli_free_null(void);
void test_pauli_from_str_case_sensitive(void);
void test_pauli_from_str_empty(void);
void test_pauli_mask_invariant(void);

/*
 * =====================================================================================================================
 * Forward declarations — pauli_phase
 * =====================================================================================================================
 */

void test_pauli_phase(void);
void test_pauli_phase_eta2(void);
void test_pauli_phase_eta3(void);

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
void test_pauli_apply_YZ(void);

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
void test_pauli_apply_YYY(void);

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
void test_pauli_apply_involution_Z_nontrivial(void);

/*
 * =====================================================================================================================
 * Forward declarations — identity / edge cases
 * =====================================================================================================================
 */

void test_pauli_apply_identity_noop(void);

/*
 * =====================================================================================================================
 * Forward declarations — pauli_exp_pure
 * =====================================================================================================================
 */

void test_pauli_exp_identity_global_phase(void);
void test_pauli_exp_pure_z(void);
void test_pauli_exp_pure_nondiaganal_X(void);
void test_pauli_exp_pure_nondiaganal_Y(void);
void test_pauli_exp_pure_nondiaganal_XX(void);
void test_pauli_exp_pure_nondiaganal_YY(void);
void test_pauli_exp_pure_nondiaganal_YYY(void);
void test_pauli_exp_pure_nondiaganal_XYZ(void);
void test_pauli_exp_theta_zero_is_identity(void);
void test_pauli_exp_theta_pi_is_neg_i_pauli(void);
void test_pauli_exp_theta_2pi_is_neg_identity(void);
void test_pauli_exp_4q(void);
void test_pauli_exp_unitarity(void);
void test_pauli_exp_inverse(void);
void test_pauli_exp_eigenstate_phases(void);
void test_pauli_exp_theta_4pi_is_identity(void);
void test_pauli_exp_consistency_with_apply(void);
void test_pauli_exp_non_unit_norm(void);
void test_pauli_exp_nondiaganal_split_bit(void);

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
    RUN_TEST(test_pauli_new_n0);
    RUN_TEST(test_pauli_new_null_labels);
    RUN_TEST(test_pauli_new_too_many_qubits);
    RUN_TEST(test_pauli_new_n64);
    RUN_TEST(test_pauli_free_null);
    RUN_TEST(test_pauli_from_str_case_sensitive);
    RUN_TEST(test_pauli_from_str_empty);
    RUN_TEST(test_pauli_mask_invariant);

    printf("\n========== pauli_phase Tests ==========\n");
    RUN_TEST(test_pauli_phase);
    RUN_TEST(test_pauli_phase_eta2);
    RUN_TEST(test_pauli_phase_eta3);

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
    RUN_TEST(test_pauli_apply_YZ);

    printf("\n========== Three-Qubit String Tests (n=3) ==========\n");
    RUN_TEST(test_pauli_apply_XYZ);
    RUN_TEST(test_pauli_apply_ZZZ);
    RUN_TEST(test_pauli_apply_IXI);
    RUN_TEST(test_pauli_apply_XII);
    RUN_TEST(test_pauli_apply_IZI);
    RUN_TEST(test_pauli_apply_ZIZ);
    RUN_TEST(test_pauli_apply_YYY);

    printf("\n========== Four-Qubit String Tests (n=4) ==========\n");
    RUN_TEST(test_pauli_apply_4q);
    RUN_TEST(test_pauli_apply_4q_split1);
    RUN_TEST(test_pauli_apply_4q_split2);
    RUN_TEST(test_pauli_apply_4q_split3);

    printf("\n========== Round-Trip and Involution Tests ==========\n");
    RUN_TEST(test_pauli_from_str_roundtrip);
    RUN_TEST(test_pauli_apply_involution);
    RUN_TEST(test_pauli_apply_involution_Z_nontrivial);

    printf("\n========== Identity / Edge Case Tests ==========\n");
    RUN_TEST(test_pauli_apply_identity_noop);

    printf("\n========== pauli_exp_pure Tests ==========\n");
    RUN_TEST(test_pauli_exp_identity_global_phase);
    RUN_TEST(test_pauli_exp_pure_z);
    RUN_TEST(test_pauli_exp_pure_nondiaganal_X);
    RUN_TEST(test_pauli_exp_pure_nondiaganal_Y);
    RUN_TEST(test_pauli_exp_pure_nondiaganal_XX);
    RUN_TEST(test_pauli_exp_pure_nondiaganal_YY);
    RUN_TEST(test_pauli_exp_pure_nondiaganal_YYY);
    RUN_TEST(test_pauli_exp_pure_nondiaganal_XYZ);
    RUN_TEST(test_pauli_exp_theta_zero_is_identity);
    RUN_TEST(test_pauli_exp_theta_pi_is_neg_i_pauli);
    RUN_TEST(test_pauli_exp_theta_2pi_is_neg_identity);
    RUN_TEST(test_pauli_exp_4q);
    RUN_TEST(test_pauli_exp_unitarity);
    RUN_TEST(test_pauli_exp_inverse);
    RUN_TEST(test_pauli_exp_eigenstate_phases);
    RUN_TEST(test_pauli_exp_theta_4pi_is_identity);
    RUN_TEST(test_pauli_exp_consistency_with_apply);
    RUN_TEST(test_pauli_exp_non_unit_norm);
    RUN_TEST(test_pauli_exp_nondiaganal_split_bit);

    return UNITY_END();
}
