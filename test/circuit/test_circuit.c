// test_circuit.c - Unity runner for circuit-module unit tests
//
// Forward declarations for all tests defined in test_pauli.c.
// setUp/tearDown are Unity required hooks; state generation is done per-test
// inside each test function using test_mk_states_pure / test_rm_states_pure.

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
void test_pauli_new_n0(void);
void test_pauli_new_n64(void);
void test_pauli_free_null(void);
void test_pauli_from_str_empty(void);
void test_pauli_mask_invariant(void);
void test_pauli_to_str(void);

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
 * Forward declarations — Pauli application (all cases in APPLY_CASES table)
 * =====================================================================================================================
 */

void test_pauli_apply_all(void);

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
    RUN_TEST(test_pauli_to_str);
    RUN_TEST(test_pauli_new_n0);
    RUN_TEST(test_pauli_new_n64);
    RUN_TEST(test_pauli_free_null);
    RUN_TEST(test_pauli_from_str_empty);
    RUN_TEST(test_pauli_mask_invariant);

    printf("\n========== pauli_phase Tests ==========\n");
    RUN_TEST(test_pauli_phase);
    RUN_TEST(test_pauli_phase_eta2);
    RUN_TEST(test_pauli_phase_eta3);

    printf("\n========== Pauli Application Tests (all Paulis) ==========\n");
    RUN_TEST(test_pauli_apply_all);

    printf("\n========== Round-Trip and Involution Tests ==========\n");
    RUN_TEST(test_pauli_from_str_roundtrip);
    RUN_TEST(test_pauli_apply_involution);
    RUN_TEST(test_pauli_apply_involution_Z_nontrivial);

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
