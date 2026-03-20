// test_channel.c - Unity runner for quantum channel tests

#include "test.h"

void setUp(void) {}
void tearDown(void) {}

/* Forward declarations: test_channel_superop.c */
void test_superop_identity(void);
void test_superop_x_gate(void);
void test_superop_bitflip_p0(void);
void test_superop_ampdamp_g0(void);
void test_superop_ampdamp_g1(void);
void test_superop_depolarizing_full(void);
void test_superop_dephasing(void);

/* Forward declarations: test_channel_packed.c — sweep */
void test_channel_1q_identity_packed(void);
void test_channel_1q_bitflip_packed(void);
void test_channel_1q_phaseflip_packed(void);
void test_channel_1q_depolar_packed(void);
void test_channel_1q_ampdamp_packed(void);

/* Forward declarations: test_channel_packed.c — analytical / boundary */
void test_ampdamp_g1_decay_packed(void);
void test_depolar_full_packed(void);
void test_phaseflip_coherence_packed(void);
void test_bitflip_p1_involution_packed(void);

/* Forward declarations: test_channel_packed.c — trace preservation */
void test_trace_preserved_bitflip_packed(void);
void test_trace_preserved_depolar_packed(void);
void test_trace_preserved_ampdamp_packed(void);

/* Forward declarations: test_channel_tiled.c — sweep */
void test_channel_1q_identity_tiled(void);
void test_channel_1q_bitflip_tiled(void);
void test_channel_1q_phaseflip_tiled(void);
void test_channel_1q_depolar_tiled(void);
void test_channel_1q_ampdamp_tiled(void);

/* Forward declarations: test_channel_tiled.c — analytical / boundary */
void test_ampdamp_g1_decay_tiled(void);
void test_depolar_full_tiled(void);
void test_phaseflip_coherence_tiled(void);
void test_bitflip_p1_involution_tiled(void);

/* Forward declarations: test_channel_tiled.c — trace preservation */
void test_trace_preserved_bitflip_tiled(void);
void test_trace_preserved_depolar_tiled(void);
void test_trace_preserved_ampdamp_tiled(void);

/* Forward declarations: test_channel_tiled.c — cross-backend consistency */
void test_packed_tiled_agree_bitflip(void);
void test_packed_tiled_agree_phaseflip(void);
void test_packed_tiled_agree_depolar(void);
void test_packed_tiled_agree_ampdamp(void);

int main(void) {
    UNITY_BEGIN();

    printf("\n========== kraus_to_superop Tests ==========\n");
    RUN_TEST(test_superop_identity);
    RUN_TEST(test_superop_x_gate);
    RUN_TEST(test_superop_bitflip_p0);
    RUN_TEST(test_superop_ampdamp_g0);
    RUN_TEST(test_superop_ampdamp_g1);
    RUN_TEST(test_superop_depolarizing_full);
    RUN_TEST(test_superop_dephasing);

    printf("\n========== MIXED_PACKED Channel Tests ==========\n");
    RUN_TEST(test_channel_1q_identity_packed);
    RUN_TEST(test_channel_1q_bitflip_packed);
    RUN_TEST(test_channel_1q_phaseflip_packed);
    RUN_TEST(test_channel_1q_depolar_packed);
    RUN_TEST(test_channel_1q_ampdamp_packed);
    RUN_TEST(test_ampdamp_g1_decay_packed);
    RUN_TEST(test_depolar_full_packed);
    RUN_TEST(test_phaseflip_coherence_packed);
    RUN_TEST(test_bitflip_p1_involution_packed);
    RUN_TEST(test_trace_preserved_bitflip_packed);
    RUN_TEST(test_trace_preserved_depolar_packed);
    RUN_TEST(test_trace_preserved_ampdamp_packed);

    printf("\n========== MIXED_TILED Channel Tests ==========\n");
    RUN_TEST(test_channel_1q_identity_tiled);
    RUN_TEST(test_channel_1q_bitflip_tiled);
    RUN_TEST(test_channel_1q_phaseflip_tiled);
    RUN_TEST(test_channel_1q_depolar_tiled);
    RUN_TEST(test_channel_1q_ampdamp_tiled);
    RUN_TEST(test_ampdamp_g1_decay_tiled);
    RUN_TEST(test_depolar_full_tiled);
    RUN_TEST(test_phaseflip_coherence_tiled);
    RUN_TEST(test_bitflip_p1_involution_tiled);
    RUN_TEST(test_trace_preserved_bitflip_tiled);
    RUN_TEST(test_trace_preserved_depolar_tiled);
    RUN_TEST(test_trace_preserved_ampdamp_tiled);

    printf("\n========== Cross-Backend Consistency Tests ==========\n");
    RUN_TEST(test_packed_tiled_agree_bitflip);
    RUN_TEST(test_packed_tiled_agree_phaseflip);
    RUN_TEST(test_packed_tiled_agree_depolar);
    RUN_TEST(test_packed_tiled_agree_ampdamp);

    return UNITY_END();
}
