// test_channel_tiled.c - Test harness for channel_1q on tiled mixed states

#include "channel.h"
#include "gatemat.h"
#include "linalg.h"
#include "test_mixed_states.h"
#include "test_mixed_utils.h"
#include "test.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#define MAXQUBITS_TILED 6

/*
 * =====================================================================================================================
 * Channel constructors
 * =====================================================================================================================
 */

static kraus_t make_bitflip(double p) {
    kraus_t k = {.n_qubits = 1, .n_terms = 2, .data = malloc(2 * 4 * sizeof(cplx_t))};
    double a = sqrt(1.0 - p), b = sqrt(p);
    k.data[0] = a;  k.data[1] = 0;  k.data[2] = 0;  k.data[3] = a;
    k.data[4] = 0;  k.data[5] = b;  k.data[6] = b;  k.data[7] = 0;
    return k;
}

static kraus_t make_phaseflip(double p) {
    kraus_t k = {.n_qubits = 1, .n_terms = 2, .data = malloc(2 * 4 * sizeof(cplx_t))};
    double a = sqrt(1.0 - p), b = sqrt(p);
    k.data[0] = a;  k.data[1] = 0;  k.data[2] = 0;   k.data[3] = a;
    k.data[4] = b;  k.data[5] = 0;  k.data[6] = 0;   k.data[7] = -b;
    return k;
}

static kraus_t make_depolarizing(double p) {
    kraus_t k = {.n_qubits = 1, .n_terms = 4, .data = malloc(4 * 4 * sizeof(cplx_t))};
    double a = sqrt(1.0 - 3.0 * p / 4.0), b = sqrt(p / 4.0);
    k.data[0] = a;  k.data[1] = 0;     k.data[2] = 0;    k.data[3] = a;
    k.data[4] = 0;  k.data[5] = b;     k.data[6] = b;    k.data[7] = 0;
    k.data[8] = 0;  k.data[9] = -b*I;  k.data[10] = b*I; k.data[11] = 0;
    k.data[12] = b; k.data[13] = 0;    k.data[14] = 0;   k.data[15] = -b;
    return k;
}

static kraus_t make_ampdamp(double gamma) {
    kraus_t k = {.n_qubits = 1, .n_terms = 2, .data = malloc(2 * 4 * sizeof(cplx_t))};
    double s = sqrt(gamma), t = sqrt(1.0 - gamma);
    k.data[0] = 1;  k.data[1] = 0;  k.data[2] = 0;  k.data[3] = t;
    k.data[4] = 0;  k.data[5] = s;  k.data[6] = 0;  k.data[7] = 0;
    return k;
}

/*
 * =====================================================================================================================
 * Tiled sweep harness
 * =====================================================================================================================
 */

void testChannel1qTiled(kraus_t *kraus) {
    state_t test_state = {0};
    test_state.type = MIXED_TILED;
    cplx_t **test_rhos = NULL;
    cplx_t *rho_full = NULL, *ref_full = NULL;
    unsigned nmixed = 0;

    cplx_t **full_kraus_mat = calloc(kraus->n_terms, sizeof(*full_kraus_mat));
    if (!full_kraus_mat) {
        fprintf(stderr, "testChannel1qTiled(): full_kraus_mat allocation failed\n");
        goto cleanup;
    }

    superop_t sop = kraus_to_superop(kraus);
    if (!sop.data) {
        fprintf(stderr, "testChannel1qTiled(): kraus_to_superop failed\n");
        goto cleanup;
    }

    for (unsigned nqubits = 1; nqubits <= MAXQUBITS_TILED; ++nqubits) {
        const unsigned dim = POW2(nqubits, unsigned);

        for (unsigned target = 0; target < nqubits; ++target) {
            if (!((test_rhos = test_mk_states_mixed(nqubits, &nmixed)))) {
                fprintf(stderr, "testChannel1qTiled(): test_rhos initialisation failed\n");
                goto cleanup;
            }

            /* Extend each 2x2 Kraus operator to full n-qubit space (column-major for BLAS) */
            for (uint64_t k = 0; k < kraus->n_terms; ++k) {
                const cplx_t *rm = kraus->data + k * 4;
                cplx_t cm[4] = {rm[0], rm[2], rm[1], rm[3]};
                if (!((full_kraus_mat[k] = mat_single_qubit_gate(nqubits, cm, target)))) {
                    goto cleanup;
                }
            }

            for (unsigned i = 0; i < nmixed; ++i) {
                rho_full = density_unpack(dim, test_rhos[i]);
                if (!rho_full) {
                    fprintf(stderr, "testChannel1qTiled(): density_unpack failed\n");
                    goto cleanup;
                }

                ref_full = zsumu(dim, kraus->n_terms, (const cplx_t **)full_kraus_mat, rho_full);
                if (!ref_full) {
                    fprintf(stderr, "testChannel1qTiled(): zsumu failed\n");
                    goto cleanup;
                }

                tiled_state_from_full(&test_state, nqubits, rho_full);
                if (!test_state.data) {
                    fprintf(stderr, "testChannel1qTiled(): tiled state creation failed\n");
                    goto cleanup;
                }

                channel_1q(&test_state, &sop, target);

                assert_tiled_equals_full(&test_state, ref_full, dim);

                free(rho_full);  rho_full = NULL;
                free(ref_full);  ref_full = NULL;
                state_free(&test_state);
            }

            for (uint64_t k = 0; k < kraus->n_terms; ++k) {
                free(full_kraus_mat[k]);  full_kraus_mat[k] = NULL;
            }

            for (unsigned i = 0; i < nmixed; ++i) {
                free(test_rhos[i]);  test_rhos[i] = NULL;
            }
            free(test_rhos);  test_rhos = NULL;
        }
    }

    free(sop.data);
    free(full_kraus_mat);
    return;

cleanup:
    if (test_rhos) {
        for (unsigned i = 0; i < nmixed; ++i) {
            free(test_rhos[i]);  test_rhos[i] = NULL;
        }
    }
    free(test_rhos);
    if (full_kraus_mat) {
        for (uint64_t k = 0; k < kraus->n_terms; ++k) {
            free(full_kraus_mat[k]);  full_kraus_mat[k] = NULL;
        }
    }
    free(full_kraus_mat);
    free(rho_full);
    free(ref_full);
    state_free(&test_state);
}

/*
 * =====================================================================================================================
 * Test wrappers
 * =====================================================================================================================
 */

void test_channel_1q_bitflip_tiled(void) {
    kraus_t k = make_bitflip(0.3);
    testChannel1qTiled(&k);
    free(k.data);
}

void test_channel_1q_phaseflip_tiled(void) {
    kraus_t k = make_phaseflip(0.5);
    testChannel1qTiled(&k);
    free(k.data);
}

void test_channel_1q_depolar_tiled(void) {
    kraus_t k = make_depolarizing(0.1);
    testChannel1qTiled(&k);
    free(k.data);
}

void test_channel_1q_ampdamp_tiled(void) {
    kraus_t k = make_ampdamp(0.3);
    testChannel1qTiled(&k);
    free(k.data);
}

void test_channel_1q_identity_tiled(void) {
    cplx_t data[4] = {1, 0, 0, 1};
    kraus_t k = {.n_qubits = 1, .n_terms = 1, .data = data};
    testChannel1qTiled(&k);
}

/*
 * =====================================================================================================================
 * Analytical / boundary tests (tiled)
 * =====================================================================================================================
 */

/* Amplitude damping gamma=1 on |1><1| → |0><0| (1 qubit) */
void test_ampdamp_g1_decay_tiled(void) {
    state_t st = {0};
    st.type = MIXED_TILED;
    tiled_state_from_full(&st, 1, (cplx_t[]){0, 0, 0, 1});  /* |1><1| col-major */

    kraus_t k = make_ampdamp(1.0);
    superop_t sop = kraus_to_superop(&k);
    channel_1q(&st, &sop, 0);

    cplx_t ref[4] = {1, 0, 0, 0};  /* |0><0| col-major */
    assert_tiled_equals_full(&st, ref, 2);

    state_free(&st);
    free(sop.data);
    free(k.data);
}

/* Depolarizing p=1 on |0><0| → I/2 */
void test_depolar_full_tiled(void) {
    state_t st = {0};
    st.type = MIXED_TILED;
    tiled_state_from_full(&st, 1, (cplx_t[]){1, 0, 0, 0});  /* |0><0| */

    kraus_t k = make_depolarizing(1.0);
    superop_t sop = kraus_to_superop(&k);
    channel_1q(&st, &sop, 0);

    cplx_t ref[4] = {0.5, 0, 0, 0.5};  /* I/2 col-major */
    assert_tiled_equals_full(&st, ref, 2);

    state_free(&st);
    free(sop.data);
    free(k.data);
}

/* Phase-flip p=0.5 on |+><+| → I/2 */
void test_phaseflip_coherence_tiled(void) {
    state_t st = {0};
    st.type = MIXED_TILED;
    tiled_state_from_full(&st, 1, (cplx_t[]){0.5, 0.5, 0.5, 0.5});  /* |+><+| */

    kraus_t k = make_phaseflip(0.5);
    superop_t sop = kraus_to_superop(&k);
    channel_1q(&st, &sop, 0);

    cplx_t ref[4] = {0.5, 0, 0, 0.5};  /* I/2 */
    assert_tiled_equals_full(&st, ref, 2);

    state_free(&st);
    free(sop.data);
    free(k.data);
}

/* Bit-flip p=1 involution: apply twice = identity */
void test_bitflip_p1_involution_tiled(void) {
    kraus_t k = make_bitflip(1.0);
    superop_t sop = kraus_to_superop(&k);

    for (unsigned nqubits = 1; nqubits <= MAXQUBITS_TILED; ++nqubits) {
        const unsigned dim = POW2(nqubits, unsigned);
        unsigned nmixed = 0;
        cplx_t **test_rhos = test_mk_states_mixed(nqubits, &nmixed);
        if (!test_rhos) goto done;

        for (unsigned target = 0; target < nqubits; ++target) {
            for (unsigned i = 0; i < nmixed; ++i) {
                cplx_t *rho_full = density_unpack(dim, test_rhos[i]);
                if (!rho_full) goto done;

                /* Save original full matrix as reference */
                cplx_t *orig = malloc(dim * dim * sizeof(cplx_t));
                memcpy(orig, rho_full, dim * dim * sizeof(cplx_t));

                state_t st = {0};
                st.type = MIXED_TILED;
                tiled_state_from_full(&st, nqubits, rho_full);

                channel_1q(&st, &sop, target);
                channel_1q(&st, &sop, target);

                assert_tiled_equals_full(&st, orig, dim);

                free(rho_full);
                free(orig);
                state_free(&st);
            }
        }
        test_rm_states_mixed(test_rhos, nmixed);
    }

done:
    free(sop.data);
    free(k.data);
}

/*
 * =====================================================================================================================
 * Trace preservation tests (tiled)
 * =====================================================================================================================
 */

static void check_trace_preserved_tiled(kraus_t *kraus) {
    superop_t sop = kraus_to_superop(kraus);

    for (unsigned nqubits = 1; nqubits <= MAXQUBITS_TILED; ++nqubits) {
        const unsigned dim = POW2(nqubits, unsigned);
        unsigned nmixed = 0;
        cplx_t **test_rhos = test_mk_states_mixed(nqubits, &nmixed);
        if (!test_rhos) { free(sop.data); return; }

        for (unsigned target = 0; target < nqubits; ++target) {
            for (unsigned i = 0; i < nmixed; ++i) {
                cplx_t *rho_full = density_unpack(dim, test_rhos[i]);
                if (!rho_full) { free(sop.data); return; }

                /* Trace before */
                double tr_before = 0.0;
                for (unsigned d = 0; d < dim; ++d)
                    tr_before += creal(rho_full[d + d * dim]);

                state_t st = {0};
                st.type = MIXED_TILED;
                tiled_state_from_full(&st, nqubits, rho_full);
                free(rho_full);

                channel_1q(&st, &sop, target);

                /* Trace after via state_get */
                double tr_after = 0.0;
                for (unsigned d = 0; d < dim; ++d)
                    tr_after += creal(state_get(&st, d, d));

                TEST_ASSERT_DOUBLE_WITHIN(PRECISION, tr_before, tr_after);

                state_free(&st);
            }
        }
        test_rm_states_mixed(test_rhos, nmixed);
    }
    free(sop.data);
}

void test_trace_preserved_bitflip_tiled(void) {
    kraus_t k = make_bitflip(0.3);
    check_trace_preserved_tiled(&k);
    free(k.data);
}

void test_trace_preserved_depolar_tiled(void) {
    kraus_t k = make_depolarizing(0.5);
    check_trace_preserved_tiled(&k);
    free(k.data);
}

void test_trace_preserved_ampdamp_tiled(void) {
    kraus_t k = make_ampdamp(0.7);
    check_trace_preserved_tiled(&k);
    free(k.data);
}

/*
 * =====================================================================================================================
 * Cross-backend consistency tests
 * =====================================================================================================================
 */

static void check_packed_tiled_agree(kraus_t *kraus) {
    superop_t sop = kraus_to_superop(kraus);

    for (unsigned nqubits = 1; nqubits <= MAXQUBITS; ++nqubits) {
        const unsigned dim = POW2(nqubits, unsigned);
        unsigned nmixed = 0;
        cplx_t **test_rhos = test_mk_states_mixed(nqubits, &nmixed);
        if (!test_rhos) { free(sop.data); return; }

        for (unsigned target = 0; target < nqubits; ++target) {
            for (unsigned i = 0; i < nmixed; ++i) {
                const unsigned packed_len = dim * (dim + 1) / 2;

                /* Packed path */
                cplx_t *copy = malloc(packed_len * sizeof(cplx_t));
                memcpy(copy, test_rhos[i], packed_len * sizeof(cplx_t));
                state_t st_packed = {0};
                st_packed.type = MIXED_PACKED;
                state_init(&st_packed, nqubits, &copy);
                channel_1q(&st_packed, &sop, target);

                /* Tiled path */
                cplx_t *rho_full = density_unpack(dim, test_rhos[i]);
                state_t st_tiled = {0};
                st_tiled.type = MIXED_TILED;
                tiled_state_from_full(&st_tiled, nqubits, rho_full);
                free(rho_full);
                channel_1q(&st_tiled, &sop, target);

                /* Compare element by element */
                for (unsigned r = 0; r < dim; ++r) {
                    for (unsigned c = 0; c <= r; ++c) {
                        cplx_t packed_val = state_get(&st_packed, r, c);
                        cplx_t tiled_val = state_get(&st_tiled, r, c);
                        TEST_ASSERT_COMPLEX_WITHIN(packed_val, tiled_val, PRECISION);
                    }
                }

                free(st_packed.data); st_packed.data = NULL;
                state_free(&st_tiled);
            }
        }
        test_rm_states_mixed(test_rhos, nmixed);
    }
    free(sop.data);
}

void test_packed_tiled_agree_bitflip(void) {
    kraus_t k = make_bitflip(0.3);
    check_packed_tiled_agree(&k);
    free(k.data);
}

void test_packed_tiled_agree_phaseflip(void) {
    kraus_t k = make_phaseflip(0.5);
    check_packed_tiled_agree(&k);
    free(k.data);
}

void test_packed_tiled_agree_depolar(void) {
    kraus_t k = make_depolarizing(0.1);
    check_packed_tiled_agree(&k);
    free(k.data);
}

void test_packed_tiled_agree_ampdamp(void) {
    kraus_t k = make_ampdamp(0.3);
    check_packed_tiled_agree(&k);
    free(k.data);
}
