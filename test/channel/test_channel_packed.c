// test_channel_packed.c - Test harness for channel_1q on packed mixed states

#include "channel.h"
#include "gatemat.h"
#include "linalg.h"
#include "test_mixed_states.h"
#include "test_mixed_utils.h"
#include "test.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

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
 * Packed sweep harness
 * =====================================================================================================================
 */

void testChannel1qPacked(kraus_t *kraus) {
    state_t test_state = {0};
    test_state.type = MIXED_PACKED;
    cplx_t **test_rhos = NULL;
    cplx_t *rho_full = NULL, *ref_full = NULL;
    unsigned nmixed = 0;

    cplx_t **full_kraus_mat = calloc(kraus->n_terms, sizeof(*full_kraus_mat));
    if (!full_kraus_mat) {
        fprintf(stderr, "testChannel1qPacked(): full_kraus_mat allocation failed\n");
        goto cleanup;
    }

    superop_t sop = kraus_to_superop(kraus);
    if (!sop.data) {
        fprintf(stderr, "testChannel1qPacked(): kraus_to_superop failed\n");
        goto cleanup;
    }

    for (unsigned nqubits = 1; nqubits <= MAXQUBITS; ++nqubits) {
        const unsigned dim = POW2(nqubits, unsigned);

        for (unsigned target = 0; target < nqubits; ++target) {
            if (!((test_rhos = test_mk_states_mixed(nqubits, &nmixed)))) {
                fprintf(stderr, "testChannel1qPacked(): test_rhos initialisation failed\n");
                goto cleanup;
            }

            /* Extend each 2x2 Kraus operator to full n-qubit space.
             * mat_single_qubit_gate expects column-major 2x2; kraus_t stores row-major.
             * Transpose: cm[0]=rm[0], cm[1]=rm[2], cm[2]=rm[1], cm[3]=rm[3]. */
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
                    fprintf(stderr, "testChannel1qPacked(): density_unpack failed\n");
                    goto cleanup;
                }

                ref_full = zsumu(dim, kraus->n_terms, (const cplx_t **)full_kraus_mat, rho_full);
                if (!ref_full) {
                    fprintf(stderr, "testChannel1qPacked(): zsumu failed\n");
                    goto cleanup;
                }

                cplx_t *ref = density_pack(dim, ref_full);
                if (!ref) {
                    fprintf(stderr, "testChannel1qPacked(): density_pack failed\n");
                    goto cleanup;
                }

                state_init(&test_state, nqubits, &test_rhos[i]);
                if (!test_state.data) {
                    fprintf(stderr, "testChannel1qPacked(): state_init failed\n");
                    goto cleanup;
                }

                channel_1q(&test_state, &sop, target);

                const unsigned packed_len = dim * (dim + 1) / 2;
                TEST_ASSERT_EQUAL_COMPLEX_ARRAY_TOL(ref, test_state.data, packed_len, PRECISION);

                free(ref);
                free(rho_full);  rho_full = NULL;
                free(ref_full);  ref_full = NULL;
                free(test_state.data);  test_state.data = NULL;
            }

            for (uint64_t k = 0; k < kraus->n_terms; ++k) {
                free(full_kraus_mat[k]);  full_kraus_mat[k] = NULL;
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

void test_channel_1q_bitflip_packed(void) {
    kraus_t k = make_bitflip(0.3);
    testChannel1qPacked(&k);
    free(k.data);
}

void test_channel_1q_phaseflip_packed(void) {
    kraus_t k = make_phaseflip(0.5);
    testChannel1qPacked(&k);
    free(k.data);
}

void test_channel_1q_depolar_packed(void) {
    kraus_t k = make_depolarizing(0.1);
    testChannel1qPacked(&k);
    free(k.data);
}

void test_channel_1q_ampdamp_packed(void) {
    kraus_t k = make_ampdamp(0.3);
    testChannel1qPacked(&k);
    free(k.data);
}

void test_channel_1q_identity_packed(void) {
    cplx_t data[4] = {1, 0, 0, 1};
    kraus_t k = {.n_qubits = 1, .n_terms = 1, .data = data};
    testChannel1qPacked(&k);
}

/*
 * =====================================================================================================================
 * Analytical / boundary tests (packed)
 * =====================================================================================================================
 */

/* Amplitude damping gamma=1 on |1><1| → |0><0| (1 qubit, hardcoded) */
void test_ampdamp_g1_decay_packed(void) {
    cplx_t *rho = calloc(3, sizeof(cplx_t));
    rho[2] = 1.0;  /* |1><1| packed: [rho(0,0), rho(1,0), rho(1,1)] = [0, 0, 1] */

    state_t st = {0};
    st.type = MIXED_PACKED;
    state_init(&st, 1, &rho);

    kraus_t k = make_ampdamp(1.0);
    superop_t sop = kraus_to_superop(&k);
    channel_1q(&st, &sop, 0);

    cplx_t expected[3] = {1.0, 0.0, 0.0};  /* |0><0| */
    TEST_ASSERT_EQUAL_COMPLEX_ARRAY_TOL(expected, st.data, 3, PRECISION);

    free(st.data); st.data = NULL;
    free(sop.data);
    free(k.data);
}

/* Depolarizing p=1 on |0><0| → I/2 (1 qubit, hardcoded) */
void test_depolar_full_packed(void) {
    cplx_t *rho = calloc(3, sizeof(cplx_t));
    rho[0] = 1.0;  /* |0><0| */

    state_t st = {0};
    st.type = MIXED_PACKED;
    state_init(&st, 1, &rho);

    kraus_t k = make_depolarizing(1.0);
    superop_t sop = kraus_to_superop(&k);
    channel_1q(&st, &sop, 0);

    cplx_t expected[3] = {0.5, 0.0, 0.5};  /* I/2 */
    TEST_ASSERT_EQUAL_COMPLEX_ARRAY_TOL(expected, st.data, 3, PRECISION);

    free(st.data); st.data = NULL;
    free(sop.data);
    free(k.data);
}

/* Phase-flip p=0.5 on |+><+| → I/2 (1 qubit, hardcoded) */
void test_phaseflip_coherence_packed(void) {
    cplx_t *rho = calloc(3, sizeof(cplx_t));
    rho[0] = 0.5; rho[1] = 0.5; rho[2] = 0.5;  /* |+><+| packed */

    state_t st = {0};
    st.type = MIXED_PACKED;
    state_init(&st, 1, &rho);

    kraus_t k = make_phaseflip(0.5);
    superop_t sop = kraus_to_superop(&k);
    channel_1q(&st, &sop, 0);

    cplx_t expected[3] = {0.5, 0.0, 0.5};  /* I/2 */
    TEST_ASSERT_EQUAL_COMPLEX_ARRAY_TOL(expected, st.data, 3, PRECISION);

    free(st.data); st.data = NULL;
    free(sop.data);
    free(k.data);
}

/* Bit-flip p=1: X channel applied twice is identity (involution) */
void test_bitflip_p1_involution_packed(void) {
    kraus_t k = make_bitflip(1.0);
    superop_t sop = kraus_to_superop(&k);

    for (unsigned nqubits = 1; nqubits <= MAXQUBITS; ++nqubits) {
        const unsigned dim = POW2(nqubits, unsigned);
        const unsigned packed_len = dim * (dim + 1) / 2;
        unsigned nmixed = 0;
        cplx_t **test_rhos = test_mk_states_mixed(nqubits, &nmixed);
        if (!test_rhos) goto done;

        for (unsigned target = 0; target < nqubits; ++target) {
            for (unsigned i = 0; i < nmixed; ++i) {
                cplx_t *orig = malloc(packed_len * sizeof(cplx_t));
                memcpy(orig, test_rhos[i], packed_len * sizeof(cplx_t));

                cplx_t *copy = malloc(packed_len * sizeof(cplx_t));
                memcpy(copy, test_rhos[i], packed_len * sizeof(cplx_t));

                state_t st = {0};
                st.type = MIXED_PACKED;
                state_init(&st, nqubits, &copy);

                channel_1q(&st, &sop, target);
                channel_1q(&st, &sop, target);  /* Apply twice */

                TEST_ASSERT_EQUAL_COMPLEX_ARRAY_TOL(orig, st.data, packed_len, PRECISION);

                free(orig);
                free(st.data); st.data = NULL;
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
 * Trace preservation tests (packed)
 * =====================================================================================================================
 */

static double trace_packed(unsigned dim, const cplx_t *data) {
    double tr = 0.0;
    for (unsigned i = 0; i < dim; ++i) {
        tr += creal(data[col_off(dim, i)]);  /* diagonal: offset 0 within column i */
    }
    return tr;
}

static void check_trace_preserved_packed(kraus_t *kraus) {
    superop_t sop = kraus_to_superop(kraus);

    for (unsigned nqubits = 1; nqubits <= MAXQUBITS; ++nqubits) {
        const unsigned dim = POW2(nqubits, unsigned);
        const unsigned packed_len = dim * (dim + 1) / 2;
        unsigned nmixed = 0;
        cplx_t **test_rhos = test_mk_states_mixed(nqubits, &nmixed);
        if (!test_rhos) { free(sop.data); return; }

        for (unsigned target = 0; target < nqubits; ++target) {
            for (unsigned i = 0; i < nmixed; ++i) {
                cplx_t *copy = malloc(packed_len * sizeof(cplx_t));
                memcpy(copy, test_rhos[i], packed_len * sizeof(cplx_t));

                double tr_before = trace_packed(dim, test_rhos[i]);

                state_t st = {0};
                st.type = MIXED_PACKED;
                state_init(&st, nqubits, &copy);
                channel_1q(&st, &sop, target);

                double tr_after = trace_packed(dim, st.data);
                TEST_ASSERT_DOUBLE_WITHIN(PRECISION, tr_before, tr_after);

                free(st.data); st.data = NULL;
            }
        }
        test_rm_states_mixed(test_rhos, nmixed);
    }
    free(sop.data);
}

void test_trace_preserved_bitflip_packed(void) {
    kraus_t k = make_bitflip(0.3);
    check_trace_preserved_packed(&k);
    free(k.data);
}

void test_trace_preserved_depolar_packed(void) {
    kraus_t k = make_depolarizing(0.5);
    check_trace_preserved_packed(&k);
    free(k.data);
}

void test_trace_preserved_ampdamp_packed(void) {
    kraus_t k = make_ampdamp(0.7);
    check_trace_preserved_packed(&k);
    free(k.data);
}
