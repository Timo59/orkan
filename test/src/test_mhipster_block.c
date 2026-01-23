/**
 * @file test_mhipster_block.c
 * @brief Test harness for blocked density matrix gates
 *
 * Tests blocked format state management and gate operations against
 * reference implementation using full matrix multiplication.
 */

#include "state_blocked.h"
#include "mhipster_block.h"
#include "test.h"
#include "test_gate.h"
#include "gatemat.h"
#include "linalg.h"

#include <stdlib.h>
#include <string.h>

/*
 * =====================================================================================================================
 * Reference implementation (from test_mhipster.c)
 * =====================================================================================================================
 */

/**
 * @brief Insert bits at specified positions
 */
static unsigned insertBits(unsigned i, const unsigned n, unsigned bits, const unsigned *pos) {
    for (unsigned j = 0; j < n; ++j) {
        const unsigned right = ((1u << pos[j]) - 1) & i;
        unsigned set = (bits >> j) & 1;
        unsigned left = (i >> pos[j]) << 1;
        left |= set;
        left <<= pos[j];
        i = left | right;
    }
    return i;
}

/**
 * @brief Apply multi-qubit gate to full Hermitian matrix in place
 */
static void multiQubitGate(const unsigned n, cplx_t *h, const unsigned k, const unsigned *a, const cplx_t *m) {
    const unsigned N = 1 << n;
    const unsigned K = 1 << k;
    const unsigned mat_size = K * K;
    const unsigned M = 1 << (n - k);

    for (unsigned col_block = 0; col_block < M; ++col_block) {
        for (unsigned row_block = 0; row_block < M; ++row_block) {
            cplx_t tmp[mat_size];

            for (unsigned j = 0; j < K; ++j) {
                unsigned col = insertBits(col_block, k, j, a);
                for (unsigned i = 0; i < K; ++i) {
                    unsigned row = insertBits(row_block, k, i, a);
                    tmp[i + j * K] = h[row + col * N];
                    h[row + col * N] = 0;
                }
            }

            for (unsigned j = 0; j < K; ++j) {
                unsigned col = insertBits(col_block, k, j, a);
                for (unsigned i = 0; i < K; ++i) {
                    unsigned row = insertBits(row_block, k, i, a);
                    for (unsigned s = 0; s < K; ++s) {
                        for (unsigned t = 0; t < K; ++t) {
                            h[row + col * N] += m[i + s * K] * tmp[s + t * K] * conj(m[j + t * K]);
                        }
                    }
                }
            }
        }
    }
}

/*
 * =====================================================================================================================
 * Conversion utilities
 * =====================================================================================================================
 */

/**
 * @brief Convert full matrix to blocked format
 */
static void full_to_blocked(const cplx_t *full, dim_t dim, state_blocked_t *blocked) {
    for (dim_t row = 0; row < dim; ++row) {
        for (dim_t col = 0; col <= row; ++col) {
            cplx_t val = full[row + col * dim];  /* Column-major full matrix */
            state_blocked_set(blocked, row, col, val);
        }
    }
}

/**
 * @brief Convert blocked format to full matrix
 */
static cplx_t* blocked_to_full(const state_blocked_t *blocked) {
    const dim_t dim = (dim_t)1 << blocked->qubits;
    cplx_t *full = malloc(dim * dim * sizeof(*full));
    if (!full) return NULL;

    for (dim_t col = 0; col < dim; ++col) {
        for (dim_t row = 0; row < dim; ++row) {
            full[row + col * dim] = state_blocked_get(blocked, row, col);
        }
    }
    return full;
}

/**
 * @brief Create full matrix for |+⟩⟨+| state
 */
static cplx_t* plus_state_full(qubit_t qubits) {
    const dim_t dim = (dim_t)1 << qubits;
    cplx_t *full = malloc(dim * dim * sizeof(*full));
    if (!full) return NULL;

    const cplx_t prefactor = 1.0 / dim + 0.0 * I;
    for (dim_t i = 0; i < dim * dim; ++i) {
        full[i] = prefactor;
    }
    return full;
}

/*
 * =====================================================================================================================
 * State tests
 * =====================================================================================================================
 */

void test_state_blocked_len(void) {
    /* 2 qubits: dim=4, n_blocks=1, n_tiles=1, storage=64*64=4096 */
    TEST_ASSERT_EQUAL_INT(BLOCK_DIM * BLOCK_DIM, state_blocked_len(2));

    /* 6 qubits: dim=64, n_blocks=1, n_tiles=1, storage=4096 */
    TEST_ASSERT_EQUAL_INT(BLOCK_DIM * BLOCK_DIM, state_blocked_len(6));

    /* 7 qubits: dim=128, n_blocks=2, n_tiles=3, storage=3*4096=12288 */
    TEST_ASSERT_EQUAL_INT(3 * BLOCK_DIM * BLOCK_DIM, state_blocked_len(7));
}

void test_state_blocked_init(void) {
    state_blocked_t state = {0};

    state_blocked_init(&state, 3);
    TEST_ASSERT_NOT_NULL(state.data);
    TEST_ASSERT_EQUAL_INT(3, state.qubits);
    TEST_ASSERT_EQUAL_INT(1, state.n_blocks);  /* dim=8 < 64, so 1 block */
    TEST_ASSERT_EQUAL_INT(BLOCK_DIM, state.padded_dim);

    state_blocked_free(&state);
    TEST_ASSERT_NULL(state.data);
    TEST_ASSERT_EQUAL_INT(0, state.qubits);
}

void test_state_blocked_plus(void) {
    state_blocked_t state = {0};
    const qubit_t qubits = 3;
    const dim_t dim = 1 << qubits;

    state_blocked_plus(&state, qubits);
    TEST_ASSERT_NOT_NULL(state.data);

    /* Check trace = 1 */
    double trace = state_blocked_trace(&state);
    TEST_ASSERT_DOUBLE_WITHIN(PRECISION, 1.0, trace);

    /* Check all elements equal 1/dim */
    const cplx_t expected = 1.0 / dim + 0.0 * I;
    for (dim_t row = 0; row < dim; ++row) {
        for (dim_t col = 0; col <= row; ++col) {
            cplx_t val = state_blocked_get(&state, row, col);
            TEST_ASSERT_COMPLEX_WITHIN(expected, val, PRECISION);
        }
    }

    state_blocked_free(&state);
}

void test_state_blocked_get_set(void) {
    state_blocked_t state = {0};
    state_blocked_init(&state, 3);

    /* Set some values */
    state_blocked_set(&state, 0, 0, 1.0 + 0.0 * I);
    state_blocked_set(&state, 1, 0, 0.5 + 0.5 * I);
    state_blocked_set(&state, 1, 1, 0.25 + 0.0 * I);

    /* Read them back */
    TEST_ASSERT_COMPLEX_WITHIN(1.0 + 0.0 * I, state_blocked_get(&state, 0, 0), PRECISION);
    TEST_ASSERT_COMPLEX_WITHIN(0.5 + 0.5 * I, state_blocked_get(&state, 1, 0), PRECISION);
    TEST_ASSERT_COMPLEX_WITHIN(0.25 + 0.0 * I, state_blocked_get(&state, 1, 1), PRECISION);

    /* Check Hermitian conjugate access */
    TEST_ASSERT_COMPLEX_WITHIN(0.5 - 0.5 * I, state_blocked_get(&state, 0, 1), PRECISION);

    state_blocked_free(&state);
}

/*
 * =====================================================================================================================
 * Gate tests
 * =====================================================================================================================
 */

/**
 * @brief Test Z gate on blocked format against reference
 */
void test_z_blocked(void) {
    for (qubit_t nqubits = 2; nqubits <= MAXQUBITS; ++nqubits) {
        const dim_t dim = (dim_t)1 << nqubits;

        for (qubit_t target = 0; target < nqubits; ++target) {
            /* Create reference full matrix (|+⟩⟨+|) */
            cplx_t *ref_full = plus_state_full(nqubits);
            TEST_ASSERT_NOT_NULL(ref_full);

            /* Apply reference gate: U * ρ * U† */
            unsigned target_uint = target;
            multiQubitGate(nqubits, ref_full, 1, &target_uint, ZMAT);

            /* Create blocked state */
            state_blocked_t blocked = {0};
            state_blocked_plus(&blocked, nqubits);

            /* Apply blocked gate */
            z_blocked(&blocked, target);

            /* Convert blocked to full for comparison */
            cplx_t *result_full = blocked_to_full(&blocked);
            TEST_ASSERT_NOT_NULL(result_full);

            /* Compare */
            TEST_ASSERT_EQUAL_COMPLEX_ARRAY_TOL(ref_full, result_full, dim * dim, PRECISION);

            /* Cleanup */
            free(ref_full);
            free(result_full);
            state_blocked_free(&blocked);
        }
    }
}

/*
 * =====================================================================================================================
 * Unity test runner
 * =====================================================================================================================
 */

void setUp(void) {}
void tearDown(void) {}

int main(void) {
    UNITY_BEGIN();

    /* State tests */
    RUN_TEST(test_state_blocked_len);
    RUN_TEST(test_state_blocked_init);
    RUN_TEST(test_state_blocked_plus);
    RUN_TEST(test_state_blocked_get_set);

    /* Gate tests */
    RUN_TEST(test_z_blocked);

    return UNITY_END();
}
