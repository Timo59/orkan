// test_channel_superop.c - Standalone unit tests for kraus_to_superop()

#include "channel.h"
#include "test.h"

#include <math.h>
#include <stdlib.h>

/*
 * =====================================================================================================================
 * Channel constructors
 * =====================================================================================================================
 */

static kraus_t make_bitflip(double p) {
    kraus_t k = {.n_qubits = 1, .n_terms = 2, .data = malloc(2 * 4 * sizeof(cplx_t))};
    double a = sqrt(1.0 - p), b = sqrt(p);
    /* K0 = sqrt(1-p)*I  (row-major) */
    k.data[0] = a;  k.data[1] = 0;  k.data[2] = 0;  k.data[3] = a;
    /* K1 = sqrt(p)*X    (row-major) */
    k.data[4] = 0;  k.data[5] = b;  k.data[6] = b;  k.data[7] = 0;
    return k;
}

static kraus_t make_phaseflip(double p) {
    kraus_t k = {.n_qubits = 1, .n_terms = 2, .data = malloc(2 * 4 * sizeof(cplx_t))};
    double a = sqrt(1.0 - p), b = sqrt(p);
    /* K0 = sqrt(1-p)*I */
    k.data[0] = a;  k.data[1] = 0;  k.data[2] = 0;   k.data[3] = a;
    /* K1 = sqrt(p)*Z */
    k.data[4] = b;  k.data[5] = 0;  k.data[6] = 0;   k.data[7] = -b;
    return k;
}

static kraus_t make_depolarizing(double p) {
    kraus_t k = {.n_qubits = 1, .n_terms = 4, .data = malloc(4 * 4 * sizeof(cplx_t))};
    double a = sqrt(1.0 - 3.0 * p / 4.0), b = sqrt(p / 4.0);
    /* K0 = sqrt(1-3p/4)*I */
    k.data[0] = a;  k.data[1] = 0;     k.data[2] = 0;    k.data[3] = a;
    /* K1 = sqrt(p/4)*X */
    k.data[4] = 0;  k.data[5] = b;     k.data[6] = b;    k.data[7] = 0;
    /* K2 = sqrt(p/4)*Y */
    k.data[8] = 0;  k.data[9] = -b*I;  k.data[10] = b*I; k.data[11] = 0;
    /* K3 = sqrt(p/4)*Z */
    k.data[12] = b; k.data[13] = 0;    k.data[14] = 0;   k.data[15] = -b;
    return k;
}

static kraus_t make_ampdamp(double gamma) {
    kraus_t k = {.n_qubits = 1, .n_terms = 2, .data = malloc(2 * 4 * sizeof(cplx_t))};
    double s = sqrt(gamma), t = sqrt(1.0 - gamma);
    /* K0 = [[1,0],[0,sqrt(1-gamma)]]  (row-major) */
    k.data[0] = 1;  k.data[1] = 0;  k.data[2] = 0;  k.data[3] = t;
    /* K1 = [[0,sqrt(gamma)],[0,0]]    (row-major) */
    k.data[4] = 0;  k.data[5] = s;  k.data[6] = 0;  k.data[7] = 0;
    return k;
}

/*
 * =====================================================================================================================
 * Tests
 * =====================================================================================================================
 */

/* Identity channel: K=I → superop is 4x4 identity */
void test_superop_identity(void) {
    cplx_t data[4] = {1, 0, 0, 1};
    kraus_t k = {.n_qubits = 1, .n_terms = 1, .data = data};
    superop_t sop = kraus_to_superop(&k);

    cplx_t eye[16] = {0};
    eye[0] = 1; eye[5] = 1; eye[10] = 1; eye[15] = 1;

    TEST_ASSERT_EQUAL_COMPLEX_ARRAY_TOL(eye, sop.data, 16, PRECISION);
    free(sop.data);
}

/* X gate: K=X → S_{(a,b),(c,d)} = X_{ac} * conj(X_{bd}) */
void test_superop_x_gate(void) {
    cplx_t data[4] = {0, 1, 1, 0};  /* X in row-major */
    kraus_t k = {.n_qubits = 1, .n_terms = 1, .data = data};
    superop_t sop = kraus_to_superop(&k);

    /* X⊗X* permutes: (00)→(11), (01)→(10), (10)→(01), (11)→(00)
     * In row-major 4x4: sop[row*4+col] = 1 when (row,col) in {(0,3),(1,2),(2,1),(3,0)} */
    cplx_t ref[16] = {0};
    ref[3] = 1; ref[6] = 1; ref[9] = 1; ref[12] = 1;

    TEST_ASSERT_EQUAL_COMPLEX_ARRAY_TOL(ref, sop.data, 16, PRECISION);
    free(sop.data);
}

/* Bit-flip p=0 → identity superop */
void test_superop_bitflip_p0(void) {
    kraus_t k = make_bitflip(0.0);
    superop_t sop = kraus_to_superop(&k);

    cplx_t eye[16] = {0};
    eye[0] = 1; eye[5] = 1; eye[10] = 1; eye[15] = 1;

    TEST_ASSERT_EQUAL_COMPLEX_ARRAY_TOL(eye, sop.data, 16, PRECISION);
    free(sop.data);
    free(k.data);
}

/* Amplitude damping gamma=0 → identity superop */
void test_superop_ampdamp_g0(void) {
    kraus_t k = make_ampdamp(0.0);
    superop_t sop = kraus_to_superop(&k);

    cplx_t eye[16] = {0};
    eye[0] = 1; eye[5] = 1; eye[10] = 1; eye[15] = 1;

    TEST_ASSERT_EQUAL_COMPLEX_ARRAY_TOL(eye, sop.data, 16, PRECISION);
    free(sop.data);
    free(k.data);
}

/* Amplitude damping gamma=1: apply to vec(|1><1|) = [0,0,0,1] → vec(|0><0|) = [1,0,0,0] */
void test_superop_ampdamp_g1(void) {
    kraus_t k = make_ampdamp(1.0);
    superop_t sop = kraus_to_superop(&k);

    /* Check: sop * [0,0,0,1]^T should give [1,0,0,0]^T */
    /* That means column 3 of sop should be [1,0,0,0] (row-major: sop[row*4+3]) */
    TEST_ASSERT_COMPLEX_WITHIN(1.0, sop.data[0 * 4 + 3], PRECISION);
    TEST_ASSERT_COMPLEX_WITHIN(0.0, sop.data[1 * 4 + 3], PRECISION);
    TEST_ASSERT_COMPLEX_WITHIN(0.0, sop.data[2 * 4 + 3], PRECISION);
    TEST_ASSERT_COMPLEX_WITHIN(0.0, sop.data[3 * 4 + 3], PRECISION);

    free(sop.data);
    free(k.data);
}

/* Depolarizing p=1: E(rho) = Tr(rho) * I/2.
 * S = [[0.5, 0, 0, 0.5],
 *      [0,   0, 0, 0  ],
 *      [0,   0, 0, 0  ],
 *      [0.5, 0, 0, 0.5]] */
void test_superop_depolarizing_full(void) {
    kraus_t k = make_depolarizing(1.0);
    superop_t sop = kraus_to_superop(&k);

    cplx_t ref[16] = {0};
    ref[0]  = 0.5;  /* (0,0) */
    ref[3]  = 0.5;  /* (0,3) */
    ref[12] = 0.5;  /* (3,0) */
    ref[15] = 0.5;  /* (3,3) */

    TEST_ASSERT_EQUAL_COMPLEX_ARRAY_TOL(ref, sop.data, 16, PRECISION);
    free(sop.data);
    free(k.data);
}

/* Dephasing channel: K0=|0><0|, K1=|1><1| → zeroes off-diagonals */
void test_superop_dephasing(void) {
    cplx_t data[8] = {
        1, 0, 0, 0,  /* K0 = |0><0| row-major */
        0, 0, 0, 1   /* K1 = |1><1| row-major */
    };
    kraus_t k = {.n_qubits = 1, .n_terms = 2, .data = data};
    superop_t sop = kraus_to_superop(&k);

    /* Dephasing preserves diagonal, zeroes off-diagonal:
     * sop * [a, b, c, d]^T = [a, 0, 0, d]^T for any input.
     * So sop = diag(1, 0, 0, 1). */
    cplx_t ref[16] = {0};
    ref[0] = 1; ref[15] = 1;

    TEST_ASSERT_EQUAL_COMPLEX_ARRAY_TOL(ref, sop.data, 16, PRECISION);
    free(sop.data);
}
