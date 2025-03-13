//
// Created by Timo Ziegler on 13.03.25.
//
/*
 * =====================================================================================================================
 *                                                      includes
 * =====================================================================================================================
 */
#ifndef PQC_H
#include "pqc.h"
#endif

#ifndef QHIPSTER_H
#include "qhipster.h"
#endif

#ifndef GATEMAT_H
#include "gateMat.h"
#endif

#ifndef UNITY_FRAMEWORK_H
#include "unity.h"
#endif

#ifndef UTILS_H
#include "utils.h"
#endif

/*
 * =====================================================================================================================
 *                                                      macros
 * =====================================================================================================================
 */
#define PRECISION               1e-8
#define MAXQUBITS               6
#include <exact.h>

/*
 * =====================================================================================================================
 *                                                  setUp and tearDown
 * =====================================================================================================================
 */
void setUp(void) {}
void tearDown(void) {}

/*
 * =====================================================================================================================
 *                                                  testMeanObs
 * =====================================================================================================================
 */
double diagObs[64] = {-0.672013803, -1.080744850, 1.989524770, 1.988930799, 1.731153443, 0.219447958, -1.496797918,
    0.153447271, 1.394199856, 1.434125830, -0.525672657, 0.082399361, -0.256029741, 0.729175258, 1.631006491,
    2.747422307, -1.333258664, 0.147877182, 0.678307326, -2.064295279, -2.955271079, 0.834284985, -0.036103279,
    -0.174833204, 2.422058764, -2.213238216, 1.295379069, -2.132723608, 2.890693833, -2.400779788, -2.791631103,
    -0.748354991, -1.296158621, 1.286226791, -1.301539827, 3.019915018, 1.863008808, 0.373589359, -1.390747670,
    -2.016940260, 0.848049481, 1.515307618, -2.638919196, 1.680740517, -2.302014922, 0.864815536, -2.727445042,
    1.404416495, 1.087979959, -0.563509955, -2.881082581, 2.242557345, 1.721171693, 2.040283920, 1.542535211,
    0.520890342, 0.069996072, -3.053401020, 0.427764464, 1.344208689, -2.416954506, -2.499271157, 2.576468900,
    0.289965927
    };

void testMeanObs(void) {
    state_t testState;
    for (qubit_t qubits = 1; qubits <= MAXQUBITS; qubits++) {
        const dim_t dim = POW2(qubits, dim_t);
        stateInitEmpty(&testState, qubits);
        cplx_t** vecs = generateTestVectors(qubits);

        for (uint8_t i = 0; i < dim + 1; i++) {
            stateInitVector(&testState, vecs[i]);
            double test = meanObs(&testState, diagObs);

            double ref = 0;
            for (dim_t j = 0; j < dim; j++) {
                ref += pow(cabs(vecs[i][j]), 2) * diagObs[j];
            }

            TEST_ASSERT_TRUE(cvectorAlmostEqual(vecs[i], testState.vec, dim, PRECISION));
        }
        freeTestVectors(vecs, qubits);
        stateFreeVector(&testState);
    }
}

/*
 * =====================================================================================================================
 *                                                 testGradientPQC
 * =====================================================================================================================
 */
double randPar[4] = {0.342377381, 1.416765874, -1.458084103, 1.767723341};

void rx2(state_t* state, const double par) {
    applyRX(state, 0, par);
    applyRX(state, 1, par);
}
void ry2(state_t* state, const double par) {
    applyRY(state, 0, par);
    applyRY(state, 1, par);
}
void rz2(state_t* state, const double par) {
    applyRZ(state, 0, par);
    applyRZ(state, 1, par);
}
void rswap2(state_t* state, const double par) {
    applyRswap(state, 0, 1, par);
}
void rx3(state_t* state, const double par) {
    applyRX(state, 0, par);
    applyRX(state, 1, par);
    applyRX(state, 2, par);
}
void ry3(state_t* state, const double par) {
    applyRY(state, 0, par);
    applyRY(state, 1, par);
    applyRY(state, 2, par);
}
void rz3(state_t* state, const double par) {
    applyRZ(state, 0, par);
    applyRZ(state, 1, par);
    applyRZ(state, 2, par);
}
void rswap3(state_t* state, const double par) {
    applyRswap(state, 0, 1, par);
}
void rx4(state_t* state, const double par) {
    applyRX(state, 0, par);
    applyRX(state, 1, par);
    applyRX(state, 2, par);
    applyRX(state, 3, par);
}
void ry4(state_t* state, const double par) {
    applyRY(state, 0, par);
    applyRY(state, 1, par);
    applyRY(state, 2, par);
    applyRY(state, 3, par);
}
void rz4(state_t* state, const double par) {
    applyRZ(state, 0, par);
    applyRZ(state, 1, par);
    applyRZ(state, 2, par);
    applyRZ(state, 3, par);
}
void rswap4(state_t* state, const double par) {
    applyRswap(state, 0, 1, par);
    applyRswap(state, 2, 3, par);
}
void rx5(state_t* state, const double par) {
    applyRX(state, 0, par);
    applyRX(state, 1, par);
    applyRX(state, 2, par);
    applyRX(state, 3, par);
    applyRX(state, 4, par);
}
void ry5(state_t* state, const double par) {
    applyRY(state, 0, par);
    applyRY(state, 1, par);
    applyRY(state, 2, par);
    applyRY(state, 3, par);
    applyRY(state, 4, par);
}
void rz5(state_t* state, const double par) {
    applyRZ(state, 0, par);
    applyRZ(state, 1, par);
    applyRZ(state, 2, par);
    applyRZ(state, 3, par);
    applyRZ(state, 4, par);
}
void rswap5(state_t* state, const double par) {
    applyRswap(state, 0, 1, par);
    applyRswap(state, 2, 3, par);
}
void rx6(state_t* state, const double par) {
    applyRX(state, 0, par);
    applyRX(state, 1, par);
    applyRX(state, 2, par);
    applyRX(state, 3, par);
    applyRX(state, 4, par);
    applyRX(state, 5, par);
}
void ry6(state_t* state, const double par) {
    applyRY(state, 0, par);
    applyRY(state, 1, par);
    applyRY(state, 2, par);
    applyRY(state, 3, par);
    applyRY(state, 4, par);
    applyRY(state, 5, par);
}
void rz6(state_t* state, const double par) {
    applyRZ(state, 0, par);
    applyRZ(state, 1, par);
    applyRZ(state, 2, par);
    applyRZ(state, 3, par);
    applyRZ(state, 4, par);
    applyRZ(state, 5, par);
}
void rswap6(state_t* state, const double par) {
    applyRswap(state, 0, 1, par);
    applyRswap(state, 2, 3, par);
    applyRswap(state, 4, 5, par);
}

applyPQB pqbs[20] = {rx2, ry2, rz2, rswap2, rx3, ry3, rz3, rswap3, rx4, ry4, rz4, rswap4, rx5, ry5, rz5, rswap5, \
    rx6, ry6, rz6, rswap6};

cplx_t* (*rMat[])(qubit_t, qubit_t, double) = {RXGateMat, RYGateMat, RZGateMat};

/*
 * =====================================================================================================================
 *                                                      main
 * =====================================================================================================================
 */
int main(void) {
    UNITY_BEGIN();
    RUN_TEST(testMeanObs);
    return UNITY_END();
}