//
// Created by Timo Ziegler on 16.03.25.
//

#ifndef TESTPQC_H
#define TESTPQC_H

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

/*
 * =====================================================================================================================
 *                                                  setUp and tearDown
 * =====================================================================================================================
 */
extern inline void setUp(void) {}
extern inline void tearDown(void) {}

/*
 * =====================================================================================================================
 *                                                  Definitions
 * =====================================================================================================================
 */
double randPar[4] = {0.342377381, 1.416765874, -1.458084103, 1.767723341};
cplx_t coeff[4] = {0.551736480 + I * 0.461109576, -2.734804444 + I * 0.567079912, 0.362535970 + I * 0.714549265,
    0.372739852 - I * -1.612172457};
const cplx_t* pauliMat[] = {XMAT, YMAT, ZMAT};
cplx_t* (*rMat[])(qubit_t, qubit_t, double) = {RXGateMat, RYGateMat, RZGateMat};

/*
 * =====================================================================================================================
 *                                                  Quantum blocks
 * =====================================================================================================================
 */
extern inline void x2(state_t* state) {
    applyX(state, 0);
    applyX(state, 1);
}
extern inline void y2(state_t* state) {
    applyY(state, 0);
    applyY(state, 1);
}
extern inline void z2(state_t* state) {
    applyZ(state, 0);
    applyZ(state, 1);

}
extern inline void swap2(state_t* state) {
    applySWAP(state, 0, 1);
}
extern inline void x3(state_t* state) {
    applyX(state, 0);
    applyX(state, 1);
    applyX(state, 2);
}
extern inline void y3(state_t* state) {
    applyY(state, 0);
    applyY(state, 1);
    applyY(state, 2);
}
extern inline void z3(state_t* state) {
    applyZ(state, 0);
    applyZ(state, 1);
    applyZ(state, 2);
}
extern inline void swap3(state_t* state) {
    applySWAP(state, 0, 1);
}
extern inline void x4(state_t* state) {
    applyX(state, 0);
    applyX(state, 1);
    applyX(state, 2);
    applyX(state, 3);
}
extern inline void y4(state_t* state) {
    applyY(state, 0);
    applyY(state, 1);
    applyY(state, 2);
    applyY(state, 3);
}
extern inline void z4(state_t* state) {
    applyZ(state, 0);
    applyZ(state, 1);
    applyZ(state, 2);
    applyZ(state, 3);
}
extern inline void swap4(state_t* state) {
    applySWAP(state, 0, 1);
    applySWAP(state, 2, 3);
}
extern inline void x5(state_t* state) {
    applyX(state, 0);
    applyX(state, 1);
    applyX(state, 2);
    applyX(state, 3);
    applyX(state, 4);
}
extern inline void y5(state_t* state) {
    applyY(state, 0);
    applyY(state, 1);
    applyY(state, 2);
    applyY(state, 3);
    applyY(state, 4);
}
extern inline void z5(state_t* state) {
    applyZ(state, 0);
    applyZ(state, 1);
    applyZ(state, 2);
    applyZ(state, 3);
    applyZ(state, 4);
}
extern inline void swap5(state_t* state) {
    applySWAP(state, 0, 1);
    applySWAP(state, 2, 3);
}
extern inline void x6(state_t* state) {
    applyX(state, 0);
    applyX(state, 1);
    applyX(state, 2);
    applyX(state, 3);
    applyX(state, 4);
    applyX(state, 5);
}
extern inline void y6(state_t* state) {
    applyY(state, 0);
    applyY(state, 1);
    applyY(state, 2);
    applyY(state, 3);
    applyY(state, 4);
    applyY(state, 5);
}
extern inline void z6(state_t* state) {
    applyZ(state, 0);
    applyZ(state, 1);
    applyZ(state, 2);
    applyZ(state, 3);
    applyZ(state, 4);
    applyZ(state, 5);
}
extern inline void swap6(state_t* state) {
    applySWAP(state, 0, 1);
    applySWAP(state, 2, 3);
    applySWAP(state, 4, 5);
}
applyQB qbs[20] = {x2, y2, z2, swap2, x3, y3, z3, swap3, x4, y4, z4, swap4, x5, y5, z5, swap5, x6, y6, z6, swap6};

// void rx2(state_t* state, const double par) {
//     applyRX(state, 0, par);
//     applyRX(state, 1, par);
// }
// void ry2(state_t* state, const double par) {
//     applyRY(state, 0, par);
//     applyRY(state, 1, par);
// }
// void rz2(state_t* state, const double par) {
//     applyRZ(state, 0, par);
//     applyRZ(state, 1, par);
// }
// void rswap2(state_t* state, const double par) {
//     applyRswap(state, 0, 1, par);
// }
// void rx3(state_t* state, const double par) {
//     applyRX(state, 0, par);
//     applyRX(state, 1, par);
//     applyRX(state, 2, par);
// }
// void ry3(state_t* state, const double par) {
//     applyRY(state, 0, par);
//     applyRY(state, 1, par);
//     applyRY(state, 2, par);
// }
// void rz3(state_t* state, const double par) {
//     applyRZ(state, 0, par);
//     applyRZ(state, 1, par);
//     applyRZ(state, 2, par);
// }
// void rswap3(state_t* state, const double par) {
//     applyRswap(state, 0, 1, par);
// }
// void rx4(state_t* state, const double par) {
//     applyRX(state, 0, par);
//     applyRX(state, 1, par);
//     applyRX(state, 2, par);
//     applyRX(state, 3, par);
// }
// void ry4(state_t* state, const double par) {
//     applyRY(state, 0, par);
//     applyRY(state, 1, par);
//     applyRY(state, 2, par);
//     applyRY(state, 3, par);
// }
// void rz4(state_t* state, const double par) {
//     applyRZ(state, 0, par);
//     applyRZ(state, 1, par);
//     applyRZ(state, 2, par);
//     applyRZ(state, 3, par);
// }
// void rswap4(state_t* state, const double par) {
//     applyRswap(state, 0, 1, par);
//     applyRswap(state, 2, 3, par);
// }
// void rx5(state_t* state, const double par) {
//     applyRX(state, 0, par);
//     applyRX(state, 1, par);
//     applyRX(state, 2, par);
//     applyRX(state, 3, par);
//     applyRX(state, 4, par);
// }
// void ry5(state_t* state, const double par) {
//     applyRY(state, 0, par);
//     applyRY(state, 1, par);
//     applyRY(state, 2, par);
//     applyRY(state, 3, par);
//     applyRY(state, 4, par);
// }
// void rz5(state_t* state, const double par) {
//     applyRZ(state, 0, par);
//     applyRZ(state, 1, par);
//     applyRZ(state, 2, par);
//     applyRZ(state, 3, par);
//     applyRZ(state, 4, par);
// }
// void rswap5(state_t* state, const double par) {
//     applyRswap(state, 0, 1, par);
//     applyRswap(state, 2, 3, par);
// }
// void rx6(state_t* state, const double par) {
//     applyRX(state, 0, par);
//     applyRX(state, 1, par);
//     applyRX(state, 2, par);
//     applyRX(state, 3, par);
//     applyRX(state, 4, par);
//     applyRX(state, 5, par);
// }
// void ry6(state_t* state, const double par) {
//     applyRY(state, 0, par);
//     applyRY(state, 1, par);
//     applyRY(state, 2, par);
//     applyRY(state, 3, par);
//     applyRY(state, 4, par);
//     applyRY(state, 5, par);
// }
// void rz6(state_t* state, const double par) {
//     applyRZ(state, 0, par);
//     applyRZ(state, 1, par);
//     applyRZ(state, 2, par);
//     applyRZ(state, 3, par);
//     applyRZ(state, 4, par);
//     applyRZ(state, 5, par);
// }
// void rswap6(state_t* state, const double par) {
//     applyRswap(state, 0, 1, par);
//     applyRswap(state, 2, 3, par);
//     applyRswap(state, 4, 5, par);
// }
//
// applyPQB pqbs[20] = {rx2, ry2, rz2, rswap2, rx3, ry3, rz3, rswap3, rx4, ry4, rz4, rswap4, rx5, ry5, rz5, rswap5, \
//     rx6, ry6, rz6, rswap6};

#endif //TESTPQC_H
