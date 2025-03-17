//
// Created by Timo Ziegler on 17.03.25.
//

#ifndef TESTEXACT_H
#define TESTEXACT_H

/*
 * =====================================================================================================================
 *                                                      includes
 * =====================================================================================================================
 */
#ifndef EXACT_H
#include <exact.h>
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
 *                                                      Macros
 * =====================================================================================================================
 */
#define APPROXPRECISION         1e-4
#define EPSILON                 1e-9
#define MAXQUBITS               6
#define PRECISION               1e-8

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

double randPar[4] = {0.342377381, 1.416765874, -1.458084103, 1.767723341};
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

extern inline void evoX2(state_t* state) {
    evoQB(state, x2, randPar[0]);
}
extern inline void evoY2(state_t* state) {
    evoQB(state, y2, randPar[1]);
}
extern inline void evoZ2(state_t* state) {
    evoQB(state, z2, randPar[2]);
}
extern inline void evoSwap2(state_t* state) {
    evoQB(state, swap2, randPar[3]);
}
extern inline void evoX3(state_t* state) {
    evoQB(state, x3, randPar[0]);
}
extern inline void evoY3(state_t* state) {
    evoQB(state, y3, randPar[1]);
}
extern inline void evoZ3(state_t* state) {
    evoQB(state, z3, randPar[2]);
}
extern inline void evoSwap3(state_t* state) {
    evoQB(state, swap3, randPar[3]);
}
extern inline void evoX4(state_t* state) {
    evoQB(state, x4, randPar[0]);
}
extern inline void evoY4(state_t* state) {
    evoQB(state, y4, randPar[1]);
}
extern inline void evoZ4(state_t* state) {
    evoQB(state, z4, randPar[2]);
}
extern inline void evoSwap4(state_t* state) {
    evoQB(state, swap4, randPar[3]);
}
extern inline void evoX5(state_t* state) {
    evoQB(state, x5, randPar[0]);
}
extern inline void evoY5(state_t* state) {
    evoQB(state, y5, randPar[1]);
}
extern inline void evoZ5(state_t* state) {
    evoQB(state, z5, randPar[2]);
}
extern inline void evoSwap5(state_t* state) {
    evoQB(state, swap5, randPar[3]);
}
extern inline void evoX6(state_t* state) {
    evoQB(state, x6, randPar[0]);
}
extern inline void evoY6(state_t* state) {
    evoQB(state, y6, randPar[1]);
}
extern inline void evoZ6(state_t* state) {
    evoQB(state, z6, randPar[2]);
}
extern inline void evoSwap6(state_t* state) {
    evoQB(state, swap6, randPar[3]);
}
applyQB pqc[20] = {evoX2, evoY2, evoZ2, evoSwap2, evoX3, evoY3, evoZ3, evoSwap3, evoX4, evoY4, evoZ4, evoSwap4, evoX5,
    evoY5, evoZ5, evoSwap5, evoX6, evoY6, evoZ6, evoSwap6};

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

#endif //TESTEXACT_H
