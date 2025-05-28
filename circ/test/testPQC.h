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
#define MAXQUBITS               7

/*
 * =====================================================================================================================
 *                                                  Type definitions
 * =====================================================================================================================
 */
/*
 * Function:  Matrix representation of a quantum block
 * ------------------------------------
 * This function returns the matrix representation of some quantum block
 */
typedef cplx_t* (*matQB)(qubit_t qubits);

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

double dcoeff[5] = {5.101260256, -3.686076046, 0.654143819, -7.576553655, -0.387642064};

cplx_t zcoeff[5] = {
        -7.930739979 * I * 8.583870403, -4.788748126 - I * 0.052266189, -1.505055930 + I * 1.203127483,
        -6.350354815 - I * 6.586643177, -8.775751469 - I * 6.272504304
};


double diagObs[64] = {
        -0.672013803, -1.080744850, 1.989524770, 1.988930799, 1.731153443, 0.219447958, -1.496797918, 0.153447271,
        1.394199856, 1.434125830, -0.525672657, 0.082399361, -0.256029741, 0.729175258, 1.631006491, 2.747422307,
        -1.333258664, 0.147877182, 0.678307326, -2.064295279, -2.955271079, 0.834284985, -0.036103279, -0.174833204,
        2.422058764, -2.213238216, 1.295379069, -2.132723608, 2.890693833, -2.400779788, -2.791631103, -0.748354991,
        -1.296158621, 1.286226791, -1.301539827, 3.019915018, 1.863008808, 0.373589359, -1.390747670, -2.016940260,
        0.848049481, 1.515307618, -2.638919196, 1.680740517, -2.302014922, 0.864815536, -2.727445042, 1.404416495,
        1.087979959, -0.563509955, -2.881082581, 2.242557345, 1.721171693, 2.040283920, 1.542535211, 0.520890342,
        0.069996072, -3.053401020, 0.427764464, 1.344208689, -2.416954506, -2.499271157, 2.576468900, 0.289965927
};

/*
 * =====================================================================================================================
 *                                                  Quantum blocks
 * =====================================================================================================================
 */
// 2 qubits
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

// 3 qubits
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

// 4 qubits
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

// 5 qubits
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

// 6 qubits
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

// Diagonal observable
extern inline void diag(state_t* state) {
    applyDiag(state, diagObs);
}

applyQB qb[25] = {x2, y2, z2, swap2, diag, x3, y3, z3, swap3, diag, x4, y4, z4, swap4, diag, x5, y5, z5, swap5, diag,
    x6, y6, z6, swap6, diag};

/*
 * =====================================================================================================================
 *                                          Reference: Quantum blocks
 * =====================================================================================================================
 */
extern inline cplx_t* xMat(qubit_t qubits) {
    dim_t dim = 1 << qubits;
    cplx_t* out = singleQubitGateMat(XMAT, qubits, qubits - 1);
    for (qubit_t i = qubits - 1; i > 0; --i) {
        cplx_t* tmp = singleQubitGateMat(XMAT, qubits, i - 1);
        cmatMulInPlace(out, tmp, dim);
        free(tmp);
    }
    return out;
}
extern inline cplx_t* yMat(qubit_t qubits) {
    dim_t dim = 1 << qubits;
    cplx_t* out = singleQubitGateMat(YMAT, qubits, qubits - 1);
    for (qubit_t i = qubits - 1; i > 0; --i) {
        cplx_t* tmp = singleQubitGateMat(YMAT, qubits, i - 1);
        cmatMulInPlace(out, tmp, dim);
        free(tmp);
    }
    return out;
}
extern inline cplx_t* zMat(qubit_t qubits) {
    dim_t dim = 1 << qubits;
    cplx_t* out = singleQubitGateMat(ZMAT, qubits, qubits - 1);
    for (qubit_t i = qubits - 1; i > 0; --i) {
        cplx_t* tmp = singleQubitGateMat(ZMAT, qubits, i - 1);
        cmatMulInPlace(out, tmp, dim);
        free(tmp);
    }
    return out;
}
extern inline cplx_t* swapMat(qubit_t qubits) {
    dim_t dim = 1 << qubits;
    qubit_t start = qubits / 2;                                     // start = k for qubits = 2k and qubits = (2k + 1)
    cplx_t* out = swapGateMat(qubits, 2 * start - 2, 2 * start - 1);    // SWAP the last qubit with odd index to left
    for (qubit_t i = start - 1; i > 0; --i) {                       // Starting at the next qubit with odd index
        cplx_t* tmp = swapGateMat(qubits, 2 * i - 2, 2 * i - 1);
        cmatMulInPlace(out, tmp, dim);
        free(tmp);
    }
    return out;
}

extern inline cplx_t* diagMat(qubit_t qubits) {
    return diagGateMat(qubits, diagObs);
}

matQB qbMat[25] = {xMat, yMat, zMat, swapMat, diagMat};

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
