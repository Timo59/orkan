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
#define MAXQUBITS               7
#define PRECISION               1e-8

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
 * Function:  Matrix representation of a parametrized quantum block
 * ------------------------------------
 * This function returns the matrix representation of some quantum block
 */
typedef cplx_t* (*matPQB)(qubit_t qubits, double par);

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

double randPar[5] = {-1.935113472, 0.170656438, -1.916424226, -0.785971871, 1.746181027};
cplx_t zcoeff[15] = {9.20389 - I * 2.07242, 0.58007 + I * 0.58377, -9.98112 - I * 5.35024, -6.38089 + I * 9.44588,
                    4.31933 - I * 7.73341, 8.49492 + I * 6.53551, -8.63475 + I * 1.27771, 2.56232 + I * 6.10205,
                    3.77543 - I * 7.55000, 3.87864 + I * 3.90134, 0.09866 - I * 6.01226, -6.48157 + I * 2.52767,
                    -1.57778 + I * 7.94151, -2.33955 - I * 4.21345, -4.49085 + I * 6.10035};



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
                  x6, y6, z6, swap6, diag
};

/*
 * =====================================================================================================================
 *                                          Quantum block matrices
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

matQB qbMat[5] = {xMat, yMat, zMat, swapMat, diagMat};

/*
 * =====================================================================================================================
 *                                          Parametrized quantum blocks
 * =====================================================================================================================
 */
// 2 qubits
extern inline void evoX2(state_t* state, const double par) {
    evoQB(state, x2, par);
}
extern inline void rx2(state_t* state, const double par) {
    applyRX(state, 0, par);
    applyRX(state, 1, par);
}

extern inline void evoY2(state_t* state, const double par) {
    evoQB(state, y2, par);
}
extern inline void ry2(state_t* state, const double par) {
    applyRY(state, 0, par);
    applyRY(state, 1, par);
}

extern inline void evoZ2(state_t* state, const double par) {
    evoQB(state, z2, par);
}
extern inline void rz2(state_t* state, const double par) {
    applyRZ(state, 0, par);
    applyRZ(state, 1, par);
}

extern inline void evoSwap2(state_t* state, const double par) {
    evoQB(state, swap2, par);
}
extern inline void rswap2(state_t* state, const double par) {
    applyRSWAP(state, 0, 1, par);
}

extern inline void evoDiag2(state_t* state, const double par) {
    evoDiag(state, diagObs, par);
}

// 3 qubits
extern inline void evoX3(state_t* state, const double par) {
    evoQB(state, x3, par);
}
extern inline void rx3(state_t* state, const double par) {
    applyRX(state, 0, par);
    applyRX(state, 1, par);
    applyRX(state, 2, par);
}

extern inline void evoY3(state_t* state, const double par) {
    evoQB(state, y3, par);
}
extern inline void ry3(state_t* state, const double par) {
    applyRY(state, 0, par);
    applyRY(state, 1, par);
    applyRY(state, 2, par);
}

extern inline void evoZ3(state_t* state, const double par) {
    evoQB(state, z3, par);
}
extern inline void rz3(state_t* state, const double par) {
    applyRZ(state, 0, par);
    applyRZ(state, 1, par);
    applyRZ(state, 2, par);
}

extern inline void evoSwap3(state_t* state, const double par) {
    evoQB(state, swap3, par);
}
extern inline void rswap3(state_t* state, const double par) {
    applyRSWAP(state, 0, 1, par);
}

extern inline void evoX4(state_t* state, const double par) {
    evoQB(state, x4, par);
}
extern inline void rx4(state_t* state, const double par) {
    applyRX(state, 0, par);
    applyRX(state, 1, par);
    applyRX(state, 2, par);
    applyRX(state, 3, par);
}

extern inline void evoY4(state_t* state, const double par) {
    evoQB(state, y4, par);
}
extern inline void ry4(state_t* state, const double par) {
    applyRY(state, 0, par);
    applyRY(state, 1, par);
    applyRY(state, 2, par);
    applyRY(state, 3, par);
}

extern inline void evoZ4(state_t* state, const double par) {
    evoQB(state, z4, par);
}
extern inline void rz4(state_t* state, const double par) {
    applyRZ(state, 0, par);
    applyRZ(state, 1, par);
    applyRZ(state, 2, par);
    applyRZ(state, 3, par);
}

extern inline void evoSwap4(state_t* state, const double par) {
    evoQB(state, swap4, par);
}
extern inline void rswap4(state_t* state, const double par) {
    applyRSWAP(state, 0, 1, par);
    applyRSWAP(state, 2, 3, par);
}

extern inline void evoX5(state_t* state, const double par) {
    evoQB(state, x5, par);
}
extern inline void rx5(state_t* state, const double par) {
    applyRX(state, 0, par);
    applyRX(state, 1, par);
    applyRX(state, 2, par);
    applyRX(state, 3, par);
    applyRX(state, 4, par);
}

extern inline void evoY5(state_t* state, const double par) {
    evoQB(state, y5, par);
}
extern inline void ry5(state_t* state, const double par) {
    applyRY(state, 0, par);
    applyRY(state, 1, par);
    applyRY(state, 2, par);
    applyRY(state, 3, par);
    applyRY(state, 4, par);
}

extern inline void evoZ5(state_t* state, const double par) {
    evoQB(state, z5, par);
}
extern inline void rz5(state_t* state, const double par) {
    applyRZ(state, 0, par);
    applyRZ(state, 1, par);
    applyRZ(state, 2, par);
    applyRZ(state, 3, par);
    applyRZ(state, 4, par);
}

extern inline void evoSwap5(state_t* state, const double par) {
    evoQB(state, swap5, par);
}
extern inline void rswap5(state_t* state, const double par) {
    applyRSWAP(state, 0, 1, par);
    applyRSWAP(state, 2, 3, par);
}

extern inline void evoX6(state_t* state, const double par) {
    evoQB(state, x6, par);
}
extern inline void rx6(state_t* state, const double par) {
    applyRX(state, 0, par);
    applyRX(state, 1, par);
    applyRX(state, 2, par);
    applyRX(state, 3, par);
    applyRX(state, 4, par);
    applyRX(state, 5, par);
}

extern inline void evoY6(state_t* state, const double par) {
    evoQB(state, y6, par);
}
extern inline void ry6(state_t* state, const double par) {
    applyRY(state, 0, par);
    applyRY(state, 1, par);
    applyRY(state, 2, par);
    applyRY(state, 3, par);
    applyRY(state, 4, par);
    applyRY(state, 5, par);
}

extern inline void evoZ6(state_t* state, const double par) {
    evoQB(state, z6, par);
}
extern inline void rz6(state_t* state, const double par) {
    applyRZ(state, 0, par);
    applyRZ(state, 1, par);
    applyRZ(state, 2, par);
    applyRZ(state, 3, par);
    applyRZ(state, 4, par);
    applyRZ(state, 5, par);
}

extern inline void evoSwap6(state_t* state, const double par) {
    evoQB(state, swap6, par);
}
extern inline void rswap6(state_t* state, const double par) {
    applyRSWAP(state, 0, 1, par);
    applyRSWAP(state, 2, 3, par);
    applyRSWAP(state, 4, 5, par);
}

applyPQB pqc[25] = {evoX2, evoY2, evoZ2, evoSwap2, evoDiag2, evoX3, evoY3, evoZ3, evoSwap3, evoDiag2, evoX4, evoY4,
                    evoZ4, evoSwap4, evoDiag2, evoX5, evoY5, evoZ5, evoSwap5, evoDiag2, evoX6, evoY6, evoZ6, evoSwap6, evoDiag2};

applyPQB rpqc[20] = {rx2, ry2, rz2, rswap2, rx3, ry3, rz3, rswap3, rx4, ry4, rz4, rswap4, rx5, ry5, rz5, rswap5, \
     rx6, ry6, rz6, rswap6};

/*
 * =====================================================================================================================
 *                                      Parametrized quantum block matrices
 * =====================================================================================================================
 */
extern inline cplx_t* evoXMat(qubit_t qubits, double par) {
    dim_t dim = 1 << qubits;
    cplx_t* gateMat = xMat(qubits);
    cplx_t* out = zexpm(gateMat, par, dim);
    free(gateMat);

    return out;
}
extern inline cplx_t* rxMat(qubit_t qubits, double par) {
    dim_t dim = 1 << qubits;
    cplx_t* out = RXGateMat(qubits, qubits - 1, par);
    for (qubit_t i = qubits - 1; i > 0; --i) {
        cplx_t* tmp = RXGateMat(qubits, i - 1, par);
        cmatMulInPlace(out, tmp, dim);
        free(tmp);
    }
    return out;
}

extern inline cplx_t* evoYMat(qubit_t qubits, double par) {
    dim_t dim = 1 << qubits;
    cplx_t* gateMat = yMat(qubits);
    cplx_t* out = zexpm(gateMat, par, dim);
    free(gateMat);

    return out;
}
extern inline cplx_t* ryMat(qubit_t qubits, double par) {
    dim_t dim = 1 << qubits;
    cplx_t* out = RYGateMat(qubits, qubits - 1, par);
    for (qubit_t i = qubits - 1; i > 0; --i) {
        cplx_t* tmp = RYGateMat(qubits, i - 1, par);
        cmatMulInPlace(out, tmp, dim);
        free(tmp);
    }
    return out;
}

extern inline cplx_t* evoZMat(qubit_t qubits, double par) {
    dim_t dim = 1 << qubits;
    cplx_t* gateMat = zMat(qubits);
    cplx_t* out = zexpm(gateMat, par, dim);
    free(gateMat);

    return out;
}
extern inline cplx_t* rzMat(qubit_t qubits, double par) {
    dim_t dim = 1 << qubits;
    cplx_t* out = RZGateMat(qubits, qubits - 1, par);
    for (qubit_t i = qubits - 1; i > 0; --i) {
        cplx_t* tmp = RZGateMat(qubits, i - 1, par);
        cmatMulInPlace(out, tmp, dim);
        free(tmp);
    }
    return out;
}

extern inline cplx_t* evoSwapMat(qubit_t qubits, double par) {
    dim_t dim = 1 << qubits;
    cplx_t* gateMat = swapMat(qubits);
    cplx_t* out = zexpm(gateMat, par, dim);
    free(gateMat);

    return out;
}
extern inline cplx_t* rSwapMat(qubit_t qubits, double par) {
    dim_t dim = 1 << qubits;
    qubit_t start = qubits / 2;                                     // start = k for qubits = 2k and qubits = (2k + 1)
    cplx_t* out = RswapGateMat(qubits, 2 * start - 2, 2 * start - 1, par);
    for (qubit_t i = start - 1; i > 0; --i) {                       // Starting at the next qubit with odd index
        cplx_t* tmp = RswapGateMat(qubits, 2 * i - 2, 2 * i - 1, par);
        cmatMulInPlace(out, tmp, dim);
        free(tmp);
    }
    return out;
}

extern inline cplx_t* evoDiagMat(qubit_t qubits, double par) {
    dim_t dim = 1 << qubits;
    cplx_t* gateMat = diagMat(qubits);
    cplx_t* out = zexpm(gateMat, par, dim);
    free(gateMat);

    return out;
}

matPQB evoMat[5] = {evoXMat, evoYMat, evoZMat, evoSwapMat, evoDiagMat};

/*
 * =====================================================================================================================
 *                                      Parametrized quantum blocks with fixed parameters
 * =====================================================================================================================
 */
extern inline void UX2(state_t* state) {
    evoX2(state, randPar[0]);
}
extern inline void URX2(state_t* state) {
    rx2(state, randPar[0]);
}

extern inline void UY2(state_t* state) {
    evoY2(state, randPar[1]);
}
extern inline void URY2(state_t* state) {
    ry2(state, randPar[1]);
}

extern inline void UZ2(state_t* state) {
    evoZ2(state, randPar[2]);
}
extern inline void URZ2(state_t* state) {
    rz2(state, randPar[2]);
}

extern inline void USwap2(state_t* state) {
    evoSwap2(state, randPar[3]);
}
extern inline void URSwap2(state_t* state) {
    rswap2(state, randPar[3]);
}

extern inline void UDiag(state_t* state) {
    evoDiag2(state, randPar[4]);
}

extern inline void UX3(state_t* state) {
    evoX3(state, randPar[0]);
}
extern inline void URX3(state_t* state) {
    rx3(state, randPar[0]);
}

extern inline void UY3(state_t* state) {
    evoY3(state, randPar[1]);
}
extern inline void URY3(state_t* state) {
    ry3(state, randPar[1]);
}

extern inline void UZ3(state_t* state) {
    evoZ3(state, randPar[2]);
}
extern inline void URZ3(state_t* state) {
    rz3(state, randPar[2]);
}

extern inline void USwap3(state_t* state) {
    evoSwap2(state, randPar[3]);
}
extern inline void URSwap3(state_t* state) {
    rswap2(state, randPar[3]);
}

extern inline void UX4(state_t* state) {
    evoX4(state, randPar[0]);
}
extern inline void URX4(state_t* state) {
    rx4(state, randPar[0]);
}

extern inline void UY4(state_t* state) {
    evoY4(state, randPar[1]);
}
extern inline void URY4(state_t* state) {
    ry4(state, randPar[1]);
}

extern inline void UZ4(state_t* state) {
    evoZ4(state, randPar[2]);
}
extern inline void URZ4(state_t* state) {
    rz4(state, randPar[2]);
}

extern inline void USwap4(state_t* state) {
    evoSwap4(state, randPar[3]);
}
extern inline void URSwap4(state_t* state) {
    rswap4(state, randPar[3]);
}

extern inline void UX5(state_t* state) {
    evoX5(state, randPar[0]);
}
extern inline void URX5(state_t* state) {
    rx5(state, randPar[0]);
}

extern inline void UY5(state_t* state) {
    evoY5(state, randPar[1]);
}
extern inline void URY5(state_t* state) {
    ry5(state, randPar[1]);
}

extern inline void UZ5(state_t* state) {
    evoZ5(state, randPar[2]);
}
extern inline void URZ5(state_t* state) {
    rz5(state, randPar[2]);
}

extern inline void USwap5(state_t* state) {
    evoSwap5(state, randPar[3]);
}
extern inline void URSwap5(state_t* state) {
    rswap5(state, randPar[3]);
}

extern inline void UX6(state_t* state) {
    evoX6(state, randPar[0]);
}
extern inline void URX6(state_t* state) {
    rx6(state, randPar[0]);
}

extern inline void UY6(state_t* state) {
    evoY6(state, randPar[1]);
}
extern inline void URY6(state_t* state) {
    ry6(state, randPar[1]);
}

extern inline void UZ6(state_t* state) {
    evoZ6(state, randPar[2]);
}
extern inline void URZ6(state_t* state) {
    rz6(state, randPar[2]);
}

extern inline void USwap6(state_t* state) {
    evoSwap6(state, randPar[3]);
}
extern inline void URSwap6(state_t* state) {
    rswap6(state, randPar[3]);
}

applyQB obs[15] = {diag, y2, x2, diag, y3, x3, diag, y4, x4, diag, y5, x5, diag, y6, x6};

applyQB channel[75] = {x2, y2, z2, swap2, diag, UX2, UY2, UZ2, USwap2, UDiag, URX2, URY2, URZ2, URSwap2, UDiag,
                       x3, y3, z3, swap3, diag, UX3, UY3, UZ3, USwap3, UDiag, URX3, URY3, URZ3, URSwap3, UDiag,
                       x4, y4, z4, swap4, diag, UX4, UY4, UZ4, USwap4, UDiag, URX4, URY4, URZ4, URSwap4, UDiag,
                       x5, y5, z5, swap5, diag, UX5, UY5, UZ5, USwap5, UDiag, URX5, URY5, URZ5, URSwap5, UDiag,
                       x6, y6, z6, swap6, diag, UX6, UY6, UZ6, USwap6, UDiag, URX6, URY6, URZ6, URSwap6, UDiag
};
/*
 * =====================================================================================================================
 *                          Parametrized quantum blocks with fixed parameters as matrices
 * =====================================================================================================================
 */
extern inline cplx_t* UXMat(qubit_t qubits) {
    return evoXMat(qubits, randPar[0]);
}
extern inline cplx_t* URXMat(qubit_t qubits) {
    return rxMat(qubits, randPar[0]);
}

extern inline cplx_t* UYMat(qubit_t qubits) {
    return evoYMat(qubits, randPar[1]);
}
extern inline cplx_t* URYMat(qubit_t qubits) {
    return ryMat(qubits, randPar[1]);
}

extern inline cplx_t* UZMat(qubit_t qubits) {
    return evoZMat(qubits, randPar[2]);
}
extern inline cplx_t* URZMat(qubit_t qubits) {
    return rzMat(qubits, randPar[2]);
}

extern inline cplx_t* USwapMat(qubit_t qubits) {
    return evoSwapMat(qubits, randPar[3]);
}
extern inline cplx_t* URSwapMat(qubit_t qubits) {
    return rSwapMat(qubits, randPar[3]);
}

extern inline cplx_t* UDiagMat(qubit_t qubits) {
    return evoDiagMat(qubits, randPar[4]);
}

matQB obsMat[3] = {diagMat, yMat, xMat};

matQB channelMat[3][5] = {{xMat, yMat, zMat, swapMat, diagMat}, {UXMat, UYMat, UZMat, USwapMat, UDiagMat},
                        {URXMat, URYMat, URZMat, URSwapMat, UDiagMat}};


#endif //TESTEXACT_H
