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
 *                                                  Type definitions
 * =====================================================================================================================
 */
/*
 * Function:  Matrix representation of a composite quantum gate
 * -------------------------------------------------------------
 * This function returns the matrix representation of some composite quantum gate
 */
typedef cplx_t* (*matCG)(qubit_t qubits);

/*
 * Function:  Matrix representation of a parametrized composite quantum gate
 * --------------------------------------------------------------------------
 * This function returns the matrix representation of some parametrized composite quantum gate
 */
typedef cplx_t* (*matPCG)(qubit_t qubits, double par);

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
double dcoeff[5] = {5.101260256, -3.686076046, 0.654143819, -7.576553655, -0.387642064};

double randPar[5] = {-1.935113472, 0.170656438, -1.916424226, -0.785971871, 1.746181027};

cplx_t zcoeff[15] = {9.20389 - I * 2.07242, 0.58007 + I * 0.58377, -9.98112 - I * 5.35024, -6.38089 + I * 9.44588,
                    4.31933 - I * 7.73341, 8.49492 + I * 6.53551, -8.63475 + I * 1.27771, 2.56232 + I * 6.10205,
                    3.77543 - I * 7.55000, 3.87864 + I * 3.90134, 0.09866 - I * 6.01226, -6.48157 + I * 2.52767,
                    -1.57778 + I * 7.94151, -2.33955 - I * 4.21345, -4.49085 + I * 6.10035};

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

/*
 * =====================================================================================================================
 *                                              Allocated memory and cleanup
 * =====================================================================================================================
 */
state_t testState = {                                       // Define the state the tested function are applied to
    .qubits = 0,
    .dim = 0,
    .vec = NULL
};

cplx_t** vecs = NULL;                                       // Define the state vectors for the function tests
cplx_t* testObsMat[2] = {NULL, NULL};                       // Define the matrix representations of the test observables
cplx_t* testOpMat = NULL;                                   // Define the matrix representation of the test operator

herm_t testHerm = {                                         // Define the hermitian operator struct
    .len = 0,
    .comp = NULL,
    .weight = NULL,
};

double* testVec[2] = {NULL};                                // Define vector valued test results
double* refVec[2] = {NULL};                                 // Define vector valued reference results

cplx_t* testLCUMat[15] = {NULL};                            // Define the matrices for the LCU sequence in testMMseq
cplx_t* testMomMat[2] = {NULL};                                 // Define the moment matrices

extern inline void cleanup(void) {
    if (testState.vec != NULL) {
        free(testState.vec);
        testState.vec = NULL;
    }
    if (vecs != NULL) {
        freeTestVectors(vecs, testState.qubits);
        vecs = NULL;
    }
    for (unsigned char i = 0; i < 2; ++i) {
        if (testObsMat[i] != NULL) {
            free(testObsMat[i]);
            testObsMat[i] = NULL;
        }
    }
    if (testOpMat != NULL) {
        free(testOpMat);
        testOpMat = NULL;
    }
    for (unsigned char i = 0; i < 2; ++i) {
        if (testVec[i] != NULL) {
            free(testVec[i]);
            testVec[i] = NULL;
        }
    }
    for (unsigned char i = 0; i < 2; ++i) {
        if (refVec[i] != NULL) {
            free(refVec[i]);
            refVec[i] = NULL;
        }
    }
    for (unsigned char i = 0; i < 15; ++i) {
        if (testLCUMat[i] != NULL) {
            free(testLCUMat[i]);
            testLCUMat[i] = NULL;
        }
    }
    for (unsigned char i = 0; i < 2; ++i) {
        if (testMomMat[i] != NULL) {
            free(testMomMat[i]);
            testMomMat[i] = NULL;
        }
    }
}

/*
 * =====================================================================================================================
 *                                                  Composite quantum gates
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

applyCG cg[25] = {x2, y2, z2, swap2, diag, x3, y3, z3, swap3, diag, x4, y4, z4, swap4, diag, x5, y5, z5, swap5, diag,
                  x6, y6, z6, swap6, diag
};

/*
 * =====================================================================================================================
 *                                  Composite quantum gates - matrix representation
 * =====================================================================================================================
 */
extern inline cplx_t* xMat(const qubit_t qubits) {
    const dim_t dim = 1 << qubits;
    cplx_t* out = singleQubitGateMat(XMAT, qubits, qubits - 1);
    for (qubit_t i = qubits - 1; i > 0; --i) {
        cplx_t* tmp = singleQubitGateMat(XMAT, qubits, i - 1);
        cmatMulInPlace(out, tmp, dim);
        free(tmp);
    }
    return out;
}
extern inline cplx_t* yMat(const qubit_t qubits) {
    const dim_t dim = 1 << qubits;
    cplx_t* out = singleQubitGateMat(YMAT, qubits, qubits - 1);
    for (qubit_t i = qubits - 1; i > 0; --i) {
        cplx_t* tmp = singleQubitGateMat(YMAT, qubits, i - 1);
        cmatMulInPlace(out, tmp, dim);
        free(tmp);
    }
    return out;
}
extern inline cplx_t* zMat(const qubit_t qubits) {
    const dim_t dim = 1 << qubits;
    cplx_t* out = singleQubitGateMat(ZMAT, qubits, qubits - 1);
    for (qubit_t i = qubits - 1; i > 0; --i) {
        cplx_t* tmp = singleQubitGateMat(ZMAT, qubits, i - 1);
        cmatMulInPlace(out, tmp, dim);
        free(tmp);
    }
    return out;
}
extern inline cplx_t* swapMat(const qubit_t qubits) {
    const dim_t dim = 1 << qubits;
    const qubit_t start = qubits / 2;                               // start = k for qubits = 2k and qubits = (2k + 1)
    cplx_t* out = swapGateMat(qubits, 2 * start - 2, 2 * start - 1);    // SWAP the last qubit with odd index to left
    for (qubit_t i = start - 1; i > 0; --i) {                       // Starting at the next qubit with odd index
        cplx_t* tmp = swapGateMat(qubits, 2 * i - 2, 2 * i - 1);
        cmatMulInPlace(out, tmp, dim);
        free(tmp);
    }
    return out;
}

extern inline cplx_t* diagMat(const qubit_t qubits) {
    return diagGateMat(qubits, diagObs);
}

matCG cgMat[5] = {xMat, yMat, zMat, swapMat, diagMat};

/*
 * =====================================================================================================================
 *                                          Parametrized composite quantum gates
 * =====================================================================================================================
 */
// 2 qubits
extern inline void evoX2(state_t* state, const double par) {
    evolve(state, x2, par);
}
extern inline void rx2(state_t* state, const double par) {
    applyRX(state, 0, par);
    applyRX(state, 1, par);
}

extern inline void evoY2(state_t* state, const double par) {
    evolve(state, y2, par);
}
extern inline void ry2(state_t* state, const double par) {
    applyRY(state, 0, par);
    applyRY(state, 1, par);
}

extern inline void evoZ2(state_t* state, const double par) {
    evolve(state, z2, par);
}
extern inline void rz2(state_t* state, const double par) {
    applyRZ(state, 0, par);
    applyRZ(state, 1, par);
}

extern inline void evoSwap2(state_t* state, const double par) {
    evolve(state, swap2, par);
}
extern inline void rswap2(state_t* state, const double par) {
    applyRSWAP(state, 0, 1, par);
}

extern inline void evoDiag2(state_t* state, const double par) {
    evoDiag(state, diagObs, par);
}

// 3 qubits
extern inline void evoX3(state_t* state, const double par) {
    evolve(state, x3, par);
}
extern inline void rx3(state_t* state, const double par) {
    applyRX(state, 0, par);
    applyRX(state, 1, par);
    applyRX(state, 2, par);
}

extern inline void evoY3(state_t* state, const double par) {
    evolve(state, y3, par);
}
extern inline void ry3(state_t* state, const double par) {
    applyRY(state, 0, par);
    applyRY(state, 1, par);
    applyRY(state, 2, par);
}

extern inline void evoZ3(state_t* state, const double par) {
    evolve(state, z3, par);
}
extern inline void rz3(state_t* state, const double par) {
    applyRZ(state, 0, par);
    applyRZ(state, 1, par);
    applyRZ(state, 2, par);
}

extern inline void evoSwap3(state_t* state, const double par) {
    evolve(state, swap3, par);
}
extern inline void rswap3(state_t* state, const double par) {
    applyRSWAP(state, 0, 1, par);
}

extern inline void evoX4(state_t* state, const double par) {
    evolve(state, x4, par);
}
extern inline void rx4(state_t* state, const double par) {
    applyRX(state, 0, par);
    applyRX(state, 1, par);
    applyRX(state, 2, par);
    applyRX(state, 3, par);
}

extern inline void evoY4(state_t* state, const double par) {
    evolve(state, y4, par);
}
extern inline void ry4(state_t* state, const double par) {
    applyRY(state, 0, par);
    applyRY(state, 1, par);
    applyRY(state, 2, par);
    applyRY(state, 3, par);
}

extern inline void evoZ4(state_t* state, const double par) {
    evolve(state, z4, par);
}
extern inline void rz4(state_t* state, const double par) {
    applyRZ(state, 0, par);
    applyRZ(state, 1, par);
    applyRZ(state, 2, par);
    applyRZ(state, 3, par);
}

extern inline void evoSwap4(state_t* state, const double par) {
    evolve(state, swap4, par);
}
extern inline void rswap4(state_t* state, const double par) {
    applyRSWAP(state, 0, 1, par);
    applyRSWAP(state, 2, 3, par);
}

extern inline void evoX5(state_t* state, const double par) {
    evolve(state, x5, par);
}
extern inline void rx5(state_t* state, const double par) {
    applyRX(state, 0, par);
    applyRX(state, 1, par);
    applyRX(state, 2, par);
    applyRX(state, 3, par);
    applyRX(state, 4, par);
}

extern inline void evoY5(state_t* state, const double par) {
    evolve(state, y5, par);
}
extern inline void ry5(state_t* state, const double par) {
    applyRY(state, 0, par);
    applyRY(state, 1, par);
    applyRY(state, 2, par);
    applyRY(state, 3, par);
    applyRY(state, 4, par);
}

extern inline void evoZ5(state_t* state, const double par) {
    evolve(state, z5, par);
}
extern inline void rz5(state_t* state, const double par) {
    applyRZ(state, 0, par);
    applyRZ(state, 1, par);
    applyRZ(state, 2, par);
    applyRZ(state, 3, par);
    applyRZ(state, 4, par);
}

extern inline void evoSwap5(state_t* state, const double par) {
    evolve(state, swap5, par);
}
extern inline void rswap5(state_t* state, const double par) {
    applyRSWAP(state, 0, 1, par);
    applyRSWAP(state, 2, 3, par);
}

extern inline void evoX6(state_t* state, const double par) {
    evolve(state, x6, par);
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
    evolve(state, y6, par);
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
    evolve(state, z6, par);
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
    evolve(state, swap6, par);
}
extern inline void rswap6(state_t* state, const double par) {
    applyRSWAP(state, 0, 1, par);
    applyRSWAP(state, 2, 3, par);
    applyRSWAP(state, 4, 5, par);
}

applyPCG pqc[25] = {evoX2, evoY2, evoZ2, evoSwap2, evoDiag2, evoX3, evoY3, evoZ3, evoSwap3, evoDiag2, evoX4, evoY4,
                    evoZ4, evoSwap4, evoDiag2, evoX5, evoY5, evoZ5, evoSwap5, evoDiag2, evoX6, evoY6, evoZ6, evoSwap6, evoDiag2};

applyPCG rpqc[20] = {rx2, ry2, rz2, rswap2, rx3, ry3, rz3, rswap3, rx4, ry4, rz4, rswap4, rx5, ry5, rz5, rswap5, \
     rx6, ry6, rz6, rswap6};

/*
 * =====================================================================================================================
 *                              Parametrized composite quantum gates - matrix representation
 * =====================================================================================================================
 */
extern inline cplx_t* evoXMat(const qubit_t qubits, const double par) {
     const dim_t dim = 1 << qubits;
     cplx_t* gateMat = xMat(qubits);
     cplx_t* out = zexpm(gateMat, par / 2., dim);
     free(gateMat);

     return out;
 }
extern inline cplx_t* rxMat(const qubit_t qubits, const double par) {
    const dim_t dim = 1 << qubits;
    cplx_t* out = RXGateMat(qubits, qubits - 1, par);
    for (qubit_t i = qubits - 1; i > 0; --i) {
        cplx_t* tmp = RXGateMat(qubits, i - 1, par);
        cmatMulInPlace(out, tmp, dim);
        free(tmp);
    }
    return out;
}

    extern inline cplx_t* evoYMat(const qubit_t qubits, const double par) {
    const dim_t dim = 1 << qubits;
    cplx_t* gateMat = yMat(qubits);
    cplx_t* out = zexpm(gateMat, par / 2., dim);
    free(gateMat);

    return out;
}
    extern inline cplx_t* ryMat(const qubit_t qubits, const double par) {
    const dim_t dim = 1 << qubits;
    cplx_t* out = RYGateMat(qubits, qubits - 1, par);
    for (qubit_t i = qubits - 1; i > 0; --i) {
        cplx_t* tmp = RYGateMat(qubits, i - 1, par);
        cmatMulInPlace(out, tmp, dim);
        free(tmp);
    }
    return out;
}

extern inline cplx_t* evoZMat(const qubit_t qubits, const double par) {
    const dim_t dim = 1 << qubits;
    cplx_t* gateMat = zMat(qubits);
    cplx_t* out = zexpm(gateMat, par / 2., dim);
    free(gateMat);

    return out;
}
extern inline cplx_t* rzMat(const qubit_t qubits, const double par) {
    const dim_t dim = 1 << qubits;
    cplx_t* out = RZGateMat(qubits, qubits - 1, par);
    for (qubit_t i = qubits - 1; i > 0; --i) {
        cplx_t* tmp = RZGateMat(qubits, i - 1, par);
        cmatMulInPlace(out, tmp, dim);
        free(tmp);
    }
    return out;
}

extern inline cplx_t* evoSwapMat(const qubit_t qubits, const double par) {
    const dim_t dim = 1 << qubits;
    cplx_t* gateMat = swapMat(qubits);
    cplx_t* out = zexpm(gateMat, par / 2., dim);
    free(gateMat);

    return out;
}
extern inline cplx_t* rSwapMat(const qubit_t qubits, const double par) {
    const dim_t dim = 1 << qubits;
    const qubit_t start = qubits / 2;                               // start = k for qubits = 2k and qubits = (2k + 1)
    cplx_t* out = RswapGateMat(qubits, 2 * start - 2, 2 * start - 1, par);
    for (qubit_t i = start - 1; i > 0; --i) {                       // Starting at the next qubit with odd index
        cplx_t* tmp = RswapGateMat(qubits, 2 * i - 2, 2 * i - 1, par);
        cmatMulInPlace(out, tmp, dim);
        free(tmp);
    }
    return out;
}

extern inline cplx_t* evoDiagMat(const qubit_t qubits, const double par) {
    const dim_t dim = 1 << qubits;
    cplx_t* gateMat = diagMat(qubits);
    cplx_t* out = zexpm(gateMat, par / 2., dim);
    free(gateMat);

    return out;
}

matPCG pqcMat[5] = {evoXMat, evoYMat, evoZMat, evoSwapMat, evoDiagMat};

matPCG rpqcMat[5] = {rxMat, ryMat, rzMat, rSwapMat};

/*
 * =====================================================================================================================
 *                                      Generators of products of parametrized gates
 * =====================================================================================================================
 */

double weight[6] = {1., 1., 1., 1., 1., 1.};

extern inline void singleX0(state_t* state) {
    applyX(state, 0);
}
extern inline void singleX1(state_t* state) {
    applyX(state, 1);
}
extern inline void singleX2(state_t* state) {
    applyX(state, 2);
}
extern inline void singleX3(state_t* state) {
    applyX(state, 3);
}
extern inline void singleX4(state_t* state) {
    applyX(state, 4);
}
extern inline void singleX5(state_t* state) {
    applyX(state, 5);
}
applyCG singleX[6] = {singleX0, singleX1, singleX2, singleX3, singleX4, singleX5};

extern inline void genRX2(state_t* state) {
    herm_t gen;
    gen.len = 2;
    gen.comp = singleX;
    gen.weight = weight;

    apply(state, &gen);
}

extern inline void genRX3(state_t* state) {
    herm_t gen;
    gen.len = 3;
    gen.comp = singleX;
    gen.weight = weight;

    apply(state, &gen);
}

extern inline void genRX4(state_t* state) {
    herm_t gen;
    gen.len = 4;
    gen.comp = singleX;
    gen.weight = weight;

    apply(state, &gen);
}

extern inline void genRX5(state_t* state) {
    herm_t gen;
    gen.len = 5;
    gen.comp = singleX;
    gen.weight = weight;

    apply(state, &gen);
}

extern inline void genRX6(state_t* state) {
    herm_t gen;
    gen.len = 6;
    gen.comp = singleX;
    gen.weight = weight;

    apply(state, &gen);
}

extern inline void singleY0(state_t* state) {
    applyY(state, 0);
}
extern inline void singleY1(state_t* state) {
    applyY(state, 1);
}
extern inline void singleY2(state_t* state) {
    applyY(state, 2);
}
extern inline void singleY3(state_t* state) {
    applyY(state, 3);
}
extern inline void singleY4(state_t* state) {
    applyY(state, 4);
}
extern inline void singleY5(state_t* state) {
    applyY(state, 5);
}
applyCG singleY[6] = {singleY0, singleY1, singleY2, singleY3, singleY4, singleY5};

extern inline void genRY2(state_t* state) {
    herm_t gen;
    gen.len = 2;
    gen.comp = singleY;
    gen.weight = weight;

    apply(state, &gen);
}

extern inline void genRY3(state_t* state) {
    herm_t gen;
    gen.len = 3;
    gen.comp = singleY;
    gen.weight = weight;

    apply(state, &gen);
}

extern inline void genRY4(state_t* state) {
    herm_t gen;
    gen.len = 4;
    gen.comp = singleY;
    gen.weight = weight;

    apply(state, &gen);
}

extern inline void genRY5(state_t* state) {
    herm_t gen;
    gen.len = 5;
    gen.comp = singleY;
    gen.weight = weight;

    apply(state, &gen);
}

extern inline void genRY6(state_t* state) {
    herm_t gen;
    gen.len = 6;
    gen.comp = singleY;
    gen.weight = weight;

    apply(state, &gen);
}

extern inline void singleZ0(state_t* state) {
    applyZ(state, 0);
}
extern inline void singleZ1(state_t* state) {
    applyZ(state, 1);
}
extern inline void singleZ2(state_t* state) {
    applyZ(state, 2);
}
extern inline void singleZ3(state_t* state) {
    applyZ(state, 3);
}
extern inline void singleZ4(state_t* state) {
    applyZ(state, 4);
}
extern inline void singleZ5(state_t* state) {
    applyZ(state, 5);
}
applyCG singleZ[6] = {singleZ0, singleZ1, singleZ2, singleZ3, singleZ4, singleZ5};

extern inline void genRZ2(state_t* state) {
    herm_t gen;
    gen.len = 2;
    gen.comp = singleZ;
    gen.weight = weight;

    apply(state, &gen);
}

extern inline void genRZ3(state_t* state) {
    herm_t gen;
    gen.len = 3;
    gen.comp = singleZ;
    gen.weight = weight;

    apply(state, &gen);
}

extern inline void genRZ4(state_t* state) {
    herm_t gen;
    gen.len = 4;
    gen.comp = singleZ;
    gen.weight = weight;

    apply(state, &gen);
}

extern inline void genRZ5(state_t* state) {
    herm_t gen;
    gen.len = 5;
    gen.comp = singleZ;
    gen.weight = weight;

    apply(state, &gen);
}

extern inline void genRZ6(state_t* state) {
    herm_t gen;
    gen.len = 6;
    gen.comp = singleZ;
    gen.weight = weight;

    apply(state, &gen);
}

extern inline void singleSWAP0(state_t* state) {
    applySWAP(state, 0, 1);
}
extern inline void singleSWAP1(state_t* state) {
    applySWAP(state, 2, 3);
}
extern inline void singleSWAP2(state_t* state) {
    applySWAP(state, 4, 5);
}
applyCG singleSWAP[3] = {singleSWAP0, singleSWAP1, singleSWAP2};

extern inline void genRSWAP2(state_t* state) {
    herm_t gen;
    gen.len = 1;
    gen.comp = singleSWAP;
    gen.weight = weight;

    apply(state, &gen);
}

extern inline void genRSWAP3(state_t* state) {
    herm_t gen;
    gen.len = 1;
    gen.comp = singleSWAP;
    gen.weight = weight;

    apply(state, &gen);
}

extern inline void genRSWAP4(state_t* state) {
    herm_t gen;
    gen.len = 2;
    gen.comp = singleSWAP;
    gen.weight = weight;

    apply(state, &gen);
}

extern inline void genRSWAP5(state_t* state) {
    herm_t gen;
    gen.len = 2;
    gen.comp = singleSWAP;
    gen.weight = weight;

    apply(state, &gen);
}

extern inline void genRSWAP6(state_t* state) {
    herm_t gen;
    gen.len = 3;
    gen.comp = singleSWAP;
    gen.weight = weight;

    apply(state, &gen);
}

applyCG gen[20] = {genRX2, genRY2, genRZ2, genRSWAP2, genRX3, genRY3, genRZ3, genRSWAP3, genRX4, genRY4, genRZ4,
    genRSWAP4, genRX5, genRY5, genRZ5, genRSWAP5, genRX6, genRY6, genRZ6, genRSWAP6
};

/*
 * =====================================================================================================================
 *                          Generators of products of parametrized gates - matrix representation
 * =====================================================================================================================
 */

extern inline cplx_t* genRXMat(const qubit_t qubits) {
    const dim_t dim = 1 << qubits;
    cplx_t* out = singleQubitGateMat(XMAT, qubits, qubits - 1);
    for (qubit_t i = qubits - 1; i > 0; --i) {
        cplx_t* tmp = singleQubitGateMat(XMAT, qubits, i - 1);
        cmatAddInPlace(out, tmp, dim);
        free(tmp);
    }
    return out;
}

extern inline cplx_t* genRYMat(const qubit_t qubits) {
    const dim_t dim = 1 << qubits;
    cplx_t* out = singleQubitGateMat(YMAT, qubits, qubits - 1);
    for (qubit_t i = qubits - 1; i > 0; --i) {
        cplx_t* tmp = singleQubitGateMat(YMAT, qubits, i - 1);
        cmatAddInPlace(out, tmp, dim);
        free(tmp);
    }
    return out;
}

extern inline cplx_t* genRZMat(const qubit_t qubits) {
    const dim_t dim = 1 << qubits;
    cplx_t* out = singleQubitGateMat(ZMAT, qubits, qubits - 1);
    for (qubit_t i = qubits - 1; i > 0; --i) {
        cplx_t* tmp = singleQubitGateMat(ZMAT, qubits, i - 1);
        cmatAddInPlace(out, tmp, dim);
        free(tmp);
    }
    return out;
}

extern inline cplx_t* genRSWAPMat(const qubit_t qubits) {
    const dim_t dim = 1 << qubits;
    const qubit_t start = qubits / 2;                               // start = k for qubits = 2k and qubits = (2k + 1)
    cplx_t* out = swapGateMat(qubits, 2 * start - 2, 2 * start - 1);    // SWAP the last qubit with odd index to left
    for (qubit_t i = start - 1; i > 0; --i) {                       // Starting at the next qubit with odd index
        cplx_t* tmp = swapGateMat(qubits, 2 * i - 2, 2 * i - 1);
        cmatAddInPlace(out, tmp, dim);
        free(tmp);
    }
    return out;
}

matCG genMat[4] = {genRXMat, genRYMat, genRZMat, genRSWAPMat};

/*
 * =====================================================================================================================
 *                                      Parametrized quantum blocks with fixed parameters
 * =====================================================================================================================
 */

inline void applyTestHerm(state_t* state) {                 // Define a function that applies testHerm to a state
    apply(state, &testHerm);
}
inline void applyDiagObs(const state_t* state) {            // Define a function that applies diagObs to a state
    apply(state, diagObs);
}
applyCG obs[2] = {applyDiagObs, applyTestHerm};             // For mmseq define the observables as function pointers

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

applyCG unitary[75] = {x2, y2, z2, swap2, diag, UX2, UY2, UZ2, USwap2, UDiag, URX2, URY2, URZ2, URSwap2, UDiag,
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

matCG unitaryMat[3][5] = {{xMat, yMat, zMat, swapMat, diagMat}, {UXMat, UYMat, UZMat, USwapMat, UDiagMat},
                        {URXMat, URYMat, URZMat, URSwapMat, UDiagMat}};


#endif //TESTEXACT_H
