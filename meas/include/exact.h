//
// Created by Timo Ziegler on 12.03.25.
//

#ifndef EXACT_H
#define EXACT_H
/*
 * =====================================================================================================================
 *                                                      Includes
 * =====================================================================================================================
 */
#ifndef PSTATE_H
#include "../../include/state.h"
#endif

#ifndef PQC_H
#include "pqc.h"
#endif

/*
 * =====================================================================================================================
 *                                                      C++ check
 * =====================================================================================================================
 */
#ifdef __cplusplus
extern "C" {
#endif

/*
 * =====================================================================================================================
 *                                              Observable mean value
 * =====================================================================================================================
 */

/*
 * @brief This function returns the mean value of a diagonal observable with respect to the quantum state
 *
 * @param[in]   state   Quantum state
 * @param[in]   obs     The diagonal entries of an observable O (in the computational basis)
 *
 * @return  The mean value of the observable with respect to the quantum state
 */
double meanDiagObs(const state_t* state, const double obs[]);

/*
 * @brief This function returns the mean value of a hermitian operator with respect to the quantum state
 *
 * @param[in]   state   Quantum state
 * @param[in]   obs     Hermitian operator; i.e., a weighted sum of composite quantum gates representing the observable
 *
 * @return  The mean value of the composite quantum gate with respect to the quantum state
 */
double meanObsHerm(const state_t* state, const herm_t* obs);

/*
 * @brief Function polymorphism of calculating the observable's mean value: Chooses the appropriate function based on
 * second input
 *
 * This macro uses C11's `_Generic` feature to perform function selection based on the type of the second parameter `Y`.
 * Supported Types and Corresponding Functions:
 * - `double*`  : Calls `meanDiagObs`, which processes diagonal observables.
 * - `applyQB`  : Calls `meanObsGC`, which processes composite quantum gates.
 * - `herm_t*`  : Calls `meanObsHerm`, which processes general hermitian operators.
 *
 * @param[in]       X   Quantum state. On exit, state after evolution with respect to exp(-i*(par/2)*H)
 * @param[in]       Y   Type of evolution; see above for supported types
 *
 * @note Ensure that all expected types are accounted for to prevent unexpected behavior or compile-time errors.
 */
#define meanObs(X, Y) _Generic((Y), \
    const double*:  meanDiagObs,    \
    herm_t*:        meanObsHerm     \
    ) (X, Y)

/*
 * =====================================================================================================================
 *                                              Operator mean value
 * =====================================================================================================================
 */
/*
 * @brief This function returns the mean value of a composite quantum gate with respect to the quantum state
 *
 * @param[in]   state   Quantum state
 * @param[in]   obs     Composite quantum gate; i.e, a function applying quantum gates to the quantum state,
 *                      representing the observable
 *
 * @return  The mean value of the composite quantum gate with respect to the quantum state
 */
cplx_t meanCG(const state_t* state, applyCG obs);

/*
 * =====================================================================================================================
 *                                              Gradient PQC
 * =====================================================================================================================
 */

/*
 * This function returns the gradient of an observable's mean value with respect to a quantum state after the
 * application of a parametrized quantum circuit where the observable is diagonal
 *
 * @param[in,out]   state   Input quantum state to the PQC. On exit, the output of the PQC with the CURRENT parameter
 *                          setting
 * @param[in]       obs[]   Observable; i.e., diagonal entries of its matrix representation
 * @param[in]       d       Number of Parametrized Quantum Blocks in the PQC
 * @param[in]       pqc[]   Parametrized quantum circuit; i.e., sequence of parametrized unitary quantum operations
 * @param[in]       par[]   Parameter value of the PQC
 * @param[in]       gen[]   Array of the hamiltonians generating the parametrized unitary transformations
 * @param[in,out]   grad[]  Array of size parc. On exit, holding the PQC's gradient
 */
void gradPQCDiag(const state_t* state, const double obs[], depth_t d, const applyPCG pqc[], const double par[],
    const applyCG gen[], double grad[]);

/*
 * This function returns the gradient of an observable's mean value with respect to a quantum state after the
 * application of a parametrized quantum circuit where the observable is any hermitian operator
 *
 * @param[in,out]   state   Input quantum state to the PQC. On exit, the output of the PQC with the CURRENT parameter
 *                          setting
 * @param[in]       obs[]   Observable; i.e., a weighted sum of elementary hermitian operations
 * @param[in]       d       Number of Parametrized Quantum Blocks in the PQC
 * @param[in]       pqc[]   Parametrized quantum circuit; i.e., sequence of parametrized unitary quantum operations
 * @pram[in]        par[]   Parameter values of the PQC
 * @param[in]       gen[]   Array of the hamiltonians generating the parametrized unitary transformations
 * @param[in,out]   grad[]  Array of size parc. On exit, holding the PQC's gradient
 *
 */
void gradPQCHerm(const state_t* state, const herm_t* obs, depth_t d, const applyPCG pqc[], const double par[],
    const applyCG gen[], double grad[]);

/*
 * @brief Function polymorphism of calculating the gradient of an observable's mean value with respect to a quantum
 * state after the application of a parametrized quantum circuit: Chooses the appropriate function based on second input
 *
 * This macro uses C11's `_Generic` feature to perform function selection based on the type of the second parameter `Y`.
 * Supported Types and Corresponding Functions:
 * - `double*`  : Calls `gradPQCDiag`, which processes diagonal observables.
 * - `applyQB`  : Calls `gradPQCHerm`, which processes general hermitian operators.
 *
 * @param[in,out]   STATE   Quantum state. On exit, state after the PQC being applied
 * @param[in]       OBS     Type of observable; see above for supported types
 * @param[in]       DEPTH   Number of Parametrized Quantum Blocks in the PQC
 * @param[in]       PQC     Parametrized quantum circuit; i.e., sequence of parametrized unitary quantum operations
 * @param[in]       PAR     Parameter values of the PQC
 * @param[in]       GEN     Array of the hamiltonians generating the parametrized unitary transformations
 * @param[in,out]   GRAD    Array of size parc. On exit, holding the PQC's gradient
 *
 * @note Ensure that all expected types are accounted for to prevent unexpected behavior or compile-time errors.
 */
#define gradPQC(STATE, OBS, DEPTH, PQC, PAR, GEN, GRAD) _Generic((OBS), \
    const double*:  gradPQCDiag,                \
    herm_t*:        gradPQCHerm                 \
    ) (STATE, OBS, DEPTH, PQC, PAR, GEN, GRAD)

/*
 * =====================================================================================================================
 *                                              Moment matrix
 * =====================================================================================================================
 */

void mmseq(state_t* state, depth_t obsc, const applyCG obs[], depth_t circdepth, const lccg_t lcu[], depth_t link,
           cplx_t* momMat[]);

#ifdef __cplusplus
}
#endif

#endif //EXACT_H
