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
#include "pstate.h"
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
 * @brief This function returns the mean value of a composite quantum gate with respect to the quantum state
 *
 * @param[in]   state   Quantum state
 * @param[in]   obs     Composite quantum gate; i.e, a function applying quantum gates to the quantum state,
 *                      representing the observable
 *
 * @return  The mean value of the composite quantum gate with respect to the quantum state
 */
double meanCG(const state_t* state, applyCG obs);

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
 *
 * @param[in]       X   Quantum state. On exit, state after evolution with respect to exp(-i*(par/2)*H)
 * @param[in]       Y   Type of evolution; see above for supported types
 *
 * @note Ensure that all expected types are accounted for to prevent unexpected behavior or compile-time errors.
 */
#define meanObs(X, Y) _Generic((Y), \
    double*: meanDiagObs,           \
    applyCG: meanCG,                \
    herm_t*: meanObsHerm            \
    ) (X, Y)

/*
 * =====================================================================================================================
 *                                              Gradient PQC
 * =====================================================================================================================
 */
    void gradPQC(state_t* state, depth_t d, const applyPCG pqc[], const double par[], const applyCG cg[],
        const double obs[], double* grad);

void mmseq(state_t* state, depth_t obsc, const applyCG obs[], depth_t circdepth, const lccg_t lcu[], depth_t link,
           cplx_t* momMat[]);

#ifdef __cplusplus
}
#endif

#endif //EXACT_H
