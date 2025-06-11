//
// Created by Timo Ziegler on 12.03.25.
//

#ifndef PQC_H
#define PQC_H
/*
 * =====================================================================================================================
 *                                                      Includes
 * =====================================================================================================================
 */
#ifndef PSTATE_H
#include "pstate.h"
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
 *                                                  Type definitions
 * =====================================================================================================================
 */

/*
 * Function: Apply Composite Gate (CG)
 * ------------------------------------
 * This function applies a fixed composition of quantum gates to a quantum state
 */
typedef void (*applyCG)(state_t* state);

/*
 * Function: Apply Parametrized Composite Gate (PCG)
 * -------------------------------------------------
 * This function applies a parametrized composition of quantum gates to a quantum state.
 */
typedef void (*applyPCG)(state_t* state, double par);

/*
 * Struct: herm
 * -------------
 * This struct represents a hermitian operator; i.e., a weighted sum of hermitian operators with DOUBLE PRECISION
 * weights
 * Contents:
 *      len     : Number of terms
 *      comp    : Elementary hermitian operators; i.e. composition of quantum gates
 *      weights : Weights (DOUBLE PRECISION)
 */
typedef struct herm {
    unsigned int    len;
    applyCG*        comp;
    double*         weight;
}herm_t;

/*
 * Struct: lccg
 * -------------
 * This struct represents a linear combination of composite gates; i.e, a weighted sum of compositions of quantum gates
 * with COMPLEX DOUBLE weights
 * Contents:
 *      len     : Number of terms
 *      comp    : Composite gates
 *      weights : Weights (COMPLEX DOUBLE)
 */
typedef struct lccg {
    depth_t     len;
    applyCG*    comp;
    cplx_t*     weight;
}lccg_t;

/*
 * Struct: pqc
 * ------------
 * This struct represents a parametrized quantum circuit
 * Contents:
 *      len     : Number of parametrized quantum blocks
 *      pcg     : Array of parametrized composite gates
 *      parIdx  : Assigns a parameter from par to a parametrized composite gate of pcg; must be of length len
 *
 */
typedef struct pqc {
    depth_t     len;
    applyPCG*   pqb;
    depth_t*    parIdx;
} pqc_t;

/*
 * =====================================================================================================================
 *                                                  Direct application
 * =====================================================================================================================
 */

/*
 * @brief This function applies the diagonal hermitian operator to the quantum state
 *
 * @param[in,out]   state   Quantum state; On exit, the statevector is the product of H.v
 * @param[in]       diag    The diagonal entries of the hermitian operator H
 *
 * @note    WARNING: This action is unphysical in the sense of a circuit model and rather belongs to the measurement
 *          part of a circuit. As long as the hermitian operator is not involutive; i.e., H**2, the output is no longer
 *          a quantum state, due to missing normalization.
 */
void applyDiag(const state_t* state, const double diag[]);

/*
 * @brief Wrapper function for the application of a composite quantum gate to a quntum state
 *
 * @param[in,out]   state   Quantum state. On exit, state after application of the composite quantum gate
 * @param[in]       cg      Composite quantum gate
 *
 */
void applyCGwrapper(state_t* state, applyCG cg);

/*
 * @brief This function applies a general hermitian operator (sum of weighted hermitian composite gates) to the quantum
 * state
 *
 * @param[in,out]   state   Quantum state; On exit, the statevector is the product of H.v
 * @param[in]       herm    The hermitian operator; i.e., a weighted sum of hermitian composite gates
 *
 * @note    WARNING: This action is unphysical in the sense of a circuit model and rather belongs to the measurement
 *          part of a circuit. As long as the hermitian operator is not unitary, the output is no longer a quantum
 *          state, due to missing normalization.
 */
void applyHerm(state_t* state, const herm_t* herm);

/*
 * @brief This function applies a linear combination of composite gates to the input state
 *
 * @param[in,out]   state   Quantum state. On exit, state after linear combination of composite gates
 * @param[in]       lccg    The linear combination of composite gates; i.e., a weighted sum of compositions of quantum
 *                          gates
 *
 * @note    WARNING: This action may be unphysical in the sense of the circuit model. As long as the linear combination
 *          is not unitary, the output is no longer a quantum state, due to missing normalization.
 */
void applyLCCG(state_t* state, const lccg_t* lccg);

/*
 * @brief Function polymorphism of the application of composite quantum gates to a quantum state: Chooses the
 * appropriate function based on second input
 *
 * This macro uses C11's `_Generic` feature to perform function selection based on the type of the second parameter `Y`.
 * Supported Types and Corresponding Functions:
 * - `double*`  : Calls `evoDiag`, which processes diagonal hermitian operators.
 * - `applyQB`  : Calls `evoQB`, which processes quantum blocks.
 * - `applyPQB` : Calls `evoPQB`, a wrapper function for parametrized quantum blocks.
 *
 * @param[in,out]   X   Quantum state. On exit, state after evolution with respect to exp(-i*(par/2)*H)
 * @param[in]       Y   Type of evolution; see above for supported types
 * @param[in]       Z   Time of the evolution
 *
 * @note Ensure that all expected types are accounted for to prevent unexpected behavior or compile-time errors.
 */
#define apply(X, Y) _Generic((Y),   \
    double*:    applyDiag,          \
    applyCG:    applyCGwrapper,     \
    herm_t*:    applyHerm,          \
    lccg_t*:    applyLCCG           \
) (X, Y)

/*
 * =====================================================================================================================
 *                                              Unitary time evolution
 * =====================================================================================================================
 */

/*
 * @brief This function unitarily evolves the quantum state for time par/2 with respect to the diagonal hermitian
 * operator
 *
 * @param[in,out]   state   Quantum state. On exit, state after evolution with respect to exp(-i*(par/2)*H)
 * @param[in]       diag    The diagonal entries of the hermitian operator H
 * @param[in]       par     Time of the evolution
 */
void evoDiag(const state_t* state, const double diag[], double par);

/*
 * @brief This function unitarily evolves the quantum state for time par/2 with respect to the specified composition of
 * quantum gates
 *
 * @param[in,out]   state   Quantum state. On exit, state after evolution with respect to exp(-i*(par/2)*cg)
 * @param[in]       cg      Composite gate; i.e, a function applying a composition of quantum gates to the quantum state
 * @param[in]       par     Time of the evolution
 *
 * @note WARNING: This function only gives appropriate results if the CG is involutory; i.e., QB*QB = Id
 */
void evoCG(state_t* state, applyCG cg, double par);

/*
 * @brief This function unitarily evolves of the quantum state for time par/2 with respect to the hermitian operator
 * (weighted sum of composite quantum gates)
 *
 * @param[in,out]   state   Quantum state. On exit, state after evolution with respect to exp(-i*(par/2)*qb)
 * @param[in]       herm    Hermitian operator; i.e, a struct inheriting the weighted sum of composite quantum gates
 * @param[in]       par     Time of the evolution
 *
 * @note WARNING:   This function only gives appropriate results if the QBs pairwise commute (Trotter-Suzuki
 *                  approximation) and each CG is involutory; i.e., CG*CG = Id.
 */
void evoHerm(state_t* state, const herm_t* herm, const double par);

/*
 * @brief Wrapper function for the unitary evolution of the quantum state for time par/2 with respect to a parametrized
 * composite gate
 *
 * @param[in,out]   state   Quantum state. On exit, state after evolution with respect to pcg(state, par)
 * @param[in]       pcg     Parametrized composite quantum gate
 * @param[in]       par     Time of the evolution
 *
 * @note In accordance with the definition of rotational Pauli gates par is internally divided by two.
 */
void evoPCGwrapper(state_t* state, applyPCG pcg, double par);

/*
 * @brief Function polymorphism of the quantum state's unitary evolution: Chooses the appropriate function based on
 * second input
 *
 * This macro uses C11's `_Generic` feature to perform function selection based on the type of the second parameter `Y`.
 * Supported Types and Corresponding Functions:
 * - `double*`  : Calls `evoDiag`, which processes diagonal hermitian operators.
 * - `applyQB`  : Calls `evoQB`, which processes quantum blocks.
 * - `applyPQB` : Calls `evoPQB`, a wrapper function for parametrized quantum blocks.
 *
 * @param[in,out]   X   Quantum state. On exit, state after evolution with respect to exp(-i*(par/2)*H)
 * @param[in]       Y   Type of evolution; see above for supported types
 * @param[in]       Z   Time of the evolution
 *
 * @note Ensure that all expected types are accounted for to prevent unexpected behavior or compile-time errors.
 */
#define evolve(X, Y, Z) _Generic((Y), \
    double*:    evoDiag,              \
    applyCG:    evoCG,                \
    herm_t*:    evoHerm,              \
    applyPCG:   evoPCGwrapper         \
    ) (X, Y, Z)

/*
 * @brief This function applies a parametrized quantum circuit to the quantum state; i.e., a sequence of parametrized
 * quantum blocks.
 *
 * @param[in,out]   state   Quantum state. On exit, states after passing the parametrized quantum circuit
 * @param[in]       pqc     Parametrized quantum circuit; struct inheriting the parametrized blocks, their total, the
 *                          parameters and an assignment of parameters to parametrized quantum blocks
 */
void applyPQC(state_t* state, pqc_t pqc, const double par[]);

#ifdef __cplusplus
}
#endif
#endif //PQC_H
