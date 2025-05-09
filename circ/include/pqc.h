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
 * Function: Apply Quantum Block (QB)
 * ------------------------------------
 * This function applies a fixed sequence of gates to a quantum state
 */
typedef void (*applyQB)(state_t* state);

/*
 * Function: Apply Parametrized Quantum Block (PQB)
 * -------------------------------------------------
 * This function applies a parametrized quantum block to a quantum state.
 */
typedef void (*applyPQB)(state_t* state, double par);

/*
 * Struct: herm
 * -------------
 * This struct represents a hermitian operator as a sum of weighted hermitian operators
 * Contents:
 *      len     : Number of terms
 *      comp    : Elementary hermitian operators
 *      coeff   : Coefficients
 */
typedef struct herm {
    unsigned int    len;
    applyQB*        comp;
    double*         coeff;
}herm_t;

/*
 * Struct: pqc
 * ------------
 * This struct represents a parametrized quantum circuit
 * Contents:
 *      len     : Number of parametrized quantum blocks
 *      pqbs    : Array of parametrized quantum blocks
 *      parIdx  : Assigns a parameter from par to a parametrized quantum block of pqbs; must be of length len
 *
 */
typedef struct pqc {
    depth_t     len;
    applyPQB*   pqb;
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
void applyDiag(state_t* state, const double diag[]);

/*
 * @brief This function applies a linear combination of quantum blocks to the input state
 *
 * @param[in,out]   state   Quantum state. On exit, state after linear combination of quantum blocks
 * @param[in]       d       Number of quantum blocks in the linear combination
 * @param[in]       qb      Quantum block; i.e, a function applying quantum gates to the quantum state
 * @param[in]       c       Coefficients of the quantum blocks in the linear combination
 */
void lcQB(state_t* state, depth_t d, const applyQB qb[], const cplx_t c[]);

/*
 * =====================================================================================================================
 *                                              Unitary time evolution
 * =====================================================================================================================
 */

/*
 * lbrief This function unitarily evolves the quantum state with respect to the diagonal hermitian operator
 *
 * @param[in,out]   state   Quantum state. On exit, state after evolution with respect to exp(-i*(par/2)*H)
 * @param[in]       diag    The diagonal entries of the hermitian operator H
 * @param[in]       par     Time of the evolution
 */
void evoDiag(state_t* state, const double diag[], double par);

/*
 * @brief This function unitarily evolves of the quantum state with respect to the specified quantum block (sequence of
 * gates)
 *
 * @param[in,out]   state   Quantum state. On exit, state after evolution with respect to exp(-i*(par/2)*qb)
 * @param[in]       qb      Quantum block; i.e, a function applying quantum gates to the quantum state
 * @param[in]       par     Time of the evolution
 *
 * @note WARNING: This function only gives appropriate results if the QB is involutory; i.e., QB*QB = Id
 */
void evoQB(state_t* state, applyQB qb, double par);

/*
 * @brief Wrapper funtion for the unitary evolution of the quantum state with respect to a parametrized quantum block
 *
 * @param[in,out]   state   Quantum state. On exit, state after evolution with respect to pqb(state, par)
 * @param[in]       pqb     Parametrized quantum block
 * @param[in]       par     Time of the evolution
 *
 * @note In accordance with the definition of rotational Pauli gates par is internally divided by two.
 */
void evoPQB(state_t* state, applyPQB pqb, const double par);

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
    applyQB:    evoQB                 \
    applyPQB:   evoPQB                \
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
