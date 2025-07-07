//
// Created by Timo Ziegler on 05.07.25.
//
/*
 * =====================================================================================================================
 *                                                      includes
 * =====================================================================================================================
 */
#ifndef __MATH__
#include <math.h>
#endif

#ifndef MSTATE_H
#include "mstate.h"
#endif

#ifndef _STDIO_H_
#include <stdio.h>
#endif

#ifndef _STDLIB_H_
#include <stdlib.h>
#endif

#ifdef MACOS
#include <vecLib/cblas_new.h>
#include <vecLib/lapack.h>
#else
#include <cblas.h>
#include <lapacke.h>
#endif

/*
 * =====================================================================================================================
 *                                                  State preparation
 * =====================================================================================================================
 */

/*
 * @brief	Defines a mixed quantum state of a specified number of qubits and allocates the required memory for the
 *			statevector
 *
 * @param[in,out]	state	Address of the state
 * @param[in]		qubits	Number of qubits
 *
 * @returns	This function has no return value; alters the state in place
 */
void stateInitEmpty(state_t* state, const qubit_t qubits) {
    state->qubits = qubits;
    state->dim = POW2(qubits, dim_t);
    if((state->vec = calloc(state->dim * (state->dim + 1) / 2, sizeof(cplx_t))) == NULL) {
        fprintf(stderr, "stateInitEmpty(): state->vec allocation failed\n");
        exit(EXIT_FAILURE);
    }
}

/*
 * @brief	Initialize a mixed quantum state to the plus state (uniform superposition of all computational basis states)
 *
 * @param[in,out]	state		Address of the state
 * @param[in]		qubits		Number of qubits
 *
 * @returns	This function has no return value; alters the state in place
 */
void stateInitPlus(state_t* state, const qubit_t qubits) {
    const double prefactor = 1. / state->dim;
    for (dim_t i = 0; i < state->dim; ++i) {
        state->vec[i] = prefactor;
    }
}

/*
 * @brief	Initialize a pure quantum state to the state given by a complex vector
 *
 * @param[in,out]	state	Address of the state
 * @param[in]		vector	Complex vector with DOUBLE PRECISION that is copied to the state's vector representation
 *
 * @returns	This function has no return value; alters the state in place
 */
void stateInitVector(state_t* state, const cplx_t vector[]) {
    const dim_t len = state->dim * (state->dim + 1) / 2;
    cblas_zcopy(len, vector, 1, state->vec, 1);
}

/*
 * =====================================================================================================================
 *                                  		        Free state
 * =====================================================================================================================
 */

void stateFreeVector(state_t* state) {
    free(state->vec);
    state->vec = NULL;
}

/*
 * =====================================================================================================================
 *                                              Binary operations
 * =====================================================================================================================
 */

unsigned char zhpev(const dim_t n, cplx_t* h[], double w[]);
void zlowertofull(dim_t n, const cplx_t* lower, cplx_t* full);

/*
 * @brief	Calculate the overlap (fidelity) between two mixed quantum states, by the Uhlmann fidelity, i.e.,
 *                          F(rho, sigma) = [Tr(sqrt(sqrt(rho) sigma sqrt(rho)))]^2
 *
 * param[in,out]	state1	rho
 * param[in,out]	state2	sigma
 */
cplx_t stateOverlap(const state_t* state1, const state_t* state2) {
    // Compute the eigenvalues and -vectors of rho and store the latter in a temporary
    const dim_t len = state1->dim * (state1->dim + 1) / 2;
    cplx_t *u, *sqrtRho, *tmp;
    if ((u = malloc(len * sizeof(cplx_t))) == NULL) {
        fprintf(stderr, "stateOverlap(): tmp allocation failed\n");
        exit(EXIT_FAILURE);
    }
    double* eigv = calloc(state1->dim, sizeof(double));
    if (eigv == NULL) {
        fprintf(stderr, "stateOverlap(): eigv allocation failed\n");
        free(u);
        exit(EXIT_FAILURE);
    }

    if (zhpev(state1->dim, &u, eigv)) {
        fprintf(stderr, "stateOverlap(): zhpev failed\n");
        exit(EXIT_FAILURE);
    }

    // Build the matrix diag(sqrt(eigv[i]))
    if ((sqrtRho = calloc(len, sizeof(cplx_t))) == NULL) {
        fprintf(stderr, "stateOverlap(): sqrtRho allocation failed\n");
        free(eigv);
        free(u);
        exit(EXIT_FAILURE);
    }
    for (dim_t j = 0; j < state1->dim; ++j) {
        if (eigv[j] > 1e-12) {
            sqrtRho[state1->dim * j - j * (j - 1) / 2] = sqrt(eigv[j]);
        } else {
            fprintf(stderr, "stateOverlap(): eigv[%ld] is negative; cannot compute srt of state1\n", j);
            free(sqrtRho);
            free(u);
            free(eigv);
            exit(EXIT_FAILURE);
        }
    }
    free(eigv);
    eigv = NULL;

    // Compute sqrt(rho) = U*sqrt(eigv)*U**H
    const cplx_t ALPHA = 1.;
    const cplx_t BETA = 0.;

    if ((tmp = malloc(len * sizeof(cplx_t))) == NULL) {
        fprintf(stderr, "stateOverlap(): tmp allocation failed\n");
        free(sqrtRho);
        free(u);
        exit(EXIT_FAILURE);
    }
    cblas_zhemm(CblasColMajor, CblasRight, CblasLower,
        state1->dim, state1->dim,
        &ALPHA, u, state1->dim,
        sqrtRho, state1->dim,
        &BETA, tmp, state1->dim);
    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasConjTrans,
        state1->dim, state1->dim, state1->dim,
        &ALPHA, tmp, state1->dim,
        u, state1->dim,
        &BETA, sqrtRho, state1->dim);
    free(tmp);
    tmp = NULL;
    free(u);
    u = NULL;

    // Expand sigma to a full matrix
    if ((tmp = malloc(state2->dim * state2->dim * sizeof(cplx_t))) == NULL) {
        fprintf(stderr, "stateOverlap(): tmp allocation failed\n");
        free(sqrtRho);
        exit(EXIT_FAILURE);
    }
    zlowertofull(state2->dim, state2->vec, tmp);

    // Compute sqrt(rho)*sigma*sqrt(rho)

}

/*
 * @brief   Computes eigenvalues and -vectors of a HERMITIAN matrix H in PACKED column major format
 *
 * @param[in]       n   Order of matrix H
 * @param[in,out]   h   Address of the Hermitian matrix H in packed column major format. On exit, the orthonormal
 *                      eigenvectors of H with the i-th column holding the eigenvector associated with w[i]
 * @param[in,out]   w   Array of size n. On exit, the eigenvalues of H in ascending order
 *
 * @returns 0 on success, 1 otherwise
 */
unsigned char zhpev(const dim_t n, cplx_t* h[], double w[]) {
    const char JOBZ = 'V';                                          // Compute eigenvalues and eigenvectors
    const char UPLO = 'L';                                          // Lower triangle of H is stored
    dim_t INFO;                                                     // Status of zhpev_/LAPACKE_zhpev
    cplx_t* z;                                                      // Orthonormal eigenvectors of H, with the i-th
    if ((z = malloc(n * n * sizeof (cplx_t))) == NULL) {        // column holding the eigenvector associated
        fprintf(stderr, "zhpev: z allocation failed\n");            // with w[i]
        exit(EXIT_FAILURE);
    }

#ifdef MACOS
    cplx_t work[2 * n - 1];
    double rwork[3 * n - 2];

    zhpev_(&JOBZ, &UPLO, &n, *h, w, z, &n, work, rwork, &INFO);
#else
    INFO = LAPACKE_zhpev(LAPACK_COL_MAJ, JOBZ, UPLO, n, *h, w, n);
#endif
    if (INFO == 0) {
        free(*h);
        *h = z;
        return 0;
    }
    if (INFO < 0) {
        fprintf(stderr, "zhpev: the %ld-th argument had an illegal value.", -INFO);
        return 1;
    }
    else {
        fprintf(stderr, "zhpev: the algorithm failed to converge; "
                        "%ld off-diagonal elements of an intermediate tridiagonal form did not converge to zero.",
                        INFO);
        return 1;
    }
}

void zlowertofull(dim_t n, const cplx_t* lower, cplx_t* full) {
    // Iterate the columns
    dim_t offset = 0;
    for (dim_t j = 0; j < n; ++j) {
        full[j + j * n] = lower[offset];
        for (dim_t i = j + 1; i < n; ++i) {
            full[i + j * n] = lower[offset + (i - j)];
            full[j + i * n] = conj(full[i + j * n]);
        }
        offset += n - j;
    }
}