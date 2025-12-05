// test_qhipster.c - Unit tests for functions from qhipster.c

#ifndef Q_TEST_QHIPSTER_H
#include "test_qhipster.h"
#endif

#ifndef _STDLIB_H_
#include <stdlib.h>
#endif

/*
 * =====================================================================================================================
 * Unity: setUp and tearDown
 * =====================================================================================================================
 */

void setUp(void) {}

void tearDown(void) {}

/*
 * =====================================================================================================================
 * Pauli gates
 * =====================================================================================================================
 */

/*
 * @brief   Unit test of the Pauli-X gate on each qubit of a multi-qubit state of up MAXQUBITS qubits
 */
void testApplyX(void) {
    for (unsigned nqubits= 1; nqubits <= MAXQUBITS; ++nqubits) {


        for (qubit_t pos = 0; pos < qubits; ++pos) {                    // For each target qubit from the right:
            cplx_t* gateMat = singleQubitGateMat(XMAT, qubits, pos);    // generate the matrix representation from the
            // Kronecker product
            for (dim_t i = 0; i < dim + 1; ++i) {                       // For each state vector:
                stateInitVector(&testState, vecs[i]);                   // Copy state vector to testState's vector,
                applyX(&testState, pos);                                // apply test function
                cplx_t* ref = cmatVecMul(gateMat, vecs[i], dim);    // Matrix multiplication as reference
                TEST_ASSERT_TRUE(cvectorAlmostEqual(ref, testState.vec, dim, PRECISION));
                free(ref);
            }
            free(gateMat);
        }
        stateFreeVector(&testState);
        freeTestVectors(vecs, qubits);
    }
}


/*
 * =====================================================================================================================
 * Linear algebra
 * =====================================================================================================================
 */

cplx_t* kron(const unsigned k, const unsigned l, const cplx_t A[],
             const unsigned m, const unsigned n, const cplx_t B[]) {
    // Allocate memory for the Kronecker product matrix
    cplx_t *out = malloc(k * l * m * n * sizeof(*out));
    if (!out) {
        fprintf(stderr, "kron(): out allocation failed\n");
        return out;
    }

    // Iterate A's columns
    for (unsigned a_col = 0; a_col < l; ++a_col) {
        // Iterate A's rows
        for (unsigned a_row = 0; a_row < k; ++a_row) {
            const cplx_t a = A[a_row + a_col * k];  // (a_row,a_col)-th element of A

            // Iterate B's columns
            for (unsigned b_col = 0; b_col < n; ++b_col) {
                // Iterate B's rows
                for (unsigned b_row = 0; b_row < m; ++b_row) {
                    // b_row is the row in the a_row-th block of rows with each block comprising m rows
                    // b_col is the column in the a_col-th block of columns with each block comprising n columns
                    out[b_row + (a_row + (b_col + a_col * n) + k) * m] = a * B[b_row + b_col * m];
                }   // b_row
            }   // b_col
        }   // a_row
    }   // a_col

    return out;
}

/*
 * =====================================================================================================================
 * Multiqubit gates
 * =====================================================================================================================
 */

cplx_t* mat_id(const unsigned nqubits) {
    // Allocate memory for identity matrix and initialize all entries to zero
    const unsigned dim = 1ULL << nqubits;
    cplx_t *out = calloc(dim * dim, sizeof(*out));
    if (!out) {
        fprintf(stderr, "mat_id(): out allocation failed\n");
        return out;
    }

    // Populate the diagonal with ones
    for (unsigned i = 0; i < dim; ++i) {
        out[i + i * dim] = 1.0 + I*0.0;
    }

    return out;
}

cplx_t* mat_single_qubit_gate(const unsigned nqubits, const cplx_t gate[], const unsigned target) {
    // Check that target is in scope
    if (target >= nqubits) {
        fprintf(stderr, "mat_single_qubit_gate(): target is out of scope; Expected %u, Was %u\n", nqubits, target);
        return NULL;
    }

    // Adjusted position from the left
    const unsigned pos = nqubits - 1 - target;
    unsigned dim_left = 1ULL << pos;    // Hilbert space dimension of qubits to the left of target

    // Identity, temporary Kronecker product and output pointers
    cplx_t *id = NULL, *tmp = NULL, *out = NULL;

    // Identity matrix for all qubits to the left of target
    if (!((id = mat_id(pos)))) goto cleanup;

    // Kronecker product of left identity with single qubit gate matrix
    if (!((tmp = kron(dim_left, dim_left, id, 2, 2, gate)))) goto cleanup;
    dim_left *= 2;  // Hilbert space dimension of qubits to the left of and including target
    free(id);
    id = NULL;

    // Identity matrix for all qubits to the right of target
    if (!((id = mat_id(target - 1)))) goto cleanup;

    // Kronecker product with right identity
    const unsigned dim_right = 1ULL << (target - 1);
    out = kron(dim_left, dim_left, tmp, dim_right, dim_right, id);

    cleanup:
        free(id);
        id = NULL;
        free(tmp);
        tmp = NULL;

    return out;
}
