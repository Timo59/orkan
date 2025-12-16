// test_qhipster.c - Unit tests for representations of quantum gates on pure states

#ifndef TEST_QHIPSTER_H
#include "test_qhipster.h"
#endif

#include <stdlib.h>

#if defined(__APPLE__)
#include <vecLib/cblas_new.h>
#elif defined(__linux__)
#include <cblas.h>
#endif

/*
 * =====================================================================================================================
 * Unity: setUp and tearDown
 * =====================================================================================================================
 */

void setUp(void) {}

void tearDown(void) {}

void test_apply_x(void) {testSingleQubitGate(applyX, XMAT);}
void test_apply_y(void) {testSingleQubitGate(applyY, YMAT);}
void test_apply_z(void) {testSingleQubitGate(applyZ, ZMAT);}

/*
 * =====================================================================================================================
 *                                                      main
 * =====================================================================================================================
 */

int main(void) {
    UNITY_BEGIN();
    RUN_TEST(test_apply_x);
    RUN_TEST(test_apply_y);
    RUN_TEST(test_apply_z);
    return UNITY_END();
}

/*
 * =====================================================================================================================
 * Single qubit gates
 * =====================================================================================================================
 */

void testSingleQubitGate(const single_qubit_gate gate, const cplx_t *mat) {
    state_t test_state = {0}; // Pure test state
    test_state.type = PURE;
    cplx_t **ref_vecs = {0}, **test_vecs = {0}; // Reference and test state vectors
    unsigned nvecs = 0;    // Number of test state vectors
    cplx_t *gateMat = NULL; // Matrix representation of the quantum gate

    // Iterate the number of qubits
    for (unsigned nqubits= 1; nqubits <= MAXQUBITS; ++nqubits) {
        const unsigned dim = POW2(nqubits, dim_t);

        // Generate reference state vectors (not touched by the functions)
        if (!((ref_vecs = test_gen_states_pure(nqubits, &nvecs)))) {
            fprintf(stderr, "testSingleQubitGate(): ref_vecs initialization failed\n");
            goto cleanup;
        }

        // Iterate target qubits
        for (unsigned pos = 0; pos < nqubits; ++pos) {
            // Generate test state vectors
            if (!((test_vecs = test_gen_states_pure(nqubits, &nvecs)))) {
                fprintf(stderr, "testSingleQubitGate(): test_vecs initialization failed\n");
                goto cleanup;
            }

            // Initialize matrix with a Pauli-X acting on the target qubit
            if (!((gateMat = mat_single_qubit_gate(nqubits, mat, pos)))) goto cleanup;

            // Iterate the test state vectors
            for (unsigned i = 0; i < nvecs; ++i) {
                // Initialize test state with i-th state vector
                state_init(&test_state, nqubits, test_vecs[i]);
                if (!test_state.data) {
                    fprintf(stderr, "testSingleQubitGate(): test state data initialization failed\n");
                    goto cleanup;
                }

                // Apply test function
                gate(&test_state, pos);
                if (!test_state.data) {
                    fprintf(stderr, "testSingleQubitGate(): gate application failed\n");
                    goto cleanup;
                }

                // Reference state vector by matrix-vector multiplication
                cplx_t* ref = mv(dim, gateMat, ref_vecs[i]);
                if (!ref) {
                    fprintf(stderr, "testSingleQubitGate(): matrix-vector multiplication failed\n");
                    goto cleanup;
                }

                // Compare state vectors
                TEST_ASSERT_EQUAL_COMPLEX_ARRAY_TOL(ref, test_state.data, dim, PRECISION);

                // Free reference vector
                free(ref);
                ref = NULL;
            }

            // Free matrix representation
            free(gateMat);
            gateMat = NULL;

            // Free test state vectors
            test_rm_states_pure(nqubits, test_vecs);
        }

        // Free reference state vectors
        test_rm_states_pure(nqubits, ref_vecs);
    }

    cleanup:
        // Free reference state vectors
        if (ref_vecs) {
            for (unsigned i = 0; i < nvecs; ++i) {
                free(ref_vecs[i]);
                ref_vecs[i] = NULL;
            }
        }
        free(ref_vecs);
        ref_vecs = NULL;

        // Free test state vectors
        if (test_vecs) {
            for (unsigned i = 0; i < nvecs; ++i) {
                free(test_vecs[i]);
                test_vecs[i] = NULL;
            }
        }
        free(test_vecs);
        test_vecs = NULL;

        // Free gate matrix representation
        free(gateMat);
        gateMat = NULL;

        // Free test state
        state_free(&test_state);
}


/*
 * =====================================================================================================================
 * Linear algebra
 * =====================================================================================================================
 */

cplx_t* mv(const unsigned n, const cplx_t *m, const cplx_t *v) {
    cplx_t *out = NULL;

    // Allocate memory for the output
    if (!((out = malloc(n * sizeof (*out))))) {
        fprintf(stderr, "mv(): out allocation failed\n");
        return out;
    }

    // Define multipliers: w = ALPHA*M.v + BETA*w
    const dim_t N = (dim_t) n;
    const cplx_t ALPHA = 1.0 + I*0.0;
    const cplx_t BETA = 0.0 + I*0.0;
    cblas_zgemv(CblasColMajor, CblasNoTrans, N, N, &ALPHA, m, N, v, 1, &BETA, out, 1);

    return out;
}


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
                    const unsigned col = b_col + a_col * n;
                    const unsigned row = b_row + a_row * m;
                    out[row + col * k * m] = a * B[b_row + b_col * m];
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

cplx_t* mat_single_qubit_gate(const unsigned nqubits, const cplx_t *gate, const unsigned target) {
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
    if (!((id = mat_id(target)))) goto cleanup;

    // Kronecker product with right identity
    const unsigned dim_right = 1ULL << target;
    out = kron(dim_left, dim_left, tmp, dim_right, dim_right, id);

    cleanup:
        free(id);
        id = NULL;
        free(tmp);
        tmp = NULL;

    return out;
}
