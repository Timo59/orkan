// test_mhipster.c  - Test harness for quantum gates on mixed states. The reference is given by matrix-matrix
//                    multiplication of the density matrix with the matrix representation of the quantum gate.

/*
 * =====================================================================================================================
 * includes
 * =====================================================================================================================
 */

#ifndef TEST_GATE_H
#include "test_gate.h"
#endif

#ifndef GATEMAT_H
#include "gatemat.h"
#endif

#ifndef LINALG_H
#include "linalg.h"
#endif

#ifndef TEST_MIXED_STATES_H
#include "test_mixed_states.h"
#endif

#ifndef UTILS_H
#include "utils.h"
#endif

#include <stdlib.h>

/*
 * =====================================================================================================================
 * Packed storage utilities
 * =====================================================================================================================
 */

static cplx_t* density_unpack(const unsigned n, const cplx_t *packed) {
    cplx_t *full = NULL;

    // Allocate memory for full n×n matrix
    if (!((full = malloc(n * n * sizeof(*full))))) {
        fprintf(stderr, "density_unpack(): full allocation failed\n");
        return full;
    }

    // Unpack lower triangle from packed storage (column-major)
    // Packed format: column 0 (n elements), column 1 (n-1 elements), ..., column n-1 (1 element)
    unsigned packed_idx = 0;
    for (unsigned col = 0; col < n; ++col) {
        for (unsigned row = col; row < n; ++row) {
            full[row + col * n] = packed[packed_idx];  // Lower triangle
            if (row != col) {
                // Hermitian property: H[col,row] = conj(H[row,col])
                full[col + row * n] = conj(packed[packed_idx]);  // Upper triangle
            }
            ++packed_idx;
        }
    }

    return full;
}


static cplx_t* density_pack(const unsigned n, const cplx_t *full) {
    cplx_t *packed = NULL;

    // Allocate memory for packed storage
    const unsigned packed_len = n * (n + 1) / 2;
    if (!((packed = malloc(packed_len * sizeof(*packed))))) {
        fprintf(stderr, "density_pack(): packed allocation failed\n");
        return packed;
    }

    // Pack lower triangle into packed storage (column-major)
    unsigned packed_idx = 0;
    for (unsigned col = 0; col < n; ++col) {
        for (unsigned row = col; row < n; ++row) {
            packed[packed_idx] = full[row + col * n];
            ++packed_idx;
        }
    }

    return packed;
}

void multiQubitGate(unsigned n, cplx_t *h, unsigned k, const unsigned *a, const cplx_t *m);

/*
 * =====================================================================================================================
 * Test: Single qubit gates
 * =====================================================================================================================
 */

void testSingleQubitGateMixed(const single_qubit_gate gate, const cplx_t *mat) {
    state_t test_state = {0}; // Mixed test state
    test_state.type = MIXED;
    cplx_t **test_rhos = NULL, *gateMat = NULL;  // Test density matrices, matrix representation of the quantum gate
    cplx_t *rho_full = NULL, *ref_full = NULL;  // Temporary matrices for reference computation
    unsigned nmixed = 0;    // Number of test density matrices

    // Iterate the number of qubits
    for (unsigned nqubits = 2; nqubits <= MAXQUBITS; ++nqubits) {
        const unsigned dim = POW2(nqubits, dim_t);

        // Iterate target qubits
        for (unsigned pos = 0; pos < nqubits; ++pos) {
            // Generate test density matrices
            if (!((test_rhos = test_mk_states_mixed(nqubits, &nmixed)))) {
                fprintf(stderr, "testSingleQubitGateMixed(): test_rhos initialization failed\n");
                goto cleanup;
            }

            // Initialize matrix representation of the gate acting on the target qubit
            if (!((gateMat = mat_single_qubit_gate(nqubits, mat, pos)))) goto cleanup;

            // Iterate the test density matrices
            for (unsigned i = 0; i < nmixed; ++i) {
                // Unpack density matrix from packed storage to full matrix
                rho_full = density_unpack(dim, test_rhos[i]);
                if (!rho_full) {
                    fprintf(stderr, "testSingleQubitGateMixed(): density_unpack failed\n");
                    goto cleanup;
                }

                // Compute reference: ref_full = U * rho * U†
                ref_full = zumu(dim, gateMat, rho_full);
                if (!ref_full) {
                    fprintf(stderr, "testSingleQubitGateMixed(): zumu failed\n");
                    goto cleanup;
                }

                // Pack reference result to packed storage
                cplx_t *ref = density_pack(dim, ref_full);
                if (!ref) {
                    fprintf(stderr, "testSingleQubitGateMixed(): density_pack failed\n");
                    goto cleanup;
                }

                // Compute reference from multiQubitGate as well
                multiQubitGate(nqubits, rho_full, 1, &pos, mat);

                // Initialize test state with the density matrix (previous density matrix gets freed inside state_init())
                state_init(&test_state, nqubits, &test_rhos[i]);
                if (!test_state.data) {
                    fprintf(stderr, "testSingleQubitGateMixed(): test state data initialization failed\n");
                    goto cleanup;
                }

                // Apply test function
                qs_error_t err = gate(&test_state, pos);
                TEST_ASSERT_EQUAL_INT(QS_OK, err);
                if (!test_state.data) {
                    fprintf(stderr, "testSingleQubitGateMixed(): gate application failed\n");
                    goto cleanup;
                }

                // Compare density matrices in packed format
                const unsigned packed_len = dim * (dim + 1) / 2;
                TEST_ASSERT_EQUAL_COMPLEX_ARRAY_TOL(ref, test_state.data, packed_len, PRECISION);
                TEST_ASSERT_EQUAL_COMPLEX_ARRAY_TOL(ref_full, rho_full, dim * dim, PRECISION);

                // Free temporary matrices
                free(ref);
                ref = NULL;
                free(rho_full);
                rho_full = NULL;
                free(ref_full);
                ref_full = NULL;

                // Free test state data
                free(test_state.data);
                test_state.data = NULL;
            }

            // Free matrix representation
            free(gateMat);
            gateMat = NULL;

            // Free test density matrix array
            free(test_rhos);
            test_rhos = NULL;
        }
    }

    return;

    cleanup:
        // Free test density matrices
        if (test_rhos) {
            for (unsigned i = 0; i < nmixed; ++i) {
                free(test_rhos[i]);
                test_rhos[i] = NULL;
            }
        }
        free(test_rhos);
        test_rhos = NULL;

        // Free gate matrix representation
        free(gateMat);
        gateMat = NULL;

        // Free temporary matrices
        free(rho_full);
        rho_full = NULL;
        free(ref_full);
        ref_full = NULL;

        // Free test state
        state_free(&test_state);
}

/*
 * =====================================================================================================================
 * TEST
 * =====================================================================================================================
 */

/*
 * @brief   Returns the integer obtained by inserting the specified bits at the corresponding positions of the input.
 *
 * @param   i       Input integer
 * @param   n       Number of bits to insert
 * @param   bits    Integer whose rightmost n bits define the bits to insert
 * @param   pos     Ordered (smallest to largest) list of positions at which the bits are inserted (little endian)
 */
static unsigned insertBits(unsigned i, const unsigned n, unsigned bits, const unsigned *pos) {
    // Iterate bits to be inserted
    for (unsigned j = 0; j < n; ++j) {
        // Mask of bits in i right to the current position
        const unsigned right = ((1u << pos[j]) - 1) & i;

        // Copy value of the j-th bit to rightmost bit of a bit set mask
        unsigned set = (bits >> j) & 1;

        // Set all bits right to current position to zero and the bit at position to target bit
        unsigned left = (i >> pos[j]) << 1;
        left |= set;
        left <<= pos[j];

        i = left | right;
    }

    return i;
}

/*
 * @brief   Applies multi-qubit operator on a full hermitian operator in place
 *
 * @param   n   Number of qubits
 * @param   h   Array of size 2**n * 2**n; concatenation of hermitain operator's columns
 * @param   k   Number of qubits the gate acts on non-trivially
 * @param   a   Array of size k; sorted list of indices of qubits acted upon non-trivially (little endian)
 * @param   m   Array of size 2**k * 2**k; of k-local operator in column-major format
 */
void multiQubitGate(const unsigned n, cplx_t *h, const unsigned k, const unsigned *a, const cplx_t *m) {
    const unsigned N = 1 << n;  // Order of the hermitian matrix
    const unsigned K = 1 << k;  // Invariant subspace dimension; number of rows and columns intertwined by the gate
    const unsigned mat_size = K * K;    // Number of entries in each invariant block matrix
    const unsigned M = 1 << (n - k);    // Number of invariant subspaces

    // Iterate the invariant column blocks
    for (unsigned col_block = 0; col_block < M; ++col_block) {

        // Iterate the invariant row blocks
        for (unsigned row_block = 0; row_block < M; ++row_block) {
            cplx_t tmp[mat_size];

            // Fetch elements from invariant matrix blocks to the stack
            for (unsigned j = 0; j < K; ++j) {
                unsigned col = insertBits(col_block, k, j, a);

                for (unsigned i = 0; i < K; ++i) {
                    unsigned row = insertBits(row_block, k, i, a);
                    tmp[i + j * K] = h[row + col * N];
                    h[row + col * N] = 0;
                }
            }

            // Add elements from of conjugation M->M*tmp*M**H of stack matrix to invariant blocks
            for (unsigned j = 0; j < K; ++j) {
                unsigned col = insertBits(col_block, k, j, a);

                for (unsigned i = 0; i < K; ++i) {
                    unsigned row = insertBits(row_block, k, i, a);

                    // H_ij = M_ik tmp_kl (M**H)_lj
                    for (unsigned s = 0; s < K; ++s) {
                        for (unsigned t = 0; t < K; ++t) {
                            h[row + col * N] += m[i + s * K] * tmp[s + t * K] * conj(m[j + t * K]);
                        }
                    }
                }
            }
        }
    }
}