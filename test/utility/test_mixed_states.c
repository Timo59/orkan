// test_mixed_states.c - Create and destroy density matrices for testing quantum circuits

/*
 * =====================================================================================================================
 * includes
 * =====================================================================================================================
 */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#ifndef Q_TEST_H
#include "test.h"
#endif

#ifndef TEST_MIXED_STATES_H
#include "test_mixed_states.h"
#endif

#ifndef TEST_PURE_STATES_H
#include "test_pure_states.h"
#endif

/*
 * =====================================================================================================================
 * Helper functions
 * =====================================================================================================================
 */

/*
 * @brief   Calculate the index in packed lower-triangular storage for element (row, col)
 *
 * @param[in]   n       Matrix dimension (order)
 * @param[in]   row     Row index (0-indexed, must be >= col)
 * @param[in]   col     Column index (0-indexed)
 *
 * @returns Index in the packed array (column-major order, LAPACK convention)
 *
 * @note    For Hermitian matrix stored in lower-triangular packed format (UPLO='L')
 *          Elements are stored column by column:
 *          Column 0: A(0,0), A(1,0), ..., A(n-1,0)
 *          Column 1: A(1,1), A(2,1), ..., A(n-1,1)
 *          etc.
 *          Element A(i,j) with i >= j is at position: i - j + j*(2n - j + 1)/2
 */
static inline unsigned packed_index(unsigned n, unsigned row, unsigned col) {
    assert(row >= col && "packed_index(): row must be >= col for lower-triangular storage");
    assert(row < n && col < n && "packed_index(): indices out of bounds");
    return row - col + (col * (2 * n - col + 1)) / 2;
}


/*
 * @brief   Convert a pure state vector to density matrix ρ = |ψ⟩⟨ψ| in packed storage
 *
 * @param[in]   psi     Pure state vector (size: dim)
 * @param[in]   dim     Hilbert space dimension (2^nqubits)
 *
 * @returns Density matrix in packed lower-triangular storage, or NULL on failure
 */
static cplx_t* pure_to_density(const cplx_t* psi, unsigned dim) {
    if (!psi) {
        fprintf(stderr, "pure_to_density(): psi is NULL\n");
        return NULL;
    }

    const unsigned packed_size = (dim * (dim + 1)) / 2;
    cplx_t* rho = NULL;

    // Allocate packed density matrix
    if (!((rho = calloc(packed_size, sizeof(*rho))))) {
        fprintf(stderr, "pure_to_density(): rho allocation failed\n");
        return NULL;
    }

    // Compute density matrix elements: ρ(i,j) = ψ(i) * conj(ψ(j))
    // Store only lower triangle: i >= j
    for (unsigned i = 0; i < dim; ++i) {
        for (unsigned j = 0; j <= i; ++j) {
            unsigned idx = packed_index(dim, i, j);
            // ρ(i,j) = ψ[i] * conj(ψ[j])
            double psi_i_re = creal(psi[i]);
            double psi_i_im = cimag(psi[i]);
            double psi_j_re = creal(psi[j]);
            double psi_j_im = cimag(psi[j]);

            // Complex conjugate: (a + ib)(c - id) = (ac + bd) + i(bc - ad)
            double re = psi_i_re * psi_j_re + psi_i_im * psi_j_im;
            double im = psi_i_im * psi_j_re - psi_i_re * psi_j_im;

            rho[idx] = re + I * im;
        }
    }

    return rho;
}


/*
 * =====================================================================================================================
 * Function definitions
 * =====================================================================================================================
 */

cplx_t** test_cb_mixed(const unsigned nqubits) {
    cplx_t** out = NULL;
    cplx_t** pure_states = NULL;

    const unsigned dim = 1 << nqubits;  // Hilbert space dimension

    // Get pure computational basis states
    if (!((pure_states = test_cb_pure(nqubits)))) {
        fprintf(stderr, "test_cb_mixed(): test_cb_pure() failed\n");
        return NULL;
    }

    // Allocate array of pointers to density matrices
    if (!((out = calloc(dim, sizeof(*out))))) {
        fprintf(stderr, "test_cb_mixed(): out allocation failed\n");
        goto cleanup_pure;
    }

    // Convert each pure state to density matrix
    for (unsigned i = 0; i < dim; ++i) {
        if (!((out[i] = pure_to_density(pure_states[i], dim)))) {
            fprintf(stderr, "test_cb_mixed(): conversion of state %u failed\n", i);
            goto cleanup;
        }
    }

    // Free pure states (no longer needed)
    for (unsigned i = 0; i < dim; ++i) {
        free(pure_states[i]);
    }
    free(pure_states);

    return out;

cleanup:
    if (out) {
        for (unsigned i = 0; i < dim; ++i) {
            free(out[i]);
        }
        free(out);
    }
cleanup_pure:
    if (pure_states) {
        for (unsigned i = 0; i < dim; ++i) {
            free(pure_states[i]);
        }
        free(pure_states);
    }
    return NULL;
}


cplx_t** test_xb_mixed(const unsigned nqubits) {
    cplx_t** out = NULL;
    cplx_t** pure_states = NULL;

    const unsigned dim = 1 << nqubits;

    // Get pure Hadamard basis states
    if (!((pure_states = test_xb_pure(nqubits)))) {
        fprintf(stderr, "test_xb_mixed(): test_xb_pure() failed\n");
        return NULL;
    }

    // Allocate array of pointers to density matrices
    if (!((out = calloc(dim, sizeof(*out))))) {
        fprintf(stderr, "test_xb_mixed(): out allocation failed\n");
        goto cleanup_pure;
    }

    // Convert each pure state to density matrix
    for (unsigned i = 0; i < dim; ++i) {
        if (!((out[i] = pure_to_density(pure_states[i], dim)))) {
            fprintf(stderr, "test_xb_mixed(): conversion of state %u failed\n", i);
            goto cleanup;
        }
    }

    // Free pure states
    for (unsigned i = 0; i < dim; ++i) {
        free(pure_states[i]);
    }
    free(pure_states);

    return out;

cleanup:
    if (out) {
        for (unsigned i = 0; i < dim; ++i) {
            free(out[i]);
        }
        free(out);
    }
cleanup_pure:
    if (pure_states) {
        for (unsigned i = 0; i < dim; ++i) {
            free(pure_states[i]);
        }
        free(pure_states);
    }
    return NULL;
}


cplx_t** test_yb_mixed(const unsigned nqubits) {
    cplx_t** out = NULL;
    cplx_t** pure_states = NULL;

    const unsigned dim = 1 << nqubits;

    // Get pure circular basis states
    if (!((pure_states = test_yb_pure(nqubits)))) {
        fprintf(stderr, "test_yb_mixed(): test_yb_pure() failed\n");
        return NULL;
    }

    // Allocate array of pointers to density matrices
    if (!((out = calloc(dim, sizeof(*out))))) {
        fprintf(stderr, "test_yb_mixed(): out allocation failed\n");
        goto cleanup_pure;
    }

    // Convert each pure state to density matrix
    for (unsigned i = 0; i < dim; ++i) {
        if (!((out[i] = pure_to_density(pure_states[i], dim)))) {
            fprintf(stderr, "test_yb_mixed(): conversion of state %u failed\n", i);
            goto cleanup;
        }
    }

    // Free pure states
    for (unsigned i = 0; i < dim; ++i) {
        free(pure_states[i]);
    }
    free(pure_states);

    return out;

cleanup:
    if (out) {
        for (unsigned i = 0; i < dim; ++i) {
            free(out[i]);
        }
        free(out);
    }
cleanup_pure:
    if (pure_states) {
        for (unsigned i = 0; i < dim; ++i) {
            free(pure_states[i]);
        }
        free(pure_states);
    }
    return NULL;
}


cplx_t** test_bell_mixed(const unsigned nqubits) {
    cplx_t** out = NULL;
    cplx_t** pure_states = NULL;

    // Bell states only defined for n=2 qubits
    if (nqubits != 2) {
        fprintf(stderr, "test_bell_mixed(): Bell states only defined for 2 qubits\n");
        return NULL;
    }

    const unsigned dim = 4;  // 2 qubits

    // Get pure Bell states
    if (!((pure_states = test_bell_pure(nqubits)))) {
        fprintf(stderr, "test_bell_mixed(): test_bell_pure() failed\n");
        return NULL;
    }

    // Allocate array of pointers to density matrices (4 Bell states)
    if (!((out = calloc(4, sizeof(*out))))) {
        fprintf(stderr, "test_bell_mixed(): out allocation failed\n");
        goto cleanup_pure;
    }

    // Convert each Bell state to density matrix
    for (unsigned i = 0; i < 4; ++i) {
        if (!((out[i] = pure_to_density(pure_states[i], dim)))) {
            fprintf(stderr, "test_bell_mixed(): conversion of state %u failed\n", i);
            goto cleanup;
        }
    }

    // Free pure states
    for (unsigned i = 0; i < 4; ++i) {
        free(pure_states[i]);
    }
    free(pure_states);

    return out;

cleanup:
    if (out) {
        for (unsigned i = 0; i < 4; ++i) {
            free(out[i]);
        }
        free(out);
    }
cleanup_pure:
    if (pure_states) {
        for (unsigned i = 0; i < 4; ++i) {
            free(pure_states[i]);
        }
        free(pure_states);
    }
    return NULL;
}


cplx_t** test_ghz_mixed(const unsigned nqubits) {
    cplx_t** out = NULL;
    cplx_t** pure_state = NULL;

    if (nqubits < 2) {
        fprintf(stderr, "test_ghz_mixed(): GHZ states require at least 2 qubits\n");
        return NULL;
    }

    const unsigned dim = 1 << nqubits;

    // Get pure GHZ state
    if (!((pure_state = test_ghz_pure(nqubits)))) {
        fprintf(stderr, "test_ghz_mixed(): test_ghz_pure() failed\n");
        return NULL;
    }

    // Allocate array for single density matrix
    if (!((out = calloc(1, sizeof(*out))))) {
        fprintf(stderr, "test_ghz_mixed(): out allocation failed\n");
        free(pure_state[0]);
        free(pure_state);
        return NULL;
    }

    // Convert to density matrix
    if (!((out[0] = pure_to_density(pure_state[0], dim)))) {
        fprintf(stderr, "test_ghz_mixed(): conversion failed\n");
        free(out);
        free(pure_state[0]);
        free(pure_state);
        return NULL;
    }

    // Free pure state
    free(pure_state[0]);
    free(pure_state);

    return out;
}


cplx_t** test_w_mixed(const unsigned nqubits) {
    cplx_t** out = NULL;
    cplx_t** pure_state = NULL;

    if (nqubits < 3) {
        fprintf(stderr, "test_w_mixed(): W states require at least 3 qubits\n");
        return NULL;
    }

    const unsigned dim = 1 << nqubits;

    // Get pure W state
    if (!((pure_state = test_w_pure(nqubits)))) {
        fprintf(stderr, "test_w_mixed(): test_w_pure() failed\n");
        return NULL;
    }

    // Allocate array for single density matrix
    if (!((out = calloc(1, sizeof(*out))))) {
        fprintf(stderr, "test_w_mixed(): out allocation failed\n");
        free(pure_state[0]);
        free(pure_state);
        return NULL;
    }

    // Convert to density matrix
    if (!((out[0] = pure_to_density(pure_state[0], dim)))) {
        fprintf(stderr, "test_w_mixed(): conversion failed\n");
        free(out);
        free(pure_state[0]);
        free(pure_state);
        return NULL;
    }

    // Free pure state
    free(pure_state[0]);
    free(pure_state);

    return out;
}


cplx_t** test_maximally_mixed(const unsigned nqubits) {
    cplx_t** out = NULL;

    const unsigned dim = 1 << nqubits;  // Hilbert space dimension: 2^n
    const unsigned packed_size = (dim * (dim + 1)) / 2;
    const double diag_val = 1.0 / (double) dim;  // I/2^n has diagonal entries 1/2^n

    // Allocate array of one density matrix pointer
    if (!((out = calloc(1, sizeof(*out))))) {
        fprintf(stderr, "test_maximally_mixed(): out allocation failed\n");
        return NULL;
    }

    // Allocate packed density matrix (initialized to zero by calloc)
    if (!((out[0] = calloc(packed_size, sizeof(*out[0]))))) {
        fprintf(stderr, "test_maximally_mixed(): out[0] allocation failed\n");
        free(out);
        return NULL;
    }

    // Set diagonal elements: ρ(i,i) = 1/2^n
    for (unsigned i = 0; i < dim; ++i) {
        unsigned idx = packed_index(dim, i, i);
        out[0][idx] = diag_val + I * 0.0;
    }

    return out;
}


cplx_t** test_random_mixture(const unsigned nqubits) {
    cplx_t** out = NULL;
    cplx_t** pure_states = NULL;
    double* probabilities = NULL;
    unsigned num_pure = 0;

    // Static counter to generate different mixtures on each call (reproducible sequence)
    static unsigned call_count = 0;

    // Hardcoded parameters for random mixture generation
    const unsigned dim = 1 << nqubits;
    const unsigned max_components = (3 * dim < 10) ? 3 * dim : 10;
    const unsigned num_components = 2 + (call_count % (max_components - 1));  // Varies with call
    const unsigned seed = 1000 + call_count;  // Deterministic seed based on call count

    call_count++;  // Increment for next call

    const unsigned packed_size = (dim * (dim + 1)) / 2;

    // Generate pure basis states to mix
    if (!((pure_states = test_mk_states_pure(nqubits, &num_pure)))) {
        fprintf(stderr, "test_random_mixture(): test_mk_states_pure failed\n");
        return NULL;
    }

    if (num_components > num_pure) {
        fprintf(stderr, "test_random_mixture(): num_components (%u) exceeds available pure states (%u)\n",
                num_components, num_pure);
        goto cleanup_pure;
    }

    // Allocate array for single density matrix
    if (!((out = calloc(1, sizeof(*out))))) {
        fprintf(stderr, "test_random_mixture(): out allocation failed\n");
        goto cleanup_pure;
    }

    // Allocate density matrix (initialized to zero)
    if (!((out[0] = calloc(packed_size, sizeof(*out[0]))))) {
        fprintf(stderr, "test_random_mixture(): rho allocation failed\n");
        free(out);
        out = NULL;
        goto cleanup_pure;
    }

    // Generate random probabilities that sum to 1
    if (!((probabilities = malloc(num_components * sizeof(*probabilities))))) {
        fprintf(stderr, "test_random_mixture(): probabilities allocation failed\n");
        goto cleanup;
    }

    // Initialize random number generator
    srand(seed);

    // Generate random positive numbers
    double sum = 0.0;
    for (unsigned k = 0; k < num_components; ++k) {
        probabilities[k] = (double)rand() / (double)RAND_MAX;
        sum += probabilities[k];
    }

    // Normalize to sum to 1
    for (unsigned k = 0; k < num_components; ++k) {
        probabilities[k] /= sum;
    }

    // Compute mixture: ρ = Σₖ pₖ |ψₖ⟩⟨ψₖ|
    // Use first num_components pure states
    for (unsigned k = 0; k < num_components; ++k) {
        double pk = probabilities[k];

        // Add pk * |ψₖ⟩⟨ψₖ| to ρ
        for (unsigned i = 0; i < dim; ++i) {
            for (unsigned j = 0; j <= i; ++j) {
                unsigned idx = packed_index(dim, i, j);

                // ρ(i,j) += pk * ψₖ[i] * conj(ψₖ[j])
                double psi_i_re = creal(pure_states[k][i]);
                double psi_i_im = cimag(pure_states[k][i]);
                double psi_j_re = creal(pure_states[k][j]);
                double psi_j_im = cimag(pure_states[k][j]);

                // Complex conjugate multiplication
                double re = pk * (psi_i_re * psi_j_re + psi_i_im * psi_j_im);
                double im = pk * (psi_i_im * psi_j_re - psi_i_re * psi_j_im);

                out[0][idx] += re + I * im;
            }
        }
    }

    free(probabilities);
    test_rm_states_pure(nqubits, pure_states);
    pure_states = NULL;

    return out;

cleanup:
    if (out) {
        free(out[0]);
        free(out);
    }
cleanup_pure:
    if (probabilities) free(probabilities);
    if (pure_states) {
        test_rm_states_pure(nqubits, pure_states);
        free(pure_states);
    }
    return NULL;
}


cplx_t** test_mk_states_mixed(const unsigned nqubits, unsigned* nmixed) {
    cplx_t** out = NULL;
    unsigned current_idx = 0;

    const unsigned dim = 1 << nqubits;  // Hilbert space dimension

    // Calculate total number of mixed states:
    // - Computational basis: dim states
    // - Hadamard basis: dim states
    // - Circular basis: dim states
    // - 1 maximally mixed state
    // - 10 random mixtures (as specified in GATES_MODULE.md)
    // - Entangled states: 4 Bell states (n=2) or 2 states (GHZ + W for n>2)
    const unsigned num_random = 10;
    *nmixed = 3 * dim + 1 + num_random;
    if (nqubits == 2) *nmixed += 4;  // Bell states
    else if (nqubits > 2) *nmixed += 2;  // GHZ + W states

    // Allocate array of pointers
    if (!((out = calloc(*nmixed, sizeof(*out))))) {
        fprintf(stderr, "test_mk_states_mixed(): out allocation failed\n");
        return NULL;
    }

    // 1. Add computational basis states as density matrices
    cplx_t** cb_mixed = test_cb_mixed(nqubits);
    if (!cb_mixed) {
        fprintf(stderr, "test_mk_states_mixed(): test_cb_mixed failed\n");
        goto cleanup;
    }
    for (unsigned i = 0; i < dim; ++i) {
        out[current_idx++] = cb_mixed[i];
    }
    free(cb_mixed);  // Free array of pointers, not the density matrices

    // 2. Add Hadamard basis states as density matrices
    cplx_t** xb_mixed = test_xb_mixed(nqubits);
    if (!xb_mixed) {
        fprintf(stderr, "test_mk_states_mixed(): test_xb_mixed failed\n");
        goto cleanup;
    }
    for (unsigned i = 0; i < dim; ++i) {
        out[current_idx++] = xb_mixed[i];
    }
    free(xb_mixed);

    // 3. Add circular basis states as density matrices
    cplx_t** yb_mixed = test_yb_mixed(nqubits);
    if (!yb_mixed) {
        fprintf(stderr, "test_mk_states_mixed(): test_yb_mixed failed\n");
        goto cleanup;
    }
    for (unsigned i = 0; i < dim; ++i) {
        out[current_idx++] = yb_mixed[i];
    }
    free(yb_mixed);

    // 4. Add maximally mixed state
    cplx_t** max_mixed = test_maximally_mixed(nqubits);
    if (!max_mixed) {
        fprintf(stderr, "test_mk_states_mixed(): test_maximally_mixed failed\n");
        goto cleanup;
    }
    out[current_idx++] = max_mixed[0];
    free(max_mixed);

    // 5. Add random mixtures
    for (unsigned i = 0; i < num_random; ++i) {
        cplx_t** random_mix = test_random_mixture(nqubits);
        if (!random_mix) {
            fprintf(stderr, "test_mk_states_mixed(): test_random_mixture %u failed\n", i);
            goto cleanup;
        }
        out[current_idx++] = random_mix[0];
        free(random_mix);
    }

    // 6. Add entangled states (Bell, GHZ, W)
    if (nqubits == 2) {
        // Add Bell states
        cplx_t** bell_mixed = test_bell_mixed(nqubits);
        if (!bell_mixed) {
            fprintf(stderr, "test_mk_states_mixed(): test_bell_mixed failed\n");
            goto cleanup;
        }
        for (unsigned i = 0; i < 4; ++i) {
            out[current_idx++] = bell_mixed[i];
        }
        free(bell_mixed);
    } else if (nqubits > 2) {
        // Add GHZ state
        cplx_t** ghz_mixed = test_ghz_mixed(nqubits);
        if (!ghz_mixed) {
            fprintf(stderr, "test_mk_states_mixed(): test_ghz_mixed failed\n");
            goto cleanup;
        }
        out[current_idx++] = ghz_mixed[0];
        free(ghz_mixed);

        // Add W state
        cplx_t** w_mixed = test_w_mixed(nqubits);
        if (!w_mixed) {
            fprintf(stderr, "test_mk_states_mixed(): test_w_mixed failed\n");
            goto cleanup;
        }
        out[current_idx++] = w_mixed[0];
        free(w_mixed);
    }

    return out;

cleanup:
    if (out) {
        for (unsigned i = 0; i < current_idx; ++i) {
            free(out[i]);
        }
        free(out);
    }
    return NULL;
}


void test_rm_states_mixed(cplx_t** states, unsigned nmixed) {
    if (states) {
        for (unsigned i = 0; i < nmixed; ++i) {
            free(states[i]);
            states[i] = NULL;
        }
    }
    free(states);
}
