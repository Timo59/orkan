// test_pure.c  - Create and destroy state vectors for testing quantum circuits

/*
 * =====================================================================================================================
 * includes
 * =====================================================================================================================
 */

#include <math.h>

#include <stdlib.h>

#ifndef Q_TEST_H
#include "test.h"
#endif

/*
 * =====================================================================================================================
 * Function definitions
 * =====================================================================================================================
 */

cplx_t** test_cb_pure(const unsigned nqubits) {
    cplx_t** out = NULL;

    // Allocate the pointers to all computational basis states
    const unsigned dim = 1 << nqubits;  // Hilbert space dimension
    if (!((out = calloc(dim, sizeof (*out))))) {
        fprintf(stderr, "test_cb_pure: out allocation failed\n");
        return out;
    }

    // Initialize the computational basis states
    for (unsigned i = 0; i < dim; ++i) {
        if (!((out[i] = calloc(dim, sizeof (*out[i]))))) {
            fprintf(stderr, "test_cb_pure: out[%u] allocation failed\n", i);
            goto cleanup;
        }

        out[i][i] = 1.0 + I*0.0;
    }

    return out;

    cleanup:
        for (unsigned i = 0; i < dim; ++i) {
            free(out[i]);
            out = NULL;
        }
        free(out);
        out = NULL;

    return out;
}


cplx_t** test_xb_pure(const unsigned nqubits) {
    cplx_t** out = NULL;
    const cplx_t ampl = (cplx_t) pow(INVSQRT2, nqubits);

    // Allocate the pointers to all Hadamard basis states
    const unsigned dim = 1 << nqubits;  // Hilbert space dimension
    if (!((out = calloc(dim, sizeof (*out))))) {
        fprintf(stderr, "test_xb_pure: out allocation failed\n");
        return out;
    }

    // Initialize the Hadamard basis states; the i-th Hadamard basis state is
    // H |i> = H |i_n-1...i_1 i_0> = (|0> + (-1)**(i_n-1) |1>) ... (|0> + (-1)**i_1 |1>) (|0> + (-1)**i_0 |1>)
    for (unsigned i = 0; i < dim; ++i) {
        if (!((out[i] = calloc(dim, sizeof (*out[i]))))) {
            fprintf(stderr, "test_xb_pure: out[%u] allocation failed\n", i);
            goto cleanup;
        }

        // The parity of the j-th computational basis state in H |i> determines the sign of the j-th entry
        for (unsigned j = 0; j < dim; ++j) {
            cplx_t sign = 1.0 + I*0.0;
            // The sign collects a minus for each bit set simultaneously in i and j
            unsigned parity = i & j;
            while (parity > 0) {
                if (parity & 1) sign *= -1.0;
                parity >>= 1;
            }

            out[i][j] = sign * ampl;
        }
    }

    return out;

    cleanup:
        for (unsigned i = 0; i < dim; ++i) {
            free(out[i]);
            out = NULL;
        }
    free(out);
    out = NULL;

    return out;
}


cplx_t** test_yb_pure(const unsigned nqubits) {
    cplx_t** out = NULL;
    const cplx_t ampl = (cplx_t) pow(INVSQRT2, nqubits);

    // Allocate the pointers to all circular basis states
    const unsigned dim = 1 << nqubits;  // Hilbert space dimension
    if (!((out = calloc(dim, sizeof (*out))))) {
        fprintf(stderr, "test_yb_pure: out allocation failed\n");
        return out;
    }

    // Initialize the circular basis states; the j-th circular basis state is
    // H_y |j> = H_y |j_n-1...j_1 j_0> = (|0> + (-1)**(j_n-1) i |1>) ... (|0> + (-1)**j_1 i |1>) (|0> + (-1)**j_0 i |1>)
    for (unsigned i = 0; i < dim; ++i) {
        if (!((out[i] = calloc(dim, sizeof (*out[i]))))) {
            fprintf(stderr, "test_yb_pure: out[%u] allocation failed\n", i);
            goto cleanup;
        }

        // The j-th amplitude gets multiplied by 'i' for each bit set in j and, additionally, by '-1' if that bit is
        // set in i as well
        for (unsigned j = 0; j < dim; ++j) {
            cplx_t sign = 1.0 + I*0.0;

            for (unsigned k = 0; k < nqubits; ++k) {
                if (j & (1 << k)) {
                    sign *= I*1.0;
                    if (i & (1 << k)) sign *= -1.0;
                }
            }

            out[i][j] = sign * ampl;
        }
    }

    return out;

    cleanup:
        for (unsigned i = 0; i < dim; ++i) {
            free(out[i]);
            out = NULL;
        }
    free(out);
    out = NULL;

    return out;
}


cplx_t** test_bell_pure(const unsigned nqubits) {
    cplx_t** out = NULL;
    const cplx_t ampl = INVSQRT2 + I*0.0;

    // Bell states only defined for n=2 qubits
    if (nqubits != 2) {
        fprintf(stderr, "test_bell_pure(): Bell states only defined for 2 qubits\n");
        return NULL;
    }

    const unsigned dim = 4;  // Hilbert space dimension for 2 qubits

    // Allocate the pointers to all 4 Bell states
    if (!((out = calloc(4, sizeof (*out))))) {
        fprintf(stderr, "test_bell_pure(): out allocation failed\n");
        return out;
    }

    // Initialize all 4 Bell states
    for (unsigned i = 0; i < 4; ++i) {
        if (!((out[i] = calloc(dim, sizeof (*out[i]))))) {
            fprintf(stderr, "test_bell_pure(): out[%u] allocation failed\n", i);
            goto cleanup;
        }
    }

    // |Φ⁺⟩ = (|00⟩ + |11⟩)/√2
    out[0][0] = ampl;   // |00⟩
    out[0][3] = ampl;   // |11⟩

    // |Φ⁻⟩ = (|00⟩ - |11⟩)/√2
    out[1][0] = ampl;   // |00⟩
    out[1][3] = -ampl;  // |11⟩

    // |Ψ⁺⟩ = (|01⟩ + |10⟩)/√2
    out[2][1] = ampl;   // |01⟩
    out[2][2] = ampl;   // |10⟩

    // |Ψ⁻⟩ = (|01⟩ - |10⟩)/√2
    out[3][1] = ampl;   // |01⟩
    out[3][2] = -ampl;  // |10⟩

    return out;

    cleanup:
        for (unsigned i = 0; i < 4; ++i) {
            free(out[i]);
            out[i] = NULL;
        }
        free(out);
        out = NULL;

    return out;
}


cplx_t** test_ghz_pure(const unsigned nqubits) {
    cplx_t** out = NULL;
    const cplx_t ampl = INVSQRT2 + I*0.0;

    // GHZ states defined for n >= 2 qubits
    if (nqubits < 2) {
        fprintf(stderr, "test_ghz_pure(): GHZ states require at least 2 qubits\n");
        return NULL;
    }

    const unsigned dim = 1 << nqubits;  // Hilbert space dimension

    // Allocate pointer for single GHZ state
    if (!((out = calloc(1, sizeof (*out))))) {
        fprintf(stderr, "test_ghz_pure(): out allocation failed\n");
        return out;
    }

    // Allocate the GHZ state vector
    if (!((out[0] = calloc(dim, sizeof (*out[0]))))) {
        fprintf(stderr, "test_ghz_pure(): out[0] allocation failed\n");
        free(out);
        return NULL;
    }

    // |GHZ⟩ = (|0...0⟩ + |1...1⟩)/√2
    out[0][0] = ampl;           // |00...0⟩
    out[0][dim - 1] = ampl;     // |11...1⟩

    return out;
}


cplx_t** test_w_pure(const unsigned nqubits) {
    cplx_t** out = NULL;
    const cplx_t ampl = (cplx_t) (1.0 / sqrt((double)nqubits)) + I*0.0;

    // W states defined for n >= 3 qubits
    if (nqubits < 3) {
        fprintf(stderr, "test_w_pure(): W states require at least 3 qubits\n");
        return NULL;
    }

    const unsigned dim = 1 << nqubits;  // Hilbert space dimension

    // Allocate pointer for single W state
    if (!((out = calloc(1, sizeof (*out))))) {
        fprintf(stderr, "test_w_pure(): out allocation failed\n");
        return out;
    }

    // Allocate the W state vector
    if (!((out[0] = calloc(dim, sizeof (*out[0]))))) {
        fprintf(stderr, "test_w_pure(): out[0] allocation failed\n");
        free(out);
        return NULL;
    }

    // |W⟩ = (|100...0⟩ + |010...0⟩ + ... + |00...01⟩)/√n
    // Each basis state with exactly one bit set gets amplitude 1/√n
    for (unsigned k = 0; k < nqubits; ++k) {
        unsigned index = 1 << k;  // Basis state with only k-th bit set
        out[0][index] = ampl;
    }

    return out;
}


cplx_t** test_mk_states_pure(const unsigned nqubits, unsigned *nvecs) {
    cplx_t** out = NULL;

    // Allocate the pointers for all test states
    const unsigned dim = 1 << nqubits;  // Hilbert space dimension
    *nvecs = 3 * dim;   // Number of basis states
    if (nqubits == 2) *nvecs += 4;   // Number of Bell states
    else if (nqubits > 2) *nvecs += 2;   // Number of GHZ + W-state
    if (!((out = calloc(*nvecs, sizeof (*out))))) {
        fprintf(stderr, "test_mk_states_pure(): out allocation failed\n");
        return out;
    }

    // Copy pointers from computational basis states
    cplx_t **tmp = test_cb_pure(nqubits);
    if (!tmp) {
        fprintf(stderr, "test_mk_states_pure(): Error in test_cb_pure\n");
        goto cleanup;
    }
    for (unsigned i = 0; i < dim; ++i) {
        out[i] = tmp[i];
    }

    free(tmp);
    tmp = NULL;

    // Copy pointers from the Hadamard basis states
    tmp = test_xb_pure(nqubits);
    if (!tmp) {
        fprintf(stderr, "test_mk_states_pure(): Initialization of Hadamard basis states failed\n");
        goto cleanup;
    }
    for (unsigned i = 0; i < dim; ++i) {
        out[dim + i] = tmp[i];
    }

    free(tmp);
    tmp = NULL;

    // Copy pointers from the circular basis states
    tmp = test_yb_pure(nqubits);
    if (!tmp) {
        fprintf(stderr, "test_mk_states_pure(): Initialization of Circular basis states failed\n");
        goto cleanup;
    }
    for (unsigned i = 0; i < dim; ++i) {
        out[2 * dim + i] = tmp[i];
    }

    free(tmp);
    tmp = NULL;

    // Copy pointers from the Bell states
    if (nqubits == 2) {
        tmp = test_bell_pure(nqubits);
        if (!tmp) {
            fprintf(stderr, "test_mk_states_pure(): Initialization of Bell states failed\n");
            goto cleanup;
        }
        for (unsigned i = 0; i < 4; ++i) {
            out[3 * dim + i] = tmp[i];
        }

        free(tmp);
        tmp = NULL;
    }

    else if (nqubits > 2) {
        // Copy pointer from the GHZ state
        tmp = test_ghz_pure(nqubits);
        if (!tmp) {
            fprintf(stderr, "test_mk_states_pure(): Initialization of GHZ state failed\n");
            goto cleanup;
        }
        out[3 * dim] = tmp[0];

        free(tmp);
        tmp = NULL;

        // Copy pointer from the W-state
        tmp = test_w_pure(nqubits);
        if (!tmp) {
            fprintf(stderr, "test_mk_states_pure(): Initialization of W-state failed\n");
            goto cleanup;
        }
        out[3 * dim + 1] = tmp[0];

        free(tmp);
        tmp = NULL;
    }

    return out;

    cleanup:
        for (unsigned i = 0; i < *nvecs; ++i) {
            free(out[i]);
            out[i] = NULL;
        }

    return out;
}


void test_rm_states_pure(const unsigned nqubits, cplx_t **states) {
    const unsigned dim = 1 << nqubits;
    unsigned nvecs = 3 * dim;
    if (nqubits == 2) nvecs += 4;
    else if (nqubits > 2) nvecs += 2;

    if (states) {
        for (unsigned i = 0; i < nvecs; ++i) {
            free(states[i]);
            states[i] = NULL;
        }
    }
}