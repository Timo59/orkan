// test.c - Functions and objects to check library

/*
 * =====================================================================================================================
 * includes
 * =====================================================================================================================
 */

#ifndef Q_TEST_H
#include "test.h"
#endif

#ifndef _STDLIB_H_
#include <stdlib.h>
#endif
#include <stdio.h>

#include "q_types.h"

/*
 * =====================================================================================================================
 * Pure quantum test states
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
            fprintf(stderr, "test_cb_pure: out[%u] allocation failed\n");
            goto cleanup;
        }

        out[i][i] = 1.0 + I*0.0;
    }

    return out;

    cleanup:
        if (out) {
            for (unsigned i = 0; i < dim; ++i) {
                free(out);
                out = NULL;
            }
        }
        free(out);
        out = NULL;

    return out;
}


cplx_t** test_xb_pure(const unsigned nqubits) {
    cplx_t** out = NULL;
    const cplx_t ampl = sqrt()

    // Allocate the pointers to all Hadamard basis states
    const unsigned dim = 1 << nqubits;  // Hilbert space dimension
    if (!((out = calloc(dim, sizeof (*out))))) {
        fprintf(stderr, "test_cb_pure: out allocation failed\n");
        return out;
    }

    // Initialize the Hadamard basis states; the i-th Hadamard basis state is
    // H |i> = H |i_n-1...i_1 i_0> = (|0> + (-1)**(i_n-1) |1>) ... (|0> + (-1)**i_1 |1>) (|0> + (-1)**i_0 |1>)
    for (unsigned i = 0; i < dim; ++i) {
        if (!((out[i] = calloc(dim, sizeof (*out[i]))))) {
            fprintf(stderr, "test_cb_pure: out[%u] allocation failed\n");
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


        }

        out[i][i] = 1.0 + I*0.0;
    }

    return out;

    cleanup:
        if (out) {
            for (unsigned i = 0; i < dim; ++i) {
                free(out);
                out = NULL;
            }
        }
    free(out);
    out = NULL;

    return out;
}


cplx_t** test_yb_pure(const unsigned nqubits);


cplx_t** test_gen_states_pure(const unsigned nqubits) {
    cplx_t** out = NULL;

    // Allocate the pointers for all test states
    unsigned nvecs = 1 << nqubits;  // Number of test states
    if (!((out = calloc(nvecs, sizeof (*out))))) {
        fprintf(stderr, "test_gen_states_pure: out allocation failed\n");
        return out;
    }

    // Copy pointers to computational basis states
    cplx_t **tmp = test_cb_pure(nqubits);
    const unsigned dim = 1 << nqubits;
    if (!tmp) {
        fprintf(stderr, "test_gen_states_pure: Error in test_cb_pure\n");
        goto cleanup;
    }
    for (unsigned i = 0; i < dim; ++i) {
        out[i] = tmp[i];
    }

    cleanup:
        if (out) {
            for (unsigned i = 0; i < nvecs; ++i) {
                free(out[i]);
                out[i] = NULL;
            }
        }

    return out;
}
