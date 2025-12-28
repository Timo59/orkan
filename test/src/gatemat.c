// gatemat.h    - Matrix representation of quantum gates

/*
 * =====================================================================================================================
 * includes
 * =====================================================================================================================
 */

#ifndef GATEMAT_H
#include "gatemat.h"
#endif

#ifndef LINALG_H
#include "linalg.h"
#endif

#include <stdio.h>
#include <stdlib.h>

/*
 * =====================================================================================================================
 * Function definitions
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
    if (!((tmp = zkron(dim_left, dim_left, id, 2, 2, gate)))) goto cleanup;
    dim_left *= 2;  // Hilbert space dimension of qubits to the left of and including target
    free(id);
    id = NULL;

    // Identity matrix for all qubits to the right of target
    if (!((id = mat_id(target)))) goto cleanup;

    // Kronecker product with right identity
    const unsigned dim_right = 1ULL << target;
    out = zkron(dim_left, dim_left, tmp, dim_right, dim_right, id);

    cleanup:
        free(id);
    id = NULL;
    free(tmp);
    tmp = NULL;

    return out;
}
