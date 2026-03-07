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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*
 * =====================================================================================================================
 * Function definitions
 * =====================================================================================================================
 */

/*
 * =====================================================================================================================
 * Pauli string matrix builder
 *
 * Builds the 2^n × 2^n matrix for a Pauli string using Kronecker products of the
 * four 2×2 Pauli matrices. This is an independent implementation from the bitmask
 * path in pauli_apply_pure — no masks are consulted.
 *
 * String convention: big-endian — s[0] = qubit n-1 (MSB), s[n-1] = qubit 0 (LSB).
 * Matrix = M(s[0]) ⊗ M(s[1]) ⊗ ... ⊗ M(s[n-1]).
 * =====================================================================================================================
 */

/* 2×2 Pauli matrices in column-major order: {M[0,0], M[1,0], M[0,1], M[1,1]} */
static const cplx_t PMAT_I[4] = { 1.0+0.0*I,  0.0+0.0*I,  0.0+0.0*I,  1.0+0.0*I};
static const cplx_t PMAT_X[4] = { 0.0+0.0*I,  1.0+0.0*I,  1.0+0.0*I,  0.0+0.0*I};
static const cplx_t PMAT_Y[4] = { 0.0+0.0*I,  0.0+1.0*I,  0.0-1.0*I,  0.0+0.0*I};
static const cplx_t PMAT_Z[4] = { 1.0+0.0*I,  0.0+0.0*I,  0.0+0.0*I, -1.0+0.0*I};

static const cplx_t *char_to_pmat(char c) {
    switch (c) {
        case 'X': return PMAT_X;
        case 'Y': return PMAT_Y;
        case 'Z': return PMAT_Z;
        default:  return PMAT_I;
    }
}

cplx_t *mat_pauli_str(const char *s, unsigned n_qubits) {
    if (n_qubits == 0) {
        /* Scalar 1: 1×1 "matrix" */
        cplx_t *out = malloc(sizeof(cplx_t));
        if (!out) {
            fprintf(stderr, "mat_pauli_str(): allocation failed for n=0\n");
            return NULL;
        }
        *out = 1.0 + 0.0*I;
        return out;
    }

    /* Start with the 2×2 matrix for qubit n-1 (s[0]) */
    cplx_t *acc = malloc(4 * sizeof(cplx_t));
    if (!acc) {
        fprintf(stderr, "mat_pauli_str(): initial allocation failed\n");
        return NULL;
    }
    memcpy(acc, char_to_pmat(s[0]), 4 * sizeof(cplx_t));
    unsigned acc_dim = 2;

    /* Kron acc with M(s[i]) for each subsequent qubit */
    for (unsigned i = 1; i < n_qubits; ++i) {
        cplx_t *tmp = zkron(acc_dim, acc_dim, acc, 2, 2, char_to_pmat(s[i]));
        free(acc);
        if (!tmp) {
            fprintf(stderr, "mat_pauli_str(): zkron failed at position %u\n", i);
            return NULL;
        }
        acc = tmp;
        acc_dim *= 2;
    }

    return acc;
}

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

cplx_t* mat_two_qubit_gate(const unsigned nqubits, const cplx_t *gate, const unsigned q1, const unsigned q2) {
    // Check that we have at least 2 qubits
    if (nqubits < 2) {
        fprintf(stderr, "mat_two_qubit_gate(): nqubits must be >= 2; Was %u\n", nqubits);
        return NULL;
    }

    // Check that both qubits are in scope
    if (q1 >= nqubits) {
        fprintf(stderr, "mat_two_qubit_gate(): q1 out of scope; Expected < %u, Was %u\n", nqubits, q1);
        return NULL;
    }
    if (q2 >= nqubits) {
        fprintf(stderr, "mat_two_qubit_gate(): q2 out of scope; Expected < %u, Was %u\n", nqubits, q2);
        return NULL;
    }

    // Check that q1 != q2
    if (q1 == q2) {
        fprintf(stderr, "mat_two_qubit_gate(): q1 and q2 must be different; Both were %u\n", q1);
        return NULL;
    }

    const unsigned dim = 1u << nqubits;
    cplx_t *out = calloc(dim * dim, sizeof(*out));
    if (!out) {
        fprintf(stderr, "mat_two_qubit_gate(): out allocation failed\n");
        return NULL;
    }

    // Build the matrix by directly computing each element
    // For each input basis state |j> and output basis state |i>, compute <i|U|j>
    // The gate acts on qubits q1 and q2, with the 4x4 matrix indexed as:
    //   gate[row + col*4] where row = 2*q1_bit + q2_bit, col = 2*q1_bit' + q2_bit'
    for (unsigned j = 0; j < dim; ++j) {
        // Extract bits at positions q1 and q2 from input state j
        unsigned j_q1 = (j >> q1) & 1;
        unsigned j_q2 = (j >> q2) & 1;
        unsigned j_gate_idx = 2 * j_q1 + j_q2;  // Index into 4x4 gate matrix column

        // The remaining bits (not q1 or q2) stay the same
        // For each possible output of the 2-qubit gate:
        for (unsigned out_q1 = 0; out_q1 < 2; ++out_q1) {
            for (unsigned out_q2 = 0; out_q2 < 2; ++out_q2) {
                unsigned i_gate_idx = 2 * out_q1 + out_q2;  // Index into 4x4 gate matrix row
                cplx_t gate_elem = gate[i_gate_idx + j_gate_idx * 4];

                if (gate_elem == 0.0) continue;

                // Construct output state i: same as j but with q1 and q2 bits replaced
                unsigned i = j;
                // Clear bits at q1 and q2
                i &= ~((1u << q1) | (1u << q2));
                // Set new bits
                i |= (out_q1 << q1) | (out_q2 << q2);

                out[i + j * dim] += gate_elem;
            }
        }
    }

    return out;
}

cplx_t* mat_three_qubit_gate(const unsigned nqubits, const cplx_t *gate,
                             const unsigned q1, const unsigned q2, const unsigned q3) {
    if (nqubits < 3) {
        fprintf(stderr, "mat_three_qubit_gate(): nqubits must be >= 3; Was %u\n", nqubits);
        return NULL;
    }
    if (q1 >= nqubits || q2 >= nqubits || q3 >= nqubits) {
        fprintf(stderr, "mat_three_qubit_gate(): qubit out of scope\n");
        return NULL;
    }
    if (q1 == q2 || q1 == q3 || q2 == q3) {
        fprintf(stderr, "mat_three_qubit_gate(): q1, q2, q3 must all differ; Were %u, %u, %u\n", q1, q2, q3);
        return NULL;
    }

    const unsigned dim = 1u << nqubits;
    cplx_t *out = calloc(dim * dim, sizeof(*out));
    if (!out) {
        fprintf(stderr, "mat_three_qubit_gate(): out allocation failed\n");
        return NULL;
    }

    for (unsigned j = 0; j < dim; ++j) {
        unsigned j_q1 = (j >> q1) & 1;
        unsigned j_q2 = (j >> q2) & 1;
        unsigned j_q3 = (j >> q3) & 1;
        unsigned j_gate_idx = 4 * j_q1 + 2 * j_q2 + j_q3;

        for (unsigned out_q1 = 0; out_q1 < 2; ++out_q1) {
            for (unsigned out_q2 = 0; out_q2 < 2; ++out_q2) {
                for (unsigned out_q3 = 0; out_q3 < 2; ++out_q3) {
                    unsigned i_gate_idx = 4 * out_q1 + 2 * out_q2 + out_q3;
                    cplx_t gate_elem = gate[i_gate_idx + j_gate_idx * 8];

                    if (gate_elem == 0.0) continue;

                    unsigned i = j;
                    i &= ~((1u << q1) | (1u << q2) | (1u << q3));
                    i |= (out_q1 << q1) | (out_q2 << q2) | (out_q3 << q3);

                    out[i + j * dim] += gate_elem;
                }
            }
        }
    }

    return out;
}

/*
 * =====================================================================================================================
 * Rotation gate matrix builders  (moved from test_gate_pure_1q.c / test_gate.c)
 * =====================================================================================================================
 */

void mat_rx(double theta, cplx_t *mat) {
    double c = cos(theta / 2.0);
    double s = sin(theta / 2.0);
    // [[cos(θ/2), -i·sin(θ/2)], [-i·sin(θ/2), cos(θ/2)]]
    // Column-major: [0,0], [1,0], [0,1], [1,1]
    mat[0] = c;
    mat[1] = -s * I;
    mat[2] = -s * I;
    mat[3] = c;
}

void mat_ry(double theta, cplx_t *mat) {
    double c = cos(theta / 2.0);
    double s = sin(theta / 2.0);
    // [[cos(θ/2), -sin(θ/2)], [sin(θ/2), cos(θ/2)]]
    // Column-major: [0,0], [1,0], [0,1], [1,1]
    mat[0] = c;
    mat[1] = s;
    mat[2] = -s;
    mat[3] = c;
}

void mat_rz(double theta, cplx_t *mat) {
    double c = cos(theta / 2.0);
    double s = sin(theta / 2.0);
    // [[e^(-iθ/2), 0], [0, e^(iθ/2)]]
    // Column-major: [0,0], [1,0], [0,1], [1,1]
    mat[0] = c - s * I;
    mat[1] = 0;
    mat[2] = 0;
    mat[3] = c + s * I;
}

void mat_p(double theta, cplx_t *mat) {
    // [[1, 0], [0, e^(iθ)]]
    // Column-major: [0,0], [1,0], [0,1], [1,1]
    mat[0] = 1.0;
    mat[1] = 0.0;
    mat[2] = 0.0;
    mat[3] = cos(theta) + sin(theta) * I;
}

void mat_rx_pi3(cplx_t *mat) {
    const double c = cos(M_PI / 6.0);  /* cos(pi/6) = sqrt(3)/2 */
    const double s = sin(M_PI / 6.0);  /* sin(pi/6) = 1/2       */
    mat[0] =  c + 0.0 * I;
    mat[1] = 0.0 -   s * I;
    mat[2] = 0.0 -   s * I;
    mat[3] =  c + 0.0 * I;
}
