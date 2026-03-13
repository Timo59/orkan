/*
 * bench_circuit.c — Circuit lifecycle, append helpers, and circuit constructors.
 *
 * Implements the circuit abstraction defined in bench_circuit.h:
 *   - Dynamic array management (alloc, free, append with doubling realloc)
 *   - Composite gate helpers (circuit_rzz)
 *   - Three circuit constructors: QAOA p=1, VQE HEA, Quantum Volume
 */

#include "bench_circuit.h"
#include <assert.h>
#include <string.h>

/* =============================================================================================================
 * Circuit lifecycle
 * ============================================================================================================= */

circuit_t circuit_alloc(int n_qubits, int initial_cap, const char *name)
{
    circuit_t c;
    c.n_qubits = n_qubits;
    c.n_ops    = 0;
    c.name     = name;

    if (initial_cap <= 0) {
        c.ops = NULL;
        c.cap = 0;
        return c;
    }

    c.ops = (circuit_op_t *)malloc((size_t)initial_cap * sizeof(circuit_op_t));
    if (!c.ops) {
        c.cap = 0;
        return c;
    }
    c.cap = initial_cap;
    return c;
}

void circuit_free(circuit_t *c)
{
    if (!c) return;
    free(c->ops);
    memset(c, 0, sizeof(*c));
}

/* =============================================================================================================
 * Internal: ensure capacity for at least one more op (doubling realloc)
 * ============================================================================================================= */

static int circuit_ensure_cap(circuit_t *c)
{
    if (c->n_ops < c->cap) return 0;

    int new_cap = (c->cap == 0) ? 16 : c->cap * 2;
    circuit_op_t *new_ops = (circuit_op_t *)realloc(c->ops,
                                                     (size_t)new_cap * sizeof(circuit_op_t));
    if (!new_ops) return -1;
    c->ops = new_ops;
    c->cap = new_cap;
    return 0;
}

/* =============================================================================================================
 * Append helpers
 * ============================================================================================================= */

int circuit_append_1q(circuit_t *c, circuit_gate_t g, qubit_t q0)
{
    if (circuit_ensure_cap(c) != 0) return -1;

    circuit_op_t *op = &c->ops[c->n_ops++];
    op->gate_type = g;
    op->q0        = q0;
    op->q1        = 0;
    op->angle     = 0.0;
    return 0;
}

int circuit_append_1q_angle(circuit_t *c, circuit_gate_t g, qubit_t q0, double angle)
{
    if (circuit_ensure_cap(c) != 0) return -1;

    circuit_op_t *op = &c->ops[c->n_ops++];
    op->gate_type = g;
    op->q0        = q0;
    op->q1        = 0;
    op->angle     = angle;
    return 0;
}

int circuit_append_2q(circuit_t *c, circuit_gate_t g, qubit_t q0, qubit_t q1)
{
    if (circuit_ensure_cap(c) != 0) return -1;

    circuit_op_t *op = &c->ops[c->n_ops++];
    op->gate_type = g;
    op->q0        = q0;
    op->q1        = q1;
    op->angle     = 0.0;
    return 0;
}

/* =============================================================================================================
 * Composite: RZZ(gamma) = CX(q0,q1) . Rz(q1, 2*gamma) . CX(q0,q1)
 * ============================================================================================================= */

int circuit_rzz(circuit_t *c, qubit_t q0, qubit_t q1, double gamma)
{
    if (circuit_append_2q(c, GATE_CX, q0, q1) != 0) return -1;
    if (circuit_append_1q_angle(c, GATE_RZ, q1, 2.0 * gamma) != 0) return -1;
    if (circuit_append_2q(c, GATE_CX, q0, q1) != 0) return -1;
    return 0;
}

/* =============================================================================================================
 * QAOA p=1 on cycle graph
 *
 * n H gates + n RZZ(gamma) decompositions (3n gates) + n Rx(2*beta) = 5n total gates.
 * ============================================================================================================= */

circuit_t circuit_qaoa_p1(int n, double gamma, double beta)
{
    int initial_cap = 5 * n;
    circuit_t c = circuit_alloc(n, initial_cap, "qaoa_p1");
    if (!c.ops) return c;

    /* Layer 1: Hadamard on all qubits */
    for (int q = 0; q < n; ++q) {
        if (circuit_append_1q(&c, GATE_H, (qubit_t)q) != 0) { c.n_ops = 0; return c; }
    }

    /* Layer 2: ZZ interactions on cycle edges (i, (i+1) mod n) */
    for (int i = 0; i < n; ++i) {
        qubit_t q0 = (qubit_t)i;
        qubit_t q1 = (qubit_t)((i + 1) % n);
        if (circuit_rzz(&c, q0, q1, gamma) != 0) { c.n_ops = 0; return c; }
    }

    /* Layer 3: Rx(2*beta) on all qubits */
    for (int q = 0; q < n; ++q) {
        if (circuit_append_1q_angle(&c, GATE_RX, (qubit_t)q, 2.0 * beta) != 0) { c.n_ops = 0; return c; }
    }

    assert(c.n_ops == 5 * n);
    return c;
}

/* =============================================================================================================
 * Hardware-Efficient VQE Ansatz
 *
 * depth layers of: n Ry(pi/4) + (n-1) CX ladder = (2n-1) ops per layer.
 * Total: depth * (2n - 1).
 * ============================================================================================================= */

circuit_t circuit_vqe_hea(int n, int depth)
{
    if (depth < 0) depth = n / 2;

    int initial_cap = depth * (2 * n - 1);
    circuit_t c = circuit_alloc(n, initial_cap, "vqe_hea");
    if (!c.ops && initial_cap > 0) return c;

    for (int d = 0; d < depth; ++d) {
        /* Ry(pi/4) on all qubits */
        for (int q = 0; q < n; ++q) {
            if (circuit_append_1q_angle(&c, GATE_RY, (qubit_t)q, M_PI / 4.0) != 0) { c.n_ops = 0; return c; }
        }
        /* CX ladder: 0->1, 1->2, ..., (n-2)->(n-1) */
        for (int q = 0; q < n - 1; ++q) {
            if (circuit_append_2q(&c, GATE_CX, (qubit_t)q, (qubit_t)(q + 1)) != 0) { c.n_ops = 0; return c; }
        }
    }

    assert(c.n_ops == depth * (2 * n - 1));
    return c;
}

/* =============================================================================================================
 * Quantum Volume circuit
 *
 * QV circuit template: random parameterized two-qubit circuit template.
 *
 * Each qubit pair receives: Ry(a0)x3 on qa + Ry(a1)x3 on qb
 *   + CX(qa,qb) . CX(qb,qa) . CX(qa,qb)
 *   + Rz(a2)x3 on qa + Rz(a3)x3 on qb = 15 ops/pair.
 * This exercises the same gate MIX as a KAK decomposition (9 single-qubit rotations
 * and 3 CNOT gates per pair) but the angles are not KAK-computed interaction parameters.
 * Do not describe this as a KAK decomposition in any external-facing text.
 *
 * n layers, each layer: random permutation of n qubits (Fisher-Yates), then n/2 pairs.
 * Total ops: n * (n/2) * 15 = 15n^2/2. Only valid for even n.
 * ============================================================================================================= */

/* 64-bit LCG: returns an angle in [0, pi] */
static double lcg_next_angle(uint64_t *state)
{
    *state = (*state) * 6364136223846793005ULL + 1442695040888963407ULL;
    return (double)(*state >> 33) / (double)(1ULL << 31) * M_PI;
}

circuit_t circuit_qv(int n, int seed)
{
    if (n % 2 != 0) {
        fprintf(stderr, "circuit_qv: n must be even, got %d\n", n);
        return circuit_alloc(n, 0, "qv");
    }

    int n_pairs  = n / 2;
    int n_layers = n;
    int total_ops = 15 * n_pairs * n_layers;
    circuit_t c = circuit_alloc(n, total_ops, "qv");
    if (!c.ops) return c;

    uint64_t lcg_state = (uint64_t)(uint32_t)seed;

    /* Permutation buffer */
    qubit_t *perm = (qubit_t *)malloc((size_t)n * sizeof(qubit_t));
    if (!perm) {
        circuit_free(&c);
        return circuit_alloc(n, 0, "qv");
    }

    for (int layer = 0; layer < n_layers; ++layer) {
        /* Initialize identity permutation */
        for (int i = 0; i < n; ++i) perm[i] = (qubit_t)i;

        /* Fisher-Yates shuffle using LCG */
        for (int j = 0; j < n - 1; ++j) {
            /* Draw a random index in [j, n-1] */
            lcg_state = lcg_state * 6364136223846793005ULL + 1442695040888963407ULL;
            int k = j + (int)((lcg_state >> 33) % (uint64_t)(n - j));
            qubit_t tmp = perm[j];
            perm[j] = perm[k];
            perm[k] = tmp;
        }

        /* Process pairs */
        for (int p = 0; p < n_pairs; ++p) {
            qubit_t qa = perm[2 * p];
            qubit_t qb = perm[2 * p + 1];

            /* 3 Ry on qa */
            for (int i = 0; i < 3; ++i) {
                double angle = lcg_next_angle(&lcg_state);
                if (circuit_append_1q_angle(&c, GATE_RY, qa, angle) != 0) { c.n_ops = 0; free(perm); return c; }
            }
            /* 3 Ry on qb */
            for (int i = 0; i < 3; ++i) {
                double angle = lcg_next_angle(&lcg_state);
                if (circuit_append_1q_angle(&c, GATE_RY, qb, angle) != 0) { c.n_ops = 0; free(perm); return c; }
            }
            /* CX(qa,qb), CX(qb,qa), CX(qa,qb) */
            if (circuit_append_2q(&c, GATE_CX, qa, qb) != 0) { c.n_ops = 0; free(perm); return c; }
            if (circuit_append_2q(&c, GATE_CX, qb, qa) != 0) { c.n_ops = 0; free(perm); return c; }
            if (circuit_append_2q(&c, GATE_CX, qa, qb) != 0) { c.n_ops = 0; free(perm); return c; }
            /* 3 Rz on qa */
            for (int i = 0; i < 3; ++i) {
                double angle = lcg_next_angle(&lcg_state);
                if (circuit_append_1q_angle(&c, GATE_RZ, qa, angle) != 0) { c.n_ops = 0; free(perm); return c; }
            }
            /* 3 Rz on qb */
            for (int i = 0; i < 3; ++i) {
                double angle = lcg_next_angle(&lcg_state);
                if (circuit_append_1q_angle(&c, GATE_RZ, qb, angle) != 0) { c.n_ops = 0; free(perm); return c; }
            }
        }
    }

    free(perm);
    assert(c.n_ops == total_ops);
    return c;
}
