/**
 * @file bench_mixed.c
 * @brief Benchmark comparing qlib's mixed state gate implementations against alternatives
 *
 * Compares:
 * 1. qlib packed sparse - Production implementation with packed lower-triangular storage
 * 2. BLAS dense         - Full matrix UρU† via zgemm
 * 3. Naive loop         - Element-by-element computation
 */

#include "bench.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <math.h>

#if defined(__APPLE__)
    #include <vecLib/cblas_new.h>
#elif defined(__linux__)
    #include <cblas.h>
#endif

/*
 * =====================================================================================================================
 * Gate matrices (2x2, column-major)
 * =====================================================================================================================
 */

#define INVSQRT2 0.7071067811865475

static const cplx_t XMAT[4] = {0.0, 1.0, 1.0, 0.0};
static const cplx_t ZMAT[4] = {1.0, 0.0, 0.0, -1.0};
static const cplx_t HMAT[4] = {INVSQRT2, INVSQRT2, INVSQRT2, -INVSQRT2};

/*
 * Two-qubit gate matrices (4x4, column-major)
 *
 * Convention: for gate(q1, q2), matrix index = b_q2 * 2 + b_q1
 * where b_qi is the bit value of qubit qi in the basis state.
 */

/* CX(control, target): flips target when control=1 */
static const cplx_t CXMAT[16] = {
    1, 0, 0, 0,   /* col 0: |c=0,t=0> -> |c=0,t=0> */
    0, 0, 0, 1,   /* col 1: |c=1,t=0> -> |c=1,t=1> */
    0, 0, 1, 0,   /* col 2: |c=0,t=1> -> |c=0,t=1> */
    0, 1, 0, 0    /* col 3: |c=1,t=1> -> |c=1,t=0> */
};

/* SWAP(q1, q2): exchanges qubit values */
static const cplx_t SWAPMAT[16] = {
    1, 0, 0, 0,   /* col 0: |0,0> -> |0,0> */
    0, 0, 1, 0,   /* col 1: |1,0> -> |0,1> */
    0, 1, 0, 0,   /* col 2: |0,1> -> |1,0> */
    0, 0, 0, 1    /* col 3: |1,1> -> |1,1> */
};

/*
 * =====================================================================================================================
 * Helper functions
 * =====================================================================================================================
 */

/** @brief Initialize a random mixed state (density matrix) */
static void init_random_mixed_state(state_t *state, qubit_t qubits, state_type_t type) {
    state->type = type;
    state_init(state, qubits, NULL);
    if (!state->data) {
        fprintf(stderr, "init_random_mixed_state: allocation failed\n");
        return;
    }

    /* Fill with random values (not a valid density matrix, but fine for benchmarking) */
    dim_t len = state_len(state);
    for (dim_t i = 0; i < len; ++i) {
        double re = (double)rand() / RAND_MAX - 0.5;
        double im = (double)rand() / RAND_MAX - 0.5;
        state->data[i] = re + im * I;
    }
}

/** @brief Kronecker product: C = A ⊗ B */
static cplx_t* kron(dim_t ma, dim_t na, const cplx_t *A,
                    dim_t mb, dim_t nb, const cplx_t *B) {
    dim_t mc = ma * mb;
    dim_t nc = na * nb;
    cplx_t *C = malloc(mc * nc * sizeof(cplx_t));
    if (!C) return NULL;

    for (dim_t ja = 0; ja < na; ++ja) {
        for (dim_t ia = 0; ia < ma; ++ia) {
            cplx_t aij = A[ia + ja * ma];
            for (dim_t jb = 0; jb < nb; ++jb) {
                for (dim_t ib = 0; ib < mb; ++ib) {
                    dim_t ic = ib + ia * mb;
                    dim_t jc = jb + ja * nb;
                    C[ic + jc * mc] = aij * B[ib + jb * mb];
                }
            }
        }
    }
    return C;
}

/** @brief Build full gate matrix I⊗...⊗U⊗...⊗I for single-qubit gate on `target` */
static cplx_t* build_full_gate_matrix(qubit_t qubits, const cplx_t gate[4], qubit_t target) {
    qubit_t pos = qubits - 1 - target;
    dim_t dim_left = (dim_t)1 << pos;

    cplx_t *id_left = calloc(dim_left * dim_left, sizeof(cplx_t));
    if (!id_left) return NULL;
    for (dim_t i = 0; i < dim_left; ++i) id_left[i + i * dim_left] = 1.0;

    cplx_t *tmp = kron(dim_left, dim_left, id_left, 2, 2, gate);
    free(id_left);
    if (!tmp) return NULL;
    dim_left *= 2;

    dim_t dim_right = (dim_t)1 << target;
    cplx_t *id_right = calloc(dim_right * dim_right, sizeof(cplx_t));
    if (!id_right) { free(tmp); return NULL; }
    for (dim_t i = 0; i < dim_right; ++i) id_right[i + i * dim_right] = 1.0;

    cplx_t *U = kron(dim_left, dim_left, tmp, dim_right, dim_right, id_right);
    free(tmp);
    free(id_right);
    return U;
}

/** @brief Build full gate matrix for two-qubit gate on qubits q1, q2 */
static cplx_t* build_full_gate_matrix_2q(qubit_t qubits, const cplx_t gate[16],
                                          qubit_t q1, qubit_t q2) {
    dim_t dim = (dim_t)1 << qubits;
    cplx_t *U = calloc(dim * dim, sizeof(cplx_t));
    if (!U) return NULL;

    for (dim_t i = 0; i < dim; ++i) {
        int b1 = (i >> q1) & 1;
        int b2 = (i >> q2) & 1;
        int g_col = b2 * 2 + b1;

        for (int g_row = 0; g_row < 4; ++g_row) {
            cplx_t val = gate[g_row + g_col * 4];
            if (cabs(val) < 1e-15) continue;

            int out_b1 = g_row & 1;
            int out_b2 = (g_row >> 1) & 1;

            dim_t j = i;
            j = (j & ~((dim_t)1 << q1)) | ((dim_t)out_b1 << q1);
            j = (j & ~((dim_t)1 << q2)) | ((dim_t)out_b2 << q2);
            U[j + i * dim] = val;
        }
    }
    return U;
}

/*
 * Build all ordered qubit pairs (q1 < q2) for an n-qubit system.
 * Returns the number of pairs written. MAX_PAIRS must be >= qubits*(qubits-1)/2.
 */
#define MAX_PAIRS 128  /* sufficient for up to ~16 qubits */

static int build_all_pairs(qubit_t qubits, qubit_t *q1_out, qubit_t *q2_out) {
    int n = 0;
    for (qubit_t q1 = 0; q1 < qubits; ++q1) {
        for (qubit_t q2 = q1 + 1; q2 < qubits; ++q2) {
            if (n >= MAX_PAIRS) break;
            q1_out[n] = q1;
            q2_out[n] = q2;
            n++;
        }
    }
    return n;
}

/*
 * Compute UρU† using BLAS zgemm.
 *
 * tmp: pre-allocated scratch buffer of dim*dim complex doubles.
 * Passing tmp from outside the hot loop eliminates per-call allocation overhead,
 * keeping the comparison focused on algorithmic cost.
 */
static void apply_uruh_blas(dim_t dim, const cplx_t *U, cplx_t *rho, cplx_t *tmp) {
    const cplx_t alpha = 1.0;
    const cplx_t beta  = 0.0;

    /* tmp = U * rho */
    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                dim, dim, dim, &alpha, U, dim, rho, dim, &beta, tmp, dim);

    /* rho = tmp * U† */
    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasConjTrans,
                dim, dim, dim, &alpha, tmp, dim, U, dim, &beta, rho, dim);
}

/*
 * Compute UρU† using naive triple-nested loops.
 *
 * tmp, out: pre-allocated scratch buffers of dim*dim complex doubles.
 * Same rationale as apply_uruh_blas: allocation lives outside the hot loop.
 */
static void apply_uruh_naive(dim_t dim, const cplx_t *U, cplx_t *rho,
                              cplx_t *tmp, cplx_t *out) {
    /* tmp = U * rho */
    for (dim_t j = 0; j < dim; ++j) {
        for (dim_t i = 0; i < dim; ++i) {
            cplx_t sum = 0.0;
            for (dim_t k = 0; k < dim; ++k)
                sum += U[i + k * dim] * rho[k + j * dim];
            tmp[i + j * dim] = sum;
        }
    }

    /* out = tmp * U† */
    for (dim_t j = 0; j < dim; ++j) {
        for (dim_t i = 0; i < dim; ++i) {
            cplx_t sum = 0.0;
            for (dim_t k = 0; k < dim; ++k)
                sum += tmp[i + k * dim] * conj(U[j + k * dim]);
            out[i + j * dim] = sum;
        }
    }

    memcpy(rho, out, dim * dim * sizeof(cplx_t));
}

/*
 * =====================================================================================================================
 * Benchmark implementations
 * =====================================================================================================================
 */

bench_result_t bench_qlib_packed(qubit_t qubits, const char *gate_name,
                                  void (*gate_fn)(state_t*, qubit_t),
                                  int iterations, int warmup) {
    bench_result_t result = {0};
    result.qubits      = qubits;
    result.dim         = (dim_t)1 << qubits;
    result.gate_name   = gate_name;
    result.method      = "qlib_packed";
    result.iterations  = iterations;

    state_t state = {0};
    init_random_mixed_state(&state, qubits, MIXED_PACKED);
    if (!state.data) return result;

    /* Theoretical size is exact: we know the allocator */
    result.memory_bytes = bench_packed_size(qubits);

    /* Warm-up: cycle through all qubit targets so every code path is hot */
    for (int i = 0; i < warmup; ++i)
        gate_fn(&state, (qubit_t)(i % qubits));

    uint64_t start = bench_time_ns();
    for (int i = 0; i < iterations; ++i)
        for (qubit_t t = 0; t < qubits; ++t)
            gate_fn(&state, t);
    uint64_t end = bench_time_ns();

    result.time_ms    = bench_ns_to_ms(end - start);
    result.ops_per_sec = (double)(iterations * qubits) / (result.time_ms / 1000.0);

    state_free(&state);
    return result;
}

bench_result_t bench_qlib_tiled(qubit_t qubits, const char *gate_name,
                                 void (*gate_fn)(state_t*, qubit_t),
                                 int iterations, int warmup) {
    bench_result_t result = {0};
    result.qubits      = qubits;
    result.dim         = (dim_t)1 << qubits;
    result.gate_name   = gate_name;
    result.method      = "qlib_tiled";
    result.iterations  = iterations;

    state_t state = {0};
    init_random_mixed_state(&state, qubits, MIXED_TILED);
    if (!state.data) return result;

    result.memory_bytes = bench_tiled_size(qubits);

    for (int i = 0; i < warmup; ++i)
        gate_fn(&state, (qubit_t)(i % qubits));

    uint64_t start = bench_time_ns();
    for (int i = 0; i < iterations; ++i)
        for (qubit_t t = 0; t < qubits; ++t)
            gate_fn(&state, t);
    uint64_t end = bench_time_ns();

    result.time_ms     = bench_ns_to_ms(end - start);
    result.ops_per_sec = (double)(iterations * qubits) / (result.time_ms / 1000.0);

    state_free(&state);
    return result;
}

bench_result_t bench_blas_dense(qubit_t qubits, const char *gate_name,
                                 const cplx_t gate_mat[4],
                                 int iterations, int warmup) {
    bench_result_t result = {0};
    result.qubits      = qubits;
    result.dim         = (dim_t)1 << qubits;
    result.gate_name   = gate_name;
    result.method      = "blas_dense";
    result.iterations  = iterations;

    dim_t dim = result.dim;

    cplx_t *rho = malloc(dim * dim * sizeof(cplx_t));
    if (!rho) return result;
    for (dim_t i = 0; i < dim * dim; ++i)
        rho[i] = (double)rand() / RAND_MAX + ((double)rand() / RAND_MAX) * I;

    cplx_t **U = malloc(qubits * sizeof(cplx_t*));
    if (!U) { free(rho); return result; }
    for (qubit_t t = 0; t < qubits; ++t)
        U[t] = build_full_gate_matrix(qubits, gate_mat, t);

    /* Single scratch buffer allocated once outside the hot loop */
    cplx_t *tmp = malloc(dim * dim * sizeof(cplx_t));
    if (!tmp) {
        for (qubit_t t = 0; t < qubits; ++t) free(U[t]);
        free(U); free(rho);
        return result;
    }

    result.memory_bytes = bench_dense_size(qubits);

    for (int i = 0; i < warmup; ++i)
        apply_uruh_blas(dim, U[i % qubits], rho, tmp);

    uint64_t start = bench_time_ns();
    for (int i = 0; i < iterations; ++i)
        for (qubit_t t = 0; t < qubits; ++t)
            apply_uruh_blas(dim, U[t], rho, tmp);
    uint64_t end = bench_time_ns();

    result.time_ms     = bench_ns_to_ms(end - start);
    result.ops_per_sec = (double)(iterations * qubits) / (result.time_ms / 1000.0);

    free(tmp);
    free(rho);
    for (qubit_t t = 0; t < qubits; ++t) free(U[t]);
    free(U);
    return result;
}

bench_result_t bench_naive_loop(qubit_t qubits, const char *gate_name,
                                 const cplx_t gate_mat[4],
                                 int iterations, int warmup) {
    bench_result_t result = {0};
    result.qubits      = qubits;
    result.dim         = (dim_t)1 << qubits;
    result.gate_name   = gate_name;
    result.method      = "naive_loop";
    result.iterations  = iterations;

    dim_t dim = result.dim;

    cplx_t *rho = malloc(dim * dim * sizeof(cplx_t));
    if (!rho) return result;
    for (dim_t i = 0; i < dim * dim; ++i)
        rho[i] = (double)rand() / RAND_MAX + ((double)rand() / RAND_MAX) * I;

    cplx_t **U = malloc(qubits * sizeof(cplx_t*));
    if (!U) { free(rho); return result; }
    for (qubit_t t = 0; t < qubits; ++t)
        U[t] = build_full_gate_matrix(qubits, gate_mat, t);

    /* Two scratch buffers allocated once outside the hot loop */
    cplx_t *tmp = malloc(dim * dim * sizeof(cplx_t));
    cplx_t *out = malloc(dim * dim * sizeof(cplx_t));
    if (!tmp || !out) {
        free(tmp); free(out);
        for (qubit_t t = 0; t < qubits; ++t) free(U[t]);
        free(U); free(rho);
        return result;
    }

    result.memory_bytes = bench_dense_size(qubits);

    for (int i = 0; i < warmup; ++i)
        apply_uruh_naive(dim, U[i % qubits], rho, tmp, out);

    uint64_t start = bench_time_ns();
    for (int i = 0; i < iterations; ++i)
        for (qubit_t t = 0; t < qubits; ++t)
            apply_uruh_naive(dim, U[t], rho, tmp, out);
    uint64_t end = bench_time_ns();

    result.time_ms     = bench_ns_to_ms(end - start);
    result.ops_per_sec = (double)(iterations * qubits) / (result.time_ms / 1000.0);

    free(tmp); free(out);
    free(rho);
    for (qubit_t t = 0; t < qubits; ++t) free(U[t]);
    free(U);
    return result;
}

/*
 * =====================================================================================================================
 * Two-qubit benchmark implementations
 * =====================================================================================================================
 *
 * Two-qubit gates iterate over ALL ordered pairs (q1, q2) with q1 < q2,
 * giving qubits*(qubits-1)/2 applications per iteration.
 * This exercises every code path: within-tile, mixed, and tile-level cases.
 */

bench_result_t bench_qlib_packed_2q(qubit_t qubits, const char *gate_name,
                                      void (*gate_fn)(state_t*, qubit_t, qubit_t),
                                      int iterations, int warmup) {
    bench_result_t result = {0};
    result.qubits      = qubits;
    result.dim         = (dim_t)1 << qubits;
    result.gate_name   = gate_name;
    result.method      = "qlib_packed";
    result.iterations  = iterations;

    state_t state = {0};
    init_random_mixed_state(&state, qubits, MIXED_PACKED);
    if (!state.data) return result;
    result.memory_bytes = bench_packed_size(qubits);

    qubit_t q1s[MAX_PAIRS], q2s[MAX_PAIRS];
    int num_pairs = build_all_pairs(qubits, q1s, q2s);
    if (num_pairs == 0) { state_free(&state); return result; }

    for (int i = 0; i < warmup; ++i)
        gate_fn(&state, q1s[i % num_pairs], q2s[i % num_pairs]);

    uint64_t start = bench_time_ns();
    for (int i = 0; i < iterations; ++i)
        for (int p = 0; p < num_pairs; ++p)
            gate_fn(&state, q1s[p], q2s[p]);
    uint64_t end = bench_time_ns();

    result.time_ms     = bench_ns_to_ms(end - start);
    result.ops_per_sec = (double)(iterations * num_pairs) / (result.time_ms / 1000.0);

    state_free(&state);
    return result;
}

bench_result_t bench_qlib_tiled_2q(qubit_t qubits, const char *gate_name,
                                     void (*gate_fn)(state_t*, qubit_t, qubit_t),
                                     int iterations, int warmup) {
    bench_result_t result = {0};
    result.qubits      = qubits;
    result.dim         = (dim_t)1 << qubits;
    result.gate_name   = gate_name;
    result.method      = "qlib_tiled";
    result.iterations  = iterations;

    state_t state = {0};
    init_random_mixed_state(&state, qubits, MIXED_TILED);
    if (!state.data) return result;
    result.memory_bytes = bench_tiled_size(qubits);

    qubit_t q1s[MAX_PAIRS], q2s[MAX_PAIRS];
    int num_pairs = build_all_pairs(qubits, q1s, q2s);
    if (num_pairs == 0) { state_free(&state); return result; }

    for (int i = 0; i < warmup; ++i)
        gate_fn(&state, q1s[i % num_pairs], q2s[i % num_pairs]);

    uint64_t start = bench_time_ns();
    for (int i = 0; i < iterations; ++i)
        for (int p = 0; p < num_pairs; ++p)
            gate_fn(&state, q1s[p], q2s[p]);
    uint64_t end = bench_time_ns();

    result.time_ms     = bench_ns_to_ms(end - start);
    result.ops_per_sec = (double)(iterations * num_pairs) / (result.time_ms / 1000.0);

    state_free(&state);
    return result;
}

bench_result_t bench_blas_dense_2q(qubit_t qubits, const char *gate_name,
                                     const cplx_t gate_mat[16],
                                     int iterations, int warmup) {
    bench_result_t result = {0};
    result.qubits      = qubits;
    result.dim         = (dim_t)1 << qubits;
    result.gate_name   = gate_name;
    result.method      = "blas_dense";
    result.iterations  = iterations;

    dim_t dim = result.dim;

    qubit_t q1s[MAX_PAIRS], q2s[MAX_PAIRS];
    int num_pairs = build_all_pairs(qubits, q1s, q2s);
    if (num_pairs == 0) return result;

    cplx_t *rho = malloc(dim * dim * sizeof(cplx_t));
    if (!rho) return result;
    for (dim_t i = 0; i < dim * dim; ++i)
        rho[i] = (double)rand() / RAND_MAX + ((double)rand() / RAND_MAX) * I;

    cplx_t **U = malloc(num_pairs * sizeof(cplx_t*));
    if (!U) { free(rho); return result; }
    for (int p = 0; p < num_pairs; ++p)
        U[p] = build_full_gate_matrix_2q(qubits, gate_mat, q1s[p], q2s[p]);

    cplx_t *tmp = malloc(dim * dim * sizeof(cplx_t));
    if (!tmp) {
        for (int p = 0; p < num_pairs; ++p) free(U[p]);
        free(U); free(rho);
        return result;
    }

    result.memory_bytes = bench_dense_size(qubits);

    for (int i = 0; i < warmup; ++i)
        apply_uruh_blas(dim, U[i % num_pairs], rho, tmp);

    uint64_t start = bench_time_ns();
    for (int i = 0; i < iterations; ++i)
        for (int p = 0; p < num_pairs; ++p)
            apply_uruh_blas(dim, U[p], rho, tmp);
    uint64_t end = bench_time_ns();

    result.time_ms     = bench_ns_to_ms(end - start);
    result.ops_per_sec = (double)(iterations * num_pairs) / (result.time_ms / 1000.0);

    free(tmp);
    free(rho);
    for (int p = 0; p < num_pairs; ++p) free(U[p]);
    free(U);
    return result;
}

/*
 * =====================================================================================================================
 * Output functions
 * =====================================================================================================================
 */

void bench_print_result(const bench_result_t *result, int verbose) {
    printf("  %-14s %10.3f ms (%'12.0f ops/sec)",
           result->method, result->time_ms, result->ops_per_sec);
    if (verbose)
        printf("  [%zu bytes]", result->memory_bytes);
    printf("\n");
}

void bench_print_csv_header(void) {
    printf("qubits,dim,gate,method,time_ms,ops_per_sec,memory_bytes,iterations\n");
}

void bench_print_csv(const bench_result_t *result) {
    printf("%u,%lld,%s,%s,%.6f,%.0f,%zu,%d\n",
           result->qubits, (long long)result->dim, result->gate_name, result->method,
           result->time_ms, result->ops_per_sec, result->memory_bytes, result->iterations);
}

/*
 * =====================================================================================================================
 * Summary struct and shared speedup/memory printing
 * =====================================================================================================================
 *
 * Collects results for one (gate, qubit-count) combination across all methods.
 * Eliminates the duplicate speedup+memory display logic between the 1Q and 2Q loops.
 */

typedef struct {
    bench_result_t packed;
    bench_result_t tiled;   int has_tiled;
    bench_result_t dense;   int has_dense;
    bench_result_t naive;   int has_naive;
    double naive_scale;     /* multiply naive speedup by this when naive ran fewer iterations */
    bench_result_t quest;   int has_quest;
    bench_result_t qulacs;  int has_qulacs;
} bench_summary_t;

static void print_speedup_and_memory(const bench_summary_t *s) {
    /* Speedup line */
    printf("  Speedup:");
    int has = 0;
    if (s->has_tiled && s->tiled.time_ms > 0) {
        printf(" %.2fx tiled vs packed", s->packed.time_ms / s->tiled.time_ms);
        has = 1;
    }
    if (s->has_dense && s->dense.time_ms > 0) {
        printf("%s%.1fx vs dense", has ? ", " : " ", s->dense.time_ms / s->packed.time_ms);
        has = 1;
    }
    if (s->has_naive && s->naive.time_ms > 0) {
        double speedup = (s->naive.time_ms / s->packed.time_ms) * s->naive_scale;
        printf("%s~%.0fx vs naive", has ? ", " : " ", speedup);
        has = 1;
    }
    if (s->has_quest && s->quest.time_ms > 0) {
        printf("%s%.1fx vs QuEST", has ? ", " : " ", s->quest.time_ms / s->packed.time_ms);
        has = 1;
    }
    if (s->has_qulacs && s->qulacs.time_ms > 0) {
        printf("%s%.1fx vs Qulacs", has ? ", " : " ", s->qulacs.time_ms / s->packed.time_ms);
        has = 1;
    }
    (void)has;
    printf("\n");

    /* Memory line */
    size_t max_mem = s->packed.memory_bytes;
    if (s->has_tiled  && s->tiled.memory_bytes  > max_mem) max_mem = s->tiled.memory_bytes;
    if (s->has_dense  && s->dense.memory_bytes  > max_mem) max_mem = s->dense.memory_bytes;
    if (s->has_naive  && s->naive.memory_bytes  > max_mem) max_mem = s->naive.memory_bytes;
    if (s->has_quest  && s->quest.memory_bytes  > max_mem) max_mem = s->quest.memory_bytes;
    if (s->has_qulacs && s->qulacs.memory_bytes > max_mem) max_mem = s->qulacs.memory_bytes;

    const char *unit  = (max_mem >= 1024 * 1024) ? "MB" : "KB";
    double      scale = (max_mem >= 1024 * 1024) ? 1024.0 * 1024.0 : 1024.0;

    printf("  Memory: packed=%.1f %s", (double)s->packed.memory_bytes / scale, unit);
    if (s->has_tiled  && s->tiled.memory_bytes  > 0)
        printf(", tiled=%.1f %s",  (double)s->tiled.memory_bytes  / scale, unit);
    if (s->has_dense  && s->dense.memory_bytes  > 0)
        printf(", dense=%.1f %s",  (double)s->dense.memory_bytes  / scale, unit);
    if (s->has_naive  && s->naive.memory_bytes  > 0)
        printf(", naive=%.1f %s",  (double)s->naive.memory_bytes  / scale, unit);
    if (s->has_quest  && s->quest.memory_bytes  > 0)
        printf(", QuEST=%.1f %s",  (double)s->quest.memory_bytes  / scale, unit);
    if (s->has_qulacs && s->qulacs.memory_bytes > 0)
        printf(", Qulacs=%.1f %s", (double)s->qulacs.memory_bytes / scale, unit);
    if (s->has_dense && s->packed.memory_bytes > 0 && s->dense.memory_bytes > 0) {
        double savings = 100.0 * (1.0 - (double)s->packed.memory_bytes / (double)s->dense.memory_bytes);
        printf(" (packed saves %.0f%%)", savings);
    }
    printf("\n");
}

/*
 * =====================================================================================================================
 * Command-line parsing
 * =====================================================================================================================
 */

static void print_usage(const char *prog) {
    printf("Usage: %s [options]\n\n", prog);
    printf("Options:\n");
    printf("  --min-qubits N   Minimum qubit count, 1..30 (default: %d)\n", BENCH_DEFAULT_MIN_QUBITS);
    printf("  --max-qubits N   Maximum qubit count, 1..30 (default: %d)\n", BENCH_DEFAULT_MAX_QUBITS);
    printf("  --step N         Qubit count increment, >=1 (default: %d)\n", BENCH_DEFAULT_STEP);
    printf("  --iterations N   Number of timed iterations, >=1 (default: %d)\n", BENCH_DEFAULT_ITERATIONS);
    printf("  --warmup N       Warm-up iterations, >=0 (default: %d)\n", BENCH_DEFAULT_WARMUP);
    printf("  --csv            Output in CSV format\n");
    printf("  --pgfplots       Output in pgfplots-compatible .dat format\n");
    printf("  --verbose        Show memory usage per result\n");
    printf("  --help           Show this help\n");
}

bench_options_t bench_parse_options(int argc, char *argv[]) {
    bench_options_t opts = {
        .min_qubits     = BENCH_DEFAULT_MIN_QUBITS,
        .max_qubits     = BENCH_DEFAULT_MAX_QUBITS,
        .step           = BENCH_DEFAULT_STEP,
        .iterations     = BENCH_DEFAULT_ITERATIONS,
        .warmup         = BENCH_DEFAULT_WARMUP,
        .csv_output     = 0,
        .pgfplots_output = 0,
        .verbose        = 0
    };

    for (int i = 1; i < argc; ++i) {
        if        (strcmp(argv[i], "--min-qubits")  == 0 && i + 1 < argc) {
            opts.min_qubits  = (qubit_t)atoi(argv[++i]);
        } else if (strcmp(argv[i], "--max-qubits")  == 0 && i + 1 < argc) {
            opts.max_qubits  = (qubit_t)atoi(argv[++i]);
        } else if (strcmp(argv[i], "--step")        == 0 && i + 1 < argc) {
            opts.step        = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--iterations")  == 0 && i + 1 < argc) {
            opts.iterations  = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--warmup")      == 0 && i + 1 < argc) {
            opts.warmup      = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--csv")         == 0) {
            opts.csv_output  = 1;
        } else if (strcmp(argv[i], "--pgfplots")    == 0) {
            opts.pgfplots_output = 1;
        } else if (strcmp(argv[i], "--verbose")     == 0) {
            opts.verbose     = 1;
        } else if (strcmp(argv[i], "--help")        == 0) {
            print_usage(argv[0]);
            exit(0);
        }
    }

    /* Input validation */
    if (opts.min_qubits < 1 || opts.min_qubits > 30) {
        fprintf(stderr, "error: --min-qubits must be between 1 and 30\n");
        exit(1);
    }
    if (opts.max_qubits < opts.min_qubits || opts.max_qubits > 30) {
        fprintf(stderr, "error: --max-qubits must be between --min-qubits (%u) and 30\n",
                opts.min_qubits);
        exit(1);
    }
    if (opts.step <= 0) {
        fprintf(stderr, "error: --step must be >= 1\n");
        exit(1);
    }
    if (opts.iterations <= 0) {
        fprintf(stderr, "error: --iterations must be >= 1\n");
        exit(1);
    }
    if (opts.warmup < 0) {
        fprintf(stderr, "error: --warmup must be >= 0\n");
        exit(1);
    }

    return opts;
}

/*
 * =====================================================================================================================
 * Gate tables
 * =====================================================================================================================
 */

typedef struct {
    const char *name;
    void (*fn)(state_t*, qubit_t);
    const cplx_t *mat;
} gate_def_t;

static const gate_def_t gates[] = {
    {"X", x, XMAT},
    {"H", h, HMAT},
    {"Z", z, ZMAT},
    {NULL, NULL, NULL}
};

typedef struct {
    const char *name;
    void (*fn)(state_t*, qubit_t, qubit_t);
    const cplx_t *mat;
    int has_tiled;
} gate_def_2q_t;

static const gate_def_2q_t gates_2q[] = {
    {"CX",   cx,        CXMAT,   1},
    {"SWAP", swap_gate, SWAPMAT, 1},
    {NULL, NULL, NULL, 0}
};

#define NUM_1Q_GATES    3
#define NUM_2Q_GATES    2
#define NUM_GATES       (NUM_1Q_GATES + NUM_2Q_GATES)
#define MAX_QUBIT_CONFIGS 32
#define NUM_METHODS     6   /* qlib_packed, qlib_tiled, blas, naive, quest, qulacs */

static const char *all_gate_names[NUM_GATES] = {"X", "H", "Z", "CX", "SWAP"};

/*
 * =====================================================================================================================
 * pgfplots output
 * =====================================================================================================================
 */

typedef struct {
    qubit_t qubits[MAX_QUBIT_CONFIGS];
    int     num_configs;
    double  time_ms[NUM_GATES][MAX_QUBIT_CONFIGS][NUM_METHODS];
    size_t  memory[NUM_GATES][MAX_QUBIT_CONFIGS][NUM_METHODS];
} pgfplots_data_t;

enum { M_QLIB = 0, M_QLIB_TILED = 1, M_BLAS = 2, M_NAIVE = 3, M_QUEST = 4, M_QULACS = 5 };

static void print_pgfplots_output(const pgfplots_data_t *pgf) {
    const char *method_names[]  = {"qlib", "qlib_tiled", "blas", "quest", "qulacs"};
    int         method_indices[] = {M_QLIB, M_QLIB_TILED, M_BLAS, M_QUEST, M_QULACS};
    int         num_methods      = 5;

    for (int g = 0; g < NUM_GATES; ++g) {
        printf("# %s gate timing (ms)\n", all_gate_names[g]);
        printf("%-8s", "qubits");
        for (int m = 0; m < num_methods; ++m) printf("  %-12s", method_names[m]);
        printf("\n");

        for (int q = 0; q < pgf->num_configs; ++q) {
            printf("%-8u", pgf->qubits[q]);
            for (int m = 0; m < num_methods; ++m) {
                double val = pgf->time_ms[g][q][method_indices[m]];
                if (val > 0)
                    printf("  %-12.6f", val);
                else
                    printf("  %-12s", "nan");
            }
            printf("\n");
        }
        printf("\n");
    }

    printf("# Memory usage (MB)\n");
    printf("%-8s", "qubits");
    for (int m = 0; m < num_methods; ++m) printf("  %-12s", method_names[m]);
    printf("\n");

    for (int q = 0; q < pgf->num_configs; ++q) {
        printf("%-8u", pgf->qubits[q]);
        for (int m = 0; m < num_methods; ++m) {
            size_t bytes = pgf->memory[0][q][method_indices[m]];
            if (bytes > 0)
                printf("  %-12.3f", (double)bytes / (1024.0 * 1024.0));
            else
                printf("  %-12s", "nan");
        }
        printf("\n");
    }
}

/*
 * =====================================================================================================================
 * Main
 * =====================================================================================================================
 */

int main(int argc, char *argv[]) {
    bench_options_t opts = bench_parse_options(argc, argv);

    /* Allocate pgfplots storage only when --pgfplots is requested */
    pgfplots_data_t *pgf = NULL;
    if (opts.pgfplots_output) {
        pgf = calloc(1, sizeof(pgfplots_data_t));
        if (!pgf) {
            fprintf(stderr, "error: failed to allocate pgfplots buffer\n");
            return 1;
        }
    }

    if (opts.csv_output) {
        bench_print_csv_header();
    } else if (!opts.pgfplots_output) {
        printf("\nMixed State Gate Benchmarks\n");
        printf("===========================\n");
        printf("Comparing: qlib packed sparse vs BLAS dense vs naive loop");
#ifdef WITH_QUEST
        printf(" vs QuEST");
#endif
#ifdef WITH_QULACS
        printf(" vs Qulacs");
#endif
        printf("\n");
        printf("Iterations: %d, Warm-up: %d\n", opts.iterations, opts.warmup);
#ifdef WITH_QUEST
        printf("QuEST: enabled\n");
#endif
#ifdef WITH_QULACS
        printf("Qulacs: enabled\n");
#endif
        printf("\n");
    }

    int qi = 0;

    for (qubit_t qubits = opts.min_qubits; qubits <= opts.max_qubits; qubits += opts.step) {
        dim_t dim    = (dim_t)1 << qubits;
        int run_dense = (qubits <= 8);  /* O(N³) BLAS/naive: practical only up to 8 qubits */

        if (pgf && qi < MAX_QUBIT_CONFIGS)
            pgf->qubits[qi] = qubits;

        if (!opts.csv_output && !opts.pgfplots_output)
            printf("qubits: %u, dim: %lld\n", qubits, (long long)dim);

        /* ── One-qubit gates ─────────────────────────────────────────────────── */
        for (int g = 0; gates[g].name != NULL; ++g) {
            if (!opts.csv_output && !opts.pgfplots_output)
                printf("\nGate: %s\n", gates[g].name);

            bench_summary_t s = {0};
            s.naive_scale = 10.0;  /* naive runs at iterations/10 due to O(N³) cost */

            srand(42); s.packed = bench_qlib_packed(qubits, gates[g].name, gates[g].fn,
                                                     opts.iterations, opts.warmup);
            srand(42); s.tiled  = bench_qlib_tiled(qubits, gates[g].name, gates[g].fn,
                                                    opts.iterations, opts.warmup);
            s.has_tiled = (s.tiled.time_ms > 0);

            if (run_dense) {
                srand(42); s.dense = bench_blas_dense(qubits, gates[g].name, gates[g].mat,
                                                       opts.iterations, opts.warmup);
                s.has_dense = 1;

                srand(42); s.naive = bench_naive_loop(qubits, gates[g].name, gates[g].mat,
                                                       opts.iterations / 10,
                                                       opts.warmup    / 10);
                s.has_naive = 1;
            }

#ifdef WITH_QUEST
            s.quest     = bench_quest(qubits, gates[g].name, opts.iterations, opts.warmup);
            s.has_quest = (s.quest.time_ms > 0);
#endif
#ifdef WITH_QULACS
            s.qulacs     = bench_qulacs(qubits, gates[g].name, opts.iterations, opts.warmup);
            s.has_qulacs = (s.qulacs.time_ms > 0);
#endif

            if (pgf && qi < MAX_QUBIT_CONFIGS) {
                pgf->time_ms[g][qi][M_QLIB]       = s.packed.time_ms;
                pgf->time_ms[g][qi][M_QLIB_TILED]  = s.tiled.time_ms;
                pgf->time_ms[g][qi][M_BLAS]        = s.dense.time_ms;
                pgf->time_ms[g][qi][M_NAIVE]       = s.naive.time_ms;
                pgf->memory[g][qi][M_QLIB]         = s.packed.memory_bytes;
                pgf->memory[g][qi][M_QLIB_TILED]   = s.tiled.memory_bytes;
                pgf->memory[g][qi][M_BLAS]         = s.dense.memory_bytes;
                pgf->memory[g][qi][M_NAIVE]        = s.naive.memory_bytes;
#ifdef WITH_QUEST
                pgf->time_ms[g][qi][M_QUEST]  = s.quest.time_ms;
                pgf->memory[g][qi][M_QUEST]   = s.quest.memory_bytes;
#endif
#ifdef WITH_QULACS
                pgf->time_ms[g][qi][M_QULACS] = s.qulacs.time_ms;
                pgf->memory[g][qi][M_QULACS]  = s.qulacs.memory_bytes;
#endif
            }

            if (opts.csv_output) {
                bench_print_csv(&s.packed);
                bench_print_csv(&s.tiled);
                if (s.has_dense) bench_print_csv(&s.dense);
                if (s.has_naive) bench_print_csv(&s.naive);
#ifdef WITH_QUEST
                if (s.has_quest) bench_print_csv(&s.quest);
#endif
#ifdef WITH_QULACS
                if (s.has_qulacs) bench_print_csv(&s.qulacs);
#endif
            } else if (!opts.pgfplots_output) {
                bench_print_result(&s.packed, opts.verbose);
                bench_print_result(&s.tiled,  opts.verbose);
                if (s.has_dense) bench_print_result(&s.dense, opts.verbose);
                if (s.has_naive) bench_print_result(&s.naive, opts.verbose);
#ifdef WITH_QUEST
                if (s.has_quest) bench_print_result(&s.quest, opts.verbose);
#endif
#ifdef WITH_QULACS
                if (s.has_qulacs) bench_print_result(&s.qulacs, opts.verbose);
#endif
                print_speedup_and_memory(&s);
            }
        }

        /* ── Two-qubit gates ─────────────────────────────────────────────────── */
        for (int g = 0; gates_2q[g].name != NULL; ++g) {
            if (!opts.csv_output && !opts.pgfplots_output)
                printf("\nGate: %s\n", gates_2q[g].name);

            bench_summary_t s = {0};
            s.naive_scale = 1.0;

            srand(42); s.packed = bench_qlib_packed_2q(qubits, gates_2q[g].name, gates_2q[g].fn,
                                                        opts.iterations, opts.warmup);

            if (gates_2q[g].has_tiled) {
                srand(42); s.tiled = bench_qlib_tiled_2q(qubits, gates_2q[g].name, gates_2q[g].fn,
                                                          opts.iterations, opts.warmup);
                s.has_tiled = (s.tiled.time_ms > 0);
            }

            if (run_dense) {
                srand(42); s.dense = bench_blas_dense_2q(qubits, gates_2q[g].name, gates_2q[g].mat,
                                                          opts.iterations, opts.warmup);
                s.has_dense = 1;
            }

#ifdef WITH_QUEST
            s.quest     = bench_quest(qubits, gates_2q[g].name, opts.iterations, opts.warmup);
            s.has_quest = (s.quest.time_ms > 0);
#endif
#ifdef WITH_QULACS
            s.qulacs     = bench_qulacs(qubits, gates_2q[g].name, opts.iterations, opts.warmup);
            s.has_qulacs = (s.qulacs.time_ms > 0);
#endif

            int pg = NUM_1Q_GATES + g;
            if (pgf && qi < MAX_QUBIT_CONFIGS) {
                pgf->time_ms[pg][qi][M_QLIB]       = s.packed.time_ms;
                pgf->time_ms[pg][qi][M_QLIB_TILED]  = s.tiled.time_ms;
                pgf->time_ms[pg][qi][M_BLAS]        = s.dense.time_ms;
                pgf->memory[pg][qi][M_QLIB]         = s.packed.memory_bytes;
                pgf->memory[pg][qi][M_QLIB_TILED]   = s.tiled.memory_bytes;
                pgf->memory[pg][qi][M_BLAS]         = s.dense.memory_bytes;
#ifdef WITH_QUEST
                pgf->time_ms[pg][qi][M_QUEST]  = s.quest.time_ms;
                pgf->memory[pg][qi][M_QUEST]   = s.quest.memory_bytes;
#endif
#ifdef WITH_QULACS
                pgf->time_ms[pg][qi][M_QULACS] = s.qulacs.time_ms;
                pgf->memory[pg][qi][M_QULACS]  = s.qulacs.memory_bytes;
#endif
            }

            if (opts.csv_output) {
                bench_print_csv(&s.packed);
                if (s.has_tiled) bench_print_csv(&s.tiled);
                if (s.has_dense) bench_print_csv(&s.dense);
#ifdef WITH_QUEST
                if (s.has_quest) bench_print_csv(&s.quest);
#endif
#ifdef WITH_QULACS
                if (s.has_qulacs) bench_print_csv(&s.qulacs);
#endif
            } else if (!opts.pgfplots_output) {
                bench_print_result(&s.packed, opts.verbose);
                if (s.has_tiled) bench_print_result(&s.tiled,  opts.verbose);
                if (s.has_dense) bench_print_result(&s.dense,  opts.verbose);
#ifdef WITH_QUEST
                if (s.has_quest) bench_print_result(&s.quest, opts.verbose);
#endif
#ifdef WITH_QULACS
                if (s.has_qulacs) bench_print_result(&s.qulacs, opts.verbose);
#endif
                print_speedup_and_memory(&s);
            }
        }

        if (!opts.csv_output && !opts.pgfplots_output)
            printf("\n" "─────────────────────────────────────────────────\n");

        if (pgf)
            qi++;
    }

    if (pgf) {
        pgf->num_configs = qi;
        print_pgfplots_output(pgf);
        free(pgf);
    }

#ifdef WITH_QUEST
    bench_quest_cleanup();
#endif

    return 0;
}
