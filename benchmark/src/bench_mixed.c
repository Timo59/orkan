/**
 * @file bench_mixed.c
 * @brief qlib packed/tiled adapter functions for mixed state benchmarks
 *
 * Contains:
 * 1. qlib packed sparse - Production implementation with packed lower-triangular storage
 * 2. qlib tiled         - Production implementation with tiled lower-triangular storage
 */

#include "bench.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

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

/*
 * Build all ordered qubit pairs (q1 < q2) for an n-qubit system.
 * Shared across all benchmark adapters (bench_baselines.c, bench_quest.c, bench_qulacs.cpp);
 * declared in bench.h. Returns the number of pairs written (at most MAX_PAIRS).
 */
int build_all_pairs(qubit_t qubits, qubit_t *q1_out, qubit_t *q2_out) {
    int n = 0;
    for (qubit_t q1 = 0; q1 < qubits && n < MAX_PAIRS; ++q1) {
        for (qubit_t q2 = q1 + 1; q2 < qubits && n < MAX_PAIRS; ++q2) {
            q1_out[n] = q1;
            q2_out[n] = q2;
            n++;
        }
    }
    return n;
}

/*
 * =====================================================================================================================
 * Internal timing helpers
 * =====================================================================================================================
 */

/**
 * @brief Descriptor for a qlib timing run (single- or two-qubit).
 *
 * Captures the parts that differ between packed and tiled variants so the
 * shared timing loop can be called from thin wrappers.
 */
typedef struct {
    state_type_t  state_type;   /**< MIXED_PACKED or MIXED_TILED */
    size_t        memory_bytes; /**< Pre-computed memory footprint */
    const char   *method;       /**< "qlib_packed" or "qlib_tiled" */
} qlib_bench_desc_t;

/**
 * @brief Shared timing loop for single-qubit qlib benchmarks.
 *
 * Initialises state, runs warmup, executes the timed loop across all qubits,
 * computes statistics, and populates result.
 */
static bench_result_t run_1q_timing(qubit_t qubits, const char *gate_name,
                                    void (*gate_fn)(state_t*, qubit_t),
                                    int iterations, int warmup, int runs,
                                    const qlib_bench_desc_t *desc) {
    bench_result_t result = {0};
    result.qubits       = qubits;
    result.dim          = (dim_t)1 << qubits;
    result.gate_name    = gate_name;
    result.method       = desc->method;
    result.iterations   = iterations;
    result.memory_bytes = desc->memory_bytes;

    state_t state = {0};
    init_random_mixed_state(&state, qubits, desc->state_type);
    if (!state.data) return result;

    for (int i = 0; i < warmup; ++i)
        gate_fn(&state, (qubit_t)(i % qubits));

    double *run_times = malloc(runs * sizeof(double));
    if (!run_times) { state_free(&state); return result; }
    for (int r = 0; r < runs; ++r) {
        uint64_t start = bench_time_ns();
        for (int i = 0; i < iterations; ++i)
            for (qubit_t t = 0; t < qubits; ++t)
                gate_fn(&state, t);
        run_times[r] = bench_ns_to_ms(bench_time_ns() - start);
    }
    bench_run_stats_t stats = bench_compute_stats(run_times, runs);
    free(run_times);
    result.time_ms        = stats.mean;
    result.time_ms_std    = stats.std_dev;
    result.time_ms_min    = stats.min;
    result.time_ms_median = stats.median;
    result.time_ms_cv     = stats.cv;
    result.runs           = runs;
    result.ops_per_sec    = (stats.mean > 0.0) ? (double)(iterations * qubits) / (stats.mean / 1000.0) : 0.0;

    state_free(&state);
    return result;
}

/**
 * @brief Shared timing loop for two-qubit qlib benchmarks.
 *
 * Iterates over all ordered pairs (q1 < q2), runs warmup, executes the timed
 * loop, computes statistics, and populates result.
 */
static bench_result_t run_2q_timing(qubit_t qubits, const char *gate_name,
                                    void (*gate_fn)(state_t*, qubit_t, qubit_t),
                                    int iterations, int warmup, int runs,
                                    const qlib_bench_desc_t *desc) {
    bench_result_t result = {0};
    result.qubits       = qubits;
    result.dim          = (dim_t)1 << qubits;
    result.gate_name    = gate_name;
    result.method       = desc->method;
    result.iterations   = iterations;
    result.memory_bytes = desc->memory_bytes;

    state_t state = {0};
    init_random_mixed_state(&state, qubits, desc->state_type);
    if (!state.data) return result;

    qubit_t q1s[MAX_PAIRS], q2s[MAX_PAIRS];
    int num_pairs = build_all_pairs(qubits, q1s, q2s);
    if (num_pairs == 0) { state_free(&state); return result; }

    for (int i = 0; i < warmup; ++i)
        gate_fn(&state, q1s[i % num_pairs], q2s[i % num_pairs]);

    double *run_times = malloc(runs * sizeof(double));
    if (!run_times) { state_free(&state); return result; }
    for (int r = 0; r < runs; ++r) {
        uint64_t start = bench_time_ns();
        for (int i = 0; i < iterations; ++i)
            for (int p = 0; p < num_pairs; ++p)
                gate_fn(&state, q1s[p], q2s[p]);
        run_times[r] = bench_ns_to_ms(bench_time_ns() - start);
    }
    bench_run_stats_t stats = bench_compute_stats(run_times, runs);
    free(run_times);
    result.time_ms        = stats.mean;
    result.time_ms_std    = stats.std_dev;
    result.time_ms_min    = stats.min;
    result.time_ms_median = stats.median;
    result.time_ms_cv     = stats.cv;
    result.runs           = runs;
    result.ops_per_sec    = (stats.mean > 0.0) ? (double)(iterations * num_pairs) / (stats.mean / 1000.0) : 0.0;

    state_free(&state);
    return result;
}

/*
 * =====================================================================================================================
 * Benchmark implementations
 * =====================================================================================================================
 */

bench_result_t bench_qlib_packed(qubit_t qubits, const char *gate_name,
                                  void (*gate_fn)(state_t*, qubit_t),
                                  int iterations, int warmup, int runs) {
    const qlib_bench_desc_t desc = {
        .state_type   = MIXED_PACKED,
        .memory_bytes = bench_packed_size(qubits),
        .method       = "qlib_packed",
    };
    return run_1q_timing(qubits, gate_name, gate_fn, iterations, warmup, runs, &desc);
}

bench_result_t bench_qlib_tiled(qubit_t qubits, const char *gate_name,
                                 void (*gate_fn)(state_t*, qubit_t),
                                 int iterations, int warmup, int runs) {
    const qlib_bench_desc_t desc = {
        .state_type   = MIXED_TILED,
        .memory_bytes = bench_tiled_size(qubits),
        .method       = "qlib_tiled",
    };
    return run_1q_timing(qubits, gate_name, gate_fn, iterations, warmup, runs, &desc);
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
                                      int iterations, int warmup, int runs) {
    const qlib_bench_desc_t desc = {
        .state_type   = MIXED_PACKED,
        .memory_bytes = bench_packed_size(qubits),
        .method       = "qlib_packed",
    };
    return run_2q_timing(qubits, gate_name, gate_fn, iterations, warmup, runs, &desc);
}

bench_result_t bench_qlib_tiled_2q(qubit_t qubits, const char *gate_name,
                                     void (*gate_fn)(state_t*, qubit_t, qubit_t),
                                     int iterations, int warmup, int runs) {
    const qlib_bench_desc_t desc = {
        .state_type   = MIXED_TILED,
        .memory_bytes = bench_tiled_size(qubits),
        .method       = "qlib_tiled",
    };
    return run_2q_timing(qubits, gate_name, gate_fn, iterations, warmup, runs, &desc);
}

/*
 * =====================================================================================================================
 * Per-qubit (_at) infrastructure
 * =====================================================================================================================
 */

/**
 * @brief Re-fill state->data with random values without re-allocating.
 *
 * Uses bench_fill_random() for OS-seeded PRNG; safe to call repeatedly
 * between timed runs to prevent the optimizer from dead-code-eliminating
 * the gate computation.
 */
static void state_reinit_random(state_t *state) {
    dim_t len = state_len(state);
    bench_fill_random(state->data, (size_t)len * sizeof(cplx_t));
}

/* --- Callback context structs and shims ---------------------------------- */

typedef struct {
    state_t  *state;
    void    (*fn)(state_t*, qubit_t);
    qubit_t   target;
} cb1q_ctx_t;

static void cb1q(void *ctx) {
    cb1q_ctx_t *c = (cb1q_ctx_t *)ctx;
    c->fn(c->state, c->target);
}

typedef struct {
    state_t  *state;
    void    (*fn)(state_t*, qubit_t, qubit_t);
    qubit_t   q1;
    qubit_t   q2;
} cb2q_ctx_t;

static void cb2q(void *ctx) {
    cb2q_ctx_t *c = (cb2q_ctx_t *)ctx;
    c->fn(c->state, c->q1, c->q2);
}

/* --- Internal helpers for the four public _at functions ------------------ */

static bench_result_perq_t run_1q_at(qubit_t qubits, const char *gate_name,
                                      void (*gate_fn)(state_t*, qubit_t),
                                      qubit_t target, int iterations, int runs,
                                      state_t *state, const char *method,
                                      size_t memory_bytes) {
    bench_result_perq_t result = {0};
    result.qubits       = qubits;
    result.dim          = (dim_t)1 << qubits;
    result.gate_name    = gate_name;
    result.method       = method;
    result.target       = target;
    result.target2      = QUBIT_NONE;
    result.is_2q        = 0;
    result.memory_bytes = memory_bytes;

    assert(runs <= BENCH_MAX_RUNS);

    state_reinit_random(state);

    cb1q_ctx_t ctx = { .state = state, .fn = gate_fn, .target = target };
    bench_harness_t h = {
        .call       = cb1q,
        .ctx        = &ctx,
        .iterations = iterations,
        .runs       = runs,
    };

    double run_times[BENCH_MAX_RUNS];
    bench_run_timed(&h, run_times, qubits);

    bench_run_stats_t stats = bench_compute_stats(run_times, runs);
    bench_fill_perq_stats(&result, &stats, iterations);
    result.runs = runs;
    return result;
}

static bench_result_perq_t run_2q_at(qubit_t qubits, const char *gate_name,
                                      void (*gate_fn)(state_t*, qubit_t, qubit_t),
                                      qubit_t q1, qubit_t q2, int iterations, int runs,
                                      state_t *state, const char *method,
                                      size_t memory_bytes) {
    bench_result_perq_t result = {0};
    result.qubits       = qubits;
    result.dim          = (dim_t)1 << qubits;
    result.gate_name    = gate_name;
    result.method       = method;
    result.target       = q1;
    result.target2      = q2;
    result.is_2q        = 1;
    result.memory_bytes = memory_bytes;

    assert(runs <= BENCH_MAX_RUNS);

    state_reinit_random(state);

    cb2q_ctx_t ctx = { .state = state, .fn = gate_fn, .q1 = q1, .q2 = q2 };
    bench_harness_t h = {
        .call       = cb2q,
        .ctx        = &ctx,
        .iterations = iterations,
        .runs       = runs,
    };

    double run_times[BENCH_MAX_RUNS];
    bench_run_timed(&h, run_times, qubits);

    bench_run_stats_t stats = bench_compute_stats(run_times, runs);
    bench_fill_perq_stats(&result, &stats, iterations);
    result.runs = runs;
    return result;
}

/*
 * =====================================================================================================================
 * Per-qubit public _at functions
 * =====================================================================================================================
 */

bench_result_perq_t bench_qlib_packed_at(qubit_t qubits, const char *gate_name,
                                          void (*gate_fn)(state_t*, qubit_t),
                                          qubit_t target, int iterations, int runs,
                                          state_t *state) {
    return run_1q_at(qubits, gate_name, gate_fn, target, iterations, runs,
                     state, "qlib_packed", bench_packed_size(qubits));
}

bench_result_perq_t bench_qlib_tiled_at(qubit_t qubits, const char *gate_name,
                                         void (*gate_fn)(state_t*, qubit_t),
                                         qubit_t target, int iterations, int runs,
                                         state_t *state) {
    return run_1q_at(qubits, gate_name, gate_fn, target, iterations, runs,
                     state, "qlib_tiled", bench_tiled_size(qubits));
}

bench_result_perq_t bench_qlib_packed_2q_at(qubit_t qubits, const char *gate_name,
                                              void (*gate_fn)(state_t*, qubit_t, qubit_t),
                                              qubit_t q1, qubit_t q2, int iterations, int runs,
                                              state_t *state) {
    return run_2q_at(qubits, gate_name, gate_fn, q1, q2, iterations, runs,
                     state, "qlib_packed", bench_packed_size(qubits));
}

bench_result_perq_t bench_qlib_tiled_2q_at(qubit_t qubits, const char *gate_name,
                                             void (*gate_fn)(state_t*, qubit_t, qubit_t),
                                             qubit_t q1, qubit_t q2, int iterations, int runs,
                                             state_t *state) {
    return run_2q_at(qubits, gate_name, gate_fn, q1, q2, iterations, runs,
                     state, "qlib_tiled", bench_tiled_size(qubits));
}

/*
 * =====================================================================================================================
 * build_all_ordered_pairs — all permutations (q1 != q2)
 * =====================================================================================================================
 */

int build_all_ordered_pairs(qubit_t qubits, qubit_t *q1s, qubit_t *q2s, int max_pairs) {
    if (qubits > MAX_QUBITS_FOR_ORDERED_PAIRS) return 0;
    int n = 0;
    for (qubit_t i = 0; i < qubits && n < max_pairs; i++) {
        for (qubit_t j = 0; j < qubits && n < max_pairs; j++) {
            if (i != j) {
                q1s[n] = i;
                q2s[n] = j;
                n++;
            }
        }
    }
    return n;
}
