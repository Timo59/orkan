/**
 * @file bench_quest.c
 * @brief QuEST framework benchmark integration
 *
 * QuEST (Quantum Exact Simulation Toolkit) is a high-performance,
 * multithreaded, distributed quantum simulator with density matrix support.
 *
 * Homepage: https://quest.qtechtheory.org/
 * GitHub: https://github.com/QuEST-Kit/QuEST
 *
 * To enable: QuEST is auto-detected at configure time. Pass -DQuEST_DIR=/path/to/quest
 *            if QuEST is not in the default extern/QuEST/ location.
 *
 * Uses an init-once pattern: bench_quest_init() calls initQuESTEnv() once at
 * benchmark start. Each bench_quest() call creates a fresh Qureg with
 * createDensityQureg, runs the benchmark, then destroys it with destroyQureg.
 * bench_quest_cleanup() calls finalizeQuESTEnv() once at program exit.
 *
 * This respects QuEST's constraint that initQuESTEnv() cannot be called again
 * after finalizeQuESTEnv() (enforced by a permanent static flag in QuEST), while
 * avoiding the overhead of forking a child process per benchmark call.
 */

#ifdef WITH_QUEST

/* quest.h must precede bench.h: bench.h's WITH_QUEST block declares
 * bench_quest_at / bench_quest_2q_at using Qureg, so the type must be
 * visible before bench.h is parsed. */
#include "quest/include/quest.h"
#include "bench.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


void bench_quest_init(void) {
    initQuESTEnv();
}

void bench_quest_cleanup(void) {
    finalizeQuESTEnv();
}

/**
 * @brief Run QuEST density matrix benchmark in-process.
 *
 * Requires bench_quest_init() to have been called. Creates and destroys a
 * Qureg per call. Performs `runs` independent timed rounds and returns the
 * full statistical summary.
 */
bench_result_t bench_quest(qubit_t qubits, const char *gate_name,
                           int iterations, int warmup, int runs) {
    bench_result_t result = {0};
    result.qubits     = qubits;
    result.dim        = (dim_t)1 << qubits;
    result.gate_name  = gate_name;
    result.method     = "quest";
    result.iterations = iterations;
    result.runs       = runs;

    Qureg rho = createDensityQureg((int)qubits);

    /*
     * Use initDebugState() to ensure all memory pages are touched.
     * initZeroState() may leave pages unmapped due to lazy allocation.
     */
    initDebugState(rho);

    /*
     * Memory: theoretical density matrix size (dim x dim complex doubles).
     * Matches the methodology used by qlib and Qulacs — all three report
     * only the state storage footprint, not framework overhead.
     */
    dim_t dim = (dim_t)1 << qubits;
    result.memory_bytes = (size_t)dim * (size_t)dim * sizeof(double) * 2; /* complex<double> */

    /* Dispatch: single-qubit gates */
    void (*gate_fn_1q)(Qureg, int) = NULL;
    if (strcmp(gate_name, "X") == 0) {
        gate_fn_1q = applyPauliX;
    } else if (strcmp(gate_name, "H") == 0) {
        gate_fn_1q = applyHadamard;
    } else if (strcmp(gate_name, "Z") == 0) {
        gate_fn_1q = applyPauliZ;
    }

    /* Dispatch: two-qubit gates */
    void (*gate_fn_2q)(Qureg, int, int) = NULL;
    if (strcmp(gate_name, "CX") == 0) {
        gate_fn_2q = applyControlledPauliX;
    } else if (strcmp(gate_name, "SWAP") == 0) {
        gate_fn_2q = applySwap;
    }

    if (!gate_fn_1q && !gate_fn_2q) {
        destroyQureg(rho);
        return result;
    }

    double run_times[200];  /* bounded by --runs validation (max 200) */
    int cnt = (runs < 200) ? runs : 200;

    if (gate_fn_1q) {
        /* Single-qubit gate: `runs` independent timed rounds */
        for (int i = 0; i < warmup; ++i)
            gate_fn_1q(rho, (int)(i % qubits));

        for (int r = 0; r < cnt; ++r) {
            uint64_t start = bench_time_ns();
            for (int i = 0; i < iterations; ++i)
                for (qubit_t t = 0; t < qubits; ++t)
                    gate_fn_1q(rho, (int)t);
            run_times[r] = bench_ns_to_ms(bench_time_ns() - start);
        }

        bench_run_stats_t stats = bench_compute_stats(run_times, cnt);
        result.time_ms        = stats.mean;
        result.time_ms_std    = stats.std_dev;
        result.time_ms_min    = stats.min;
        result.time_ms_median = stats.median;
        result.time_ms_cv     = stats.cv;
        result.ops_per_sec    = (stats.mean > 0.0) ? (double)(iterations * qubits) / (stats.mean / 1000.0) : 0.0;
    } else {
        /* Two-qubit gate: all ordered pairs (q1 < q2) */
        qubit_t q1s[MAX_PAIRS], q2s[MAX_PAIRS];
        int num_pairs = build_all_pairs(qubits, q1s, q2s);

        if (num_pairs == 0) {
            destroyQureg(rho);
            return result;
        }

        for (int i = 0; i < warmup; ++i)
            gate_fn_2q(rho, (int)q1s[i % num_pairs], (int)q2s[i % num_pairs]);

        for (int r = 0; r < cnt; ++r) {
            uint64_t start = bench_time_ns();
            for (int i = 0; i < iterations; ++i)
                for (int p = 0; p < num_pairs; ++p)
                    gate_fn_2q(rho, (int)q1s[p], (int)q2s[p]);
            run_times[r] = bench_ns_to_ms(bench_time_ns() - start);
        }

        bench_run_stats_t stats = bench_compute_stats(run_times, cnt);
        result.time_ms        = stats.mean;
        result.time_ms_std    = stats.std_dev;
        result.time_ms_min    = stats.min;
        result.time_ms_median = stats.median;
        result.time_ms_cv     = stats.cv;
        result.ops_per_sec    = (stats.mean > 0.0) ? (double)(iterations * num_pairs) / (stats.mean / 1000.0) : 0.0;
    }

    destroyQureg(rho);
    return result;
}

/* ---------------------------------------------------------------------------
 * Callback context types for the per-qubit harness
 * ---------------------------------------------------------------------------*/

typedef struct {
    Qureg qureg;
    void (*fn)(Qureg, int);
    int target;
} quest_cb1q_ctx_t;

static void quest_cb1q(void *ctx) {
    quest_cb1q_ctx_t *c = ctx;
    c->fn(c->qureg, c->target);
}

typedef struct {
    Qureg qureg;
    void (*fn)(Qureg, int, int);
    int q1;
    int q2;
} quest_cb2q_ctx_t;

static void quest_cb2q(void *ctx) {
    quest_cb2q_ctx_t *c = ctx;
    c->fn(c->qureg, c->q1, c->q2);
}

/**
 * @brief Per-qubit QuEST 1Q benchmark: time a single gate on a fixed target qubit.
 *
 * The caller owns and pre-allocates *qureg. State is re-initialised with
 * initDebugState before warmup to ensure all pages are touched. Uses the
 * bench_harness_t / bench_run_timed infrastructure for warmup and timing.
 */
bench_result_perq_t bench_quest_at(qubit_t qubits, const char *gate_name,
    qubit_t target, int iterations, int runs, Qureg *qureg)
{
    bench_result_perq_t result = {0};
    result.qubits    = qubits;
    result.dim       = (dim_t)1 << qubits;
    result.gate_name = gate_name;
    result.method    = "quest";
    result.target    = target;
    result.target2   = QUBIT_NONE;
    result.is_2q     = 0;

    /* Resolve gate name to function pointer */
    void (*fn)(Qureg, int) = NULL;
    if (strcmp(gate_name, "X") == 0) {
        fn = applyPauliX;
    } else if (strcmp(gate_name, "H") == 0) {
        fn = applyHadamard;
    } else if (strcmp(gate_name, "Z") == 0) {
        fn = applyPauliZ;
    }

    if (!fn) {
        return result;
    }

    initDebugState(*qureg);

    quest_cb1q_ctx_t ctx = { *qureg, fn, (int)target };
    bench_harness_t h = { quest_cb1q, &ctx, iterations, runs };

    double run_times[BENCH_MAX_RUNS];
    assert(runs <= BENCH_MAX_RUNS);
    bench_run_timed(&h, run_times, (int)qubits);

    bench_run_stats_t stats = bench_compute_stats(run_times, runs);
    bench_fill_perq_stats(&result, &stats, iterations);
    result.runs = runs;

    dim_t dim = (dim_t)1 << qubits;
    result.memory_bytes = (size_t)dim * (size_t)dim * sizeof(double) * 2;

    return result;
}

/**
 * @brief Per-qubit QuEST 2Q benchmark: time a single gate on a fixed (q1, q2) pair.
 *
 * For CX: q1=control, q2=target, matching applyControlledPauliX(qureg, control, target).
 * The caller owns and pre-allocates *qureg.
 */
bench_result_perq_t bench_quest_2q_at(qubit_t qubits, const char *gate_name,
    qubit_t q1, qubit_t q2, int iterations, int runs, Qureg *qureg)
{
    bench_result_perq_t result = {0};
    result.qubits    = qubits;
    result.dim       = (dim_t)1 << qubits;
    result.gate_name = gate_name;
    result.method    = "quest";
    result.target    = q1;
    result.target2   = q2;
    result.is_2q     = 1;

    /* Resolve gate name to function pointer */
    void (*fn)(Qureg, int, int) = NULL;
    if (strcmp(gate_name, "CX") == 0) {
        fn = applyControlledPauliX;
    } else if (strcmp(gate_name, "SWAP") == 0) {
        fn = applySwap;
    }

    if (!fn) {
        return result;
    }

    initDebugState(*qureg);

    quest_cb2q_ctx_t ctx = { *qureg, fn, (int)q1, (int)q2 };
    bench_harness_t h = { quest_cb2q, &ctx, iterations, runs };

    double run_times[BENCH_MAX_RUNS];
    assert(runs <= BENCH_MAX_RUNS);
    bench_run_timed(&h, run_times, (int)qubits);

    bench_run_stats_t stats = bench_compute_stats(run_times, runs);
    bench_fill_perq_stats(&result, &stats, iterations);
    result.runs = runs;

    dim_t dim = (dim_t)1 << qubits;
    result.memory_bytes = (size_t)dim * (size_t)dim * sizeof(double) * 2;

    return result;
}

#endif /* WITH_QUEST */
