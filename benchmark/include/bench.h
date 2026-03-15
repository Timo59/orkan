#ifndef BENCH_H
#define BENCH_H

/*
 * =====================================================================================================================
 * Includes
 * =====================================================================================================================
 */

#include "q_types.h"
#include "state.h"
#include "gate.h"
#include <math.h>
#include <time.h>
#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
 * =====================================================================================================================
 * Configuration
 * =====================================================================================================================
 */

#define BENCH_DEFAULT_MIN_QUBITS    2
#define BENCH_DEFAULT_MAX_QUBITS    12
#define BENCH_DEFAULT_STEP          1
#define BENCH_DEFAULT_ITERATIONS    1000
#define BENCH_DEFAULT_WARMUP        100
#define BENCH_DEFAULT_RUNS          10   /**< Independent timing runs per (gate, qubits, method) */
#define BENCH_MAX_RUNS              200  /**< Hard upper bound on --runs; shared with bench_compute_stats */

/*
 * =====================================================================================================================
 * Timing serialization barrier
 * =====================================================================================================================
 */

/** @brief Emit a CPU/compiler barrier to prevent reordering across timing points. */
void bench_timing_barrier(void);

/*
 * =====================================================================================================================
 * Random fill helper
 * =====================================================================================================================
 */

/**
 * @brief Fill `n` bytes at `buf` with pseudo-random data seeded from OS entropy.
 *
 * Uses xoshiro256** expansion from a splitmix64-seeded state. Safe to call from
 * both C and C++ translation units (operates on raw uint8_t* bytes only).
 */
void bench_fill_random(void *buf, size_t n);

/*
 * =====================================================================================================================
 * Per-qubit benchmark constants
 * =====================================================================================================================
 */

/* Per-qubit benchmark constants */
#define BENCH_PQ_DEFAULT_RUNS              15      /**< 15 gives honest CI */
#define BENCH_PQ_DEFAULT_ITERATIONS       100      /**< overridden by calibration */
#define BENCH_PQ_PROBE_ITERS              100      /**< iterations in each probe run */
#define BENCH_PQ_PROBE_RUNS               5        /**< runs for calibration probe; take minimum */
#define BENCH_PQ_WARMUP_MS               100.0     /**< minimum warmup duration in ms */
#define MIN_MEASUREMENT_MS                1.0      /**< minimum timed run duration to target */
#define BENCH_PQ_FALLBACK_ITERATIONS      10000    /**< when probe rounds to 0 ns */
#define BENCH_PQ_MAX_CALIBRATED_ITERATIONS 1000000
#define QUBIT_NONE                        ((qubit_t)-1)    /**< max value of qubit_t */
#define MAX_QUBITS_FOR_ORDERED_PAIRS      16
#define BENCH_HUGEPAGE_THRESHOLD          ((size_t)1 << 30)

/*
 * =====================================================================================================================
 * Types
 * =====================================================================================================================
 */

/**
 * @brief Benchmark result for a single (gate, qubits, method) configuration.
 *
 * Timing fields are statistics computed across `runs` independent timed runs,
 * each consisting of `iterations` gate applications.
 */
typedef struct bench_result {
    qubit_t    qubits;          /**< Number of qubits */
    dim_t      dim;             /**< Hilbert space dimension (2^qubits) */
    const char *gate_name;      /**< Name of the gate tested */
    const char *method;         /**< Implementation method name */
    /* Timing statistics across `runs` independent runs */
    double time_ms;             /**< Mean total time per run (ms) */
    double time_ms_std;         /**< Sample std dev across runs, Bessel-corrected (ms) */
    double time_ms_min;         /**< Minimum run time — best-case hardware throughput (ms) */
    double time_ms_median;      /**< Median run time — robust to scheduler spikes (ms) */
    double time_ms_cv;          /**< Coefficient of variation = std/mean × 100 (%) */
    double ops_per_sec;         /**< Gate applications per second, from mean time */
    double time_per_gate_ms;    /**< Mean time normalised per gate application: time_ms / (iterations × sweep_size) */
    size_t memory_bytes;        /**< Memory used for state representation */
    int    iterations;          /**< Gate calls per timed run */
    int    sweep_size;          /**< Gate applications per timed iteration (qubits for 1Q, pairs for 2Q) */
    int    runs;                /**< Number of independent timing runs */
} bench_result_t;

/** @brief Opaque gate callback type for the per-qubit harness */
typedef void (*bench_gate_cb_t)(void *ctx);

/**
 * @brief Timed-run harness: encapsulates a callback, its context, and loop counts.
 *
 * `call` is declared volatile to prevent the compiler from devirtualizing or
 * collapsing repeated invocations across loop iterations.
 */
typedef struct {
    volatile bench_gate_cb_t call;  /**< Opaque callback — volatile forces a fresh load of the
                                     *   function pointer on every iteration, preventing the
                                     *   compiler from hoisting or caching the value across calls */
    void    *ctx;                   /**< Backend context passed to call() */
    int      iterations;            /**< Gate applications per timed run (> 0) */
    int      runs;                  /**< Number of independent timed runs (1..BENCH_MAX_RUNS) */
} bench_harness_t;

/**
 * @brief Per-qubit benchmark result for a single (gate, qubits, target, method) configuration.
 *
 * Extends bench_result_t with per-gate microsecond timings and explicit target qubit fields,
 * enabling per-position analysis of gate application cost.
 */
typedef struct bench_result_perq {
    qubit_t    qubits;
    dim_t      dim;
    const char *gate_name;
    const char *method;
    qubit_t    target;           /**< Target qubit (q1 for 2Q gates) */
    qubit_t    target2;          /**< Second qubit for 2Q gates; QUBIT_NONE for 1Q */
    int        is_2q;
    double     time_ms_mean;
    double     time_ms_std;      /**< Sample stddev (Bessel-corrected); honest n=15 */
    double     time_ms_min;
    double     time_ms_median;
    double     time_ms_cv;
    double     ops_per_sec;
    double     time_per_gate_us;        /**< mean time per gate application (microseconds) */
    double     time_per_gate_us_median;
    double     time_per_gate_us_min;
    size_t     memory_bytes;
    int        iterations;
    int        runs;
} bench_result_perq_t;

/** @brief Command-line options */
typedef struct bench_options {
    qubit_t     min_qubits;
    qubit_t     max_qubits;
    int         step;
    int         iterations;
    int         warmup;
    int         runs;
    int         csv_output;
    int         pgfplots_output;
    int         verbose;
    int         per_qubit;          /**< 1 if --per-qubit is set; enables single-target throughput mode */
    int         runs_explicit;      /**< 1 if --runs was explicitly provided (suppresses per-qubit default override) */
    int         iterations_explicit;/**< 1 if --iterations was explicitly provided (suppresses per-qubit default override) */
    const char *gate_filter;        /**< Gate name filter from --gate (NULL = all gates); only used in --per-qubit mode */
} bench_options_t;

/*
 * =====================================================================================================================
 * Timing utilities
 * =====================================================================================================================
 */

/** @brief Get current time in nanoseconds */
static inline uint64_t bench_time_ns(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (uint64_t)ts.tv_sec * 1000000000ULL + (uint64_t)ts.tv_nsec;
}

/** @brief Convert nanoseconds to milliseconds */
static inline double bench_ns_to_ms(uint64_t ns) {
    return (double)ns / 1000000.0;
}

/*
 * =====================================================================================================================
 * Statistical helpers
 * =====================================================================================================================
 */

/** @brief Statistical summary computed from K independent timing runs */
typedef struct {
    double mean;     /**< Arithmetic mean (ms) */
    double std_dev;  /**< Sample standard deviation, Bessel-corrected (ms); 0 if runs==1 */
    double min;      /**< Minimum run time — best-case hardware throughput (ms) */
    double median;   /**< Median run time — robust to right-skewed timing distributions (ms) */
    double cv;       /**< Coefficient of variation = std_dev/mean × 100 (%); stability indicator */
} bench_run_stats_t;

/**
 * @brief Compute timing statistics from an array of K run times.
 *
 * Computes mean, sample standard deviation (Bessel's correction), minimum,
 * median, and CV from K independent timing measurements.
 *
 * @param times  Array of K run times in milliseconds
 * @param n      Number of runs (K); must be >= 1
 * @return       Populated bench_run_stats_t
 */
bench_run_stats_t bench_compute_stats(const double *times, int n);

/**
 * @brief Execute a warmed-up timed benchmark using a bench_harness_t.
 *
 * Warms up by invoking the callback until BENCH_PQ_WARMUP_MS has elapsed,
 * then records `h->runs` independent timing samples into `out` (ms each).
 *
 * @param h      Populated harness (call, ctx, iterations, runs).
 * @param out    Output array of at least h->runs doubles (run times in ms).
 * @param qubits Number of qubits — controls warmup batch size.
 * @return 0 on success.
 */
int bench_run_timed(const bench_harness_t *h, double *out, int qubits);

/**
 * @brief Populate per-qubit timing stats in `r` from a bench_run_stats_t.
 *
 * @param r          Destination result struct (all timing fields and runs populated on return).
 * @param s          Statistics computed by bench_compute_stats().
 * @param iterations Gate applications per timed run (used to derive per-gate times).
 * @param runs       Number of timed runs; written to r->runs.
 */
void bench_fill_perq_stats(bench_result_perq_t *r,
                            const bench_run_stats_t *s,
                            int iterations, int runs);

/*
 * =====================================================================================================================
 * Memory utilities
 * =====================================================================================================================
 */

/** @brief Calculate packed storage size for n-qubit mixed state */
static inline size_t bench_packed_size(qubit_t qubits) {
    dim_t dim = (dim_t)1 << qubits;
    return (size_t)(dim * (dim + 1) / 2) * sizeof(cplx_t);
}

/** @brief Calculate dense storage size for n-qubit mixed state */
static inline size_t bench_dense_size(qubit_t qubits) {
    dim_t dim = (dim_t)1 << qubits;
    return (size_t)(dim * dim) * sizeof(cplx_t);
}

/** @brief Calculate tiled storage size for n-qubit mixed state */
static inline size_t bench_tiled_size(qubit_t qubits) {
    dim_t dim = (dim_t)1 << qubits;
    dim_t n_tiles = (dim + TILE_DIM - 1) / TILE_DIM;
    dim_t n_tile_pairs = n_tiles * (n_tiles + 1) / 2;
    return (size_t)(n_tile_pairs * TILE_SIZE) * sizeof(cplx_t);
}

/**
 * @brief Allocate `n` bytes, using mmap + huge-page hints for large allocations.
 *
 * On Linux, allocations >= BENCH_HUGEPAGE_THRESHOLD use MAP_ANONYMOUS|MAP_PRIVATE
 * with MADV_HUGEPAGE. On macOS the same threshold triggers MAP_ANONYMOUS|MAP_PRIVATE
 * (no explicit huge-page hint available). Smaller allocations fall back to malloc.
 */
void *bench_alloc_huge(size_t n);

/**
 * @brief Free memory previously allocated with bench_alloc_huge().
 *
 * Must be called with the same `n` passed to bench_alloc_huge(); large allocations
 * that used mmap are released with munmap, others with free.
 */
void bench_free_huge(void *p, size_t n);

/*
 * =====================================================================================================================
 * Shared pair-building utility (implemented in bench_mixed.c)
 * =====================================================================================================================
 */

/**
 * @brief Maximum number of qubit pairs supported by build_all_pairs() and build_all_ordered_pairs().
 *
 * Unordered pairs (build_all_pairs): n*(n-1)/2 — at 16 qubits that is 120 pairs.
 * Ordered pairs (build_all_ordered_pairs): n*(n-1) — at 16 qubits that is 240 pairs.
 * MAX_PAIRS=256 covers both cases up to and including 16 qubits.
 */
#define MAX_PAIRS 256
#define MAX_TARGETS                       (2 * MAX_PAIRS)  /**< 512; covers all 1Q targets and 2Q pairs up to 16 qubits */

/**
 * @brief Build all ordered qubit pairs (q1 < q2) for an n-qubit system.
 *
 * Writes up to MAX_PAIRS pairs into q1_out/q2_out. Returns the number of pairs written.
 * Fills pairs until max_pairs is reached and returns the count; callers must not pass qubits > 16.
 */
int build_all_pairs(qubit_t qubits, qubit_t *q1_out, qubit_t *q2_out);

/*
 * =====================================================================================================================
 * Benchmark functions — single-qubit gates (qlib: bench_mixed.c; baseline: bench_baselines.c)
 * =====================================================================================================================
 */

bench_result_t bench_qlib_packed(qubit_t qubits, const char *gate_name,
                                  void (*gate_fn)(state_t*, qubit_t),
                                  int iterations, int warmup, int runs);

bench_result_t bench_qlib_tiled(qubit_t qubits, const char *gate_name,
                                 void (*gate_fn)(state_t*, qubit_t),
                                 int iterations, int warmup, int runs);

bench_result_t bench_blas_dense(qubit_t qubits, const char *gate_name,
                                 const cplx_t gate_mat[4],
                                 int iterations, int warmup, int runs);

bench_result_t bench_naive_loop(qubit_t qubits, const char *gate_name,
                                 const cplx_t gate_mat[4],
                                 int iterations, int warmup, int runs);

/*
 * =====================================================================================================================
 * Benchmark functions — two-qubit gates (qlib: bench_mixed.c; baseline: bench_baselines.c)
 * =====================================================================================================================
 */

bench_result_t bench_qlib_packed_2q(qubit_t qubits, const char *gate_name,
                                      void (*gate_fn)(state_t*, qubit_t, qubit_t),
                                      int iterations, int warmup, int runs);

bench_result_t bench_qlib_tiled_2q(qubit_t qubits, const char *gate_name,
                                     void (*gate_fn)(state_t*, qubit_t, qubit_t),
                                     int iterations, int warmup, int runs);

bench_result_t bench_blas_dense_2q(qubit_t qubits, const char *gate_name,
                                     const cplx_t gate_mat[16],
                                     int iterations, int warmup, int runs);

/*
 * =====================================================================================================================
 * Output functions and CLI parsing (implemented in bench_main.c)
 * =====================================================================================================================
 */

/** @brief Print benchmark result to console */
void bench_print_result(const bench_result_t *result, int verbose);

/** @brief Print benchmark result as CSV row */
void bench_print_csv(const bench_result_t *result);

/** @brief Print CSV header */
void bench_print_csv_header(void);

/** @brief Parse command-line options */
bench_options_t bench_parse_options(int argc, char *argv[]);

/*
 * =====================================================================================================================
 * External framework benchmarks (conditionally compiled)
 * =====================================================================================================================
 */

#ifdef WITH_QUEST
/** @brief Initialise the QuEST environment. Must be called once before any bench_quest() calls. */
void bench_quest_init(void);
/** @brief Finalise the QuEST environment. Must be called once after all bench_quest() calls. */
void bench_quest_cleanup(void);

/**
 * @brief Run QuEST density matrix benchmark in-process.
 *
 * Requires bench_quest_init() to have been called. Creates and destroys a
 * Qureg per call; performs `runs` independent timed rounds and returns the
 * full statistical summary.
 */
bench_result_t bench_quest(qubit_t qubits, const char *gate_name,
                           int iterations, int warmup, int runs);

/**
 * @brief Per-qubit QuEST 1Q benchmark: probe-calibrate-measure all target qubits.
 *
 * Creates and destroys a Qureg internally. Writes qubits results to out[0..qubits-1].
 * Requires bench_quest_init() to have been called.
 */
void bench_quest_perq_1q(qubit_t qubits, const char *gate_name,
                          int iters_explicit, int explicit_iters, int runs,
                          bench_result_perq_t *out);

/**
 * @brief Per-qubit QuEST 2Q benchmark: probe-calibrate-measure all qubit pairs.
 *
 * Creates and destroys a Qureg internally. Writes n_pairs results to out[0..n_pairs-1].
 * Requires bench_quest_init() to have been called.
 */
void bench_quest_perq_2q(qubit_t qubits, const char *gate_name,
                          const qubit_t *q1s, const qubit_t *q2s, int n_pairs,
                          int iters_explicit, int explicit_iters, int runs,
                          bench_result_perq_t *out);
#endif

#ifdef WITH_QULACS
/**
 * @brief Run Qulacs density matrix benchmark in-process.
 *
 * Qulacs csim has no global state; this is a direct malloc/compute/free call.
 * Performs `runs` independent timed rounds and returns the full statistical summary.
 */
bench_result_t bench_qulacs(qubit_t qubits, const char *gate_name,
                            int iterations, int warmup, int runs);
#endif

#ifdef WITH_AER_DM
/**
 * @brief Run Qiskit-Aer density matrix benchmark in-process.
 *
 * AER::QV::DensityMatrix<double> has no global state; this is a direct
 * construct/compute/destruct call. Performs `runs` independent timed rounds
 * and returns the full statistical summary. Exceptions from C++ code are
 * caught internally; a zero time_ms result indicates failure.
 */
bench_result_t bench_aer_dm(qubit_t qubits, const char *gate_name,
                            int iterations, int warmup, int runs);
#endif

/*
 * =====================================================================================================================
 * Per-qubit benchmark functions — _at variants (bench_mixed.c / bench_baselines.c)
 * =====================================================================================================================
 */

/* qlib backends (bench_mixed.c) */
bench_result_perq_t bench_qlib_packed_at(qubit_t qubits, const char *gate_name,
    void (*gate_fn)(state_t*, qubit_t),
    qubit_t target, int iterations, int runs, state_t *state);

bench_result_perq_t bench_qlib_tiled_at(qubit_t qubits, const char *gate_name,
    void (*gate_fn)(state_t*, qubit_t),
    qubit_t target, int iterations, int runs, state_t *state);

bench_result_perq_t bench_qlib_packed_2q_at(qubit_t qubits, const char *gate_name,
    void (*gate_fn)(state_t*, qubit_t, qubit_t),
    qubit_t q1, qubit_t q2, int iterations, int runs, state_t *state);

bench_result_perq_t bench_qlib_tiled_2q_at(qubit_t qubits, const char *gate_name,
    void (*gate_fn)(state_t*, qubit_t, qubit_t),
    qubit_t q1, qubit_t q2, int iterations, int runs, state_t *state);

/* BLAS baseline (bench_baselines.c) */
bench_result_perq_t bench_blas_dense_at(qubit_t qubits, const char *gate_name,
    const cplx_t gate_mat[4],
    qubit_t target, int iterations, int runs, cplx_t *rho);

bench_result_perq_t bench_blas_dense_2q_at(qubit_t qubits, const char *gate_name,
    const cplx_t gate_mat[16],
    qubit_t q1, qubit_t q2, int iterations, int runs, cplx_t *rho);

/* bench_quest_at / bench_quest_2q_at are intentionally NOT declared here.
 * They take a Qureg* parameter and are only callable from bench_quest.c
 * (which includes quest.h first).  External callers use the higher-level
 * bench_quest_perq_1q / bench_quest_perq_2q wrappers declared in the
 * WITH_QUEST block above. */

#ifdef WITH_QULACS
bench_result_perq_t bench_qulacs_at(qubit_t qubits, const char *gate_name,
    qubit_t target, int iterations, int runs, void *rho);

bench_result_perq_t bench_qulacs_2q_at(qubit_t qubits, const char *gate_name,
    qubit_t q1, qubit_t q2, int iterations, int runs, void *rho);
#endif

#ifdef WITH_AER_DM
bench_result_perq_t bench_aer_dm_at(qubit_t qubits, const char *gate_name,
    qubit_t target, int iterations, int runs, void *dm);

bench_result_perq_t bench_aer_dm_2q_at(qubit_t qubits, const char *gate_name,
    qubit_t q1, qubit_t q2, int iterations, int runs, void *dm);

void *bench_aer_dm_alloc(qubit_t qubits);
void  bench_aer_dm_reinit(void *dm, qubit_t qubits);
void  bench_aer_dm_free(void *dm);
#endif

/*
 * =====================================================================================================================
 * Per-qubit probe–calibrate–measure helpers (implemented in bench_main.c)
 * =====================================================================================================================
 */

/** @brief Callback type for 1Q per-qubit benchmark functions. */
typedef bench_result_perq_t (*perq_bench_1q_fn)(
    qubit_t qubits, const char *gate_name,
    qubit_t target, int iterations, int runs, void *state);

/** @brief Callback type for 2Q per-qubit benchmark functions. */
typedef bench_result_perq_t (*perq_bench_2q_fn)(
    qubit_t qubits, const char *gate_name,
    qubit_t q1, qubit_t q2, int iterations, int runs, void *state);

/**
 * @brief Probe, calibrate, and measure all targets for one 1Q backend.
 *
 * Performs a quick single-shot probe, a multi-run calibration probe, then
 * measures all targets in [0, qubits) using the calibrated iteration count.
 * Results are written to out[0..qubits-1].
 */
void perq_run_1q(qubit_t qubits, const char *gate_name,
                 perq_bench_1q_fn fn, void *state,
                 int iters_explicit, int explicit_iters, int runs,
                 bench_result_perq_t *out);

/**
 * @brief Probe, calibrate, and measure all pairs for one 2Q backend.
 *
 * Performs a quick single-shot probe on pairs[0], a multi-run calibration probe,
 * then measures all n_pairs using the calibrated iteration count.
 * Results are written to out[0..n_pairs-1].
 */
void perq_run_2q(qubit_t qubits, const char *gate_name,
                 perq_bench_2q_fn fn, void *state,
                 const qubit_t *q1s, const qubit_t *q2s, int n_pairs,
                 int iters_explicit, int explicit_iters, int runs,
                 bench_result_perq_t *out);

/*
 * =====================================================================================================================
 * Per-qubit pair-building utility
 * =====================================================================================================================
 */

/**
 * @brief Build all ordered qubit pairs (q1 != q2, all permutations) for an n-qubit system.
 *
 * Writes up to max_pairs pairs into q1s/q2s. Returns the number of pairs written.
 * Callers must not pass qubits > MAX_QUBITS_FOR_ORDERED_PAIRS.
 */
int build_all_ordered_pairs(qubit_t qubits, qubit_t *q1s, qubit_t *q2s, int max_pairs);

/*
 * =====================================================================================================================
 * Per-qubit output functions (implemented in bench_main.c)
 * =====================================================================================================================
 */

/** @brief Print CSV header for per-qubit benchmark results */
void bench_print_perq_csv_header(FILE *f);

/** @brief Print per-qubit benchmark results as CSV rows */
void bench_print_perq_csv(FILE *f, const bench_result_perq_t *results, int n);

/** @brief Print per-qubit benchmark results to console */
void bench_print_perq_console(const bench_result_perq_t *results, int n);

#ifdef __cplusplus
}
#endif

#endif /* BENCH_H */
