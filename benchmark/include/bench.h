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

/** @brief Command-line options */
typedef struct bench_options {
    qubit_t min_qubits;
    qubit_t max_qubits;
    int step;
    int iterations;
    int warmup;
    int runs;
    int csv_output;
    int pgfplots_output;
    int verbose;
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
static inline bench_run_stats_t bench_compute_stats(const double *times, int n) {
    /* Mean and minimum */
    double sum = 0.0;
    double mn  = times[0];
    for (int i = 0; i < n; ++i) {
        sum += times[i];
        if (times[i] < mn) mn = times[i];
    }
    double mean = sum / n;

    /* Sample standard deviation (Bessel's correction: divide by n-1) */
    double var = 0.0;
    for (int i = 0; i < n; ++i) {
        double d = times[i] - mean;
        var += d * d;
    }
    double std_dev = (n > 1) ? sqrt(var / (n - 1)) : 0.0;
    double cv = (mean > 0.0) ? (std_dev / mean) * 100.0 : 0.0;

    /*
     * Median: insertion-sort a stack copy.
     * n is bounded by the --runs validation (max 200), so stack allocation is safe.
     */
    double sorted[200];
    int cnt = (n < 200) ? n : 200;
    for (int i = 0; i < cnt; ++i) sorted[i] = times[i];
    for (int i = 1; i < cnt; ++i) {
        double key = sorted[i];
        int j = i - 1;
        while (j >= 0 && sorted[j] > key) { sorted[j + 1] = sorted[j]; --j; }
        sorted[j + 1] = key;
    }
    double median = (cnt % 2 == 0)
        ? (sorted[cnt / 2 - 1] + sorted[cnt / 2]) / 2.0
        : sorted[cnt / 2];

    bench_run_stats_t s;
    s.mean    = mean;
    s.std_dev = std_dev;
    s.min     = mn;
    s.median  = median;
    s.cv      = cv;
    return s;
}

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

/*
 * =====================================================================================================================
 * Shared pair-building utility (implemented in bench_mixed.c)
 * =====================================================================================================================
 */

/** @brief Maximum number of qubit pairs supported by build_all_pairs(). Sufficient for up to 16 qubits. */
#define MAX_PAIRS 128

/**
 * @brief Build all ordered qubit pairs (q1 < q2) for an n-qubit system.
 *
 * Writes up to MAX_PAIRS pairs into q1_out/q2_out. Returns the number of pairs written.
 * Asserts that qubits*(qubits-1)/2 <= MAX_PAIRS; callers must not pass qubits > 16.
 */
int build_all_pairs(qubit_t qubits, qubit_t *q1_out, qubit_t *q2_out);

/*
 * =====================================================================================================================
 * Benchmark functions — single-qubit gates (qlib: bench_mixed.c; baselines: bench_baselines.c)
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
 * Benchmark functions — two-qubit gates (qlib: bench_mixed.c; baselines: bench_baselines.c)
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

#ifdef __cplusplus
}
#endif

#endif /* BENCH_H */
