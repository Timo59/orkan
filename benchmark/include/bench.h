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
#define BENCH_DEFAULT_ITERATIONS    1000
#define BENCH_DEFAULT_WARMUP        100

/*
 * =====================================================================================================================
 * Types
 * =====================================================================================================================
 */

/** @brief Benchmark result for a single configuration */
typedef struct bench_result {
    qubit_t qubits;             /**< Number of qubits */
    dim_t dim;                  /**< Hilbert space dimension (2^qubits) */
    const char *gate_name;      /**< Name of the gate tested */
    const char *method;         /**< Implementation method name */
    double time_ms;             /**< Total time in milliseconds */
    double ops_per_sec;         /**< Operations per second */
    size_t memory_bytes;        /**< Memory used for state representation */
    int iterations;             /**< Number of iterations */
} bench_result_t;

/** @brief Command-line options */
typedef struct bench_options {
    qubit_t min_qubits;
    qubit_t max_qubits;
    int iterations;
    int warmup;
    int csv_output;
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

/*
 * =====================================================================================================================
 * Benchmark functions (implemented in bench_mixed.c)
 * =====================================================================================================================
 */

/** @brief Run qlib packed implementation benchmark */
bench_result_t bench_qlib_packed(qubit_t qubits, const char *gate_name,
                                  qs_error_t (*gate_fn)(state_t*, qubit_t),
                                  int iterations, int warmup);

/** @brief Run dense BLAS reference benchmark */
bench_result_t bench_blas_dense(qubit_t qubits, const char *gate_name,
                                 const cplx_t gate_mat[4],
                                 int iterations, int warmup);

/** @brief Run naive loop implementation benchmark */
bench_result_t bench_naive_loop(qubit_t qubits, const char *gate_name,
                                 const cplx_t gate_mat[4],
                                 int iterations, int warmup);

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
/**
 * @brief Initialize QuEST environment (call once at startup)
 */
void bench_quest_init(void);

/**
 * @brief Cleanup QuEST environment (call once at shutdown)
 */
void bench_quest_cleanup(void);

/**
 * @brief Run QuEST density matrix benchmark
 * @see https://quest.qtechtheory.org/
 */
bench_result_t bench_quest(qubit_t qubits, const char *gate_name,
                           int iterations, int warmup);
#endif

#ifdef WITH_QPP
/**
 * @brief Run Quantum++ density matrix benchmark
 * @see https://github.com/softwareQinc/qpp
 */
bench_result_t bench_qpp(qubit_t qubits, const char *gate_name,
                         int iterations, int warmup);
#endif

#ifdef __cplusplus
}
#endif

#endif /* BENCH_H */
