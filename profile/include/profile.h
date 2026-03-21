#ifndef PROFILE_H
#define PROFILE_H

/**
 * @file profile.h
 * @brief Micro-profiling framework for qlib gate functions
 *
 * Provides high-resolution timing, descriptive statistics, memory bandwidth
 * estimation, and CSV output for detailed performance analysis of individual
 * gate implementations.
 *
 * Usage pattern:
 *   1. Allocate and initialize state
 *   2. Warm up (run gate N times, discard timings)
 *   3. Measure: for each iteration, time a single gate call
 *   4. Compute statistics over the sample array
 *   5. Print human-readable summary and/or CSV row
 */

#include "q_types.h"
#include "state.h"
#include "gate.h"
#include <time.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
 * =====================================================================================================================
 * Configuration
 * =====================================================================================================================
 */

#define PROF_DEFAULT_WARMUP     200
#define PROF_DEFAULT_ITERS      2000
#define PROF_DEFAULT_MIN_QUBITS 2
#define PROF_DEFAULT_MAX_QUBITS 12

/*
 * =====================================================================================================================
 * High-resolution timing
 * =====================================================================================================================
 */

/** @brief Read monotonic clock in nanoseconds */
static inline uint64_t prof_time_ns(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (uint64_t)ts.tv_sec * 1000000000ULL + (uint64_t)ts.tv_nsec;
}

/*
 * =====================================================================================================================
 * Descriptive statistics over a sample array of doubles (nanoseconds)
 * =====================================================================================================================
 */

typedef struct prof_stats {
    double mean_ns;
    double median_ns;
    double min_ns;
    double max_ns;
    double stddev_ns;
    double p5_ns;       /* 5th percentile */
    double p95_ns;      /* 95th percentile */
    int    n;           /* sample count */
} prof_stats_t;

/** @brief Comparison function for qsort on doubles */
static int prof_cmp_double(const void *a, const void *b) {
    double da = *(const double *)a, db = *(const double *)b;
    return (da > db) - (da < db);
}

/**
 * @brief Compute descriptive statistics from an array of timing samples (ns).
 *
 * Sorts the input array in place for percentile computation.
 */
static inline prof_stats_t prof_compute_stats(double *samples, int n) {
    prof_stats_t s = {0};
    s.n = n;
    if (n == 0) return s;

    qsort(samples, (size_t)n, sizeof(double), prof_cmp_double);

    s.min_ns = samples[0];
    s.max_ns = samples[n - 1];
    s.median_ns = (n % 2 == 1) ? samples[n / 2]
                                : 0.5 * (samples[n / 2 - 1] + samples[n / 2]);
    s.p5_ns  = samples[(int)(0.05 * (n - 1))];
    s.p95_ns = samples[(int)(0.95 * (n - 1))];

    double sum = 0;
    for (int i = 0; i < n; i++) sum += samples[i];
    s.mean_ns = sum / n;

    double var = 0;
    for (int i = 0; i < n; i++) {
        double d = samples[i] - s.mean_ns;
        var += d * d;
    }
    s.stddev_ns = (n > 1) ? sqrt(var / (n - 1)) : 0;

    return s;
}

/*
 * =====================================================================================================================
 * Memory bandwidth estimation
 * =====================================================================================================================
 */

/**
 * @brief Estimate bytes touched by a single-qubit tiled gate application.
 *
 * For X gate (swap-only, no arithmetic): each element is read once and written
 * once. The density matrix has n_tiles*(n_tiles+1)/2 tiles, each TILE_SIZE
 * elements of sizeof(cplx_t) bytes.
 *
 * Returns total bytes read + written (i.e., 2 * data_size for a full sweep).
 */
static inline size_t prof_tiled_bytes(qubit_t qubits) {
    uint64_t dim = (uint64_t)1 << qubits;
    uint64_t n_tiles = (dim + TILE_DIM - 1) / TILE_DIM;
    uint64_t n_tile_pairs = n_tiles * (n_tiles + 1) / 2;
    size_t data_bytes = (size_t)(n_tile_pairs * TILE_SIZE) * sizeof(cplx_t);
    return 2 * data_bytes;  /* read + write */
}

/** @brief Same for packed storage */
static inline size_t prof_packed_bytes(qubit_t qubits) {
    uint64_t dim = (uint64_t)1 << qubits;
    size_t data_bytes = (size_t)(dim * (dim + 1) / 2) * sizeof(cplx_t);
    return 2 * data_bytes;
}

/*
 * =====================================================================================================================
 * Result structure for a single profiling run
 * =====================================================================================================================
 */

typedef struct prof_result {
    const char   *func_name;    /* e.g. "x_tiled" */
    qubit_t       qubits;
    qubit_t       target;
    const char   *path;         /* "within-tile" or "cross-tile" */
    prof_stats_t  stats;
    size_t        data_bytes;   /* state data size */
    size_t        touch_bytes;  /* estimated bytes read+written per call */
    double        bw_gbs;       /* effective bandwidth in GB/s (median) */
    double        elem_per_ns;  /* elements processed per nanosecond (median) */
} prof_result_t;

/*
 * =====================================================================================================================
 * Output helpers
 * =====================================================================================================================
 */

static inline void prof_print_header(void) {
    printf("%-14s %3s %3s %-12s %10s %10s %10s %10s %10s %10s %8s %10s\n",
           "function", "n", "tgt", "path",
           "median(ns)", "mean(ns)", "min(ns)", "p5(ns)", "p95(ns)", "stddev",
           "BW(GB/s)", "elem/ns");
    printf("%-14s %3s %3s %-12s %10s %10s %10s %10s %10s %10s %8s %10s\n",
           "--------------", "---", "---", "------------",
           "----------", "----------", "----------", "----------",
           "----------", "----------", "--------", "----------");
}

static inline void prof_print_result(const prof_result_t *r) {
    printf("%-14s %3d %3d %-12s %10.0f %10.0f %10.0f %10.0f %10.0f %10.0f %8.2f %10.2f\n",
           r->func_name, r->qubits, r->target, r->path,
           r->stats.median_ns, r->stats.mean_ns, r->stats.min_ns,
           r->stats.p5_ns, r->stats.p95_ns, r->stats.stddev_ns,
           r->bw_gbs, r->elem_per_ns);
}

static inline void prof_print_csv_header(FILE *f) {
    fprintf(f, "function,qubits,target,path,median_ns,mean_ns,min_ns,p5_ns,p95_ns,"
               "stddev_ns,max_ns,data_bytes,touch_bytes,bw_gbs,elem_per_ns,samples\n");
}

static inline void prof_print_csv_row(FILE *f, const prof_result_t *r) {
    fprintf(f, "%s,%d,%d,%s,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%zu,%zu,%.3f,%.3f,%d\n",
            r->func_name, r->qubits, r->target, r->path,
            r->stats.median_ns, r->stats.mean_ns, r->stats.min_ns,
            r->stats.p5_ns, r->stats.p95_ns, r->stats.stddev_ns,
            r->stats.max_ns,
            r->data_bytes, r->touch_bytes,
            r->bw_gbs, r->elem_per_ns, r->stats.n);
}

/*
 * =====================================================================================================================
 * Profiling runner
 * =====================================================================================================================
 *
 * Generic profiling loop for a single-qubit gate function.
 * Allocates state, warms up, measures, computes stats.
 */

typedef void (*gate_1q_fn)(state_t *, qubit_t);

static inline prof_result_t prof_run_1q(
    gate_1q_fn      fn,
    const char     *func_name,
    state_type_t    stype,
    qubit_t         qubits,
    qubit_t         target,
    int             warmup,
    int             iters)
{
    prof_result_t r = {0};
    r.func_name = func_name;
    r.qubits = qubits;
    r.target = target;

    if (stype == MIXED_TILED) {
        r.path = (target < LOG_TILE_DIM) ? "within-tile" : "cross-tile";
        r.touch_bytes = prof_tiled_bytes(qubits);
    } else {
        r.path = "packed";
        r.touch_bytes = prof_packed_bytes(qubits);
    }

    /* Allocate and initialize state to |+>^n (nonzero everywhere) */
    state_t state = {0};
    state.type = stype;
    state_plus(&state, qubits);
    r.data_bytes = (size_t)state_len(&state) * sizeof(cplx_t);

    /* Warm up */
    for (int i = 0; i < warmup; i++) {
        fn(&state, target);
    }

    /* Measure */
    double *samples = (double *)malloc((size_t)iters * sizeof(double));
    for (int i = 0; i < iters; i++) {
        uint64_t t0 = prof_time_ns();
        fn(&state, target);
        uint64_t t1 = prof_time_ns();
        samples[i] = (double)(t1 - t0);
    }

    r.stats = prof_compute_stats(samples, iters);
    free(samples);
    state_free(&state);

    /* Derived metrics from median */
    if (r.stats.median_ns > 0) {
        r.bw_gbs = (double)r.touch_bytes / r.stats.median_ns;  /* bytes/ns = GB/s */

        uint64_t dim = (uint64_t)1 << qubits;
        uint64_t total_elems;
        if (stype == MIXED_TILED) {
            uint64_t n_tiles = (dim + TILE_DIM - 1) / TILE_DIM;
            total_elems = n_tiles * (n_tiles + 1) / 2 * TILE_SIZE;
        } else {
            total_elems = dim * (dim + 1) / 2;
        }
        r.elem_per_ns = (double)total_elems / r.stats.median_ns;
    }

    return r;
}

#ifdef __cplusplus
}
#endif

#endif /* PROFILE_H */
