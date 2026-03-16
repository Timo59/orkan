/**
 * @file bench_result.h
 * @brief Result API for the gate benchmark.
 *
 * Provides the bench_result_t type and the function to compute timing
 * statistics (mean, median, min, percentile-based 95% CI) and effective
 * bandwidth from an array of nanosecond-resolution samples.
 */

#ifndef BENCH_RESULT_H
#define BENCH_RESULT_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/** Maximum number of timed samples (outer loop of two-level timing). */
#define BENCH_MAX_SAMPLES 1024

/**
 * @brief Statistical summary of a benchmark measurement.
 *
 * All timing fields are in microseconds.  Bandwidth is in GB/s.
 */
typedef struct {
    double mean_us;    /**< Arithmetic mean of samples */
    double median_us;  /**< 50th percentile */
    double min_us;     /**< Minimum sample */
    double ci_lo_us;   /**< 2.5th percentile (lower bound of 95% CI) */
    double ci_hi_us;   /**< 97.5th percentile (upper bound of 95% CI) */
    double bw_gbs;     /**< Effective bandwidth: touch_bytes / median_ns */
} bench_result_t;

/**
 * @brief Compute timing result and bandwidth from raw nanosecond samples.
 *
 * Sorts the samples array in place (caller's array is modified), then
 * derives mean, median, min, percentile-based 95% CI, and effective
 * bandwidth.  Values are converted from nanoseconds to microseconds.
 *
 * @param samples_ns  Array of per-gate-application times in nanoseconds.
 *                    Modified in place: sorted and converted to microseconds.
 * @param n           Number of samples (must be >= 1, at most BENCH_MAX_SAMPLES).
 * @param touch_bytes Bytes read+written per gate application (for bandwidth).
 * @return            Populated bench_result_t.
 */
bench_result_t bench_compute_result(double *samples_ns, int n,
                                    size_t touch_bytes);

/** @brief Return a bench_result_t with all fields set to NAN. */
bench_result_t bench_result_nan(void);

#ifdef __cplusplus
}
#endif

#endif /* BENCH_RESULT_H */
