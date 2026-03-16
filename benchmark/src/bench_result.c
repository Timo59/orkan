/**
 * @file bench_result.c
 * @brief Timing result computation: mean, median, min, percentile CI, bandwidth.
 */

#include "bench_result.h"
#include <math.h>

bench_result_t bench_result_nan(void) {
    bench_result_t r;
    r.mean_us   = NAN;
    r.median_us = NAN;
    r.min_us    = NAN;
    r.ci_lo_us  = NAN;
    r.ci_hi_us  = NAN;
    r.bw_gbs    = NAN;
    return r;
}

bench_result_t bench_compute_result(double *samples_ns, int n,
                                    size_t touch_bytes) {
    if (n <= 0) return bench_result_nan();

    /* Insertion sort in nanoseconds */
    for (int i = 1; i < n; ++i) {
        double key = samples_ns[i];
        int j = i - 1;
        while (j >= 0 && samples_ns[j] > key) {
            samples_ns[j + 1] = samples_ns[j];
            --j;
        }
        samples_ns[j + 1] = key;
    }

    /* Bandwidth from median in nanoseconds (bytes/ns = GB/s) */
    double median_ns;
    if (n % 2 == 0)
        median_ns = (samples_ns[n / 2 - 1] + samples_ns[n / 2]) / 2.0;
    else
        median_ns = samples_ns[n / 2];
    double bw = (median_ns > 0.0) ? (double)touch_bytes / median_ns : 0.0;

    /* Convert sorted array to microseconds in place */
    for (int i = 0; i < n; ++i)
        samples_ns[i] /= 1000.0;

    /* Statistics in microseconds */
    double sum = 0.0;
    for (int i = 0; i < n; ++i)
        sum += samples_ns[i];

    /* Percentile-based 95% CI */
    int lo_idx = (int)(0.025 * n);
    int hi_idx = (int)ceil(0.975 * n) - 1;
    if (lo_idx < 0) lo_idx = 0;
    if (hi_idx >= n) hi_idx = n - 1;
    if (hi_idx < 0) hi_idx = 0;

    bench_result_t r;
    r.mean_us   = sum / n;
    r.median_us = median_ns / 1000.0;
    r.min_us    = samples_ns[0];
    r.ci_lo_us  = samples_ns[lo_idx];
    r.ci_hi_us  = samples_ns[hi_idx];
    r.bw_gbs    = bw;
    return r;
}
