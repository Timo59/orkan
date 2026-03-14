/**
 * @file bench_util.c
 * @brief Non-inline benchmark utility implementations.
 *
 * Contains implementations for functions declared in bench.h that are too
 * large to inline: bench_timing_barrier, bench_fill_random, bench_run_timed,
 * bench_fill_perq_stats, bench_alloc_huge, and bench_free_huge.
 */

#include "bench.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/mman.h>
#if !defined(__APPLE__) && !defined(__FreeBSD__)
#include <sys/random.h>
#endif

/* =====================================================================================================================
 * Timing serialization barrier
 * =====================================================================================================================
 */

void bench_timing_barrier(void) {
#if defined(__aarch64__)
    __asm__ volatile("isb" ::: "memory");
#elif defined(__x86_64__)
    __asm__ volatile("lfence" ::: "memory");
#else
    __asm__ volatile("" ::: "memory");  /* compiler barrier only */
#endif
}

/* =====================================================================================================================
 * Random fill helper
 * =====================================================================================================================
 */

void bench_fill_random(void *buf, size_t n) {
    /* Seed one uint64_t from OS entropy */
    uint64_t seed;
#if defined(__APPLE__) || defined(__FreeBSD__)
    arc4random_buf(&seed, sizeof(seed));
#else
    /* Linux: requires <sys/random.h> (glibc >= 2.25).
     * getrandom() can return -1 on error (EINTR, ENOSYS) or fewer bytes than
     * requested (though for <= 256 bytes it is guaranteed atomic). Fall back
     * to a time-based seed if the call fails. */
    if (getrandom(&seed, sizeof(seed), 0) != (ssize_t)sizeof(seed)) {
        struct timespec ts;
        clock_gettime(CLOCK_MONOTONIC, &ts);
        seed = (uint64_t)ts.tv_sec * 1000000000ULL + (uint64_t)ts.tv_nsec;
    }
#endif

    /* Expand seed via splitmix64 into 4 x uint64 xoshiro256** state */
    uint64_t z;
    z = (seed += 0x9e3779b97f4a7c15ULL); z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL; uint64_t s0 = z ^ (z >> 31);
    z = (seed += 0x9e3779b97f4a7c15ULL); z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL; uint64_t s1 = z ^ (z >> 31);
    z = (seed += 0x9e3779b97f4a7c15ULL); z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL; uint64_t s2 = z ^ (z >> 31);
    z = (seed += 0x9e3779b97f4a7c15ULL); z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL; uint64_t s3 = z ^ (z >> 31);

    /* xoshiro256** bulk fill — 8 bytes per step */
    uint8_t *p = (uint8_t *)buf;
    size_t remaining = n;
    while (remaining >= 8) {
        uint64_t r = s1 * 5; r = ((r << 7) | (r >> 57)) * 9;
        memcpy(p, &r, 8); p += 8; remaining -= 8;
        uint64_t t = s1 << 17;
        s2 ^= s0; s3 ^= s1; s1 ^= s2; s0 ^= s3; s2 ^= t;
        s3 = (s3 << 45) | (s3 >> 19);
    }
    if (remaining > 0) {
        uint64_t r = s1 * 5; r = ((r << 7) | (r >> 57)) * 9;
        memcpy(p, &r, remaining);
    }
}

/* =====================================================================================================================
 * Timed run harness
 * =====================================================================================================================
 */

int bench_run_timed(const bench_harness_t *h, double *out, int qubits) {
    int warmup_batch = (qubits <= 4) ? 128 : (qubits <= 7) ? 16 : (qubits <= 10) ? 4 : 1;

    uint64_t wt0 = bench_time_ns();
    do {
        for (int wb = 0; wb < warmup_batch; ++wb)
            h->call(h->ctx);
    } while (bench_ns_to_ms(bench_time_ns() - wt0) < BENCH_PQ_WARMUP_MS);

    for (int r = 0; r < h->runs; ++r) {
        bench_timing_barrier();
        uint64_t t0 = bench_time_ns();
        bench_timing_barrier();
        for (int i = 0; i < h->iterations; ++i) {
            h->call(h->ctx);
        }
        bench_timing_barrier();
        out[r] = bench_ns_to_ms(bench_time_ns() - t0);
    }
    return 0;
}

/* =====================================================================================================================
 * Per-qubit stats fill
 * =====================================================================================================================
 */

void bench_fill_perq_stats(bench_result_perq_t *r,
                            const bench_run_stats_t *s,
                            int iterations) {
    r->time_ms_mean   = s->mean;
    r->time_ms_std    = s->std_dev;
    r->time_ms_min    = s->min;
    r->time_ms_median = s->median;
    r->time_ms_cv     = s->cv;
    r->ops_per_sec    = (s->mean > 0.0) ? ((double)iterations / s->mean) * 1000.0 : 0.0;
    r->time_per_gate_us        = (iterations > 0) ? (s->mean   / iterations) * 1000.0 : 0.0;
    r->time_per_gate_us_median = (iterations > 0) ? (s->median / iterations) * 1000.0 : 0.0;
    r->time_per_gate_us_min    = (iterations > 0) ? (s->min    / iterations) * 1000.0 : 0.0;
    r->iterations = iterations;
    r->runs       = 0;  /* bench_run_stats_t carries no n field; caller must set r->runs */
}

/* =====================================================================================================================
 * Huge-page memory utilities
 * =====================================================================================================================
 */

void *bench_alloc_huge(size_t n) {
#if defined(__linux__)
    if (n >= BENCH_HUGEPAGE_THRESHOLD) {
        void *p = mmap(NULL, n, PROT_READ|PROT_WRITE, MAP_ANONYMOUS|MAP_PRIVATE, -1, 0);
        if (p != MAP_FAILED) {
            madvise(p, n, MADV_HUGEPAGE);
            return p;
        }
    }
    return malloc(n);
#elif defined(__APPLE__)
    if (n >= BENCH_HUGEPAGE_THRESHOLD) {
        void *p = mmap(NULL, n, PROT_READ|PROT_WRITE, MAP_ANONYMOUS|MAP_PRIVATE, -1, 0);
        if (p != MAP_FAILED) return p;
    }
    return malloc(n);
#else
    return malloc(n);
#endif
}

void bench_free_huge(void *p, size_t n) {
#if defined(__linux__) || defined(__APPLE__)
    if (n >= BENCH_HUGEPAGE_THRESHOLD && p != NULL) {
        munmap(p, n);
        return;
    }
#else
    (void)n;
#endif
    free(p);
}
