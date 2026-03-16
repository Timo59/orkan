/**
 * @file bench_console.c
 * @brief Console output formatters for bench_gate
 *
 * Extracted from bench_main.c — pure output, no orchestration logic.
 */

#include "bench.h"
#include <stdio.h>

void bench_print_result(const bench_result_t *result, int verbose) {
    printf("  %-14s %8.3f ± %6.3f ms (CV %4.1f%%)  %12.0f ops/sec",
           result->method,
           result->time_ms, result->time_ms_std, result->time_ms_cv,
           result->ops_per_sec);
    if (verbose)
        printf("  [min %8.3f ms  med %8.3f ms  mem %zu B]",
               result->time_ms_min, result->time_ms_median, result->memory_bytes);
    printf("\n");
}

void bench_print_perq_console(const bench_result_perq_t *results, int n) {
    if (n <= 0 || results == NULL) return;
    const bench_result_perq_t *r0 = &results[0];
    int is_2q = r0->is_2q;

    printf("\n--- %s | %s | %u qubits ---\n",
        r0->gate_name ? r0->gate_name : "?",
        r0->method    ? r0->method    : "?",
        (unsigned)r0->qubits);

    if (!is_2q) {
        /* 1Q: one row per target */
        printf("  %-6s  %-12s  %-12s  %-12s  %-8s  %-14s\n",
            "target", "mean(us/g)", "median(us/g)", "min(us/g)", "CV(%)", "ops/sec");
        for (int i = 0; i < n; i++) {
            const bench_result_perq_t *r = &results[i];
            if (r->time_ms_mean <= 0.0) continue;
            const char *noisy = (r->time_ms_cv > 50.0) ? " (noisy)" : "";
            printf("  %-6u  %-12.3f  %-12.3f  %-12.3f  %-8.1f  %-14.0f%s\n",
                (unsigned)r->target,
                r->time_per_gate_us, r->time_per_gate_us_median, r->time_per_gate_us_min,
                r->time_ms_cv, r->ops_per_sec, noisy);
        }
    } else {
        /* 2Q: summary line showing min/mean/max with pair coordinates */
        double min_us = 1e18, max_us = -1.0, sum_us = 0.0;
        int count = 0;
        int min_i = -1, max_i = -1;
        for (int i = 0; i < n; i++) {
            const bench_result_perq_t *r = &results[i];
            if (r->time_ms_mean <= 0.0) continue;
            double v = r->time_per_gate_us;
            if (v < min_us) { min_us = v; min_i = i; }
            if (v > max_us) { max_us = v; max_i = i; }
            sum_us += v;
            count++;
        }
        if (count > 0) {
            double mean_us = sum_us / count;
            printf("  min=%.3f us (q%u,q%u)  mean=%.3f us  max=%.3f us (q%u,q%u)\n",
                min_us,
                min_i >= 0 ? (unsigned)results[min_i].target : 0,
                min_i >= 0 ? (unsigned)results[min_i].target2 : 0,
                mean_us,
                max_us,
                max_i >= 0 ? (unsigned)results[max_i].target : 0,
                max_i >= 0 ? (unsigned)results[max_i].target2 : 0);
            printf("  (full per-pair data in CSV)\n");
        }
    }
}

void print_speedup_and_memory(const bench_summary_t *s) {
    if (s->packed.time_ms <= 0.0) return;  /* no valid baseline to compare against */

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
    if (s->has_aer_dm && s->aer_dm.time_ms > 0) {
        printf("%s%.1fx vs Aer-DM", has ? ", " : " ", s->aer_dm.time_ms / s->packed.time_ms);
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
    if (s->has_aer_dm && s->aer_dm.memory_bytes > max_mem) max_mem = s->aer_dm.memory_bytes;

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
    if (s->has_aer_dm && s->aer_dm.memory_bytes > 0)
        printf(", Aer-DM=%.1f %s", (double)s->aer_dm.memory_bytes / scale, unit);
    if (s->has_dense && s->packed.memory_bytes > 0 && s->dense.memory_bytes > 0) {
        double savings = 100.0 * (1.0 - (double)s->packed.memory_bytes / (double)s->dense.memory_bytes);
        printf(" (packed saves %.0f%%)", savings);
    }
    printf("\n");
}
