/**
 * @file bench_csv.c
 * @brief CSV output formatters for bench_gate
 *
 * Extracted from bench_main.c — pure output, no orchestration logic.
 */

#include "bench.h"
#include <stdio.h>

void bench_print_csv_header(void) {
    printf("qubits,dim,gate,method,runs,iterations,sweep_size,"
           "time_ms_mean,time_ms_std,time_ms_min,time_ms_median,time_ms_cv,"
           "time_per_gate_ms,"
           "ops_per_sec,memory_bytes\n");
}

void bench_print_csv(const bench_result_t *result) {
    printf("%u,%lld,%s,%s,%d,%d,%d,%.6f,%.6f,%.6f,%.6f,%.3f,%.9f,%.0f,%zu\n",
           result->qubits, (long long)result->dim, result->gate_name, result->method,
           result->runs, result->iterations, result->sweep_size,
           result->time_ms, result->time_ms_std, result->time_ms_min,
           result->time_ms_median, result->time_ms_cv,
           result->time_per_gate_ms,
           result->ops_per_sec, result->memory_bytes);
}

void bench_print_perq_csv_header(FILE *f) {
    fprintf(f, "qubits,dim,gate,method,runs,iterations,target,target2,is_2q,"
               "time_ms_mean,time_ms_std,time_ms_min,time_ms_median,time_ms_cv,"
               "time_per_gate_us,time_per_gate_us_median,time_per_gate_us_min,"
               "ops_per_sec,memory_bytes\n");
}

void bench_print_perq_csv(FILE *f, const bench_result_perq_t *results, int n) {
    for (int i = 0; i < n; i++) {
        const bench_result_perq_t *r = &results[i];
        if (r->time_ms_mean <= 0.0) continue;  /* skip failed measurements */
        fprintf(f, "%u,%lld,%s,%s,%d,%d,%u",
            (unsigned)r->qubits, (long long)r->dim,
            r->gate_name, r->method,
            r->runs, r->iterations, (unsigned)r->target);
        /* target2: emit -1 for 1Q gates */
        if (r->target2 == QUBIT_NONE) fprintf(f, ",-1");
        else                          fprintf(f, ",%u", (unsigned)r->target2);
        fprintf(f, ",%d", r->is_2q);
        fprintf(f, ",%.6f,%.6f,%.6f,%.6f,%.4f",
            r->time_ms_mean, r->time_ms_std, r->time_ms_min,
            r->time_ms_median, r->time_ms_cv);
        fprintf(f, ",%.6f,%.6f,%.6f",
            r->time_per_gate_us, r->time_per_gate_us_median, r->time_per_gate_us_min);
        fprintf(f, ",%.2f,%zu\n",
            r->ops_per_sec, r->memory_bytes);
    }
}
