/**
 * @file bench_console.c
 * @brief Console output formatter for the gate benchmark.
 */

#include "bench.h"
#include <math.h>

/* =====================================================================
 * Header line (shared)
 * ===================================================================== */

static void print_header(qubit_t qubits, const bench_gate_info_t *gate,
                         const bench_options_t *opts) {
    printf("Gate: %s   Qubits: %d   Samples: %d   Iterations: %d   Warm-up: %d\n\n",
           gate->name, (int)qubits, opts->samples, opts->iterations, opts->warm_up);
}

/* =====================================================================
 * Aggregate mode
 * ===================================================================== */

void bench_print_console_agg(qubit_t qubits, const bench_gate_info_t *gate,
                             const bench_options_t *opts,
                             const bench_result_t results[BACKEND_COUNT]) {
    print_header(qubits, gate, opts);

    printf("%-16s %13s %13s %13s %30s %13s\n",
           "Backend", "Mean (us)", "Median (us)", "Min (us)", "95% CI (us)", "BW (GB/s)");
    printf("-------------------------------------------------------------------------------------------------------------\n");

    for (int b = 0; b < BACKEND_COUNT; ++b) {
        const bench_result_t *r = &results[b];
        if (isnan(r->mean_us)) {
            printf("%-16s %13s %13s %13s %30s %13s\n",
                   bench_backends[b].name, "nan", "nan", "nan", "nan", "nan");
        } else {
            printf("%-16s %13.6f %13.6f %13.6f   [%12.6f, %12.6f] %13.6f\n",
                   bench_backends[b].name,
                   r->mean_us, r->median_us, r->min_us,
                   r->ci_lo_us, r->ci_hi_us, r->bw_gbs);
        }
    }
    printf("\n");
}

/* =====================================================================
 * Per-qubit mode
 * ===================================================================== */

static void print_pos_header(int gate_qubits) {
    if (gate_qubits == 1)
        printf("%-16s %6s", "Backend", "Target");
    else if (gate_qubits == 2)
        printf("%-16s %7s %6s", "Backend", "Control", "Target");
    else
        printf("%-16s %4s %4s %6s", "Backend", "C1", "C2", "Target");

    printf(" %13s %13s %13s %30s %13s\n",
           "Mean (us)", "Median (us)", "Min (us)", "95% CI (us)", "BW (GB/s)");

    int pos_width = (gate_qubits == 1) ? 7 : (gate_qubits == 2) ? 15 : 16;
    int total = 16 + pos_width + 13 + 13 + 13 + 30 + 13 + 6;
    for (int i = 0; i < total; ++i) putchar('-');
    putchar('\n');
}

static void print_pos_cols(const qubit_t *pos, int gate_qubits) {
    if (gate_qubits == 1)
        printf(" %6d", (int)pos[0]);
    else if (gate_qubits == 2)
        printf(" %7d %6d", (int)pos[0], (int)pos[1]);
    else
        printf(" %4d %4d %6d", (int)pos[0], (int)pos[1], (int)pos[2]);
}

void bench_print_console_perq(qubit_t qubits, const bench_gate_info_t *gate,
                              const bench_options_t *opts,
                              const bench_pos_result_t *per_backend[BACKEND_COUNT],
                              int n_pos, int gate_qubits) {
    print_header(qubits, gate, opts);
    print_pos_header(gate_qubits);

    for (int b = 0; b < BACKEND_COUNT; ++b) {
        if (!per_backend[b]) {
            printf("%-16s", bench_backends[b].name);
            print_pos_cols((qubit_t[]){0, 0, 0}, gate_qubits);
            printf(" %13s %13s %13s %30s %13s\n",
                   "nan", "nan", "nan", "nan", "nan");
            continue;
        }
        for (int p = 0; p < n_pos; ++p) {
            const bench_pos_result_t *pr = &per_backend[b][p];
            const bench_result_t *r = &pr->result;
            printf("%-16s", bench_backends[b].name);
            print_pos_cols(pr->pos, gate_qubits);
            if (isnan(r->mean_us)) {
                printf(" %13s %13s %13s %30s %13s\n",
                       "nan", "nan", "nan", "nan", "nan");
            } else {
                printf(" %13.6f %13.6f %13.6f   [%12.6f, %12.6f] %13.6f\n",
                       r->mean_us, r->median_us, r->min_us,
                       r->ci_lo_us, r->ci_hi_us, r->bw_gbs);
            }
        }
    }
    printf("\n");
}
