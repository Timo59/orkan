/**
 * @file bench_csv.c
 * @brief CSV output formatter for the gate benchmark.
 */

#include "bench.h"
#include <math.h>

/* =====================================================================
 * Portable NaN-safe printing
 * ===================================================================== */

static void print_val(double v) {
    if (isnan(v)) printf("nan");
    else          printf("%.6f", v);
}

static void print_result_csv(const bench_result_t *r) {
    print_val(r->mean_us);   printf(",");
    print_val(r->median_us); printf(",");
    print_val(r->min_us);    printf(",");
    print_val(r->ci_lo_us);  printf(",");
    print_val(r->ci_hi_us);  printf(",");
    print_val(r->bw_gbs);    printf("\n");
}

/* =====================================================================
 * Header
 * ===================================================================== */

void bench_print_csv_header(const bench_gate_info_t *gate, int per_qubit) {
    if (!per_qubit) {
        printf("qubits,backend,mean_us,median_us,min_us,ci_lo_us,ci_hi_us,bw_gbs\n");
        return;
    }
    if (gate->n_qubits == 1)
        printf("qubits,backend,target,");
    else if (gate->n_qubits == 2)
        printf("qubits,backend,control,target,");
    else
        printf("qubits,backend,c1,c2,target,");
    printf("mean_us,median_us,min_us,ci_lo_us,ci_hi_us,bw_gbs\n");
}

/* =====================================================================
 * Aggregate mode
 * ===================================================================== */

void bench_print_csv_agg(qubit_t qubits,
                         const bench_result_t results[BACKEND_COUNT]) {
    for (int b = 0; b < BACKEND_COUNT; ++b) {
        printf("%d,%s,", (int)qubits, bench_backends[b].name);
        print_result_csv(&results[b]);
    }
}

/* =====================================================================
 * Per-qubit mode
 * ===================================================================== */

void bench_print_csv_perq(qubit_t qubits, bench_backend_id_t backend,
                          const bench_pos_result_t *results, int n_pos,
                          int gate_qubits) {
    for (int p = 0; p < n_pos; ++p) {
        const bench_pos_result_t *pr = &results[p];
        printf("%d,%s,", (int)qubits, bench_backends[backend].name);
        if (gate_qubits == 1)
            printf("%d,", (int)pr->pos[0]);
        else if (gate_qubits == 2)
            printf("%d,%d,", (int)pr->pos[0], (int)pr->pos[1]);
        else
            printf("%d,%d,%d,", (int)pr->pos[0], (int)pr->pos[1], (int)pr->pos[2]);
        print_result_csv(&pr->result);
    }
}
