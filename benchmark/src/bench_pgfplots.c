/**
 * @file bench_pgfplots.c
 * @brief PGFplots output formatter for the gate benchmark.
 */

#include "bench.h"
#include <math.h>

/* =====================================================================
 * Helpers
 * ===================================================================== */

static const char *quantity_names[] = {
    "mean_us", "median_us", "min_us", "ci_lo_us", "ci_hi_us", "bw_gbs"
};

static double get_field(const bench_result_t *r, int field) {
    switch (field) {
    case 0: return r->mean_us;
    case 1: return r->median_us;
    case 2: return r->min_us;
    case 3: return r->ci_lo_us;
    case 4: return r->ci_hi_us;
    case 5: return r->bw_gbs;
    default: return NAN;
    }
}

/** @brief Print a value in a 16-char wide field, NaN-safe. */
static void print_val_w(double v) {
    if (isnan(v)) printf(" %16s", "nan");
    else          printf(" %16.6f", v);
}

static void print_comment(const bench_gate_info_t *gate,
                           const bench_options_t *opts) {
    printf("%% Gate: %s   Samples: %d   Iterations: %d   Warm-up: %d\n",
           gate->name, opts->samples, opts->iterations, opts->warm_up);
}

static void print_backend_columns(void) {
    for (int b = 0; b < BACKEND_COUNT; ++b)
        printf(" %16s", bench_backends[b].name);
    printf("\n");
}

/* =====================================================================
 * Aggregate mode: one table per quantity
 * ===================================================================== */

void bench_print_pgfplots_agg(const bench_gate_info_t *gate,
                              const bench_options_t *opts,
                              const qubit_t *qubit_list, int n_qubits,
                              const bench_result_t (*results)[BACKEND_COUNT]) {
    for (int f = 0; f < 6; ++f) {
        print_comment(gate, opts);
        printf("%% Quantity: %s\n", quantity_names[f]);
        printf("%8s", "qubits");
        print_backend_columns();

        for (int qi = 0; qi < n_qubits; ++qi) {
            printf("%8d", (int)qubit_list[qi]);
            for (int b = 0; b < BACKEND_COUNT; ++b)
                print_val_w(get_field(&results[qi][b], f));
            printf("\n");
        }
        printf("\n");
    }
}

/* =====================================================================
 * Per-qubit mode: one table per (qubit count, quantity)
 * ===================================================================== */

void bench_print_pgfplots_perq(const bench_gate_info_t *gate,
                               const bench_options_t *opts,
                               qubit_t qubits,
                               const bench_pos_result_t *per_backend[BACKEND_COUNT],
                               int n_pos, int gate_qubits) {
    for (int f = 0; f < 6; ++f) {
        printf("%% Gate: %s   Qubits: %d   Samples: %d   Iterations: %d   Warm-up: %d\n",
               gate->name, (int)qubits, opts->samples, opts->iterations, opts->warm_up);
        printf("%% Quantity: %s\n", quantity_names[f]);

        if (gate_qubits == 1)
            printf("%8s", "target");
        else if (gate_qubits == 2)
            printf("%8s %8s", "control", "target");
        else
            printf("%8s %8s %8s", "c1", "c2", "target");
        print_backend_columns();

        for (int p = 0; p < n_pos; ++p) {
            const qubit_t *pos = NULL;
            for (int b = 0; b < BACKEND_COUNT; ++b) {
                if (per_backend[b]) { pos = per_backend[b][p].pos; break; }
            }
            if (!pos) continue;

            if (gate_qubits == 1)
                printf("%8d", (int)pos[0]);
            else if (gate_qubits == 2)
                printf("%8d %8d", (int)pos[0], (int)pos[1]);
            else
                printf("%8d %8d %8d", (int)pos[0], (int)pos[1], (int)pos[2]);

            for (int b = 0; b < BACKEND_COUNT; ++b) {
                if (per_backend[b])
                    print_val_w(get_field(&per_backend[b][p].result, f));
                else
                    printf(" %16s", "nan");
            }
            printf("\n");
        }
        printf("\n");
    }
}
