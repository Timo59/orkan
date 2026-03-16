/**
 * @file bench_pgfplots.c
 * @brief pgfplots .dat output formatters for bench_gate
 *
 * Extracted from bench_main.c — pure output, no orchestration logic.
 */

#include "bench.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*
 * t-distribution critical values for 95% two-tailed confidence interval.
 * Indexed by degrees of freedom df = runs - 1.
 * t_crit_table[0] is df=1, t_crit_table[48] is df=49.
 * For df >= 50 use 1.96 (normal approximation).
 */
static const double t_crit_table[49] = {
    12.706, /* df=1  */
     4.303, /* df=2  */
     3.182, /* df=3  */
     2.776, /* df=4  */
     2.571, /* df=5  */
     2.447, /* df=6  */
     2.365, /* df=7  */
     2.306, /* df=8  */
     2.262, /* df=9  */
     2.228, /* df=10 */
     2.201, /* df=11 */
     2.179, /* df=12 */
     2.160, /* df=13 */
     2.145, /* df=14 */
     2.131, /* df=15 */
     2.120, /* df=16 */
     2.110, /* df=17 */
     2.101, /* df=18 */
     2.093, /* df=19 */
     2.086, /* df=20 */
     2.080, /* df=21 */
     2.074, /* df=22 */
     2.069, /* df=23 */
     2.064, /* df=24 */
     2.060, /* df=25 */
     2.056, /* df=26 */
     2.052, /* df=27 */
     2.048, /* df=28 */
     2.045, /* df=29 */
     2.042, /* df=30 */
     2.040, /* df=31 */
     2.037, /* df=32 */
     2.035, /* df=33 */
     2.032, /* df=34 */
     2.030, /* df=35 */
     2.028, /* df=36 */
     2.026, /* df=37 */
     2.024, /* df=38 */
     2.023, /* df=39 */
     2.021, /* df=40 */
     2.020, /* df=41 */
     2.018, /* df=42 */
     2.017, /* df=43 */
     2.015, /* df=44 */
     2.014, /* df=45 */
     2.013, /* df=46 */
     2.012, /* df=47 */
     2.011, /* df=48 */
     2.010, /* df=49 */
};

/** @brief Look up t critical value for given degrees of freedom (df = runs - 1). */
double t_crit(int df) {
    if (df < 1)  return t_crit_table[0];  /* clamp to df=1 */
    if (df < 50) return t_crit_table[df - 1];
    return 1.96;
}

/**
 * @brief Compute sweep_size for a gate at a given qubit count.
 *
 * 1Q gates: sweep over all qubits (sweep_size = qubits).
 * 2Q gates: sweep over all ordered pairs (sweep_size = qubits*(qubits-1)/2).
 * Gate index g: 0..NUM_1Q_GATES-1 are 1Q; NUM_1Q_GATES..NUM_GATES-1 are 2Q.
 */
int gate_sweep_size(int g, qubit_t qubits) {
    if (g < NUM_1Q_GATES)
        return (int)qubits;
    else
        return (int)(qubits * (qubits - 1) / 2);
}

void print_pgfplots_output(const pgfplots_data_t *pgf) {
    const char *method_names[]   = {"qlib_packed", "qlib_tiled", "blas", "quest", "qulacs", "qiskit_aer_dm"};
    int         method_indices[] = {M_QLIB, M_QLIB_TILED, M_BLAS, M_QUEST, M_QULACS, M_AER_DM};
    int         num_methods      = 6;

    /* ── Per-gate timing tables (5 tables each: mean, ci95, min, median, cv) ─── */
    for (int g = 0; g < NUM_GATES; ++g) {

        /* time_per_gate_us — mean normalised per gate application */
        printf("# %s time_per_gate_us (mean µs per gate application)\n", all_gate_names[g]);
        printf("%-8s", "qubits");
        for (int m = 0; m < num_methods; ++m) printf("  %-14s", method_names[m]);
        printf("\n");
        for (int q = 0; q < pgf->num_configs; ++q) {
            qubit_t qubits    = pgf->qubits[q];
            double  norm      = (double)pgf->iterations * (double)gate_sweep_size(g, qubits);
            printf("%-8u", qubits);
            for (int m = 0; m < num_methods; ++m) {
                double val = pgf->time_ms[g][q][method_indices[m]];
                if (val > 0 && norm > 0) printf("  %-14.6f", val * 1000.0 / norm);
                else                     printf("  %-14s", "nan");
            }
            printf("\n");
        }
        printf("\n");

        /* ci95_per_gate — 95% confidence interval half-width, normalised per gate */
        printf("# %s ci95_per_gate (µs per gate application)\n", all_gate_names[g]);
        printf("%-8s", "qubits");
        for (int m = 0; m < num_methods; ++m) printf("  %-14s", method_names[m]);
        printf("\n");
        for (int q = 0; q < pgf->num_configs; ++q) {
            qubit_t qubits    = pgf->qubits[q];
            double  norm      = (double)pgf->iterations * (double)gate_sweep_size(g, qubits);
            int     runs_q    = pgf->runs[q];
            double  tc        = t_crit(runs_q - 1);
            printf("%-8u", qubits);
            for (int m = 0; m < num_methods; ++m) {
                double std_val = pgf->time_ms_std[g][q][method_indices[m]];
                double mean_val = pgf->time_ms[g][q][method_indices[m]];
                if (mean_val > 0 && norm > 0 && runs_q >= 2) {
                    double ci95 = tc * (std_val / sqrt((double)runs_q)) * 1000.0 / norm;
                    printf("  %-14.6f", ci95);
                } else {
                    printf("  %-14s", "nan");
                }
            }
            printf("\n");
        }
        printf("\n");

        /* min_per_gate — minimum run time normalised per gate application */
        printf("# %s min_per_gate (µs per gate application)\n", all_gate_names[g]);
        printf("%-8s", "qubits");
        for (int m = 0; m < num_methods; ++m) printf("  %-14s", method_names[m]);
        printf("\n");
        for (int q = 0; q < pgf->num_configs; ++q) {
            qubit_t qubits    = pgf->qubits[q];
            double  norm      = (double)pgf->iterations * (double)gate_sweep_size(g, qubits);
            printf("%-8u", qubits);
            for (int m = 0; m < num_methods; ++m) {
                double val = pgf->time_ms_min[g][q][method_indices[m]];
                if (val > 0 && norm > 0) printf("  %-14.6f", val * 1000.0 / norm);
                else                     printf("  %-14s", "nan");
            }
            printf("\n");
        }
        printf("\n");

        /* median_per_gate — median run time normalised per gate application */
        printf("# %s median_per_gate (µs per gate application)\n", all_gate_names[g]);
        printf("%-8s", "qubits");
        for (int m = 0; m < num_methods; ++m) printf("  %-14s", method_names[m]);
        printf("\n");
        for (int q = 0; q < pgf->num_configs; ++q) {
            qubit_t qubits    = pgf->qubits[q];
            double  norm      = (double)pgf->iterations * (double)gate_sweep_size(g, qubits);
            printf("%-8u", qubits);
            for (int m = 0; m < num_methods; ++m) {
                double val = pgf->time_ms_median[g][q][method_indices[m]];
                if (val > 0 && norm > 0) printf("  %-14.6f", val * 1000.0 / norm);
                else                     printf("  %-14s", "nan");
            }
            printf("\n");
        }
        printf("\n");

        /* cv — coefficient of variation (dimensionless; invariant under normalization) */
        printf("# %s cv (coefficient of variation, %%)\n", all_gate_names[g]);
        printf("%-8s", "qubits");
        for (int m = 0; m < num_methods; ++m) printf("  %-14s", method_names[m]);
        printf("\n");
        for (int q = 0; q < pgf->num_configs; ++q) {
            printf("%-8u", pgf->qubits[q]);
            for (int m = 0; m < num_methods; ++m) {
                double mean_val = pgf->time_ms[g][q][method_indices[m]];
                if (mean_val > 0) printf("  %-14.6f", pgf->time_ms_cv[g][q][method_indices[m]]);
                else              printf("  %-14s", "nan");
            }
            printf("\n");
        }
        printf("\n");
    }

    /* ── Memory table ──────────────────────────────────────────────────────── */
    printf("# Memory usage (MB)\n");
    printf("%-8s", "qubits");
    for (int m = 0; m < num_methods; ++m) printf("  %-14s", method_names[m]);
    printf("\n");

    for (int q = 0; q < pgf->num_configs; ++q) {
        printf("%-8u", pgf->qubits[q]);
        for (int m = 0; m < num_methods; ++m) {
            size_t bytes = pgf->memory[0][q][method_indices[m]];
            if (bytes > 0) printf("  %-14.3f", (double)bytes / (1024.0 * 1024.0));
            else           printf("  %-14s", "nan");
        }
        printf("\n");
    }
}

/*
 * =====================================================================================================================
 * Per-qubit pgfplots output
 * =====================================================================================================================
 */

/**
 * @brief Reconstruct the pair (q1, q2) for pair index t in a 2Q gate sweep.
 *
 * Mirrors the ordering produced by build_all_pairs (is_symmetric=1) and
 * build_all_ordered_pairs (is_symmetric=0).
 */
static int reconstruct_pair(qubit_t qubits, int t, int is_symmetric,
                             qubit_t *q1_out, qubit_t *q2_out)
{
    int idx = 0;
    if (is_symmetric) {
        /* build_all_pairs: q1 < q2 */
        for (qubit_t q1 = 0; q1 < qubits; q1++) {
            for (qubit_t q2 = q1 + 1; q2 < qubits; q2++) {
                if (idx == t) { *q1_out = q1; *q2_out = q2; return 1; }
                idx++;
            }
        }
    } else {
        /* build_all_ordered_pairs: all i != j */
        for (qubit_t i = 0; i < qubits; i++) {
            for (qubit_t j = 0; j < qubits; j++) {
                if (i == j) continue;
                if (idx == t) { *q1_out = i; *q2_out = j; return 1; }
                idx++;
            }
        }
    }
    return 0;
}

/**
 * @brief Check whether gate g has any populated data in the accumulator.
 */
static int gate_has_data(const pgfplots_perq_data_t *pgf, int g)
{
    for (int q = 0; q < pgf->num_configs; ++q) {
        for (int m = 0; m < NUM_METHODS; ++m) {
            int nt = pgf->n_targets[g][q][m];
            for (int t = 0; t < nt; ++t) {
                if (pgf->cells[g][q][m][t].time_per_gate_us > 0.0)
                    return 1;
            }
        }
    }
    return 0;
}

/* ── Per-qubit pgfplots table field selector ──────────────────────────────── */
typedef enum {
    PERQ_FIELD_TIME,    /* time_per_gate_us (mean) */
    PERQ_FIELD_CI95,    /* 95% CI half-width */
    PERQ_FIELD_MIN,     /* time_per_gate_us_min */
    PERQ_FIELD_MEDIAN,  /* time_per_gate_us_median */
    PERQ_FIELD_CV       /* time_ms_cv */
} perq_field_t;

/*
 * emit_perq_table — emit one .dat block for a single (gate, qubit_count, field).
 */
static void emit_perq_table(const pgfplots_perq_data_t *pgf,
                             const char *gname, int g, int q,
                             int nt, int is_2q, int is_sym,
                             perq_field_t field,
                             const char *header_suffix,
                             const char * const *method_names,
                             const int *method_indices, int num_methods)
{
    qubit_t qubits = pgf->qubits[q];

    printf("# %s %uq %s\n", gname, (unsigned)qubits, header_suffix);
    printf("%-14s", is_2q ? "pair" : "target");
    for (int m = 0; m < num_methods; ++m) printf("  %-14s", method_names[m]);
    printf("\n");

    for (int t = 0; t < nt; ++t) {
        /* Print row label */
        if (is_2q) {
            qubit_t q1 = 0, q2 = 0;
            if (!reconstruct_pair(qubits, t, is_sym, &q1, &q2)) abort();
            char buf[32];
            snprintf(buf, sizeof(buf), "q%u_q%u", (unsigned)q1, (unsigned)q2);
            printf("%-14s", buf);
        } else {
            printf("%-14u", (unsigned)t);
        }

        /* Print one value per method */
        for (int m = 0; m < num_methods; ++m) {
            int mi = method_indices[m];
            const perq_cell_t *c = &pgf->cells[g][q][mi][t];
            double val = 0.0;
            int    have = 0;

            switch (field) {
            case PERQ_FIELD_TIME:
                val  = c->time_per_gate_us;
                have = (val > 0.0);
                break;
            case PERQ_FIELD_CI95:
                if (c->time_per_gate_us > 0.0 && c->runs >= 2 && c->iterations > 0) {
                    double tc = t_crit(c->runs - 1);
                    val  = tc * (c->time_ms_std / sqrt((double)c->runs))
                           / (double)c->iterations * 1000.0;
                    have = 1;
                }
                break;
            case PERQ_FIELD_MIN:
                val  = c->time_per_gate_us_min;
                have = (val > 0.0);
                break;
            case PERQ_FIELD_MEDIAN:
                val  = c->time_per_gate_us_median;
                have = (val > 0.0);
                break;
            case PERQ_FIELD_CV:
                val  = c->time_ms_cv;
                have = (c->time_per_gate_us > 0.0);
                break;
            }

            if (have) printf("  %-14.6f", val);
            else      printf("  %-14s", "nan");
        }
        printf("\n");
    }
    printf("\n");
}

void print_pgfplots_perq_output(const pgfplots_perq_data_t *pgf)
{
    /* Method index mapping — identical to print_pgfplots_output(); M_NAIVE excluded */
    static const char *method_names[]   = {"qlib_packed", "qlib_tiled", "blas",
                                            "quest",       "qulacs",     "qiskit_aer_dm"};
    static const int   method_indices[] = {M_QLIB, M_QLIB_TILED, M_BLAS,
                                            M_QUEST, M_QULACS, M_AER_DM};
    static const int   num_methods      = (int)(sizeof(method_indices) / sizeof(method_indices[0]));

    for (int g = 0; g < NUM_GATES; ++g) {

        /* Bug 2 fix: skip gates with no data (e.g. when --gate X filters output) */
        if (!gate_has_data(pgf, g))
            continue;

        const char *gname  = all_gate_names[g];
        int         is_2q  = (g >= NUM_1Q_GATES);
        int         g2     = g - NUM_1Q_GATES;   /* 2Q gate index; valid only when is_2q */
        /* Read .is_symmetric directly from the extern table. */
        int         is_sym = is_2q ? gates_2q_perq[g2].is_symmetric : 0;

        /*
         * Bug 1 fix: emit one block per (gate, qubit_count).
         * Within each block, rows are indexed by target position (x-axis),
         * and columns are methods.  This preserves per-target variation instead
         * of averaging it away.
         */
        for (int q = 0; q < pgf->num_configs; ++q) {
            /* Determine n_targets for this (gate, qubit_config) from the first
             * method that has data.  All methods should agree on n_targets. */
            int nt = 0;
            for (int m = 0; m < NUM_METHODS && nt == 0; ++m)
                nt = pgf->n_targets[g][q][m];
            if (nt == 0)
                continue;  /* no data for this qubit count — skip block */

            emit_perq_table(pgf, gname, g, q, nt, is_2q, is_sym,
                            PERQ_FIELD_TIME,
                            "time_per_gate_us (\xc2\xb5s/gate application per target)",
                            method_names, method_indices, num_methods);
            emit_perq_table(pgf, gname, g, q, nt, is_2q, is_sym,
                            PERQ_FIELD_CI95,
                            "ci95_per_gate (\xc2\xb5s/gate, 95% CI half-width per target)",
                            method_names, method_indices, num_methods);
            emit_perq_table(pgf, gname, g, q, nt, is_2q, is_sym,
                            PERQ_FIELD_MIN,
                            "min_per_gate (\xc2\xb5s/gate, minimum run per target)",
                            method_names, method_indices, num_methods);
            emit_perq_table(pgf, gname, g, q, nt, is_2q, is_sym,
                            PERQ_FIELD_MEDIAN,
                            "median_per_gate (\xc2\xb5s/gate, median run per target)",
                            method_names, method_indices, num_methods);
            emit_perq_table(pgf, gname, g, q, nt, is_2q, is_sym,
                            PERQ_FIELD_CV,
                            "cv (coefficient of variation, % per target)",
                            method_names, method_indices, num_methods);

        } /* qubit-config loop */
    } /* gate loop */

    /* ── Memory table ────────────────────────────────────────────────────── */
    /* Memory is constant across targets for a given (method, qubit count).  */
    /* Bug 3 fix: scan all gate indices for the first one with data rather   */
    /* than hardcoding g=0, so --gate CX works correctly.                    */
    printf("# Memory usage (MB)\n");
    printf("%-8s", "qubits");
    for (int m = 0; m < num_methods; ++m) printf("  %-14s", method_names[m]);
    printf("\n");
    for (int q = 0; q < pgf->num_configs; ++q) {
        printf("%-8u", (unsigned)pgf->qubits[q]);
        for (int m = 0; m < num_methods; ++m) {
            int mi = method_indices[m];
            size_t bytes = 0;
            /* Scan all gates for the first populated cell with memory data */
            for (int g = 0; g < NUM_GATES && bytes == 0; ++g) {
                int nt = pgf->n_targets[g][q][mi];
                for (int t = 0; t < nt && bytes == 0; ++t) {
                    if (pgf->cells[g][q][mi][t].memory_bytes > 0)
                        bytes = pgf->cells[g][q][mi][t].memory_bytes;
                }
            }
            if (bytes > 0) printf("  %-14.3f", (double)bytes / (1024.0 * 1024.0));
            else           printf("  %-14s", "nan");
        }
        printf("\n");
    }
}
