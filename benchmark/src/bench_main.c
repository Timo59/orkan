/**
 * @file bench_main.c
 * @brief Orchestration, output formatting, CLI parsing, gate tables, and main() for bench_mixed
 *
 * Calls benchmark functions from bench_mixed.c, bench_baselines.c, bench_quest.c,
 * and bench_qulacs.cpp — all declared in bench.h.
 */

#include "bench.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*
 * =====================================================================================================================
 * Output functions
 * =====================================================================================================================
 */

void bench_print_result(const bench_result_t *result, int verbose) {
    printf("  %-14s %8.3f ± %6.3f ms (CV %4.1f%%)  %'12.0f ops/sec",
           result->method,
           result->time_ms, result->time_ms_std, result->time_ms_cv,
           result->ops_per_sec);
    if (verbose)
        printf("  [min %8.3f ms  med %8.3f ms  mem %zu B]",
               result->time_ms_min, result->time_ms_median, result->memory_bytes);
    printf("\n");
}

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

/*
 * =====================================================================================================================
 * Summary struct and shared speedup/memory printing
 * =====================================================================================================================
 *
 * Collects results for one (gate, qubit-count) combination across all methods.
 * Eliminates the duplicate speedup+memory display logic between the 1Q and 2Q loops.
 */

typedef struct {
    bench_result_t packed;
    bench_result_t tiled;   int has_tiled;
    bench_result_t dense;   int has_dense;
    bench_result_t naive;   int has_naive;
    double naive_scale;     /* multiply naive speedup by this when naive ran fewer iterations */
    bench_result_t quest;   int has_quest;
    bench_result_t qulacs;  int has_qulacs;
} bench_summary_t;

static void print_speedup_and_memory(const bench_summary_t *s) {
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
    (void)has;
    printf("\n");

    /* Memory line */
    size_t max_mem = s->packed.memory_bytes;
    if (s->has_tiled  && s->tiled.memory_bytes  > max_mem) max_mem = s->tiled.memory_bytes;
    if (s->has_dense  && s->dense.memory_bytes  > max_mem) max_mem = s->dense.memory_bytes;
    if (s->has_naive  && s->naive.memory_bytes  > max_mem) max_mem = s->naive.memory_bytes;
    if (s->has_quest  && s->quest.memory_bytes  > max_mem) max_mem = s->quest.memory_bytes;
    if (s->has_qulacs && s->qulacs.memory_bytes > max_mem) max_mem = s->qulacs.memory_bytes;

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
    if (s->has_dense && s->packed.memory_bytes > 0 && s->dense.memory_bytes > 0) {
        double savings = 100.0 * (1.0 - (double)s->packed.memory_bytes / (double)s->dense.memory_bytes);
        printf(" (packed saves %.0f%%)", savings);
    }
    printf("\n");
}

/*
 * =====================================================================================================================
 * Command-line parsing
 * =====================================================================================================================
 */

static void print_usage(const char *prog) {
    printf("Usage: %s [options]\n\n", prog);
    printf("Options:\n");
    printf("  --min-qubits N   Minimum qubit count, 1..30 (default: %d)\n", BENCH_DEFAULT_MIN_QUBITS);
    printf("  --max-qubits N   Maximum qubit count, 1..30 (default: %d)\n", BENCH_DEFAULT_MAX_QUBITS);
    printf("  --step N         Qubit count increment, >=1 (default: %d)\n", BENCH_DEFAULT_STEP);
    printf("  --iterations N   Gate calls per timed run, >=1 (default: %d)\n", BENCH_DEFAULT_ITERATIONS);
    printf("  --warmup N       Warm-up iterations before timing, >=0 (default: %d)\n", BENCH_DEFAULT_WARMUP);
    printf("  --runs N         Independent timing runs for statistics, 1..200 (default: %d)\n", BENCH_DEFAULT_RUNS);
    printf("  --csv            Output in CSV format\n");
    printf("  --pgfplots       Output in pgfplots-compatible .dat format (includes std dev)\n");
    printf("  --verbose        Show min, median, and memory per result\n");
    printf("  --help           Show this help\n");
}

bench_options_t bench_parse_options(int argc, char *argv[]) {
    bench_options_t opts = {
        .min_qubits      = BENCH_DEFAULT_MIN_QUBITS,
        .max_qubits      = BENCH_DEFAULT_MAX_QUBITS,
        .step            = BENCH_DEFAULT_STEP,
        .iterations      = BENCH_DEFAULT_ITERATIONS,
        .warmup          = BENCH_DEFAULT_WARMUP,
        .runs            = BENCH_DEFAULT_RUNS,
        .csv_output      = 0,
        .pgfplots_output = 0,
        .verbose         = 0
    };

    for (int i = 1; i < argc; ++i) {
        if        (strcmp(argv[i], "--min-qubits")  == 0 && i + 1 < argc) {
            opts.min_qubits  = (qubit_t)atoi(argv[++i]);
        } else if (strcmp(argv[i], "--max-qubits")  == 0 && i + 1 < argc) {
            opts.max_qubits  = (qubit_t)atoi(argv[++i]);
        } else if (strcmp(argv[i], "--step")        == 0 && i + 1 < argc) {
            opts.step        = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--iterations")  == 0 && i + 1 < argc) {
            opts.iterations  = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--warmup")      == 0 && i + 1 < argc) {
            opts.warmup      = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--runs")        == 0 && i + 1 < argc) {
            opts.runs        = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--csv")         == 0) {
            opts.csv_output  = 1;
        } else if (strcmp(argv[i], "--pgfplots")    == 0) {
            opts.pgfplots_output = 1;
        } else if (strcmp(argv[i], "--verbose")     == 0) {
            opts.verbose     = 1;
        } else if (strcmp(argv[i], "--help")        == 0) {
            print_usage(argv[0]);
            exit(0);
        }
    }

    /* Input validation */
    if (opts.min_qubits < 1 || opts.min_qubits > 30) {
        fprintf(stderr, "error: --min-qubits must be between 1 and 30\n");
        exit(1);
    }
    if (opts.max_qubits < opts.min_qubits || opts.max_qubits > 30) {
        fprintf(stderr, "error: --max-qubits must be between --min-qubits (%u) and 30\n",
                opts.min_qubits);
        exit(1);
    }
    if (opts.step <= 0) {
        fprintf(stderr, "error: --step must be >= 1\n");
        exit(1);
    }
    if (opts.iterations <= 0) {
        fprintf(stderr, "error: --iterations must be >= 1\n");
        exit(1);
    }
    if (opts.warmup < 0) {
        fprintf(stderr, "error: --warmup must be >= 0\n");
        exit(1);
    }
    if (opts.runs < 1 || opts.runs > 200) {
        fprintf(stderr, "error: --runs must be between 1 and 200\n");
        exit(1);
    }

    return opts;
}

/*
 * =====================================================================================================================
 * Gate tables
 * =====================================================================================================================
 */

/* Gate matrix constants defined in bench_baselines.c */
extern const cplx_t XMAT[4];
extern const cplx_t ZMAT[4];
extern const cplx_t HMAT[4];
extern const cplx_t CXMAT[16];
extern const cplx_t SWAPMAT[16];

typedef struct {
    const char *name;
    void (*fn)(state_t*, qubit_t);
    const cplx_t *mat;
} gate_def_t;

static const gate_def_t gates[] = {
    {"X", x, XMAT},
    {"H", h, HMAT},
    {"Z", z, ZMAT},
    {NULL, NULL, NULL}
};

typedef struct {
    const char *name;
    void (*fn)(state_t*, qubit_t, qubit_t);
    const cplx_t *mat;
    int has_tiled;
} gate_def_2q_t;

static const gate_def_2q_t gates_2q[] = {
    {"CX",   cx,        CXMAT,   1},
    {"SWAP", swap_gate, SWAPMAT, 1},
    {NULL, NULL, NULL, 0}
};

#define NUM_1Q_GATES    3
#define NUM_2Q_GATES    2
#define NUM_GATES       (NUM_1Q_GATES + NUM_2Q_GATES)
#define MAX_QUBIT_CONFIGS 32
#define NUM_METHODS     6   /* qlib_packed, qlib_tiled, blas, naive, quest, qulacs */

static const char *all_gate_names[NUM_GATES] = {"X", "H", "Z", "CX", "SWAP"};

/*
 * =====================================================================================================================
 * pgfplots output
 * =====================================================================================================================
 */

typedef struct {
    qubit_t qubits[MAX_QUBIT_CONFIGS];
    int     num_configs;
    int     iterations;                                                /**< Gate calls per timed run (shared across methods) */
    double  time_ms[NUM_GATES][MAX_QUBIT_CONFIGS][NUM_METHODS];        /**< Mean time (ms) */
    double  time_ms_std[NUM_GATES][MAX_QUBIT_CONFIGS][NUM_METHODS];    /**< Std dev (ms) for error bars */
    double  time_ms_min[NUM_GATES][MAX_QUBIT_CONFIGS][NUM_METHODS];    /**< Min time (ms) */
    double  time_ms_median[NUM_GATES][MAX_QUBIT_CONFIGS][NUM_METHODS]; /**< Median time (ms) */
    double  time_ms_cv[NUM_GATES][MAX_QUBIT_CONFIGS][NUM_METHODS];     /**< CV = std/mean × 100 (%) */
    int     runs[MAX_QUBIT_CONFIGS];                                   /**< Number of timing runs (may vary; stored per qubit-config) */
    size_t  memory[NUM_GATES][MAX_QUBIT_CONFIGS][NUM_METHODS];
} pgfplots_data_t;

enum { M_QLIB = 0, M_QLIB_TILED = 1, M_BLAS = 2, M_NAIVE = 3, M_QUEST = 4, M_QULACS = 5 };

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
static double t_crit(int df) {
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
static int gate_sweep_size(int g, qubit_t qubits) {
    if (g < NUM_1Q_GATES)
        return (int)qubits;
    else
        return (int)(qubits * (qubits - 1) / 2);
}

static void print_pgfplots_output(const pgfplots_data_t *pgf) {
    const char *method_names[]   = {"qlib_packed", "qlib_tiled", "blas", "quest", "qulacs"};
    int         method_indices[] = {M_QLIB, M_QLIB_TILED, M_BLAS, M_QUEST, M_QULACS};
    int         num_methods      = 5;

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
                else              printf("  %-14s", "0.000000");
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
 * Main
 * =====================================================================================================================
 */

int main(int argc, char *argv[]) {
    bench_options_t opts = bench_parse_options(argc, argv);

    /* Allocate pgfplots storage only when --pgfplots is requested */
    pgfplots_data_t *pgf = NULL;
    if (opts.pgfplots_output) {
        pgf = calloc(1, sizeof(pgfplots_data_t));
        if (!pgf) {
            fprintf(stderr, "error: failed to allocate pgfplots buffer\n");
            return 1;
        }
        pgf->iterations = opts.iterations;
    }

    if (opts.csv_output) {
        bench_print_csv_header();
    } else if (!opts.pgfplots_output) {
        printf("\nMixed State Gate Benchmarks\n");
        printf("===========================\n");
        printf("Comparing: qlib packed sparse vs BLAS dense vs naive loop");
#ifdef WITH_QUEST
        printf(" vs QuEST");
#endif
#ifdef WITH_QULACS
        printf(" vs Qulacs");
#endif
        printf("\n");
        printf("Iterations: %d, Warm-up: %d, Runs: %d\n",
               opts.iterations, opts.warmup, opts.runs);
#ifdef WITH_QUEST
        printf("QuEST: enabled\n");
#endif
#ifdef WITH_QULACS
        printf("Qulacs: enabled\n");
#endif
        printf("\n");
    }

    int qi = 0;

    /*
     * Annotate a result with sweep_size and derived time_per_gate_ms.
     * Called just before bench_print_csv() so the CSV row carries normalized values.
     * sweep_size: number of gate applications per timed iteration.
     *   1Q gates: qubits (one application per qubit per iteration)
     *   2Q gates: qubits*(qubits-1)/2 (one application per ordered pair per iteration)
     */
#define SET_SWEEP(res, sw) do {                                                  \
    (res).sweep_size       = (sw);                                               \
    double _norm = (double)(res).iterations * (double)(sw);                      \
    (res).time_per_gate_ms = (_norm > 0.0) ? (res).time_ms / _norm : 0.0;       \
} while (0)

    for (qubit_t qubits = opts.min_qubits; qubits <= opts.max_qubits; qubits += opts.step) {
        dim_t dim    = (dim_t)1 << qubits;
        int run_dense = (qubits <= 8);  /* O(N³) BLAS/naive: practical only up to 8 qubits */

        if (pgf && qi < MAX_QUBIT_CONFIGS)
            pgf->qubits[qi] = qubits;

        if (!opts.csv_output && !opts.pgfplots_output)
            printf("qubits: %u, dim: %lld\n", qubits, (long long)dim);

        /* ── One-qubit gates ─────────────────────────────────────────────────── */
        for (int g = 0; gates[g].name != NULL; ++g) {
            if (!opts.csv_output && !opts.pgfplots_output)
                printf("\nGate: %s\n", gates[g].name);

            bench_summary_t s = {0};
            s.naive_scale = 10.0;  /* naive runs at iterations/10 due to O(N³) cost */

            srand(42); s.packed = bench_qlib_packed(qubits, gates[g].name, gates[g].fn,
                                                     opts.iterations, opts.warmup, opts.runs);
            srand(42); s.tiled  = bench_qlib_tiled(qubits, gates[g].name, gates[g].fn,
                                                    opts.iterations, opts.warmup, opts.runs);
            s.has_tiled = (s.tiled.time_ms > 0);

            if (run_dense) {
                srand(42); s.dense = bench_blas_dense(qubits, gates[g].name, gates[g].mat,
                                                       opts.iterations, opts.warmup, opts.runs);
                s.has_dense = 1;

                srand(42); s.naive = bench_naive_loop(qubits, gates[g].name, gates[g].mat,
                                                       opts.iterations / 10,
                                                       opts.warmup    / 10,
                                                       opts.runs);
                s.has_naive = 1;
            }

#ifdef WITH_QUEST
            s.quest     = bench_quest(qubits, gates[g].name, opts.iterations, opts.warmup, opts.runs);
            s.has_quest = (s.quest.time_ms > 0);
#endif
#ifdef WITH_QULACS
            s.qulacs     = bench_qulacs(qubits, gates[g].name, opts.iterations, opts.warmup, opts.runs);
            s.has_qulacs = (s.qulacs.time_ms > 0);
#endif

            if (pgf && qi < MAX_QUBIT_CONFIGS) {
                pgf->runs[qi]                          = opts.runs;
                pgf->time_ms[g][qi][M_QLIB]           = s.packed.time_ms;
                pgf->time_ms[g][qi][M_QLIB_TILED]     = s.tiled.time_ms;
                pgf->time_ms[g][qi][M_BLAS]            = s.dense.time_ms;
                pgf->time_ms[g][qi][M_NAIVE]           = s.naive.time_ms;
                pgf->time_ms_std[g][qi][M_QLIB]        = s.packed.time_ms_std;
                pgf->time_ms_std[g][qi][M_QLIB_TILED]  = s.tiled.time_ms_std;
                pgf->time_ms_std[g][qi][M_BLAS]        = s.dense.time_ms_std;
                pgf->time_ms_std[g][qi][M_NAIVE]       = s.naive.time_ms_std;
                pgf->time_ms_min[g][qi][M_QLIB]        = s.packed.time_ms_min;
                pgf->time_ms_min[g][qi][M_QLIB_TILED]  = s.tiled.time_ms_min;
                pgf->time_ms_min[g][qi][M_BLAS]        = s.dense.time_ms_min;
                pgf->time_ms_min[g][qi][M_NAIVE]       = s.naive.time_ms_min;
                pgf->time_ms_median[g][qi][M_QLIB]        = s.packed.time_ms_median;
                pgf->time_ms_median[g][qi][M_QLIB_TILED]  = s.tiled.time_ms_median;
                pgf->time_ms_median[g][qi][M_BLAS]        = s.dense.time_ms_median;
                pgf->time_ms_median[g][qi][M_NAIVE]       = s.naive.time_ms_median;
                pgf->time_ms_cv[g][qi][M_QLIB]         = s.packed.time_ms_cv;
                pgf->time_ms_cv[g][qi][M_QLIB_TILED]   = s.tiled.time_ms_cv;
                pgf->time_ms_cv[g][qi][M_BLAS]         = s.dense.time_ms_cv;
                pgf->time_ms_cv[g][qi][M_NAIVE]        = s.naive.time_ms_cv;
                pgf->memory[g][qi][M_QLIB]             = s.packed.memory_bytes;
                pgf->memory[g][qi][M_QLIB_TILED]       = s.tiled.memory_bytes;
                pgf->memory[g][qi][M_BLAS]             = s.dense.memory_bytes;
                pgf->memory[g][qi][M_NAIVE]            = s.naive.memory_bytes;
#ifdef WITH_QUEST
                pgf->time_ms[g][qi][M_QUEST]           = s.quest.time_ms;
                pgf->time_ms_std[g][qi][M_QUEST]       = s.quest.time_ms_std;
                pgf->time_ms_min[g][qi][M_QUEST]       = s.quest.time_ms_min;
                pgf->time_ms_median[g][qi][M_QUEST]    = s.quest.time_ms_median;
                pgf->time_ms_cv[g][qi][M_QUEST]        = s.quest.time_ms_cv;
                pgf->memory[g][qi][M_QUEST]            = s.quest.memory_bytes;
#endif
#ifdef WITH_QULACS
                pgf->time_ms[g][qi][M_QULACS]          = s.qulacs.time_ms;
                pgf->time_ms_std[g][qi][M_QULACS]      = s.qulacs.time_ms_std;
                pgf->time_ms_min[g][qi][M_QULACS]      = s.qulacs.time_ms_min;
                pgf->time_ms_median[g][qi][M_QULACS]   = s.qulacs.time_ms_median;
                pgf->time_ms_cv[g][qi][M_QULACS]       = s.qulacs.time_ms_cv;
                pgf->memory[g][qi][M_QULACS]           = s.qulacs.memory_bytes;
#endif
            }

            if (opts.csv_output) {
                int sw1q = (int)qubits;
                SET_SWEEP(s.packed, sw1q);
                SET_SWEEP(s.tiled,  sw1q);
                bench_print_csv(&s.packed);
                bench_print_csv(&s.tiled);
                if (s.has_dense) { SET_SWEEP(s.dense, sw1q); bench_print_csv(&s.dense); }
                if (s.has_naive) { SET_SWEEP(s.naive, sw1q); bench_print_csv(&s.naive); }
#ifdef WITH_QUEST
                if (s.has_quest) { SET_SWEEP(s.quest, sw1q); bench_print_csv(&s.quest); }
#endif
#ifdef WITH_QULACS
                if (s.has_qulacs) { SET_SWEEP(s.qulacs, sw1q); bench_print_csv(&s.qulacs); }
#endif
            } else if (!opts.pgfplots_output) {
                bench_print_result(&s.packed, opts.verbose);
                bench_print_result(&s.tiled,  opts.verbose);
                if (s.has_dense) bench_print_result(&s.dense, opts.verbose);
                if (s.has_naive) bench_print_result(&s.naive, opts.verbose);
#ifdef WITH_QUEST
                if (s.has_quest) bench_print_result(&s.quest, opts.verbose);
#endif
#ifdef WITH_QULACS
                if (s.has_qulacs) bench_print_result(&s.qulacs, opts.verbose);
#endif
                print_speedup_and_memory(&s);
            }
        }

        /* ── Two-qubit gates ─────────────────────────────────────────────────── */
        for (int g = 0; gates_2q[g].name != NULL; ++g) {
            if (!opts.csv_output && !opts.pgfplots_output)
                printf("\nGate: %s\n", gates_2q[g].name);

            bench_summary_t s = {0};
            s.naive_scale = 1.0;

            srand(42); s.packed = bench_qlib_packed_2q(qubits, gates_2q[g].name, gates_2q[g].fn,
                                                        opts.iterations, opts.warmup, opts.runs);

            if (gates_2q[g].has_tiled) {
                srand(42); s.tiled = bench_qlib_tiled_2q(qubits, gates_2q[g].name, gates_2q[g].fn,
                                                          opts.iterations, opts.warmup, opts.runs);
                s.has_tiled = (s.tiled.time_ms > 0);
            }

            if (run_dense) {
                srand(42); s.dense = bench_blas_dense_2q(qubits, gates_2q[g].name, gates_2q[g].mat,
                                                          opts.iterations, opts.warmup, opts.runs);
                s.has_dense = 1;
            }

#ifdef WITH_QUEST
            s.quest     = bench_quest(qubits, gates_2q[g].name, opts.iterations, opts.warmup, opts.runs);
            s.has_quest = (s.quest.time_ms > 0);
#endif
#ifdef WITH_QULACS
            s.qulacs     = bench_qulacs(qubits, gates_2q[g].name, opts.iterations, opts.warmup, opts.runs);
            s.has_qulacs = (s.qulacs.time_ms > 0);
#endif

            int pg = NUM_1Q_GATES + g;
            if (pgf && qi < MAX_QUBIT_CONFIGS) {
                /* pgf->runs[qi] already written in the 1Q loop for this qubit config */
                pgf->time_ms[pg][qi][M_QLIB]           = s.packed.time_ms;
                pgf->time_ms[pg][qi][M_QLIB_TILED]     = s.tiled.time_ms;
                pgf->time_ms[pg][qi][M_BLAS]            = s.dense.time_ms;
                pgf->time_ms_std[pg][qi][M_QLIB]        = s.packed.time_ms_std;
                pgf->time_ms_std[pg][qi][M_QLIB_TILED]  = s.tiled.time_ms_std;
                pgf->time_ms_std[pg][qi][M_BLAS]        = s.dense.time_ms_std;
                pgf->time_ms_min[pg][qi][M_QLIB]        = s.packed.time_ms_min;
                pgf->time_ms_min[pg][qi][M_QLIB_TILED]  = s.tiled.time_ms_min;
                pgf->time_ms_min[pg][qi][M_BLAS]        = s.dense.time_ms_min;
                pgf->time_ms_median[pg][qi][M_QLIB]        = s.packed.time_ms_median;
                pgf->time_ms_median[pg][qi][M_QLIB_TILED]  = s.tiled.time_ms_median;
                pgf->time_ms_median[pg][qi][M_BLAS]        = s.dense.time_ms_median;
                pgf->time_ms_cv[pg][qi][M_QLIB]         = s.packed.time_ms_cv;
                pgf->time_ms_cv[pg][qi][M_QLIB_TILED]   = s.tiled.time_ms_cv;
                pgf->time_ms_cv[pg][qi][M_BLAS]         = s.dense.time_ms_cv;
                pgf->memory[pg][qi][M_QLIB]             = s.packed.memory_bytes;
                pgf->memory[pg][qi][M_QLIB_TILED]       = s.tiled.memory_bytes;
                pgf->memory[pg][qi][M_BLAS]             = s.dense.memory_bytes;
#ifdef WITH_QUEST
                pgf->time_ms[pg][qi][M_QUEST]           = s.quest.time_ms;
                pgf->time_ms_std[pg][qi][M_QUEST]       = s.quest.time_ms_std;
                pgf->time_ms_min[pg][qi][M_QUEST]       = s.quest.time_ms_min;
                pgf->time_ms_median[pg][qi][M_QUEST]    = s.quest.time_ms_median;
                pgf->time_ms_cv[pg][qi][M_QUEST]        = s.quest.time_ms_cv;
                pgf->memory[pg][qi][M_QUEST]            = s.quest.memory_bytes;
#endif
#ifdef WITH_QULACS
                pgf->time_ms[pg][qi][M_QULACS]          = s.qulacs.time_ms;
                pgf->time_ms_std[pg][qi][M_QULACS]      = s.qulacs.time_ms_std;
                pgf->time_ms_min[pg][qi][M_QULACS]      = s.qulacs.time_ms_min;
                pgf->time_ms_median[pg][qi][M_QULACS]   = s.qulacs.time_ms_median;
                pgf->time_ms_cv[pg][qi][M_QULACS]       = s.qulacs.time_ms_cv;
                pgf->memory[pg][qi][M_QULACS]           = s.qulacs.memory_bytes;
#endif
            }

            if (opts.csv_output) {
                int sw2q = (int)(qubits * (qubits - 1) / 2);
                SET_SWEEP(s.packed, sw2q);
                bench_print_csv(&s.packed);
                if (s.has_tiled) { SET_SWEEP(s.tiled, sw2q); bench_print_csv(&s.tiled); }
                if (s.has_dense) { SET_SWEEP(s.dense, sw2q); bench_print_csv(&s.dense); }
#ifdef WITH_QUEST
                if (s.has_quest) { SET_SWEEP(s.quest, sw2q); bench_print_csv(&s.quest); }
#endif
#ifdef WITH_QULACS
                if (s.has_qulacs) { SET_SWEEP(s.qulacs, sw2q); bench_print_csv(&s.qulacs); }
#endif
            } else if (!opts.pgfplots_output) {
                bench_print_result(&s.packed, opts.verbose);
                if (s.has_tiled) bench_print_result(&s.tiled,  opts.verbose);
                if (s.has_dense) bench_print_result(&s.dense,  opts.verbose);
#ifdef WITH_QUEST
                if (s.has_quest) bench_print_result(&s.quest, opts.verbose);
#endif
#ifdef WITH_QULACS
                if (s.has_qulacs) bench_print_result(&s.qulacs, opts.verbose);
#endif
                print_speedup_and_memory(&s);
            }
        }

        if (!opts.csv_output && !opts.pgfplots_output)
            printf("\n" "─────────────────────────────────────────────────\n");

        if (pgf)
            qi++;
    }

    if (pgf) {
        pgf->num_configs = qi;
        print_pgfplots_output(pgf);
        free(pgf);
    }

#ifdef WITH_QUEST
    bench_quest_cleanup();
#endif

    return 0;
}
