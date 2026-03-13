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
#include <strings.h>   /* strcasecmp — POSIX; available on macOS and Linux */
#include <math.h>      /* ceil */

#if defined(BENCH_AER_ONLY_MODE) && !defined(WITH_AER_DM)
#error "BENCH_AER_ONLY_MODE requires WITH_AER_DM to be defined"
#endif

/*
 * =====================================================================================================================
 * Output functions
 * =====================================================================================================================
 */

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
 * Per-qubit output functions
 * =====================================================================================================================
 */

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
    bench_result_t aer_dm;  int has_aer_dm;   /* NEW */
} bench_summary_t;

static void print_speedup_and_memory(const bench_summary_t *s) {
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
    printf("  --pgfplots       Output pgfplots-compatible .dat tables to stdout.\n"
           "                   With --per-qubit: per-qubit stats aggregated over targets.\n");
    printf("  --verbose        Show min, median, and memory per result\n");
    printf("  --per-qubit      Single-target throughput mode: benchmark each qubit target\n"
           "                   independently (steady-state, not circuit-representative).\n");
    printf("  --gate NAME      Restrict --per-qubit run to a single named gate (e.g. X, CX).\n"
           "                   Has no effect without --per-qubit.\n");
    printf("  --help           Show this help\n");
}

bench_options_t bench_parse_options(int argc, char *argv[]) {
    bench_options_t opts = {
        .min_qubits          = BENCH_DEFAULT_MIN_QUBITS,
        .max_qubits          = BENCH_DEFAULT_MAX_QUBITS,
        .step                = BENCH_DEFAULT_STEP,
        .iterations          = BENCH_DEFAULT_ITERATIONS,
        .warmup              = BENCH_DEFAULT_WARMUP,
        .runs                = BENCH_DEFAULT_RUNS,
        .csv_output          = 0,
        .pgfplots_output     = 0,
        .verbose             = 0,
        .per_qubit           = 0,
        .runs_explicit       = 0,
        .iterations_explicit = 0,
        .gate_filter         = NULL
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
            opts.iterations_explicit = 1;
            if (opts.iterations == 0) {
                fprintf(stderr, "error: --iterations must be > 0\n");
                exit(EXIT_FAILURE);
            }
        } else if (strcmp(argv[i], "--warmup")      == 0 && i + 1 < argc) {
            opts.warmup      = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--runs")        == 0 && i + 1 < argc) {
            opts.runs        = atoi(argv[++i]);
            opts.runs_explicit = 1;
        } else if (strcmp(argv[i], "--per-qubit")   == 0) {
            opts.per_qubit   = 1;
        } else if (strcmp(argv[i], "--gate")        == 0 && i + 1 < argc) {
            opts.gate_filter = argv[++i];
        } else if (strcmp(argv[i], "--csv")         == 0) {
            opts.csv_output  = 1;
        } else if (strcmp(argv[i], "--pgfplots")    == 0) {
            opts.pgfplots_output = 1;
        } else if (strcmp(argv[i], "--verbose")     == 0) {
            opts.verbose     = 1;
        } else if (strcmp(argv[i], "--help")        == 0) {
            print_usage(argv[0]);
            exit(0);
        } else {
            fprintf(stderr, "error: unrecognized option '%s'\n", argv[i]);
            exit(1);
        }
    }

    /* --gate without --per-qubit has no effect in sweep mode; warn the user */
    if (opts.gate_filter != NULL && !opts.per_qubit) {
        fprintf(stderr, "warning: --gate has no effect without --per-qubit\n");
    }

    /* When --per-qubit and --pgfplots are combined, --csv is silently ignored.
     * Warn so the user does not expect CSV output. */
    if (opts.per_qubit && opts.pgfplots_output && opts.csv_output) {
        fprintf(stderr, "warning: --csv is ignored when --per-qubit and --pgfplots "
                "are both specified\n");
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

    /* Per-qubit mode: apply default runs/iterations if not explicitly set */
    if (opts.per_qubit) {
        if (!opts.runs_explicit)       opts.runs       = BENCH_PQ_DEFAULT_RUNS;
        if (!opts.iterations_explicit) opts.iterations = BENCH_PQ_DEFAULT_ITERATIONS;
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
#define NUM_METHODS     7   /* qlib_packed, qlib_tiled, blas, naive, quest, qulacs, aer_dm */
/* Note: M_NAIVE (index 3) is excluded from pgfplots output; print_pgfplots_output
 * uses its own local method list of 6 entries that skips M_NAIVE. */

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

/* Maximum number of target positions (1Q: qubits; 2Q: ordered pairs).
 * Derived from MAX_PAIRS defined in bench.h. */
#define MAX_TARGETS (2 * MAX_PAIRS)  /* 256 */

/**
 * Per-target timing result for --per-qubit mode.
 * Defined here as a forward declaration for perq_cell_from_result().
 * The full per-qubit benchmark infrastructure (perq_main, bench_*_at functions)
 * will be added in a subsequent implementation step.
 */
typedef struct bench_result_perq {
    qubit_t    qubits;
    dim_t      dim;
    const char *gate_name;
    const char *method;
    qubit_t    target;                  /**< Target qubit (q1 for 2Q gates) */
    qubit_t    target2;                 /**< Second qubit for 2Q gates; QUBIT_NONE for 1Q */
    int        is_2q;
    double     time_ms_mean;
    double     time_ms_std;             /**< Sample stddev (Bessel-corrected) */
    double     time_ms_min;
    double     time_ms_median;
    double     time_ms_cv;
    double     ops_per_sec;
    double     time_per_gate_us;        /**< Mean time per gate application (µs) */
    double     time_per_gate_us_median;
    double     time_per_gate_us_min;
    size_t     memory_bytes;
    int        iterations;
    int        runs;
} bench_result_perq_t;

/**
 * Compact per-target timing record for the pgfplots per-qubit accumulator.
 *
 * Stores the 7 values needed by print_pgfplots_perq_output(). No pointer
 * fields: gate_name/method are implicit in the surrounding index.
 * memory_bytes is a plain size_t — no dangling-pointer risk.
 *
 * CV is scale-invariant (CV(c·X) = CV(X)), so time_ms_cv copied from
 * bench_result_perq_t equals the CV of per-gate times. No recomputation needed.
 *
 * sizeof = 56 bytes on LP64: 5×double(40) + size_t(8) + int(4) + 4 pad.
 */
typedef struct {
    double time_per_gate_us;        /**< Mean per-gate time (µs)   */
    double time_per_gate_us_median; /**< Median per-gate time (µs) */
    double time_per_gate_us_min;    /**< Min per-gate time (µs)    */
    double time_ms_cv;              /**< CV = std/mean × 100 (%)   */
    double ops_per_sec;             /**< Gate applications/second  */
    size_t memory_bytes;            /**< State memory footprint    */
    int    iterations;              /**< Calibrated iters per run  */
} perq_cell_t;

/**
 * Accumulator for --per-qubit --pgfplots combined output.
 *
 * Dimensions: NUM_GATES × MAX_QUBIT_CONFIGS × NUM_METHODS × MAX_TARGETS.
 * Heap size: 5×32×7×256×56 ≈ 15.3 MB (cells only; pair arrays add ~6 MB
 * more if qubit_t is 4 bytes). Allocated via calloc() — all cells start at
 * zero, which the output function treats as "missing data" (emits nan).
 *
 * Gate indexing: gates 0..NUM_1Q_GATES-1 are 1Q (X,H,Z); gates
 * NUM_1Q_GATES..NUM_GATES-1 are 2Q (CX,SWAP). Matches perq_main() loop order.
 *
 * M_NAIVE (index 3): always zero-initialized, never written in per-qubit mode
 * (naive_loop has no _at variant). print_pgfplots_perq_output() skips it via
 * method_indices[], matching the pattern in print_pgfplots_output().
 *
 * q1s_2q / q2s_2q: pair coordinates for 2Q gates. Needed because perq_cell_t
 * strips target/target2; required to annotate pair rows in the output tables.
 * For 1Q gates these arrays remain zero (calloc).
 */
typedef struct {
    qubit_t     qubits[MAX_QUBIT_CONFIGS];   /**< Qubit count per config slot   */
    int         num_configs;                  /**< Populated qubit-config slots  */
    int         n_targets[NUM_GATES][MAX_QUBIT_CONFIGS][NUM_METHODS];
    perq_cell_t cells[NUM_GATES][MAX_QUBIT_CONFIGS][NUM_METHODS][MAX_TARGETS];
    qubit_t     q1s_2q[NUM_GATES][MAX_QUBIT_CONFIGS][NUM_METHODS][MAX_TARGETS];
    qubit_t     q2s_2q[NUM_GATES][MAX_QUBIT_CONFIGS][NUM_METHODS][MAX_TARGETS];
} pgfplots_perq_data_t;

/** Copy timing fields from a bench_result_perq_t into a perq_cell_t. */
static inline void perq_cell_from_result(perq_cell_t *cell,
                                          const bench_result_perq_t *r)
{
    cell->time_per_gate_us        = r->time_per_gate_us;
    cell->time_per_gate_us_median = r->time_per_gate_us_median;
    cell->time_per_gate_us_min    = r->time_per_gate_us_min;
    cell->time_ms_cv              = r->time_ms_cv;
    cell->ops_per_sec             = r->ops_per_sec;
    cell->memory_bytes            = r->memory_bytes;
    cell->iterations              = r->iterations;
}

enum { M_QLIB = 0, M_QLIB_TILED = 1, M_BLAS = 2, M_NAIVE = 3,
       M_QUEST = 4, M_QULACS = 5, M_AER_DM = 6 };

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

/**
 * Emit pgfplots-compatible .dat tables for --per-qubit --pgfplots mode.
 *
 * For each of NUM_GATES gates: five tables with qubits as the x-axis,
 * one column per method (6 methods; M_NAIVE excluded). Each row aggregates
 * over all target positions for that qubit count:
 *   time_per_gate_us  — arithmetic mean of per-target means
 *   ci95_per_gate     — all nan (std/runs not stored; saves ~10 MB)
 *   min_per_gate      — global minimum of per-target time_per_gate_us_min
 *   median_per_gate   — arithmetic mean of per-target medians
 *   cv                — arithmetic mean of per-target CVs (scale-invariant)
 *
 * Memory table appended at end: one row per qubit-config, MB units.
 * Missing data (zero cell from calloc) emits literal "nan" string.
 *
 * t_crit() is defined earlier in this file; no forward declaration needed.
 */
static void print_pgfplots_perq_output(const pgfplots_perq_data_t *pgf)
{
    /* Method index mapping — identical to print_pgfplots_output(); M_NAIVE excluded */
    static const char *method_names[]   = {"qlib_packed", "qlib_tiled", "blas",
                                            "quest",       "qulacs",     "qiskit_aer_dm"};
    static const int   method_indices[] = {M_QLIB, M_QLIB_TILED, M_BLAS,
                                            M_QUEST, M_QULACS, M_AER_DM};
    static const int   num_methods      = 6;

    /* Gate names must match the order used when populating pgf->cells[g][...].
     * Look up the actual gate name arrays used in perq_main() (gates_1q_perq,
     * gates_2q_perq) and use their .name fields here. */

    for (int g = 0; g < NUM_GATES; ++g) {
        /* Gate name: use the same static name array as the rest of bench_main.c */
        static const char *pgf_gate_names[NUM_GATES] = {"X", "H", "Z", "CX", "SWAP"};
        const char *gname = pgf_gate_names[g];

        /* ── time_per_gate_us: mean over targets ───────────────────────── */
        printf("# %s time_per_gate_us"
               " (mean µs/gate application, averaged over target positions)\n", gname);
        printf("%-8s", "qubits");
        for (int m = 0; m < num_methods; ++m) printf("  %-14s", method_names[m]);
        printf("\n");
        for (int q = 0; q < pgf->num_configs; ++q) {
            printf("%-8u", (unsigned)pgf->qubits[q]);
            for (int m = 0; m < num_methods; ++m) {
                int mi = method_indices[m];
                int nt = pgf->n_targets[g][q][mi];
                double sum = 0.0; int valid = 0;
                for (int t = 0; t < nt; ++t) {
                    double v = pgf->cells[g][q][mi][t].time_per_gate_us;
                    if (v > 0.0) { sum += v; valid++; }
                }
                if (valid > 0) printf("  %-14.6f", sum / valid);
                else           printf("  %-14s", "nan");
            }
            printf("\n");
        }
        printf("\n");

        /* ── ci95_per_gate: not available (std not stored) ─────────────── */
        printf("# %s ci95_per_gate"
               " (µs/gate) — not available in per-qubit pgfplots mode\n", gname);
        printf("%-8s", "qubits");
        for (int m = 0; m < num_methods; ++m) printf("  %-14s", method_names[m]);
        printf("\n");
        for (int q = 0; q < pgf->num_configs; ++q) {
            printf("%-8u", (unsigned)pgf->qubits[q]);
            for (int m = 0; m < num_methods; ++m) printf("  %-14s", "nan");
            printf("\n");
        }
        printf("\n");

        /* ── min_per_gate: global min of per-target mins ───────────────── */
        printf("# %s min_per_gate"
               " (µs/gate, minimum over all target positions)\n", gname);
        printf("%-8s", "qubits");
        for (int m = 0; m < num_methods; ++m) printf("  %-14s", method_names[m]);
        printf("\n");
        for (int q = 0; q < pgf->num_configs; ++q) {
            printf("%-8u", (unsigned)pgf->qubits[q]);
            for (int m = 0; m < num_methods; ++m) {
                int mi = method_indices[m];
                int nt = pgf->n_targets[g][q][mi];
                double mn = -1.0;
                for (int t = 0; t < nt; ++t) {
                    double v = pgf->cells[g][q][mi][t].time_per_gate_us_min;
                    if (v > 0.0 && (mn < 0.0 || v < mn)) mn = v;
                }
                if (mn > 0.0) printf("  %-14.6f", mn);
                else          printf("  %-14s", "nan");
            }
            printf("\n");
        }
        printf("\n");

        /* ── median_per_gate: mean of per-target medians ───────────────── */
        printf("# %s median_per_gate"
               " (µs/gate, mean of per-target medians)\n", gname);
        printf("%-8s", "qubits");
        for (int m = 0; m < num_methods; ++m) printf("  %-14s", method_names[m]);
        printf("\n");
        for (int q = 0; q < pgf->num_configs; ++q) {
            printf("%-8u", (unsigned)pgf->qubits[q]);
            for (int m = 0; m < num_methods; ++m) {
                int mi = method_indices[m];
                int nt = pgf->n_targets[g][q][mi];
                double sum = 0.0; int valid = 0;
                for (int t = 0; t < nt; ++t) {
                    double v = pgf->cells[g][q][mi][t].time_per_gate_us_median;
                    if (v > 0.0) { sum += v; valid++; }
                }
                if (valid > 0) printf("  %-14.6f", sum / valid);
                else           printf("  %-14s", "nan");
            }
            printf("\n");
        }
        printf("\n");

        /* ── cv: mean of per-target CVs ────────────────────────────────── */
        /* CV is scale-invariant: CV(total_time) == CV(per_gate_time).     */
        printf("# %s cv (coefficient of variation, %%, mean of per-target CVs)\n", gname);
        printf("%-8s", "qubits");
        for (int m = 0; m < num_methods; ++m) printf("  %-14s", method_names[m]);
        printf("\n");
        for (int q = 0; q < pgf->num_configs; ++q) {
            printf("%-8u", (unsigned)pgf->qubits[q]);
            for (int m = 0; m < num_methods; ++m) {
                int mi = method_indices[m];
                int nt = pgf->n_targets[g][q][mi];
                double sum = 0.0; int valid = 0;
                for (int t = 0; t < nt; ++t) {
                    if (pgf->cells[g][q][mi][t].time_per_gate_us > 0.0) {
                        sum += pgf->cells[g][q][mi][t].time_ms_cv;
                        valid++;
                    }
                }
                if (valid > 0) printf("  %-14.6f", sum / valid);
                else           printf("  %-14s", "nan");
            }
            printf("\n");
        }
        printf("\n");

    } /* gate loop */

    /* ── Memory table ────────────────────────────────────────────────────── */
    /* Memory is constant across targets for a given (method, qubit count).  */
    /* Read from the first populated target for gate 0 (X gate).             */
    printf("# Memory usage (MB)\n");
    printf("%-8s", "qubits");
    for (int m = 0; m < num_methods; ++m) printf("  %-14s", method_names[m]);
    printf("\n");
    for (int q = 0; q < pgf->num_configs; ++q) {
        printf("%-8u", (unsigned)pgf->qubits[q]);
        for (int m = 0; m < num_methods; ++m) {
            int mi = method_indices[m];
            int nt = pgf->n_targets[0][q][mi];
            size_t bytes = 0;
            for (int t = 0; t < nt; ++t) {
                if (pgf->cells[0][q][mi][t].memory_bytes > 0) {
                    bytes = pgf->cells[0][q][mi][t].memory_bytes;
                    break;
                }
            }
            if (bytes > 0) printf("  %-14.3f", (double)bytes / (1024.0 * 1024.0));
            else           printf("  %-14s", "nan");
        }
        printf("\n");
    }
}

/*
 * =====================================================================================================================
 * Per-qubit benchmark
 * =====================================================================================================================
 */

/**
 * @brief Gate definition for 2Q gates in per-qubit mode.
 *
 * Separate from gate_def_2q_t to avoid modifying the existing sweep table.
 * is_symmetric: 1 = SWAP (unordered pairs via build_all_pairs);
 *               0 = CX   (ordered pairs via build_all_ordered_pairs).
 */
typedef struct {
    const char *name;
    void (*fn)(state_t*, qubit_t, qubit_t);
    const cplx_t *mat;
    int has_tiled;
    int is_symmetric;
} gate_def_2q_perq_t;

/**
 * @brief Gate definition for 1Q gates in per-qubit mode.
 *
 * Adds has_tiled flag; all current 1Q gates have a tiled implementation.
 */
typedef struct {
    const char *name;
    void (*fn)(state_t*, qubit_t);
    const cplx_t *mat;
    int has_tiled;
} gate_def_1q_perq_t;

/**
 * @brief Calibrate the number of iterations for one method/gate/qubit combination.
 *
 * Uses the minimum run time from the probe result to estimate how many iterations
 * achieve MIN_MEASUREMENT_MS per run.  The divisor is probe_result->iterations
 * (the actual number of iterations the probe ran) rather than the nominal
 * BENCH_PQ_PROBE_ITERS constant, so the calculation is correct when the adaptive
 * probe cap reduces the probe to 1 iteration for slow gates.
 * Falls back to BENCH_PQ_FALLBACK_ITERATIONS if the probe rounds to zero.
 * Caller may override by setting opts->iterations_explicit.
 *
 * @param probe_result  Result of a probe run.
 * @param iters_explicit If non-zero, ignore calibration and use explicit value.
 * @param explicit_iters Explicit iteration count when iters_explicit != 0.
 * @return Calibrated iteration count.
 */
static int calibrate_iterations(const bench_result_perq_t *probe_result,
                                 int iters_explicit, int explicit_iters) {
    if (iters_explicit) return explicit_iters;

    double probe_per_iter_ms = probe_result->time_ms_min / (double)probe_result->iterations;
    int calibrated;
    if (probe_per_iter_ms <= 0.0) {
        calibrated = BENCH_PQ_FALLBACK_ITERATIONS;
    } else {
        calibrated = (int)ceil(MIN_MEASUREMENT_MS / probe_per_iter_ms);
        if (calibrated < 1) calibrated = 1;
        if (calibrated > BENCH_PQ_MAX_CALIBRATED_ITERATIONS)
            calibrated = BENCH_PQ_MAX_CALIBRATED_ITERATIONS;
    }
    return calibrated;
}

/**
 * @brief Per-qubit benchmark: time each gate at every target qubit/pair, for every method.
 *
 * Iterates over qubits from opts->min_qubits to opts->max_qubits, runs all
 * 1Q and 2Q gates, and records per-target timing for each enabled method.
 * Calibrates iteration count per (gate, qubits, method) combination.
 */
static int perq_main(const bench_options_t *opts) {
    /* Method index assignments — must match the M_* enum in the sweep benchmark */
    const int m_packed = M_QLIB;
    const int m_tiled  = M_QLIB_TILED;
    const int m_blas   = M_BLAS;
#ifdef WITH_QUEST
    const int m_quest  = M_QUEST;
#endif
#ifdef WITH_QULACS
    const int m_qulacs = M_QULACS;
#endif
#ifdef WITH_AER_DM
    const int m_aer    = M_AER_DM;
#endif

    /* Result storage: [NUM_METHODS][MAX_TARGETS] */
    bench_result_perq_t (*perq)[MAX_TARGETS] = calloc(NUM_METHODS, sizeof(*perq));
    if (!perq) {
        fprintf(stderr, "error: calloc failed\n");
        return EXIT_FAILURE;
    }

    /* CSV output goes to stdout, same as the sweep benchmark */
    FILE *csv_f = stdout;

    if (opts->csv_output) {
        bench_print_perq_csv_header(csv_f);
    } else {
        printf("\nPer-Qubit Gate Benchmarks\n");
        printf("=========================\n");
        printf("Iterations: calibrated, Runs: %d\n", opts->runs);
#ifdef WITH_QUEST
        printf("QuEST: enabled\n");
#endif
#ifdef WITH_QULACS
        printf("Qulacs: enabled\n");
#endif
#ifdef WITH_AER_DM
        printf("Qiskit-Aer-DM: enabled\n");
#endif
        printf("\n");
    }

#ifdef WITH_QUEST
    bench_quest_init();
#endif

    int results_count = 0;

    /* 1Q gate table for the per-qubit path */
    static const gate_def_1q_perq_t gates_1q_perq[] = {
        {"X", x, XMAT, 1},
        {"H", h, HMAT, 1},
        {"Z", z, ZMAT, 1},
        {NULL, NULL, NULL, 0}
    };

    /* 2Q gate table for the per-qubit path */
    static const gate_def_2q_perq_t gates_2q_perq[] = {
        {"CX",   cx,        CXMAT,   1, 0},   /* ordered pairs (CX is not symmetric) */
        {"SWAP", swap_gate, SWAPMAT, 1, 1},   /* unordered pairs (SWAP is symmetric) */
        {NULL, NULL, NULL, 0, 0}
    };

    for (qubit_t qubits = opts->min_qubits; qubits <= opts->max_qubits;
         qubits += (qubit_t)opts->step) {

        dim_t  dim       = (dim_t)1 << qubits;
        int    run_dense = (qubits <= 8);
        size_t state_bytes = (size_t)dim * (size_t)dim * sizeof(cplx_t);

        /* ── 1Q GATES ──────────────────────────────────────────────────────── */

        for (int g = 0; gates_1q_perq[g].name != NULL; ++g) {
            const gate_def_1q_perq_t *gd = &gates_1q_perq[g];

            if (opts->gate_filter && strcasecmp(opts->gate_filter, gd->name) != 0)
                continue;

            results_count++;
            int n_targets = (int)qubits;
            memset(perq, 0, NUM_METHODS * sizeof(*perq));

            /* --- qlib packed -------------------------------------------- */
            {
                state_t state = {0};
                state.type = MIXED_PACKED;
                state_init(&state, qubits, NULL);
                if (state.data) {
                    bench_fill_random(state.data,
                        (size_t)state_len(&state) * sizeof(cplx_t));

                    /* Adaptive probe: single-shot first to detect slow gates */
                    bench_result_perq_t quick =
                        bench_qlib_packed_at(qubits, gd->name, gd->fn, 0,
                                             1, 1, &state);
                    int probe_iters = (quick.time_per_gate_us / 1000.0 >= MIN_MEASUREMENT_MS)
                                      ? 1 : BENCH_PQ_PROBE_ITERS;
                    bench_result_perq_t probe =
                        bench_qlib_packed_at(qubits, gd->name, gd->fn, 0,
                                             probe_iters,
                                             BENCH_PQ_PROBE_RUNS, &state);
                    int cal_iters = calibrate_iterations(&probe,
                                        opts->iterations_explicit,
                                        opts->iterations);

                    for (qubit_t tgt = 0; tgt < qubits; tgt++) {
                        perq[m_packed][tgt] = bench_qlib_packed_at(
                            qubits, gd->name, gd->fn, tgt,
                            cal_iters, opts->runs, &state);
                    }
                    state_free(&state);
                }
            }

            /* --- qlib tiled --------------------------------------------- */
            if (gd->has_tiled) {
                state_t state = {0};
                state.type = MIXED_TILED;
                state_init(&state, qubits, NULL);
                if (state.data) {
                    bench_fill_random(state.data,
                        (size_t)state_len(&state) * sizeof(cplx_t));

                    bench_result_perq_t quick =
                        bench_qlib_tiled_at(qubits, gd->name, gd->fn, 0,
                                            1, 1, &state);
                    int probe_iters = (quick.time_per_gate_us / 1000.0 >= MIN_MEASUREMENT_MS)
                                      ? 1 : BENCH_PQ_PROBE_ITERS;
                    bench_result_perq_t probe =
                        bench_qlib_tiled_at(qubits, gd->name, gd->fn, 0,
                                            probe_iters,
                                            BENCH_PQ_PROBE_RUNS, &state);
                    int cal_iters = calibrate_iterations(&probe,
                                        opts->iterations_explicit,
                                        opts->iterations);

                    for (qubit_t tgt = 0; tgt < qubits; tgt++) {
                        perq[m_tiled][tgt] = bench_qlib_tiled_at(
                            qubits, gd->name, gd->fn, tgt,
                            cal_iters, opts->runs, &state);
                    }
                    state_free(&state);
                }
            }

            /* --- BLAS dense (capped at 8 qubits) ------------------------ */
            if (run_dense) {
                size_t rho_bytes = state_bytes;
                cplx_t *rho = bench_alloc_huge(rho_bytes);
                if (rho) {
                    bench_fill_random(rho, rho_bytes);

                    /* Adaptive probe: single-shot first to detect slow gates */
                    bench_result_perq_t quick =
                        bench_blas_dense_at(qubits, gd->name, gd->mat, 0,
                                            1, 1, rho);
                    int probe_iters = (quick.time_per_gate_us / 1000.0 >= MIN_MEASUREMENT_MS)
                                      ? 1 : BENCH_PQ_PROBE_ITERS;
                    bench_result_perq_t probe =
                        bench_blas_dense_at(qubits, gd->name, gd->mat, 0,
                                            probe_iters,
                                            BENCH_PQ_PROBE_RUNS, rho);
                    int cal_iters = calibrate_iterations(&probe,
                                        opts->iterations_explicit,
                                        opts->iterations);

                    for (qubit_t tgt = 0; tgt < qubits; tgt++) {
                        perq[m_blas][tgt] = bench_blas_dense_at(
                            qubits, gd->name, gd->mat, tgt,
                            cal_iters, opts->runs, rho);
                    }
                    bench_free_huge(rho, rho_bytes);
                }
            }

#if defined(WITH_QUEST) && defined(QUEST_H)
            /* --- QuEST -------------------------------------------------- */
            {
                Qureg qureg = createDensityQureg((int)qubits);
                initDebugState(qureg);

                bench_result_perq_t quick =
                    bench_quest_at(qubits, gd->name, 0, 1, 1, &qureg);
                int probe_iters = (quick.time_per_gate_us / 1000.0 >= MIN_MEASUREMENT_MS)
                                  ? 1 : BENCH_PQ_PROBE_ITERS;
                bench_result_perq_t probe =
                    bench_quest_at(qubits, gd->name, 0,
                                   probe_iters,
                                   BENCH_PQ_PROBE_RUNS, &qureg);
                int cal_iters = calibrate_iterations(&probe,
                                    opts->iterations_explicit,
                                    opts->iterations);

                for (qubit_t tgt = 0; tgt < qubits; tgt++) {
                    perq[m_quest][tgt] = bench_quest_at(
                        qubits, gd->name, tgt,
                        cal_iters, opts->runs, &qureg);
                }
                destroyQureg(qureg);
            }
#endif /* WITH_QUEST && QUEST_H */

#ifdef WITH_QULACS
            /* --- Qulacs ------------------------------------------------- */
            {
                size_t rho_bytes = state_bytes;
                void *rho = bench_alloc_huge(rho_bytes);
                if (rho) {
                    bench_fill_random(rho, rho_bytes);

                    bench_result_perq_t quick =
                        bench_qulacs_at(qubits, gd->name, 0, 1, 1, rho);
                    int probe_iters = (quick.time_per_gate_us / 1000.0 >= MIN_MEASUREMENT_MS)
                                      ? 1 : BENCH_PQ_PROBE_ITERS;
                    bench_result_perq_t probe =
                        bench_qulacs_at(qubits, gd->name, 0,
                                        probe_iters,
                                        BENCH_PQ_PROBE_RUNS, rho);
                    int cal_iters = calibrate_iterations(&probe,
                                        opts->iterations_explicit,
                                        opts->iterations);

                    for (qubit_t tgt = 0; tgt < qubits; tgt++) {
                        perq[m_qulacs][tgt] = bench_qulacs_at(
                            qubits, gd->name, tgt,
                            cal_iters, opts->runs, rho);
                    }
                    bench_free_huge(rho, rho_bytes);
                }
            }
#endif /* WITH_QULACS */

#ifdef WITH_AER_DM
            /* --- Qiskit-Aer DM ------------------------------------------ */
            {
                void *dm = bench_aer_dm_alloc(qubits);
                if (dm) {
                    bench_aer_dm_reinit(dm, qubits);

                    bench_result_perq_t quick =
                        bench_aer_dm_at(qubits, gd->name, 0, 1, 1, dm);
                    int probe_iters = (quick.time_per_gate_us / 1000.0 >= MIN_MEASUREMENT_MS)
                                      ? 1 : BENCH_PQ_PROBE_ITERS;
                    bench_result_perq_t probe =
                        bench_aer_dm_at(qubits, gd->name, 0,
                                        probe_iters,
                                        BENCH_PQ_PROBE_RUNS, dm);
                    int cal_iters = calibrate_iterations(&probe,
                                        opts->iterations_explicit,
                                        opts->iterations);

                    for (qubit_t tgt = 0; tgt < qubits; tgt++) {
                        perq[m_aer][tgt] = bench_aer_dm_at(
                            qubits, gd->name, tgt,
                            cal_iters, opts->runs, dm);
                    }
                    bench_aer_dm_free(dm);
                }
            }
#endif /* WITH_AER_DM */

            /* Output results for this gate × qubit count */
            if (opts->csv_output) {
                for (int m = 0; m < NUM_METHODS; m++)
                    bench_print_perq_csv(csv_f, perq[m], n_targets);
            } else {
                for (int m = 0; m < NUM_METHODS; m++) {
                    /* Skip methods that produced no results */
                    int has_results = 0;
                    for (int t = 0; t < n_targets; t++) {
                        if (perq[m][t].time_ms_mean > 0.0) { has_results = 1; break; }
                    }
                    if (!has_results) continue;
                    bench_print_perq_console(perq[m], n_targets);
                }
            }
        } /* 1Q gate loop */

        /* ── 2Q GATES ──────────────────────────────────────────────────────── */

        for (int g = 0; gates_2q_perq[g].name != NULL; ++g) {
            const gate_def_2q_perq_t *gd = &gates_2q_perq[g];

            if (opts->gate_filter && strcasecmp(opts->gate_filter, gd->name) != 0)
                continue;

            /* Build the pair list */
            qubit_t q1s[MAX_PAIRS], q2s[MAX_PAIRS];
            int n_pairs;
            if (gd->is_symmetric) {
                /* SWAP: unordered pairs (q1 < q2) */
                n_pairs = build_all_pairs(qubits, q1s, q2s);
            } else {
                /* CX: ordered pairs (q1 != q2, all permutations) */
                n_pairs = build_all_ordered_pairs(qubits, q1s, q2s, MAX_PAIRS);
            }

            if (n_pairs <= 0) {
                fprintf(stderr, "warning: %s on %u qubits: no pairs, skipping\n",
                        gd->name, (unsigned)qubits);
                continue;
            }
            if (n_pairs > MAX_TARGETS) {
                fprintf(stderr, "warning: %s on %u qubits: %d pairs exceed MAX_TARGETS=%d, skipping\n",
                        gd->name, (unsigned)qubits, n_pairs, MAX_TARGETS);
                continue;
            }

            results_count++;
            memset(perq, 0, NUM_METHODS * sizeof(*perq));

            /* --- qlib packed -------------------------------------------- */
            {
                state_t state = {0};
                state.type = MIXED_PACKED;
                state_init(&state, qubits, NULL);
                if (state.data) {
                    bench_fill_random(state.data,
                        (size_t)state_len(&state) * sizeof(cplx_t));

                    bench_result_perq_t quick =
                        bench_qlib_packed_2q_at(qubits, gd->name, gd->fn,
                                                q1s[0], q2s[0],
                                                1, 1, &state);
                    int probe_iters = (quick.time_per_gate_us / 1000.0 >= MIN_MEASUREMENT_MS)
                                      ? 1 : BENCH_PQ_PROBE_ITERS;
                    bench_result_perq_t probe =
                        bench_qlib_packed_2q_at(qubits, gd->name, gd->fn,
                                                q1s[0], q2s[0],
                                                probe_iters,
                                                BENCH_PQ_PROBE_RUNS, &state);
                    int cal_iters = calibrate_iterations(&probe,
                                        opts->iterations_explicit,
                                        opts->iterations);

                    for (int p = 0; p < n_pairs; p++) {
                        perq[m_packed][p] = bench_qlib_packed_2q_at(
                            qubits, gd->name, gd->fn,
                            q1s[p], q2s[p],
                            cal_iters, opts->runs, &state);
                    }
                    state_free(&state);
                }
            }

            /* --- qlib tiled --------------------------------------------- */
            if (gd->has_tiled) {
                state_t state = {0};
                state.type = MIXED_TILED;
                state_init(&state, qubits, NULL);
                if (state.data) {
                    bench_fill_random(state.data,
                        (size_t)state_len(&state) * sizeof(cplx_t));

                    bench_result_perq_t quick =
                        bench_qlib_tiled_2q_at(qubits, gd->name, gd->fn,
                                               q1s[0], q2s[0],
                                               1, 1, &state);
                    int probe_iters = (quick.time_per_gate_us / 1000.0 >= MIN_MEASUREMENT_MS)
                                      ? 1 : BENCH_PQ_PROBE_ITERS;
                    bench_result_perq_t probe =
                        bench_qlib_tiled_2q_at(qubits, gd->name, gd->fn,
                                               q1s[0], q2s[0],
                                               probe_iters,
                                               BENCH_PQ_PROBE_RUNS, &state);
                    int cal_iters = calibrate_iterations(&probe,
                                        opts->iterations_explicit,
                                        opts->iterations);

                    for (int p = 0; p < n_pairs; p++) {
                        perq[m_tiled][p] = bench_qlib_tiled_2q_at(
                            qubits, gd->name, gd->fn,
                            q1s[p], q2s[p],
                            cal_iters, opts->runs, &state);
                    }
                    state_free(&state);
                }
            }

            /* --- BLAS dense (capped at 8 qubits) ------------------------ */
            if (run_dense) {
                size_t rho_bytes = state_bytes;
                cplx_t *rho = bench_alloc_huge(rho_bytes);
                if (rho) {
                    bench_fill_random(rho, rho_bytes);

                    bench_result_perq_t quick =
                        bench_blas_dense_2q_at(qubits, gd->name, gd->mat,
                                               q1s[0], q2s[0],
                                               1, 1, rho);
                    int probe_iters = (quick.time_per_gate_us / 1000.0 >= MIN_MEASUREMENT_MS)
                                      ? 1 : BENCH_PQ_PROBE_ITERS;
                    bench_result_perq_t probe =
                        bench_blas_dense_2q_at(qubits, gd->name, gd->mat,
                                               q1s[0], q2s[0],
                                               probe_iters,
                                               BENCH_PQ_PROBE_RUNS, rho);
                    int cal_iters = calibrate_iterations(&probe,
                                        opts->iterations_explicit,
                                        opts->iterations);

                    for (int p = 0; p < n_pairs; p++) {
                        perq[m_blas][p] = bench_blas_dense_2q_at(
                            qubits, gd->name, gd->mat,
                            q1s[p], q2s[p],
                            cal_iters, opts->runs, rho);
                    }
                    bench_free_huge(rho, rho_bytes);
                }
            }

#if defined(WITH_QUEST) && defined(QUEST_H)
            /* --- QuEST -------------------------------------------------- */
            {
                Qureg qureg = createDensityQureg((int)qubits);
                initDebugState(qureg);

                bench_result_perq_t quick =
                    bench_quest_2q_at(qubits, gd->name,
                                      q1s[0], q2s[0],
                                      1, 1, &qureg);
                int probe_iters = (quick.time_per_gate_us / 1000.0 >= MIN_MEASUREMENT_MS)
                                  ? 1 : BENCH_PQ_PROBE_ITERS;
                bench_result_perq_t probe =
                    bench_quest_2q_at(qubits, gd->name,
                                      q1s[0], q2s[0],
                                      probe_iters,
                                      BENCH_PQ_PROBE_RUNS, &qureg);
                int cal_iters = calibrate_iterations(&probe,
                                    opts->iterations_explicit,
                                    opts->iterations);

                for (int p = 0; p < n_pairs; p++) {
                    perq[m_quest][p] = bench_quest_2q_at(
                        qubits, gd->name,
                        q1s[p], q2s[p],
                        cal_iters, opts->runs, &qureg);
                }
                destroyQureg(qureg);
            }
#endif /* WITH_QUEST && QUEST_H */

#ifdef WITH_QULACS
            /* --- Qulacs ------------------------------------------------- */
            {
                size_t rho_bytes = state_bytes;
                void *rho = bench_alloc_huge(rho_bytes);
                if (rho) {
                    bench_fill_random(rho, rho_bytes);

                    bench_result_perq_t quick =
                        bench_qulacs_2q_at(qubits, gd->name,
                                           q1s[0], q2s[0],
                                           1, 1, rho);
                    int probe_iters = (quick.time_per_gate_us / 1000.0 >= MIN_MEASUREMENT_MS)
                                      ? 1 : BENCH_PQ_PROBE_ITERS;
                    bench_result_perq_t probe =
                        bench_qulacs_2q_at(qubits, gd->name,
                                           q1s[0], q2s[0],
                                           probe_iters,
                                           BENCH_PQ_PROBE_RUNS, rho);
                    int cal_iters = calibrate_iterations(&probe,
                                        opts->iterations_explicit,
                                        opts->iterations);

                    for (int p = 0; p < n_pairs; p++) {
                        perq[m_qulacs][p] = bench_qulacs_2q_at(
                            qubits, gd->name,
                            q1s[p], q2s[p],
                            cal_iters, opts->runs, rho);
                    }
                    bench_free_huge(rho, rho_bytes);
                }
            }
#endif /* WITH_QULACS */

#ifdef WITH_AER_DM
            /* --- Qiskit-Aer DM ------------------------------------------ */
            {
                void *dm = bench_aer_dm_alloc(qubits);
                if (dm) {
                    bench_aer_dm_reinit(dm, qubits);

                    bench_result_perq_t quick =
                        bench_aer_dm_2q_at(qubits, gd->name,
                                           q1s[0], q2s[0],
                                           1, 1, dm);
                    int probe_iters = (quick.time_per_gate_us / 1000.0 >= MIN_MEASUREMENT_MS)
                                      ? 1 : BENCH_PQ_PROBE_ITERS;
                    bench_result_perq_t probe =
                        bench_aer_dm_2q_at(qubits, gd->name,
                                           q1s[0], q2s[0],
                                           probe_iters,
                                           BENCH_PQ_PROBE_RUNS, dm);
                    int cal_iters = calibrate_iterations(&probe,
                                        opts->iterations_explicit,
                                        opts->iterations);

                    for (int p = 0; p < n_pairs; p++) {
                        perq[m_aer][p] = bench_aer_dm_2q_at(
                            qubits, gd->name,
                            q1s[p], q2s[p],
                            cal_iters, opts->runs, dm);
                    }
                    bench_aer_dm_free(dm);
                }
            }
#endif /* WITH_AER_DM */

            /* Output results for this gate × qubit count */
            if (opts->csv_output) {
                for (int m = 0; m < NUM_METHODS; m++)
                    bench_print_perq_csv(csv_f, perq[m], n_pairs);
            } else {
                for (int m = 0; m < NUM_METHODS; m++) {
                    /* Skip methods that produced no results */
                    int has_results = 0;
                    for (int t = 0; t < n_pairs; t++) {
                        if (perq[m][t].time_ms_mean > 0.0) { has_results = 1; break; }
                    }
                    if (!has_results) continue;
                    bench_print_perq_console(perq[m], n_pairs);
                }
            }
        } /* 2Q gate loop */

    } /* qubits loop */

#ifdef WITH_QUEST
    bench_quest_cleanup();
#endif

    if (results_count == 0 && opts->gate_filter) {
        fprintf(stderr, "warning: no results — gate '%s' not found\n",
                opts->gate_filter);
        free(perq);
        return EXIT_FAILURE;
    }

    free(perq);
    return 0;
}

/*
 * =====================================================================================================================
 * Main
 * =====================================================================================================================
 */

int main(int argc, char *argv[]) {
    bench_options_t opts = bench_parse_options(argc, argv);

    if (opts.per_qubit) {
        return perq_main(&opts);
    }

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
#ifndef BENCH_AER_ONLY_MODE
        printf("Comparing: qlib packed sparse vs BLAS dense vs naive loop");
#ifdef WITH_QUEST
        printf(" vs QuEST");
#endif
#ifdef WITH_QULACS
        printf(" vs Qulacs");
#endif
#ifdef WITH_AER_DM
        printf(" vs Qiskit-Aer-DM");
#endif
        printf("\n");
#else /* BENCH_AER_ONLY_MODE */
#ifdef WITH_AER_DM
        printf("Comparing: Qiskit-Aer-DM only (bench_aer_only mode)\n");
#endif
#endif /* !BENCH_AER_ONLY_MODE */
        printf("Iterations: %d, Warm-up: %d, Runs: %d\n",
               opts.iterations, opts.warmup, opts.runs);
#ifndef BENCH_AER_ONLY_MODE
#ifdef WITH_QUEST
        printf("QuEST: enabled\n");
#endif
#ifdef WITH_QULACS
        printf("Qulacs: enabled\n");
#endif
#endif
#ifdef WITH_AER_DM
        printf("Qiskit-Aer-DM: enabled\n");
#endif
        printf("\n");
    }

    int qi = 0;

#ifdef WITH_QUEST
    bench_quest_init();
#endif

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

#ifndef BENCH_AER_ONLY_MODE
            s.packed = bench_qlib_packed(qubits, gates[g].name, gates[g].fn,
                                         opts.iterations, opts.warmup, opts.runs);
            s.tiled  = bench_qlib_tiled(qubits, gates[g].name, gates[g].fn,
                                        opts.iterations, opts.warmup, opts.runs);
            s.has_tiled = (s.tiled.time_ms > 0);

            if (run_dense) {
                s.dense = bench_blas_dense(qubits, gates[g].name, gates[g].mat,
                                           opts.iterations, opts.warmup, opts.runs);
                s.has_dense = 1;

                if (!opts.pgfplots_output) {
                    s.naive = bench_naive_loop(qubits, gates[g].name, gates[g].mat,
                                               opts.iterations / 10,
                                               opts.warmup    / 10,
                                               opts.runs);
                    s.has_naive = 1;
                }
            }

#ifdef WITH_QUEST
            s.quest     = bench_quest(qubits, gates[g].name, opts.iterations, opts.warmup, opts.runs);
            s.has_quest = (s.quest.time_ms > 0);
#endif
#ifdef WITH_QULACS
            s.qulacs     = bench_qulacs(qubits, gates[g].name, opts.iterations, opts.warmup, opts.runs);
            s.has_qulacs = (s.qulacs.time_ms > 0);
#endif
#endif /* !BENCH_AER_ONLY_MODE */

#ifdef WITH_AER_DM
            s.aer_dm     = bench_aer_dm(qubits, gates[g].name, opts.iterations, opts.warmup, opts.runs);
            s.has_aer_dm = (s.aer_dm.time_ms > 0);
#endif

            if (pgf && qi < MAX_QUBIT_CONFIGS) {
                pgf->runs[qi] = opts.runs;
#ifndef BENCH_AER_ONLY_MODE
                pgf->time_ms[g][qi][M_QLIB]           = s.packed.time_ms;
                pgf->time_ms[g][qi][M_QLIB_TILED]     = s.tiled.time_ms;
                pgf->time_ms[g][qi][M_BLAS]            = s.dense.time_ms;
                pgf->time_ms_std[g][qi][M_QLIB]        = s.packed.time_ms_std;
                pgf->time_ms_std[g][qi][M_QLIB_TILED]  = s.tiled.time_ms_std;
                pgf->time_ms_std[g][qi][M_BLAS]        = s.dense.time_ms_std;
                pgf->time_ms_min[g][qi][M_QLIB]        = s.packed.time_ms_min;
                pgf->time_ms_min[g][qi][M_QLIB_TILED]  = s.tiled.time_ms_min;
                pgf->time_ms_min[g][qi][M_BLAS]        = s.dense.time_ms_min;
                pgf->time_ms_median[g][qi][M_QLIB]        = s.packed.time_ms_median;
                pgf->time_ms_median[g][qi][M_QLIB_TILED]  = s.tiled.time_ms_median;
                pgf->time_ms_median[g][qi][M_BLAS]        = s.dense.time_ms_median;
                pgf->time_ms_cv[g][qi][M_QLIB]         = s.packed.time_ms_cv;
                pgf->time_ms_cv[g][qi][M_QLIB_TILED]   = s.tiled.time_ms_cv;
                pgf->time_ms_cv[g][qi][M_BLAS]         = s.dense.time_ms_cv;
                pgf->memory[g][qi][M_QLIB]             = s.packed.memory_bytes;
                pgf->memory[g][qi][M_QLIB_TILED]       = s.tiled.memory_bytes;
                pgf->memory[g][qi][M_BLAS]             = s.dense.memory_bytes;
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
#endif /* !BENCH_AER_ONLY_MODE */
#ifdef WITH_AER_DM
                pgf->time_ms[g][qi][M_AER_DM]          = s.aer_dm.time_ms;
                pgf->time_ms_std[g][qi][M_AER_DM]      = s.aer_dm.time_ms_std;
                pgf->time_ms_min[g][qi][M_AER_DM]      = s.aer_dm.time_ms_min;
                pgf->time_ms_median[g][qi][M_AER_DM]   = s.aer_dm.time_ms_median;
                pgf->time_ms_cv[g][qi][M_AER_DM]       = s.aer_dm.time_ms_cv;
                pgf->memory[g][qi][M_AER_DM]           = s.aer_dm.memory_bytes;
#endif
            }

            if (opts.csv_output) {
                int sw1q = (int)qubits;
#ifndef BENCH_AER_ONLY_MODE
                SET_SWEEP(s.packed, sw1q); bench_print_csv(&s.packed);
                SET_SWEEP(s.tiled,  sw1q); bench_print_csv(&s.tiled);
                if (s.has_dense) { SET_SWEEP(s.dense, sw1q); bench_print_csv(&s.dense); }
                if (s.has_naive) { SET_SWEEP(s.naive, sw1q); bench_print_csv(&s.naive); }
#ifdef WITH_QUEST
                if (s.has_quest) { SET_SWEEP(s.quest, sw1q); bench_print_csv(&s.quest); }
#endif
#ifdef WITH_QULACS
                if (s.has_qulacs) { SET_SWEEP(s.qulacs, sw1q); bench_print_csv(&s.qulacs); }
#endif
#endif /* !BENCH_AER_ONLY_MODE */
#ifdef WITH_AER_DM
                if (s.has_aer_dm) { SET_SWEEP(s.aer_dm, sw1q); bench_print_csv(&s.aer_dm); }
#endif
            } else if (!opts.pgfplots_output) {
#ifndef BENCH_AER_ONLY_MODE
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
#endif /* !BENCH_AER_ONLY_MODE */
#ifdef WITH_AER_DM
                if (s.has_aer_dm) bench_print_result(&s.aer_dm, opts.verbose);
#endif
#ifndef BENCH_AER_ONLY_MODE
                print_speedup_and_memory(&s);
#endif
            }
        }

        /* ── Two-qubit gates ─────────────────────────────────────────────────── */
        for (int g = 0; gates_2q[g].name != NULL; ++g) {
            if (!opts.csv_output && !opts.pgfplots_output)
                printf("\nGate: %s\n", gates_2q[g].name);

            bench_summary_t s = {0};
            s.naive_scale = 1.0;

#ifndef BENCH_AER_ONLY_MODE
            s.packed = bench_qlib_packed_2q(qubits, gates_2q[g].name, gates_2q[g].fn,
                                            opts.iterations, opts.warmup, opts.runs);

            if (gates_2q[g].has_tiled) {
                s.tiled = bench_qlib_tiled_2q(qubits, gates_2q[g].name, gates_2q[g].fn,
                                              opts.iterations, opts.warmup, opts.runs);
                s.has_tiled = (s.tiled.time_ms > 0);
            }

            if (run_dense) {
                s.dense = bench_blas_dense_2q(qubits, gates_2q[g].name, gates_2q[g].mat,
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
#endif /* !BENCH_AER_ONLY_MODE */

#ifdef WITH_AER_DM
            s.aer_dm     = bench_aer_dm(qubits, gates_2q[g].name, opts.iterations, opts.warmup, opts.runs);
            s.has_aer_dm = (s.aer_dm.time_ms > 0);
#endif

            int pg = NUM_1Q_GATES + g;
            if (pgf && qi < MAX_QUBIT_CONFIGS) {
                /* pgf->runs[qi] already written in the 1Q loop for this qubit config */
#ifndef BENCH_AER_ONLY_MODE
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
#endif /* !BENCH_AER_ONLY_MODE */
#ifdef WITH_AER_DM
                pgf->time_ms[pg][qi][M_AER_DM]          = s.aer_dm.time_ms;
                pgf->time_ms_std[pg][qi][M_AER_DM]      = s.aer_dm.time_ms_std;
                pgf->time_ms_min[pg][qi][M_AER_DM]      = s.aer_dm.time_ms_min;
                pgf->time_ms_median[pg][qi][M_AER_DM]   = s.aer_dm.time_ms_median;
                pgf->time_ms_cv[pg][qi][M_AER_DM]       = s.aer_dm.time_ms_cv;
                pgf->memory[pg][qi][M_AER_DM]           = s.aer_dm.memory_bytes;
#endif
            }

            if (opts.csv_output) {
                int sw2q = (int)(qubits * (qubits - 1) / 2);
#ifndef BENCH_AER_ONLY_MODE
                SET_SWEEP(s.packed, sw2q); bench_print_csv(&s.packed);
                if (s.has_tiled) { SET_SWEEP(s.tiled, sw2q); bench_print_csv(&s.tiled); }
                if (s.has_dense) { SET_SWEEP(s.dense, sw2q); bench_print_csv(&s.dense); }
#ifdef WITH_QUEST
                if (s.has_quest) { SET_SWEEP(s.quest, sw2q); bench_print_csv(&s.quest); }
#endif
#ifdef WITH_QULACS
                if (s.has_qulacs) { SET_SWEEP(s.qulacs, sw2q); bench_print_csv(&s.qulacs); }
#endif
#endif /* !BENCH_AER_ONLY_MODE */
#ifdef WITH_AER_DM
                if (s.has_aer_dm) { SET_SWEEP(s.aer_dm, sw2q); bench_print_csv(&s.aer_dm); }
#endif
            } else if (!opts.pgfplots_output) {
#ifndef BENCH_AER_ONLY_MODE
                bench_print_result(&s.packed, opts.verbose);
                if (s.has_tiled) bench_print_result(&s.tiled,  opts.verbose);
                if (s.has_dense) bench_print_result(&s.dense,  opts.verbose);
#ifdef WITH_QUEST
                if (s.has_quest) bench_print_result(&s.quest, opts.verbose);
#endif
#ifdef WITH_QULACS
                if (s.has_qulacs) bench_print_result(&s.qulacs, opts.verbose);
#endif
#endif /* !BENCH_AER_ONLY_MODE */
#ifdef WITH_AER_DM
                if (s.has_aer_dm) bench_print_result(&s.aer_dm, opts.verbose);
#endif
#ifndef BENCH_AER_ONLY_MODE
                print_speedup_and_memory(&s);
#endif
            }
        }

        if (!opts.csv_output && !opts.pgfplots_output)
            printf("\n" "─────────────────────────────────────────────────\n");

        if (pgf && qi < MAX_QUBIT_CONFIGS)
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
