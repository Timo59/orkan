/**
 * @file bench_main.c
 * @brief Orchestration, CLI parsing, gate tables, and main() for bench_gate
 *
 * Output formatting extracted to bench_csv.c, bench_console.c, bench_pgfplots.c.
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
           "                   With --per-qubit: one block per (gate, qubit-count) with\n"
           "                   rows indexed by target qubit (1Q) or pair label (2Q).\n");
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

    /* When --per-qubit and --pgfplots are combined, --csv is ambiguous: the
     * caller must choose one output format explicitly. */
    if (opts.per_qubit && opts.pgfplots_output && opts.csv_output) {
        fprintf(stderr, "error: --csv and --pgfplots are mutually exclusive in "
                "--per-qubit mode; choose one output format\n");
        exit(1);
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

const char *all_gate_names[NUM_GATES] = {"X", "H", "Z", "CX", "SWAP"};

/* Per-qubit gate tables — shared with bench_pgfplots.c via extern in bench.h */

const gate_def_1q_perq_t gates_1q_perq[] = {
    {"X", x, XMAT, 1},
    {"H", h, HMAT, 1},
    {"Z", z, ZMAT, 1},
    {NULL, NULL, NULL, 0}
};

const gate_def_2q_perq_t gates_2q_perq[] = {
    {"CX",   cx,        CXMAT,   1, 0},   /* ordered pairs (CX is not symmetric) */
    {"SWAP", swap_gate, SWAPMAT, 1, 1},   /* unordered pairs (SWAP is symmetric) */
    {NULL, NULL, NULL, 0, 0}
};

/*
 * =====================================================================================================================
 * Per-qubit benchmark
 * =====================================================================================================================
 */

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

/*
 * =====================================================================================================================
 * Per-qubit probe–calibrate–measure helpers
 * =====================================================================================================================
 *
 * Each backend in perq_main() follows an identical three-phase pattern:
 *   1. Quick single-shot probe to detect slow gates (avoids 100× over-measurement).
 *   2. Multi-run probe to estimate per-iteration cost for calibration.
 *   3. Calibrated measurement loop over all targets (1Q) or pairs (2Q).
 *
 * perq_run_1q / perq_run_2q extract that pattern, taking a function pointer for the
 * actual bench_XXX_at call so backends can differ only in their callback signature
 * while sharing identical probe/calibrate/measure logic.
 */

/* perq_bench_1q_fn and perq_bench_2q_fn typedefs are defined in bench.h */

/**
 * @brief Probe, calibrate, and measure all targets for one 1Q backend.
 *
 * @param qubits     Number of qubits.
 * @param gate_name  Gate name string for results.
 * @param fn         Backend benchmark callback (matches perq_bench_1q_fn).
 * @param state      Opaque backend state pointer passed through to fn.
 * @param iters_explicit Non-zero if --iterations was provided explicitly.
 * @param explicit_iters Explicit iteration count when iters_explicit != 0.
 * @param runs       Number of timing runs for the final measurement.
 * @param out        Output array: out[tgt] is populated for tgt in [0, qubits).
 */
void perq_run_1q(qubit_t qubits, const char *gate_name,
                 perq_bench_1q_fn fn, void *state,
                 int iters_explicit, int explicit_iters, int runs,
                 bench_result_perq_t *out)
{
    bench_result_perq_t quick = fn(qubits, gate_name, 0, 1, 1, state);
    int probe_iters = (quick.time_per_gate_us / 1000.0 >= MIN_MEASUREMENT_MS)
                      ? 1 : BENCH_PQ_PROBE_ITERS;
    bench_result_perq_t probe = fn(qubits, gate_name, 0,
                                   probe_iters, BENCH_PQ_PROBE_RUNS, state);
    int cal_iters = calibrate_iterations(&probe, iters_explicit, explicit_iters);

    for (qubit_t tgt = 0; tgt < qubits; tgt++)
        out[tgt] = fn(qubits, gate_name, tgt, cal_iters, runs, state);
}

/**
 * @brief Probe, calibrate, and measure all pairs for one 2Q backend.
 *
 * @param qubits     Number of qubits.
 * @param gate_name  Gate name string for results.
 * @param fn         Backend benchmark callback (matches perq_bench_2q_fn).
 * @param state      Opaque backend state pointer passed through to fn.
 * @param q1s        Array of first-qubit indices for each pair.
 * @param q2s        Array of second-qubit indices for each pair.
 * @param n_pairs    Number of pairs in q1s/q2s.
 * @param iters_explicit Non-zero if --iterations was provided explicitly.
 * @param explicit_iters Explicit iteration count when iters_explicit != 0.
 * @param runs       Number of timing runs for the final measurement.
 * @param out        Output array: out[p] is populated for p in [0, n_pairs).
 */
void perq_run_2q(qubit_t qubits, const char *gate_name,
                 perq_bench_2q_fn fn, void *state,
                 const qubit_t *q1s, const qubit_t *q2s, int n_pairs,
                 int iters_explicit, int explicit_iters, int runs,
                 bench_result_perq_t *out)
{
    bench_result_perq_t quick = fn(qubits, gate_name,
                                   q1s[0], q2s[0], 1, 1, state);
    int probe_iters = (quick.time_per_gate_us / 1000.0 >= MIN_MEASUREMENT_MS)
                      ? 1 : BENCH_PQ_PROBE_ITERS;
    bench_result_perq_t probe = fn(qubits, gate_name,
                                   q1s[0], q2s[0],
                                   probe_iters, BENCH_PQ_PROBE_RUNS, state);
    int cal_iters = calibrate_iterations(&probe, iters_explicit, explicit_iters);

    for (int p = 0; p < n_pairs; p++)
        out[p] = fn(qubits, gate_name, q1s[p], q2s[p], cal_iters, runs, state);
}

/*
 * =====================================================================================================================
 * Per-backend wrapper contexts and trampolines for perq_run_1q / perq_run_2q
 * =====================================================================================================================
 *
 * Each backend bench_XXX_at function has a unique typed signature (different gate fn/mat
 * argument). To use perq_run_1q / perq_run_2q we need thin wrappers that match the
 * perq_bench_1q_fn / perq_bench_2q_fn callback signature:
 *   bench_result_perq_t cb(qubits, name, target, iters, runs, void *ctx)
 *
 * The ctx pointer points to a small stack-allocated context struct that carries the
 * extra per-backend arguments (gate function pointer, gate matrix, or state pointer).
 *
 * For backends where the state is part of the context (QuEST, Qulacs, Aer) the ctx
 * IS the state pointer; no extra struct needed.
 */

/* Context for qlib packed/tiled 1Q wrappers */
typedef struct { void (*gate_fn)(state_t*, qubit_t); state_t *state; } perq_ctx_qlib_1q;
/* Context for qlib packed/tiled 2Q wrappers */
typedef struct { void (*gate_fn)(state_t*, qubit_t, qubit_t); state_t *state; } perq_ctx_qlib_2q;
/* Context for BLAS dense 1Q wrapper */
typedef struct { const cplx_t *mat; cplx_t *rho; } perq_ctx_blas_1q;
/* Context for BLAS dense 2Q wrapper */
typedef struct { const cplx_t *mat; cplx_t *rho; } perq_ctx_blas_2q;

/* ─── qlib packed 1Q ─────────────────────────────────────────────────────── */
static bench_result_perq_t perq_wrap_qlib_packed_1q(
    qubit_t qubits, const char *name, qubit_t tgt, int iters, int runs, void *ctx)
{
    perq_ctx_qlib_1q *c = (perq_ctx_qlib_1q *)ctx;
    return bench_qlib_packed_at(qubits, name, c->gate_fn, tgt, iters, runs, c->state);
}

/* ─── qlib tiled 1Q ──────────────────────────────────────────────────────── */
static bench_result_perq_t perq_wrap_qlib_tiled_1q(
    qubit_t qubits, const char *name, qubit_t tgt, int iters, int runs, void *ctx)
{
    perq_ctx_qlib_1q *c = (perq_ctx_qlib_1q *)ctx;
    return bench_qlib_tiled_at(qubits, name, c->gate_fn, tgt, iters, runs, c->state);
}

/* ─── qlib packed 2Q ─────────────────────────────────────────────────────── */
static bench_result_perq_t perq_wrap_qlib_packed_2q(
    qubit_t qubits, const char *name, qubit_t q1, qubit_t q2, int iters, int runs, void *ctx)
{
    perq_ctx_qlib_2q *c = (perq_ctx_qlib_2q *)ctx;
    return bench_qlib_packed_2q_at(qubits, name, c->gate_fn, q1, q2, iters, runs, c->state);
}

/* ─── qlib tiled 2Q ──────────────────────────────────────────────────────── */
static bench_result_perq_t perq_wrap_qlib_tiled_2q(
    qubit_t qubits, const char *name, qubit_t q1, qubit_t q2, int iters, int runs, void *ctx)
{
    perq_ctx_qlib_2q *c = (perq_ctx_qlib_2q *)ctx;
    return bench_qlib_tiled_2q_at(qubits, name, c->gate_fn, q1, q2, iters, runs, c->state);
}

/* ─── BLAS dense 1Q ──────────────────────────────────────────────────────── */
static bench_result_perq_t perq_wrap_blas_1q(
    qubit_t qubits, const char *name, qubit_t tgt, int iters, int runs, void *ctx)
{
    perq_ctx_blas_1q *c = (perq_ctx_blas_1q *)ctx;
    return bench_blas_dense_at(qubits, name, c->mat, tgt, iters, runs, c->rho);
}

/* ─── BLAS dense 2Q ──────────────────────────────────────────────────────── */
static bench_result_perq_t perq_wrap_blas_2q(
    qubit_t qubits, const char *name, qubit_t q1, qubit_t q2, int iters, int runs, void *ctx)
{
    perq_ctx_blas_2q *c = (perq_ctx_blas_2q *)ctx;
    return bench_blas_dense_2q_at(qubits, name, c->mat, q1, q2, iters, runs, c->rho);
}

/* QuEST trampolines live in bench_quest.c — no Qureg type needed here. */

#ifdef WITH_QULACS
/* ─── Qulacs 1Q — ctx IS the rho pointer ─────────────────────────────────── */
static bench_result_perq_t perq_wrap_qulacs_1q(
    qubit_t qubits, const char *name, qubit_t tgt, int iters, int runs, void *ctx)
{
    return bench_qulacs_at(qubits, name, tgt, iters, runs, ctx);
}
/* ─── Qulacs 2Q ──────────────────────────────────────────────────────────── */
static bench_result_perq_t perq_wrap_qulacs_2q(
    qubit_t qubits, const char *name, qubit_t q1, qubit_t q2, int iters, int runs, void *ctx)
{
    return bench_qulacs_2q_at(qubits, name, q1, q2, iters, runs, ctx);
}
#endif /* WITH_QULACS */

#ifdef WITH_AER_DM
/* ─── Aer DM 1Q — ctx IS the dm pointer ──────────────────────────────────── */
static bench_result_perq_t perq_wrap_aer_1q(
    qubit_t qubits, const char *name, qubit_t tgt, int iters, int runs, void *ctx)
{
    return bench_aer_dm_at(qubits, name, tgt, iters, runs, ctx);
}
/* ─── Aer DM 2Q ──────────────────────────────────────────────────────────── */
static bench_result_perq_t perq_wrap_aer_2q(
    qubit_t qubits, const char *name, qubit_t q1, qubit_t q2, int iters, int runs, void *ctx)
{
    return bench_aer_dm_2q_at(qubits, name, q1, q2, iters, runs, ctx);
}
#endif /* WITH_AER_DM */

/*
 * =====================================================================================================================
 * Per-qubit benchmark orchestration
 * =====================================================================================================================
 */

/**
 * @brief Per-qubit benchmark: time each gate at every target qubit/pair, for every method.
 *
 * Iterates over qubits from opts->min_qubits to opts->max_qubits, runs all
 * 1Q and 2Q gates, and records per-target timing for each enabled method.
 * Calibrates iteration count per (gate, qubits, method) combination.
 * When opts->pgfplots_output is set, accumulates results into a
 * pgfplots_perq_data_t and emits them via print_pgfplots_perq_output()
 * at the end instead of printing per-gate console/CSV output.
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

    /* Allocate pgfplots accumulator when --pgfplots is requested */
    pgfplots_perq_data_t *pgfq = NULL;
    if (opts->pgfplots_output) {
        pgfq = calloc(1, sizeof(pgfplots_perq_data_t));
        if (!pgfq) {
            fprintf(stderr, "error: calloc failed for pgfplots_perq_data_t\n");
            free(perq);
            return EXIT_FAILURE;
        }
    }

    /* CSV output goes to stdout, same as the sweep benchmark */
    FILE *csv_f = stdout;

    /* Emit header — suppressed when accumulating for pgfplots */
    if (opts->csv_output && !opts->pgfplots_output) {
        bench_print_perq_csv_header(csv_f);
    } else if (!opts->pgfplots_output) {
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
    /* pgfplots mode: emit nothing until print_pgfplots_perq_output() */

#ifdef WITH_QUEST
    bench_quest_init();
#endif

    int results_count = 0;
    int qi = 0;  /* qubit-config index for pgfq accumulation */

    for (qubit_t qubits = opts->min_qubits; qubits <= opts->max_qubits;
         qubits += (qubit_t)opts->step) {

        dim_t  dim       = (dim_t)1 << qubits;
        int    run_dense = (qubits <= 8);
        size_t state_bytes = (size_t)dim * (size_t)dim * sizeof(cplx_t);

        /* Register qubit count in pgfq accumulator */
        if (pgfq && qi < MAX_QUBIT_CONFIGS) {
            pgfq->qubits[qi] = qubits;
            pgfq->num_configs = qi + 1;
        }

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
                    perq_ctx_qlib_1q ctx = {gd->fn, &state};
                    perq_run_1q(qubits, gd->name, perq_wrap_qlib_packed_1q, &ctx,
                                opts->iterations_explicit, opts->iterations,
                                opts->runs, perq[m_packed]);
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
                    perq_ctx_qlib_1q ctx = {gd->fn, &state};
                    perq_run_1q(qubits, gd->name, perq_wrap_qlib_tiled_1q, &ctx,
                                opts->iterations_explicit, opts->iterations,
                                opts->runs, perq[m_tiled]);
                    state_free(&state);
                }
            }

            /* --- BLAS dense (capped at 8 qubits) ------------------------ */
            if (run_dense) {
                size_t rho_bytes = state_bytes;
                cplx_t *rho = bench_alloc_huge(rho_bytes);
                if (rho) {
                    bench_fill_random(rho, rho_bytes);
                    perq_ctx_blas_1q ctx = {gd->mat, rho};
                    perq_run_1q(qubits, gd->name, perq_wrap_blas_1q, &ctx,
                                opts->iterations_explicit, opts->iterations,
                                opts->runs, perq[m_blas]);
                    bench_free_huge(rho, rho_bytes);
                }
            }

#ifdef WITH_QUEST
            /* --- QuEST -------------------------------------------------- */
            bench_quest_perq_1q(qubits, gd->name,
                                opts->iterations_explicit, opts->iterations,
                                opts->runs, perq[m_quest]);
#endif /* WITH_QUEST */

#ifdef WITH_QULACS
            /* --- Qulacs ------------------------------------------------- */
            {
                size_t rho_bytes = state_bytes;
                void *rho = bench_alloc_huge(rho_bytes);
                if (rho) {
                    bench_fill_random(rho, rho_bytes);
                    perq_run_1q(qubits, gd->name, perq_wrap_qulacs_1q, rho,
                                opts->iterations_explicit, opts->iterations,
                                opts->runs, perq[m_qulacs]);
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
                    perq_run_1q(qubits, gd->name, perq_wrap_aer_1q, dm,
                                opts->iterations_explicit, opts->iterations,
                                opts->runs, perq[m_aer]);
                    bench_aer_dm_free(dm);
                }
            }
#endif /* WITH_AER_DM */

            /* Output / accumulate results for this 1Q gate × qubit count */
            if (opts->pgfplots_output) {
                /* Accumulate into pgfq; emit deferred via print_pgfplots_perq_output() */
                if (pgfq && qi < MAX_QUBIT_CONFIGS) {
                    for (int m = 0; m < NUM_METHODS; m++) {
                        pgfq->n_targets[g][qi][m] = n_targets;
                        for (int t = 0; t < n_targets && t < MAX_TARGETS; t++) {
                            if (perq[m][t].time_ms_mean > 0.0)
                                perq_cell_from_result(&pgfq->cells[g][qi][m][t], &perq[m][t]);
                        }
                    }
                }
            } else if (opts->csv_output) {
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
                    perq_ctx_qlib_2q ctx = {gd->fn, &state};
                    perq_run_2q(qubits, gd->name, perq_wrap_qlib_packed_2q, &ctx,
                                q1s, q2s, n_pairs,
                                opts->iterations_explicit, opts->iterations,
                                opts->runs, perq[m_packed]);
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
                    perq_ctx_qlib_2q ctx = {gd->fn, &state};
                    perq_run_2q(qubits, gd->name, perq_wrap_qlib_tiled_2q, &ctx,
                                q1s, q2s, n_pairs,
                                opts->iterations_explicit, opts->iterations,
                                opts->runs, perq[m_tiled]);
                    state_free(&state);
                }
            }

            /* --- BLAS dense (capped at 8 qubits) ------------------------ */
            if (run_dense) {
                size_t rho_bytes = state_bytes;
                cplx_t *rho = bench_alloc_huge(rho_bytes);
                if (rho) {
                    bench_fill_random(rho, rho_bytes);
                    perq_ctx_blas_2q ctx = {gd->mat, rho};
                    perq_run_2q(qubits, gd->name, perq_wrap_blas_2q, &ctx,
                                q1s, q2s, n_pairs,
                                opts->iterations_explicit, opts->iterations,
                                opts->runs, perq[m_blas]);
                    bench_free_huge(rho, rho_bytes);
                }
            }

#ifdef WITH_QUEST
            /* --- QuEST -------------------------------------------------- */
            bench_quest_perq_2q(qubits, gd->name,
                                q1s, q2s, n_pairs,
                                opts->iterations_explicit, opts->iterations,
                                opts->runs, perq[m_quest]);
#endif /* WITH_QUEST */

#ifdef WITH_QULACS
            /* --- Qulacs ------------------------------------------------- */
            {
                size_t rho_bytes = state_bytes;
                void *rho = bench_alloc_huge(rho_bytes);
                if (rho) {
                    bench_fill_random(rho, rho_bytes);
                    perq_run_2q(qubits, gd->name, perq_wrap_qulacs_2q, rho,
                                q1s, q2s, n_pairs,
                                opts->iterations_explicit, opts->iterations,
                                opts->runs, perq[m_qulacs]);
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
                    perq_run_2q(qubits, gd->name, perq_wrap_aer_2q, dm,
                                q1s, q2s, n_pairs,
                                opts->iterations_explicit, opts->iterations,
                                opts->runs, perq[m_aer]);
                    bench_aer_dm_free(dm);
                }
            }
#endif /* WITH_AER_DM */

            /* Output / accumulate results for this 2Q gate × qubit count */
            if (opts->pgfplots_output) {
                /* Accumulate into pgfq; emit deferred via print_pgfplots_perq_output() */
                if (pgfq && qi < MAX_QUBIT_CONFIGS) {
                    int pg = NUM_1Q_GATES + g;  /* 2Q gate index offset */
                    for (int m = 0; m < NUM_METHODS; m++) {
                        pgfq->n_targets[pg][qi][m] = n_pairs;
                        for (int t = 0; t < n_pairs && t < MAX_TARGETS; t++) {
                            if (perq[m][t].time_ms_mean > 0.0) {
                                perq_cell_from_result(&pgfq->cells[pg][qi][m][t], &perq[m][t]);
                            }
                        }
                    }
                }
            } else if (opts->csv_output) {
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

        /* Advance qubit-config index for pgfq */
        if (pgfq && qi < MAX_QUBIT_CONFIGS)
            qi++;

    } /* qubits loop */

#ifdef WITH_QUEST
    bench_quest_cleanup();
#endif

    if (results_count == 0 && opts->gate_filter) {
        fprintf(stderr, "warning: no results — gate '%s' not found\n",
                opts->gate_filter);
        free(pgfq);
        free(perq);
        return EXIT_FAILURE;
    }

    /* Emit pgfplots output now that all data is accumulated */
    if (pgfq) {
        print_pgfplots_perq_output(pgfq);
        free(pgfq);
        pgfq = NULL;
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
