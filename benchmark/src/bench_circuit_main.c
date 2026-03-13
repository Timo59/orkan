/*
 * bench_circuit_main.c -- Circuit-level benchmark orchestration
 *
 * Measurement methodology:
 *   - All competitors are initialized to the maximally mixed state rho = I/dim before
 *     each of the `runs` independent timed runs.
 *   - Within a single timed run, the state evolves across `iterations` circuit applications.
 *     State is NOT re-initialized between iterations within a run.
 *   - This measures SUSTAINED THROUGHPUT -- the steady-state gate-application rate --
 *     rather than single-shot latency from a fixed initial state.
 *   - Timing covers only the gate application loop. State allocation and re-initialization
 *     are outside the timed region.
 *   - All competitors receive the identical circuit_op_t sequence: same gate types,
 *     same qubit indices, same angles.
 *   - All competitors run in density matrix mode.
 *   - OpenMP: benchmarks run with OpenMP enabled (matching qlib's production configuration).
 *     Thread count is controlled by OMP_NUM_THREADS. For thesis results, record the
 *     thread count alongside benchmark output.
 *
 * State re-initialized before each run, not between iterations -- measures sustained throughput.
 */

#include "bench_circuit.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/*
 * =====================================================================================================================
 * Options struct
 * =====================================================================================================================
 */

typedef struct {
    int qubit_min;
    int qubit_max;
    int qubit_step;
    int runs;
    int seed;
    int run_qaoa;
    int run_vqe;
    int run_qv;
    double gamma;
    double beta;
    int iterations_explicit;   /**< 0 = calibrate, >0 = use this value */
    int csv_output;
    int pgfplots_output;
    int verbose;
    int csv_header_printed;  /**< tracks whether the CSV header row has been emitted */
} bench_circuit_options_t;

/*
 * =====================================================================================================================
 * CLI parsing
 * =====================================================================================================================
 */

static void print_usage(const char *prog) {
    printf("Usage: %s [options]\n\n", prog);
    printf("Circuit-level benchmark: QAOA p=1, VQE HEA, Quantum Volume\n\n");
    printf("Options:\n");
    printf("  --qubits-min N    Minimum qubit count (default: 4)\n");
    printf("  --qubits-max N    Maximum qubit count (default: 12)\n");
    printf("  --qubits-step N   Qubit count increment (default: 2)\n");
    printf("  --runs N          Independent timing runs, 1..%d (default: 10)\n", BENCH_MAX_RUNS);
    printf("  --seed N          RNG seed for QV circuit (default: 42)\n");
    printf("  --qaoa            Run only QAOA circuit\n");
    printf("  --vqe             Run only VQE HEA circuit\n");
    printf("  --qv              Run only Quantum Volume circuit\n");
    printf("  --all             Run all circuits (default)\n");
    printf("  --gamma FLOAT     QAOA gamma parameter (default: pi/4)\n");
    printf("  --beta FLOAT      QAOA beta parameter (default: pi/4)\n");
    printf("  --iterations N    Circuit applications per timed run; 0 = calibrate (default: 0)\n");
    printf("  --csv             Output in CSV format\n");
    printf("  --pgfplots        Output in pgfplots-compatible format\n");
    printf("  --verbose         Extended per-result info\n");
    printf("  --help            Show this help\n");
}

static bench_circuit_options_t parse_options(int argc, char *argv[]) {
    bench_circuit_options_t opts;
    opts.qubit_min           = 4;
    opts.qubit_max           = 12;
    opts.qubit_step          = 2;
    opts.runs                = 10;
    opts.seed                = 42;
    opts.run_qaoa            = 0;
    opts.run_vqe             = 0;
    opts.run_qv              = 0;
    opts.gamma               = M_PI / 4.0;
    opts.beta                = M_PI / 4.0;
    opts.iterations_explicit = 0;
    opts.csv_output          = 0;
    opts.pgfplots_output     = 0;
    opts.verbose             = 0;
    opts.csv_header_printed  = 0;

    int circuit_selected = 0;

    for (int i = 1; i < argc; ++i) {
        if        (strcmp(argv[i], "--qubits-min") == 0 && i + 1 < argc) {
            opts.qubit_min = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--qubits-max") == 0 && i + 1 < argc) {
            opts.qubit_max = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--qubits-step") == 0 && i + 1 < argc) {
            opts.qubit_step = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--runs") == 0 && i + 1 < argc) {
            opts.runs = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--seed") == 0 && i + 1 < argc) {
            opts.seed = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--qaoa") == 0) {
            opts.run_qaoa = 1;
            circuit_selected = 1;
        } else if (strcmp(argv[i], "--vqe") == 0) {
            opts.run_vqe = 1;
            circuit_selected = 1;
        } else if (strcmp(argv[i], "--qv") == 0) {
            opts.run_qv = 1;
            circuit_selected = 1;
        } else if (strcmp(argv[i], "--all") == 0) {
            opts.run_qaoa = 1;
            opts.run_vqe  = 1;
            opts.run_qv   = 1;
            circuit_selected = 1;
        } else if (strcmp(argv[i], "--gamma") == 0 && i + 1 < argc) {
            opts.gamma = atof(argv[++i]);
        } else if (strcmp(argv[i], "--beta") == 0 && i + 1 < argc) {
            opts.beta = atof(argv[++i]);
        } else if (strcmp(argv[i], "--iterations") == 0 && i + 1 < argc) {
            opts.iterations_explicit = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--csv") == 0) {
            opts.csv_output = 1;
        } else if (strcmp(argv[i], "--pgfplots") == 0) {
            opts.pgfplots_output = 1;
        } else if (strcmp(argv[i], "--verbose") == 0) {
            opts.verbose = 1;
        } else if (strcmp(argv[i], "--help") == 0) {
            print_usage(argv[0]);
            exit(0);
        } else {
            fprintf(stderr, "error: unrecognized option '%s'\n", argv[i]);
            print_usage(argv[0]);
            exit(1);
        }
    }

    /* Default: run all circuits if none selected */
    if (!circuit_selected) {
        opts.run_qaoa = 1;
        opts.run_vqe  = 1;
        opts.run_qv   = 1;
    }

    /* Validation */
    if (opts.qubit_min < 2)  opts.qubit_min = 2;  /* n=1 causes QAOA rzz(0,0): CX control==target, UB */
    if (opts.qubit_max < opts.qubit_min) opts.qubit_max = opts.qubit_min;
    if (opts.qubit_step < 1) opts.qubit_step = 1;
    if (opts.runs < 1)       opts.runs = 1;
    if (opts.runs > BENCH_MAX_RUNS) opts.runs = BENCH_MAX_RUNS;
    if (opts.seed < 0) {
        fprintf(stderr, "error: --seed must be >= 0\n");
        exit(1);
    }
    if (!isfinite(opts.gamma)) {
        fprintf(stderr, "error: --gamma must be a finite number\n");
        exit(1);
    }
    if (!isfinite(opts.beta)) {
        fprintf(stderr, "error: --beta must be a finite number\n");
        exit(1);
    }

    return opts;
}

/*
 * =====================================================================================================================
 * Runner wrappers for calibrate_iterations (runner_fn signature)
 * =====================================================================================================================
 */

static void qlib_packed_runner(const circuit_t *c, void *ctx, int iterations) {
    (void)ctx;  /* bench_circuit_qlib_packed allocates its own state internally */
    bench_circuit_qlib_packed(c, iterations, 0, 1);
}

static void qlib_tiled_runner(const circuit_t *c, void *ctx, int iterations) {
    (void)ctx;  /* bench_circuit_qlib_tiled allocates its own state internally */
    bench_circuit_qlib_tiled(c, iterations, 0, 1);
}

/*
 * =====================================================================================================================
 * Iteration calibration
 * =====================================================================================================================
 */

/**
 * @brief Calibrate the number of iterations so each timed run lasts at least MIN_MEASUREMENT_MS.
 *
 * Runs BENCH_PQ_PROBE_ITERS iterations for BENCH_PQ_PROBE_RUNS probe rounds, takes the minimum
 * probe time, then scales to target MIN_MEASUREMENT_MS. Clamps result to
 * [1, BENCH_PQ_MAX_CALIBRATED_ITERATIONS].
 */
static int calibrate_iterations(const circuit_t *c, runner_fn fn, void *ctx) {
    double min_probe_ms = 1e30;

    for (int r = 0; r < BENCH_PQ_PROBE_RUNS; ++r) {
        bench_timing_barrier();
        uint64_t t0 = bench_time_ns();
        bench_timing_barrier();

        fn(c, ctx, BENCH_PQ_PROBE_ITERS);

        bench_timing_barrier();
        double elapsed_ms = bench_ns_to_ms(bench_time_ns() - t0);
        if (elapsed_ms < min_probe_ms) min_probe_ms = elapsed_ms;
    }

    /* If probe rounds to effectively zero, use fallback */
    if (min_probe_ms < 1e-6) {
        return BENCH_PQ_FALLBACK_ITERATIONS;
    }

    /* Scale: we want total_ms >= MIN_MEASUREMENT_MS
     * probe ran BENCH_PQ_PROBE_ITERS iterations in min_probe_ms
     * target_iters = ceil(MIN_MEASUREMENT_MS / min_probe_ms * BENCH_PQ_PROBE_ITERS) */
    double iters_needed = (MIN_MEASUREMENT_MS / min_probe_ms) * (double)BENCH_PQ_PROBE_ITERS;
    int result = (int)ceil(iters_needed);

    if (result < 1) result = 1;
    if (result > BENCH_PQ_MAX_CALIBRATED_ITERATIONS) result = BENCH_PQ_MAX_CALIBRATED_ITERATIONS;

    return result;
}

/*
 * =====================================================================================================================
 * CSV and console output
 * =====================================================================================================================
 */

static void print_csv_header(bench_circuit_options_t *opts) {
    if (opts->csv_header_printed) return;
    opts->csv_header_printed = 1;
    printf("circuit,method,n_qubits,n_ops,iterations,runs,mean_ms,ci95_ms,min_ms,median_ms,cv\n");
}

static void print_csv_row(const bench_result_t *r) {
    /* CI95 = 1.96 * std / sqrt(runs) */
    double ci95 = (r->runs > 1)
        ? 1.96 * r->time_ms_std / sqrt((double)r->runs)
        : 0.0;
    printf("%s,%s,%u,%d,%d,%d,%.6f,%.6f,%.6f,%.6f,%.3f\n",
           r->gate_name, r->method,
           r->qubits, r->sweep_size,
           r->iterations, r->runs,
           r->time_ms, ci95,
           r->time_ms_min, r->time_ms_median,
           r->time_ms_cv);
}

static void print_pgfplots_header(const char *circuit_name) {
    printf("# CIRCUIT: %s\n", circuit_name);
    printf("qubits method time_us_per_gate_mean time_us_per_gate_min "
           "time_us_per_gate_median time_us_per_gate_cv_pct memory_bytes\n");
}

static void print_pgfplots_row(const bench_result_t *r) {
    /* Convert time_per_gate from ms to us */
    double us_mean   = r->time_per_gate_ms * 1000.0;
    double total_ops = (double)r->iterations * (double)r->sweep_size;
    double us_min    = (total_ops > 0.0) ? (r->time_ms_min / total_ops) * 1000.0 : 0.0;
    double us_median = (total_ops > 0.0) ? (r->time_ms_median / total_ops) * 1000.0 : 0.0;
    printf("%u %s %.9f %.9f %.9f %.3f %zu\n",
           r->qubits, r->method,
           us_mean, us_min, us_median,
           r->time_ms_cv, r->memory_bytes);
}

static void print_console_result(const bench_result_t *r, int verbose) {
    printf("  %-14s %8.3f +/- %6.3f ms (CV %4.1f%%)  %12.0f ops/sec",
           r->method,
           r->time_ms, r->time_ms_std, r->time_ms_cv,
           r->ops_per_sec);
    if (verbose)
        printf("  [min %8.3f ms  med %8.3f ms  iters %d  mem %zu B]",
               r->time_ms_min, r->time_ms_median, r->iterations, r->memory_bytes);
    printf("\n");
}

/*
 * =====================================================================================================================
 * Correctness verification (debug builds only)
 * =====================================================================================================================
 */

#ifndef NDEBUG

#ifdef WITH_QUEST
/* Pull in QuEST headers for cpuAmps access under debug verification only.
 * This include is guarded so it does not affect non-QuEST builds. */
#include "quest/include/quest.h"
#endif

/**
 * @brief Compare final density matrices from all compiled-in backends (n <= 4 only).
 *
 * Each backend starts from rho = I/dim (maximally mixed) and applies the circuit once.
 * Full density matrices are extracted for qlib_packed, qlib_tiled, and QuEST.
 * For Qulacs, the CTYPE* buffer is ABI-compatible with double _Complex* (row-major,
 * same interleaved layout) so we can read it directly from C.
 * For Aer-DM: the DensityMatrix<double> object is opaque from C; only the diagonal
 * is accessible indirectly. Limitation: Aer-DM is excluded from full-matrix comparison.
 *
 * Comparison: pairwise Frobenius norm of the full dim x dim complex matrices.
 * Tolerance: 1e-6. Prints to stderr so it does not pollute stdout CSV/pgfplots output.
 */
static void verify_correctness(const circuit_t *c) {
    if (c->n_qubits > 4) return;

    const int n = c->n_qubits;
    const dim_t dim = (dim_t)1 << n;

    #define MAX_METHODS 5
    const char *method_names[MAX_METHODS];
    size_t dense_size = (size_t)dim * (size_t)dim * sizeof(cplx_t);
    cplx_t *matrices[MAX_METHODS];
    int n_methods = 0;

    /* --- qlib packed --- */
    {
        state_t s = {0};
        s.type = MIXED_PACKED;
        state_init(&s, (qubit_t)n, NULL);
        if (s.data) {
            size_t nbytes = bench_packed_size((qubit_t)n);
            memset(s.data, 0, nbytes);
            cplx_t dval = 1.0 / (double)dim;
            for (dim_t i = 0; i < dim; ++i) state_set(&s, i, i, dval);

            for (int op = 0; op < c->n_ops; ++op) {
                const circuit_op_t *o = &c->ops[op];
                switch (o->gate_type) {
                    case GATE_H:  h(&s, o->q0);              break;
                    case GATE_X:  x(&s, o->q0);              break;
                    case GATE_Z:  z(&s, o->q0);              break;
                    case GATE_RX: rx(&s, o->q0, o->angle);   break;
                    case GATE_RY: ry(&s, o->q0, o->angle);   break;
                    case GATE_RZ: rz(&s, o->q0, o->angle);   break;
                    case GATE_CX: cx(&s, o->q0, o->q1);      break;
                    default: break;
                }
            }

            cplx_t *m = (cplx_t *)malloc(dense_size);
            if (m) {
                for (dim_t r = 0; r < dim; ++r)
                    for (dim_t col = 0; col < dim; ++col)
                        m[r * dim + col] = state_get(&s, r, col);
                matrices[n_methods] = m;
                method_names[n_methods] = "qlib_packed";
                ++n_methods;
            }
            state_free(&s);
        }
    }

    /* --- qlib tiled --- */
    {
        state_t s = {0};
        s.type = MIXED_TILED;
        state_init(&s, (qubit_t)n, NULL);
        if (s.data) {
            size_t nbytes = bench_tiled_size((qubit_t)n);
            memset(s.data, 0, nbytes);
            cplx_t dval = 1.0 / (double)dim;
            for (dim_t i = 0; i < dim; ++i) state_set(&s, i, i, dval);

            for (int op = 0; op < c->n_ops; ++op) {
                const circuit_op_t *o = &c->ops[op];
                switch (o->gate_type) {
                    case GATE_H:  h(&s, o->q0);              break;
                    case GATE_X:  x(&s, o->q0);              break;
                    case GATE_Z:  z(&s, o->q0);              break;
                    case GATE_RX: rx(&s, o->q0, o->angle);   break;
                    case GATE_RY: ry(&s, o->q0, o->angle);   break;
                    case GATE_RZ: rz(&s, o->q0, o->angle);   break;
                    case GATE_CX: cx(&s, o->q0, o->q1);      break;
                    default: break;
                }
            }

            cplx_t *m = (cplx_t *)malloc(dense_size);
            if (m) {
                for (dim_t r = 0; r < dim; ++r)
                    for (dim_t col = 0; col < dim; ++col)
                        m[r * dim + col] = state_get(&s, r, col);
                matrices[n_methods] = m;
                method_names[n_methods] = "qlib_tiled";
                ++n_methods;
            }
            state_free(&s);
        }
    }

#ifdef WITH_QUEST
    /* --- QuEST --- */
    /* cpuAmps is qcomp* (double _Complex* in C), column-major (col * dim + row), dim^2 elements.
     * Each element is one complex amplitude; direct cast to cplx_t (also double _Complex). */
    {
        Qureg rho = createDensityQureg(n);
        initBlankState(rho);
        {
            qindex qdim = (qindex)1 << n;
            qreal inv_dim = 1.0 / (qreal)qdim;
            for (qindex i = 0; i < qdim; ++i)
                rho.cpuAmps[i * qdim + i] = inv_dim;
        }

        for (int op = 0; op < c->n_ops; ++op) {
            const circuit_op_t *o = &c->ops[op];
            switch (o->gate_type) {
                case GATE_H:  applyHadamard(rho, (int)o->q0);                          break;
                case GATE_X:  applyPauliX(rho, (int)o->q0);                            break;
                case GATE_Z:  applyPauliZ(rho, (int)o->q0);                            break;
                case GATE_RX: applyRotateX(rho, (int)o->q0, o->angle);                 break;
                case GATE_RY: applyRotateY(rho, (int)o->q0, o->angle);                 break;
                case GATE_RZ: applyRotateZ(rho, (int)o->q0, o->angle);                 break;
                case GATE_CX: applyControlledPauliX(rho, (int)o->q0, (int)o->q1);      break;
                default: break;
            }
        }

        /* cpuAmps is qcomp* (= double _Complex* in C), column-major (col * dim + row), dim^2 elements.
         * m[] is row-major (row * dim + col), so transpose during extraction. */
        cplx_t *m = (cplx_t *)malloc(dense_size);
        if (m) {
            qindex qdim = (qindex)1 << n;
            for (qindex r = 0; r < qdim; ++r)
                for (qindex col = 0; col < qdim; ++col)
                    m[r * qdim + col] = (cplx_t)rho.cpuAmps[col * qdim + r];
            matrices[n_methods] = m;
            method_names[n_methods] = "quest";
            ++n_methods;
        }
        destroyQureg(rho);
    }
#endif /* WITH_QUEST */

#ifdef WITH_QULACS
    /* --- Qulacs --- */
    /* Limitation: bench_circuit_qulacs allocates its own internal CTYPE buffer and does
     * not expose it through a C-callable interface. The csim gate functions (dm_H_gate,
     * dm_CNOT_gate, etc.) are declared in C++ headers and cannot be called from this
     * C translation unit. Full Frobenius-norm cross-check requires a C helper such as
     * qulacs_apply_circuit_once(ctx, c) to be added to bench_circuit_qulacs.cpp.
     * TODO: implement that helper to enable Qulacs cross-check here. */
    fprintf(stderr, "  [NOTE] qulacs: full-matrix cross-check skipped"
            " (csim headers are C++ only; add qulacs_apply_circuit_once() helper)\n");
#endif /* WITH_QULACS */

#ifdef WITH_AER_DM
    /* --- Aer-DM --- */
    /* AER::QV::DensityMatrix<double> is an opaque C++ object. Its internal data()
     * pointer is not accessible from C without including C++ Aer headers.
     * Limitation: Aer-DM is excluded from full Frobenius-norm cross-check.
     * To add it, expose an aer_dm_extract_matrix(ctx, out, dim) C helper in
     * bench_circuit_aer_dm.cpp. */
    fprintf(stderr, "  [NOTE] qiskit_aer_dm: full-matrix cross-check skipped"
            " (DensityMatrix<double> is opaque from C; add aer_dm_extract_matrix() helper)\n");
#endif /* WITH_AER_DM */

    /* Pairwise Frobenius norm comparison */
    const double tol = 1e-6;
    for (int a = 0; a < n_methods; ++a) {
        for (int b = a + 1; b < n_methods; ++b) {
            double frob_sq = 0.0;
            for (size_t k = 0; k < (size_t)dim * (size_t)dim; ++k) {
                cplx_t diff = matrices[a][k] - matrices[b][k];
                frob_sq += creal(diff) * creal(diff) + cimag(diff) * cimag(diff);
            }
            double frob = sqrt(frob_sq);
            if (frob < tol) {
                fprintf(stderr, "  [PASS] %s vs %s: ||diff||_F = %.2e\n",
                        method_names[a], method_names[b], frob);
            } else {
                fprintf(stderr, "  [FAIL] %s vs %s: err = %.2e\n",
                        method_names[a], method_names[b], frob);
            }
        }
    }

    for (int i = 0; i < n_methods; ++i) free(matrices[i]);
    #undef MAX_METHODS
}
#endif /* NDEBUG */

/*
 * =====================================================================================================================
 * run_all_methods: benchmark one circuit across all available backends
 * =====================================================================================================================
 */

static void run_all_methods(const circuit_t *c, const bench_circuit_options_t *opts) {
    const int runs = opts->runs;
    const int warmup = 0;  /* time-based warmup, not count-based */

    /* Collect results for console output */
    #define MAX_RESULTS 8
    bench_result_t results[MAX_RESULTS];
    int n_results = 0;

    /* --- qlib packed --- */
    {
        int iters = opts->iterations_explicit;
        if (iters == 0) {
            qlib_circuit_ctx_t *ctx = qlib_circuit_ctx_alloc(c->n_qubits, MIXED_PACKED);
            if (ctx) {
                iters = calibrate_iterations(c, qlib_packed_runner, ctx);
                qlib_circuit_ctx_free(ctx);
            } else {
                iters = BENCH_PQ_FALLBACK_ITERATIONS;
            }
        }
        bench_result_t r = bench_circuit_qlib_packed(c, iters, warmup, runs);
        if (opts->csv_output) print_csv_row(&r);
        else if (opts->pgfplots_output) print_pgfplots_row(&r);
        else print_console_result(&r, opts->verbose);
        if (n_results < MAX_RESULTS) results[n_results++] = r;
    }

    /* --- qlib tiled --- */
    {
        int iters = opts->iterations_explicit;
        if (iters == 0) {
            qlib_circuit_ctx_t *ctx = qlib_circuit_ctx_alloc(c->n_qubits, MIXED_TILED);
            if (ctx) {
                iters = calibrate_iterations(c, qlib_tiled_runner, ctx);
                qlib_circuit_ctx_free(ctx);
            } else {
                iters = BENCH_PQ_FALLBACK_ITERATIONS;
            }
        }
        bench_result_t r = bench_circuit_qlib_tiled(c, iters, warmup, runs);
        if (opts->csv_output) print_csv_row(&r);
        else if (opts->pgfplots_output) print_pgfplots_row(&r);
        else print_console_result(&r, opts->verbose);
        if (n_results < MAX_RESULTS) results[n_results++] = r;
    }

#ifdef WITH_QUEST
    {
        int iters = opts->iterations_explicit;
        if (iters == 0) {
            /* QuEST headers are not available in this C TU; calibrate via
             * the bench function itself (runs 1 timed run with probe iters). */
            double min_probe_ms = 1e30;
            for (int p = 0; p < BENCH_PQ_PROBE_RUNS; ++p) {
                bench_result_t probe = bench_circuit_quest(c, BENCH_PQ_PROBE_ITERS, 0, 1);
                if (probe.time_ms < min_probe_ms) min_probe_ms = probe.time_ms;
            }
            if (min_probe_ms < 1e-6) {
                iters = BENCH_PQ_FALLBACK_ITERATIONS;
            } else {
                double needed = (MIN_MEASUREMENT_MS / min_probe_ms) * (double)BENCH_PQ_PROBE_ITERS;
                iters = (int)ceil(needed);
                if (iters < 1) iters = 1;
                if (iters > BENCH_PQ_MAX_CALIBRATED_ITERATIONS) iters = BENCH_PQ_MAX_CALIBRATED_ITERATIONS;
            }
        }
        bench_result_t r = bench_circuit_quest(c, iters, warmup, runs);
        if (opts->csv_output) print_csv_row(&r);
        else if (opts->pgfplots_output) print_pgfplots_row(&r);
        else print_console_result(&r, opts->verbose);
        if (n_results < MAX_RESULTS) results[n_results++] = r;
    }
#endif

#ifdef WITH_QULACS
    {
        int iters = opts->iterations_explicit;
        if (iters == 0) {
            /* Calibrate via the bench function's own timing */
            double min_probe_ms = 1e30;
            for (int p = 0; p < BENCH_PQ_PROBE_RUNS; ++p) {
                bench_result_t probe = bench_circuit_qulacs(c, BENCH_PQ_PROBE_ITERS, 0, 1);
                if (probe.time_ms < min_probe_ms) min_probe_ms = probe.time_ms;
            }
            if (min_probe_ms < 1e-6) {
                iters = BENCH_PQ_FALLBACK_ITERATIONS;
            } else {
                double needed = (MIN_MEASUREMENT_MS / min_probe_ms) * (double)BENCH_PQ_PROBE_ITERS;
                iters = (int)ceil(needed);
                if (iters < 1) iters = 1;
                if (iters > BENCH_PQ_MAX_CALIBRATED_ITERATIONS) iters = BENCH_PQ_MAX_CALIBRATED_ITERATIONS;
            }
        }
        bench_result_t r = bench_circuit_qulacs(c, iters, warmup, runs);
        if (opts->csv_output) print_csv_row(&r);
        else if (opts->pgfplots_output) print_pgfplots_row(&r);
        else print_console_result(&r, opts->verbose);
        if (n_results < MAX_RESULTS) results[n_results++] = r;
    }
#endif

#ifdef WITH_AER_DM
    {
        int iters = opts->iterations_explicit;
        if (iters == 0) {
            /* Calibrate via the bench function's own timing */
            double min_probe_ms = 1e30;
            for (int p = 0; p < BENCH_PQ_PROBE_RUNS; ++p) {
                bench_result_t probe = bench_circuit_aer_dm(c, BENCH_PQ_PROBE_ITERS, 0, 1);
                if (probe.time_ms < min_probe_ms) min_probe_ms = probe.time_ms;
            }
            if (min_probe_ms < 1e-6) {
                iters = BENCH_PQ_FALLBACK_ITERATIONS;
            } else {
                double needed = (MIN_MEASUREMENT_MS / min_probe_ms) * (double)BENCH_PQ_PROBE_ITERS;
                iters = (int)ceil(needed);
                if (iters < 1) iters = 1;
                if (iters > BENCH_PQ_MAX_CALIBRATED_ITERATIONS) iters = BENCH_PQ_MAX_CALIBRATED_ITERATIONS;
            }
        }
        bench_result_t r = bench_circuit_aer_dm(c, iters, warmup, runs);
        if (opts->csv_output) print_csv_row(&r);
        else if (opts->pgfplots_output) print_pgfplots_row(&r);
        else print_console_result(&r, opts->verbose);
        if (n_results < MAX_RESULTS) results[n_results++] = r;
    }
#endif

    /* Console: print speedup line (packed vs tiled) */
    if (!opts->csv_output && !opts->pgfplots_output && n_results >= 2) {
        if (results[0].time_ms > 0.0 && results[1].time_ms > 0.0) {
            printf("  Speedup: %.2fx tiled vs packed\n",
                   results[0].time_ms / results[1].time_ms);
        }
    }

#ifndef NDEBUG
    verify_correctness(c);
#endif

    #undef MAX_RESULTS
}

/*
 * =====================================================================================================================
 * main
 * =====================================================================================================================
 */

int main(int argc, char *argv[]) {
    bench_circuit_options_t opts = parse_options(argc, argv);

    /* Print header comment with date, build config, OpenMP thread count */
    if (!opts.csv_output && !opts.pgfplots_output) {
        printf("bench_circuit -- Circuit-level density matrix benchmark\n");
        printf("Date: %s %s\n", __DATE__, __TIME__);
#ifdef NDEBUG
        printf("Build: Release\n");
#else
        printf("Build: Debug\n");
#endif
#ifdef _OPENMP
        printf("OpenMP threads: %d\n", omp_get_max_threads());
#else
        printf("OpenMP: disabled\n");
#endif
        printf("Qubits: %d..%d step %d, runs=%d\n",
               opts.qubit_min, opts.qubit_max, opts.qubit_step, opts.runs);
        if (opts.iterations_explicit > 0)
            printf("Iterations: %d (fixed)\n", opts.iterations_explicit);
        else
            printf("Iterations: calibrated (target >= %.1f ms per timed run)\n", MIN_MEASUREMENT_MS);
        printf("Backends: qlib_packed, qlib_tiled");
#ifdef WITH_QUEST
        printf(", quest");
#endif
#ifdef WITH_QULACS
        printf(", qulacs");
#endif
#ifdef WITH_AER_DM
        printf(", qiskit_aer_dm");
#endif
        printf("\n");
        printf("===========================================================\n\n");
    }

#ifdef WITH_QUEST
    bench_quest_init();
#endif

    if (opts.csv_output) print_csv_header(&opts);

    for (int qubits = opts.qubit_min;
         qubits <= opts.qubit_max;
         qubits += opts.qubit_step)
    {
        if (opts.run_qaoa) {
            circuit_t c = circuit_qaoa_p1(qubits, opts.gamma, opts.beta);
            if (!opts.csv_output && !opts.pgfplots_output)
                printf("Circuit: qaoa_p1  qubits: %d  n_ops: %d\n", qubits, c.n_ops);
            if (opts.pgfplots_output) print_pgfplots_header("qaoa_p1");
            run_all_methods(&c, &opts);
            circuit_free(&c);
            if (!opts.csv_output && !opts.pgfplots_output)
                printf("-----------------------------------------------------------\n");
        }

        if (opts.run_vqe) {
            circuit_t c = circuit_vqe_hea(qubits, -1);
            if (!opts.csv_output && !opts.pgfplots_output)
                printf("Circuit: vqe_hea  qubits: %d  n_ops: %d\n", qubits, c.n_ops);
            if (opts.pgfplots_output) print_pgfplots_header("vqe_hea");
            run_all_methods(&c, &opts);
            circuit_free(&c);
            if (!opts.csv_output && !opts.pgfplots_output)
                printf("-----------------------------------------------------------\n");
        }

        if (opts.run_qv) {
            if (qubits % 2 == 0) {
                circuit_t c = circuit_qv(qubits, opts.seed);
                if (!opts.csv_output && !opts.pgfplots_output)
                    printf("Circuit: qv  qubits: %d  n_ops: %d\n", qubits, c.n_ops);
                if (opts.pgfplots_output) print_pgfplots_header("qv");
                run_all_methods(&c, &opts);
                circuit_free(&c);
                if (!opts.csv_output && !opts.pgfplots_output)
                    printf("-----------------------------------------------------------\n");
            } else {
                /* QV requires even qubit count; skip odd with warning (no continue) */
                fprintf(stderr, "warn: QV skipped for odd qubit count %d\n", qubits);
            }
        }
    }

#ifdef WITH_QUEST
    bench_quest_cleanup();
#endif

    return 0;
}
