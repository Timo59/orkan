/**
 * @file bench_main.c
 * @brief CLI parsing, gate table, backend table, and orchestration.
 */

#include "bench.h"
#include <strings.h>   /* strcasecmp */
#include <math.h>

/* =====================================================================
 * Gate table
 * ===================================================================== */

const bench_gate_info_t bench_gate_table[BG_COUNT] = {
    [BG_X]   = { BG_X,   "X",     1, 0 },
    [BG_Y]   = { BG_Y,   "Y",     1, 0 },
    [BG_Z]   = { BG_Z,   "Z",     1, 0 },
    [BG_H]   = { BG_H,   "H",     1, 0 },
    [BG_S]   = { BG_S,   "S",     1, 0 },
    [BG_T]   = { BG_T,   "T",     1, 0 },
    [BG_RX]  = { BG_RX,  "RX",    1, 1 },
    [BG_RY]  = { BG_RY,  "RY",    1, 1 },
    [BG_RZ]  = { BG_RZ,  "RZ",    1, 1 },
    [BG_CX]  = { BG_CX,  "CX",    2, 0 },
    [BG_CY]  = { BG_CY,  "CY",    2, 0 },
    [BG_CZ]  = { BG_CZ,  "CZ",    2, 0 },
    [BG_CCX] = { BG_CCX, "CCX",   3, 0 },
};

/* =====================================================================
 * Gate resolver — case-insensitive with aliases
 * ===================================================================== */

bench_gate_id_t bench_resolve_gate(const char *name) {
    for (int i = 0; i < BG_COUNT; ++i)
        if (strcasecmp(name, bench_gate_table[i].name) == 0)
            return (bench_gate_id_t)i;

    if (strcasecmp(name, "CNOT") == 0)    return BG_CX;
    if (strcasecmp(name, "TOFFOLI") == 0) return BG_CCX;

    return BG_COUNT;
}

/* =====================================================================
 * Stub backend — placeholder for backends not yet compiled in
 * ===================================================================== */

static void *stub_init(qubit_t q, bench_gate_id_t g, double p) {
    (void)q; (void)g; (void)p;
    return NULL;
}
static void   stub_apply(void *c, const qubit_t *p) { (void)c; (void)p; }
static void   stub_cleanup(void *c) { (void)c; }
static size_t stub_touch(qubit_t q) { (void)q; return 0; }
static int    stub_available(qubit_t q) { (void)q; return 0; }

/* =====================================================================
 * Backend table — populated at startup
 * ===================================================================== */

bench_backend_t bench_backends[BACKEND_COUNT];

static void init_backends(void) {
    bench_backend_t stub = {
        .name = NULL, .init = stub_init, .apply = stub_apply,
        .cleanup = stub_cleanup, .touch_bytes = stub_touch,
        .available = stub_available
    };

    bench_backends[BACKEND_QLIB_TILED]  = qlib_tiled_backend;
    bench_backends[BACKEND_QLIB_PACKED] = qlib_packed_backend;

    bench_backends[BACKEND_BLAS] = blas_backend;
#ifdef HAVE_QUEST
    bench_backends[BACKEND_QUEST] = quest_backend;
#else
    bench_backends[BACKEND_QUEST]    = stub;
    bench_backends[BACKEND_QUEST].name = "QuEST";
#endif
#ifdef HAVE_QULACS
    bench_backends[BACKEND_QULACS] = qulacs_backend;
#else
    bench_backends[BACKEND_QULACS]    = stub;
    bench_backends[BACKEND_QULACS].name = "Qulacs";
#endif
#ifdef HAVE_AER
    bench_backends[BACKEND_AER] = aer_backend;
#else
    bench_backends[BACKEND_AER]      = stub;
    bench_backends[BACKEND_AER].name = "Qiskit Aer";
#endif
}

/* =====================================================================
 * CLI parsing
 * ===================================================================== */

static void usage(const char *prog) {
    fprintf(stderr,
        "Usage: %s <GATE> [OPTIONS]\n\n"
        "Gates: X Y Z H S T RX RY RZ CX/CNOT CY CZ CCX/TOFFOLI\n\n"
        "Options:\n"
        "  --par <float>     Gate parameter (mandatory for RX/RY/RZ)\n"
        "  --min-qubits <N>  Minimum qubits (default: 2)\n"
        "  --max-qubits <N>  Maximum qubits (default: 8)\n"
        "  --samples <N>     Timed blocks (default: 10)\n"
        "  --iterations <N>  Repetitions per block (default: 100)\n"
        "  --warm-up <N>     Warm-up iterations (default: 3)\n"
        "  --per-qubit       Report per-position statistics\n"
        "  --output <FMT>    console, csv, pgfplots (default: console)\n"
        "  --backends <LIST> Comma-separated backends to run (default: all)\n"
        "                    Names: qlib_tiled,qlib_packed,BLAS,QuEST,Qulacs,\"Qiskit Aer\"\n",
        prog);
    exit(1);
}

bench_options_t bench_parse_options(int argc, char *argv[]) {
    bench_options_t opts = {
        .gate = BG_COUNT,
        .par = 0.0,
        .par_provided = 0,
        .min_qubits = 2,
        .max_qubits = 8,
        .samples = 10,
        .iterations = 100,
        .iterations_explicit = 0,
        .warm_up = 3,
        .per_qubit = 0,
        .output = FMT_CONSOLE,
        .backend_mask = ~0u
    };

    if (argc < 2) usage(argv[0]);

    opts.gate = bench_resolve_gate(argv[1]);
    if (opts.gate == BG_COUNT) {
        fprintf(stderr, "Error: unrecognised gate '%s'\n", argv[1]);
        exit(1);
    }

    for (int i = 2; i < argc; ++i) {
        if (strcmp(argv[i], "--par") == 0 && i + 1 < argc) {
            opts.par = atof(argv[++i]);
            opts.par_provided = 1;
        } else if (strcmp(argv[i], "--min-qubits") == 0 && i + 1 < argc) {
            int val = atoi(argv[++i]);
            if (val < 0 || val > 30) {
                fprintf(stderr, "Error: --min-qubits %d out of range [0, 30]\n", val);
                exit(1);
            }
            opts.min_qubits = (qubit_t)val;
        } else if (strcmp(argv[i], "--max-qubits") == 0 && i + 1 < argc) {
            int val = atoi(argv[++i]);
            if (val < 0 || val > 30) {
                fprintf(stderr, "Error: --max-qubits %d out of range [0, 30]\n", val);
                exit(1);
            }
            opts.max_qubits = (qubit_t)val;
        } else if (strcmp(argv[i], "--samples") == 0 && i + 1 < argc) {
            opts.samples = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--iterations") == 0 && i + 1 < argc) {
            opts.iterations = atoi(argv[++i]);
            opts.iterations_explicit = 1;
        } else if (strcmp(argv[i], "--warm-up") == 0 && i + 1 < argc) {
            opts.warm_up = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--per-qubit") == 0) {
            opts.per_qubit = 1;
        } else if (strcmp(argv[i], "--output") == 0 && i + 1 < argc) {
            ++i;
            if (strcasecmp(argv[i], "console") == 0)       opts.output = FMT_CONSOLE;
            else if (strcasecmp(argv[i], "csv") == 0)       opts.output = FMT_CSV;
            else if (strcasecmp(argv[i], "pgfplots") == 0)  opts.output = FMT_PGFPLOTS;
            else { fprintf(stderr, "Error: unknown output format '%s'\n", argv[i]); exit(1); }
        } else if (strcmp(argv[i], "--backends") == 0 && i + 1 < argc) {
            opts.backend_mask = 0;
            char buf[256];
            strncpy(buf, argv[++i], sizeof(buf) - 1);
            buf[sizeof(buf) - 1] = '\0';
            for (char *tok = strtok(buf, ","); tok; tok = strtok(NULL, ",")) {
                int found = 0;
                for (int b = 0; b < BACKEND_COUNT; ++b) {
                    if (bench_backends[b].name &&
                        strcasecmp(tok, bench_backends[b].name) == 0) {
                        opts.backend_mask |= (1u << b);
                        found = 1;
                        break;
                    }
                }
                if (!found) {
                    fprintf(stderr, "Error: unknown backend '%s'\n", tok);
                    exit(1);
                }
            }
        } else {
            fprintf(stderr, "Error: unknown option '%s'\n", argv[i]);
            exit(1);
        }
    }

    if (opts.samples < 1) {
        fprintf(stderr, "Error: --samples must be >= 1\n");
        exit(1);
    }
    if (opts.iterations < 1) {
        fprintf(stderr, "Error: --iterations must be >= 1\n");
        exit(1);
    }
    if (opts.warm_up < 0) {
        fprintf(stderr, "Error: --warm-up must be >= 0\n");
        exit(1);
    }
    if (opts.min_qubits > opts.max_qubits) {
        fprintf(stderr, "Error: --min-qubits (%d) > --max-qubits (%d)\n",
                (int)opts.min_qubits, (int)opts.max_qubits);
        exit(1);
    }
    if (opts.max_qubits > 30) {
        fprintf(stderr, "Error: --max-qubits %d exceeds limit of 30\n",
                (int)opts.max_qubits);
        exit(1);
    }
    if (bench_gate_table[opts.gate].parametrised && !opts.par_provided) {
        fprintf(stderr, "Error: gate %s requires --par <angle>\n",
                bench_gate_table[opts.gate].name);
        exit(1);
    }
    if (opts.min_qubits < bench_gate_table[opts.gate].n_qubits) {
        fprintf(stderr, "Error: gate %s requires at least %d qubits (--min-qubits %d)\n",
                bench_gate_table[opts.gate].name,
                bench_gate_table[opts.gate].n_qubits,
                (int)opts.min_qubits);
        exit(1);
    }

    return opts;
}

/* =====================================================================
 * Auto-tuned iteration count per qubit count
 *
 * Targets a timed-block duration of ~1 ms (1000x CLOCK_MONOTONIC
 * resolution).  Calibrated on Apple M2 Pro; other machines may
 * overshoot or undershoot, but the order of magnitude is right.
 * ===================================================================== */

static int auto_iterations(qubit_t q) {
    static const int table[] = {
     /* 0      1      2      3      4      5      6     7    8    9   10   11 */
        10000, 10000, 10000, 10000, 10000, 10000, 1000, 500, 200, 100, 50, 20
    };
    if (q < (qubit_t)(sizeof table / sizeof table[0])) return table[q];
    return 10;
}

/* =====================================================================
 * Main
 * ===================================================================== */

int main(int argc, char *argv[]) {
    init_backends();

    bench_options_t opts = bench_parse_options(argc, argv);
    const bench_gate_info_t *gate = &bench_gate_table[opts.gate];
    int gq = gate->n_qubits;

    qubit_t positions[BENCH_MAX_POSITIONS * 3];

    if (opts.output == FMT_CSV)
        bench_print_csv_header(gate, opts.per_qubit);

    /* Pgfplots aggregate accumulator (max 31 configs: qubits 0..30) */
#define BENCH_MAX_QUBIT_CONFIGS 32
    int n_qubit_configs = opts.max_qubits - opts.min_qubits + 1;
    qubit_t qubit_list[BENCH_MAX_QUBIT_CONFIGS];
    bench_result_t agg_all[BENCH_MAX_QUBIT_CONFIGS][BACKEND_COUNT];

    for (qubit_t q = opts.min_qubits; q <= opts.max_qubits; ++q) {
        if (!opts.iterations_explicit)
            opts.iterations = auto_iterations(q);

        int qi = q - opts.min_qubits;
        qubit_list[qi] = q;
        int n_pos = bench_gen_positions(q, gq, positions);

        if (!opts.per_qubit) {
            bench_result_t results[BACKEND_COUNT];

            for (int b = 0; b < BACKEND_COUNT; ++b) {
                if (!(opts.backend_mask & (1u << b)) ||
                    !bench_backends[b].available(q)) {
                    results[b] = bench_result_nan();
                    continue;
                }
                void *ctx = bench_backends[b].init(q, opts.gate, opts.par);
                if (!ctx) { results[b] = bench_result_nan(); continue; }
                bench_run_aggregate(&bench_backends[b], ctx, q,
                                    positions, n_pos, gq, &opts, &results[b]);
                bench_backends[b].cleanup(ctx);
            }

            if (opts.output == FMT_CONSOLE)
                bench_print_console_agg(q, gate, &opts, results);
            else if (opts.output == FMT_CSV)
                bench_print_csv_agg(q, results);

            memcpy(agg_all[qi], results, sizeof(results));

        } else {
            bench_pos_result_t *alloc_results[BACKEND_COUNT];
            const bench_pos_result_t *per_backend[BACKEND_COUNT];

            for (int b = 0; b < BACKEND_COUNT; ++b) {
                if (!(opts.backend_mask & (1u << b)) ||
                    !bench_backends[b].available(q)) {
                    alloc_results[b] = NULL;
                    per_backend[b] = NULL;
                    continue;
                }
                alloc_results[b] = calloc((size_t)n_pos, sizeof(bench_pos_result_t));
                if (!alloc_results[b]) {
                    per_backend[b] = NULL;
                    continue;
                }
                void *ctx = bench_backends[b].init(q, opts.gate, opts.par);
                if (!ctx) {
                    free(alloc_results[b]);
                    alloc_results[b] = NULL;
                    per_backend[b] = NULL;
                    continue;
                }
                bench_run_per_qubit(&bench_backends[b], ctx, q,
                                    positions, n_pos, gq, &opts,
                                    alloc_results[b]);
                bench_backends[b].cleanup(ctx);
                per_backend[b] = alloc_results[b];
            }

            if (opts.output == FMT_CONSOLE) {
                bench_print_console_perq(q, gate, &opts, per_backend, n_pos, gq);
            } else if (opts.output == FMT_CSV) {
                for (int b = 0; b < BACKEND_COUNT; ++b) {
                    if (per_backend[b]) {
                        bench_print_csv_perq(q, (bench_backend_id_t)b,
                                             per_backend[b], n_pos, gq);
                    } else {
                        /* Emit nan rows with correct positions */
                        for (int p = 0; p < n_pos; ++p) {
                            const qubit_t *pos = &positions[p * gq];
                            printf("%d,%s,", (int)q, bench_backends[b].name);
                            if (gq == 1)      printf("%d,", (int)pos[0]);
                            else if (gq == 2) printf("%d,%d,", (int)pos[0], (int)pos[1]);
                            else              printf("%d,%d,%d,", (int)pos[0], (int)pos[1], (int)pos[2]);
                            printf("nan,nan,nan,nan,nan,nan\n");
                        }
                    }
                }
            } else if (opts.output == FMT_PGFPLOTS) {
                bench_print_pgfplots_perq(gate, &opts, q,
                                          per_backend, n_pos, gq);
            }

            for (int b = 0; b < BACKEND_COUNT; ++b)
                free(alloc_results[b]);
        }
    }

    if (!opts.per_qubit && opts.output == FMT_PGFPLOTS)
        bench_print_pgfplots_agg(gate, &opts, qubit_list, n_qubit_configs,
                                 (const bench_result_t (*)[BACKEND_COUNT])agg_all);

    return 0;
}
