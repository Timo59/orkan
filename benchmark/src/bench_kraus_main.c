/**
 * @file bench_kraus_main.c
 * @brief CLI parsing, backend table, and orchestration for the Kraus benchmark.
 *
 * Bridges to the gate benchmark formatters by populating bench_backends[]
 * with names and wrapping results into the expected types.
 */

#include "bench_kraus.h"
#include <strings.h>
#include <math.h>

/* =====================================================================
 * Kraus backend table
 * ===================================================================== */

static const bench_kraus_backend_t *kraus_backends[KRAUS_BACKEND_COUNT];

static void init_kraus_backends(void) {
    kraus_backends[0] = &kraus_orkan_tiled_backend;
    kraus_backends[1] = &kraus_orkan_packed_backend;
    kraus_backends[2] = &kraus_blas_backend;
#ifdef HAVE_QUEST
    kraus_backends[3] = &kraus_quest_backend;
#else
    kraus_backends[3] = NULL;
#endif
#ifdef HAVE_QULACS
    kraus_backends[4] = &kraus_qulacs_backend;
#else
    kraus_backends[4] = NULL;
#endif
#ifdef HAVE_AER
    kraus_backends[5] = &kraus_aer_backend;
#else
    kraus_backends[5] = NULL;
#endif
}

/* =====================================================================
 * Populate bench_backends[] for the formatters
 *
 * The formatter code reads bench_backends[b].name and BACKEND_COUNT.
 * We fill all BACKEND_COUNT slots: Kraus backends get their real names,
 * unused gate-only slots (orkan_tiled, orkan_packed) get stub names with
 * available() returning 0.
 * ===================================================================== */

static void *stub_init(qubit_t q, bench_gate_id_t g, double p) {
    (void)q; (void)g; (void)p; return NULL;
}
static void   stub_apply(void *c, const qubit_t *p) { (void)c; (void)p; }
static void   stub_cleanup(void *c) { (void)c; }
static size_t stub_touch(qubit_t q) { (void)q; return 0; }
static int    stub_available(qubit_t q) { (void)q; return 0; }

bench_backend_t bench_backends[BACKEND_COUNT];

static const char *kraus_backend_names[BACKEND_COUNT] = {
    [BACKEND_ORKAN_TILED]  = "orkan_tiled",
    [BACKEND_ORKAN_PACKED] = "orkan_packed",
    [BACKEND_BLAS]        = "BLAS",
    [BACKEND_QUEST]       = "QuEST",
    [BACKEND_QULACS]      = "Qulacs",
    [BACKEND_AER]         = "Qiskit Aer",
};

static void init_formatter_backends(void) {
    bench_backend_t stub = {
        .name = NULL, .init = stub_init, .apply = stub_apply,
        .cleanup = stub_cleanup, .touch_bytes = stub_touch,
        .available = stub_available
    };
    for (int b = 0; b < BACKEND_COUNT; ++b) {
        bench_backends[b] = stub;
        bench_backends[b].name = kraus_backend_names[b];
    }
}

/* =====================================================================
 * Mapping: BACKEND_COUNT slots -> KRAUS_BACKEND_COUNT slots
 *
 * gate backend enum          kraus index
 * BACKEND_ORKAN_TILED  (0)     0
 * BACKEND_ORKAN_PACKED (1)     1
 * BACKEND_BLAS        (2)     2
 * BACKEND_QUEST       (3)     3
 * BACKEND_QULACS      (4)     4
 * BACKEND_AER         (5)     5
 * ===================================================================== */

static int backend_to_kraus(int b) {
    switch (b) {
    case BACKEND_ORKAN_TILED:  return 0;
    case BACKEND_ORKAN_PACKED: return 1;
    case BACKEND_BLAS:        return 2;
    case BACKEND_QUEST:       return 3;
    case BACKEND_QULACS:      return 4;
    case BACKEND_AER:         return 5;
    default:                  return -1;
    }
}

/* =====================================================================
 * CLI parsing
 * ===================================================================== */

static void usage(const char *prog) {
    fprintf(stderr,
        "Usage: %s [OPTIONS]\n\n"
        "Options:\n"
        "  --channel <NAME>  Channel to benchmark (default: depolarize)\n"
        "  --par <float>     Channel parameter (default: 0.1)\n"
        "  --min-qubits <N>  Minimum qubits (default: 2)\n"
        "  --max-qubits <N>  Maximum qubits (default: 10)\n"
        "  --samples <N>     Timed blocks (default: 10)\n"
        "  --iterations <N>  Repetitions per block (default: auto)\n"
        "  --warm-up <N>     Warm-up iterations (default: 3)\n"
        "  --per-qubit       Report per-position statistics\n"
        "  --output <FMT>    console, csv, pgfplots (default: console)\n"
        "  --backends <LIST> Comma-separated backends to run (default: all)\n"
        "                    Names: orkan_tiled,orkan_packed,BLAS,QuEST,Qulacs,\"Qiskit Aer\"\n",
        prog);
    exit(1);
}

/* Gate table stub — required by bench.h (extern), unused by Kraus benchmark */
const bench_gate_info_t bench_gate_table[BG_COUNT] = {{0}};

static bench_channel_id_t parsed_channel;

static bench_options_t parse_kraus_options(int argc, char *argv[]) {
    bench_options_t opts = {
        .gate = BG_X,  /* unused, but must be valid */
        .par = 0.1,
        .par_provided = 0,
        .min_qubits = 2,
        .max_qubits = 10,
        .samples = 10,
        .iterations = 100,
        .iterations_explicit = 0,
        .warm_up = 3,
        .per_qubit = 0,
        .output = FMT_CONSOLE,
        .backend_mask = ~0u
    };
    parsed_channel = BC_DEPOLARIZE;

    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "--channel") == 0 && i + 1 < argc) {
            parsed_channel = bench_resolve_channel(argv[++i]);
            if (parsed_channel == BC_COUNT) {
                fprintf(stderr, "Error: unknown channel '%s'\n", argv[i]);
                exit(1);
            }
        } else if (strcmp(argv[i], "--par") == 0 && i + 1 < argc) {
            opts.par = atof(argv[++i]);
            opts.par_provided = 1;
        } else if (strcmp(argv[i], "--min-qubits") == 0 && i + 1 < argc) {
            int val = atoi(argv[++i]);
            if (val < 0 || val > 30) {
                fprintf(stderr, "Error: --min-qubits %d out of range\n", val);
                exit(1);
            }
            opts.min_qubits = (qubit_t)val;
        } else if (strcmp(argv[i], "--max-qubits") == 0 && i + 1 < argc) {
            int val = atoi(argv[++i]);
            if (val < 0 || val > 30) {
                fprintf(stderr, "Error: --max-qubits %d out of range\n", val);
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
            else { fprintf(stderr, "Error: unknown format '%s'\n", argv[i]); exit(1); }
        } else if (strcmp(argv[i], "--backends") == 0 && i + 1 < argc) {
            opts.backend_mask = 0;
            char buf[256];
            strncpy(buf, argv[++i], sizeof(buf) - 1);
            buf[sizeof(buf) - 1] = '\0';
            for (char *tok = strtok(buf, ","); tok; tok = strtok(NULL, ",")) {
                int found = 0;
                for (int b = 0; b < BACKEND_COUNT; ++b) {
                    if (kraus_backend_names[b] &&
                        strcasecmp(tok, kraus_backend_names[b]) == 0) {
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
        } else if (strcmp(argv[i], "--help") == 0 || strcmp(argv[i], "-h") == 0) {
            usage(argv[0]);
        } else {
            fprintf(stderr, "Error: unknown option '%s'\n", argv[i]);
            exit(1);
        }
    }

    if (opts.min_qubits > opts.max_qubits) {
        fprintf(stderr, "Error: --min-qubits > --max-qubits\n");
        exit(1);
    }

    const bench_channel_info_t *ch = &bench_channel_table[parsed_channel];
    if (ch->parametrised && !opts.par_provided)
        opts.par = 0.1;  /* default for depolarizing */

    if (opts.min_qubits < ch->tgt_qubits) {
        opts.min_qubits = (qubit_t)ch->tgt_qubits;
    }

    return opts;
}

/* =====================================================================
 * Auto-tuned iteration count
 * ===================================================================== */

static int auto_iterations(qubit_t q) {
    /* Kraus maps are ~4x more expensive than single gates (4 operators),
       so scale down iterations accordingly */
    static const int table[] = {
     /* 0     1     2     3     4     5     6    7    8    9   10   11 */
        5000, 5000, 5000, 5000, 5000, 5000, 500, 200, 100, 50, 20, 10
    };
    if (q < (qubit_t)(sizeof table / sizeof table[0])) return table[q];
    return 5;
}

/* =====================================================================
 * Main
 * ===================================================================== */

int main(int argc, char *argv[]) {
    init_kraus_backends();
    init_formatter_backends();

    bench_options_t opts = parse_kraus_options(argc, argv);
    const bench_channel_info_t *ch = &bench_channel_table[parsed_channel];
    int tq = ch->tgt_qubits;

    /* Build Kraus operators */
    idx_t d = (idx_t)1 << tq;
    cplx_t ops[ch->n_ops * d * d];
    ch->build(opts.par, ops);

    /* Fake gate info for formatters — channel name, tgt_qubits */
    bench_gate_info_t fake_gate = {
        .id = BG_X, .name = ch->name,
        .n_qubits = tq, .parametrised = ch->parametrised
    };

    qubit_t positions[BENCH_MAX_POSITIONS * 3];

    if (opts.output == FMT_CSV)
        bench_print_csv_header(&fake_gate, opts.per_qubit);

#define BENCH_MAX_QUBIT_CONFIGS 32
    int n_qubit_configs = opts.max_qubits - opts.min_qubits + 1;
    qubit_t qubit_list[BENCH_MAX_QUBIT_CONFIGS];
    bench_result_t agg_all[BENCH_MAX_QUBIT_CONFIGS][BACKEND_COUNT];

    for (qubit_t q = opts.min_qubits; q <= opts.max_qubits; ++q) {
        if (!opts.iterations_explicit)
            opts.iterations = auto_iterations(q);

        int qi = q - opts.min_qubits;
        qubit_list[qi] = q;
        int n_pos = kraus_gen_positions(q, tq, positions);

        if (!opts.per_qubit) {
            bench_result_t results[BACKEND_COUNT];

            for (int b = 0; b < BACKEND_COUNT; ++b) {
                int ki = backend_to_kraus(b);
                if (!(opts.backend_mask & (1u << b)) ||
                    ki < 0 || !kraus_backends[ki] ||
                    !kraus_backends[ki]->available(q)) {
                    results[b] = bench_result_nan();
                    continue;
                }
                void *ctx = kraus_backends[ki]->init(q, ch->n_ops, tq, ops);
                if (!ctx) { results[b] = bench_result_nan(); continue; }
                kraus_run_aggregate(kraus_backends[ki], ctx, q,
                                    positions, n_pos, tq, &opts, &results[b]);
                kraus_backends[ki]->cleanup(ctx);
            }

            if (opts.output == FMT_CONSOLE)
                bench_print_console_agg(q, &fake_gate, &opts, results);
            else if (opts.output == FMT_CSV)
                bench_print_csv_agg(q, results);

            memcpy(agg_all[qi], results, sizeof(results));

        } else {
            bench_pos_result_t *alloc_results[BACKEND_COUNT];
            const bench_pos_result_t *per_backend[BACKEND_COUNT];

            for (int b = 0; b < BACKEND_COUNT; ++b) {
                int ki = backend_to_kraus(b);
                if (!(opts.backend_mask & (1u << b)) ||
                    ki < 0 || !kraus_backends[ki] ||
                    !kraus_backends[ki]->available(q)) {
                    alloc_results[b] = NULL;
                    per_backend[b] = NULL;
                    continue;
                }
                alloc_results[b] = calloc((size_t)n_pos, sizeof(bench_pos_result_t));
                if (!alloc_results[b]) {
                    per_backend[b] = NULL;
                    continue;
                }
                void *ctx = kraus_backends[ki]->init(q, ch->n_ops, tq, ops);
                if (!ctx) {
                    free(alloc_results[b]);
                    alloc_results[b] = NULL;
                    per_backend[b] = NULL;
                    continue;
                }
                kraus_run_per_qubit(kraus_backends[ki], ctx, q,
                                    positions, n_pos, tq, &opts,
                                    alloc_results[b]);
                kraus_backends[ki]->cleanup(ctx);
                per_backend[b] = alloc_results[b];
            }

            if (opts.output == FMT_CONSOLE) {
                bench_print_console_perq(q, &fake_gate, &opts,
                                         per_backend, n_pos, tq);
            } else if (opts.output == FMT_CSV) {
                for (int b = 0; b < BACKEND_COUNT; ++b) {
                    if (per_backend[b]) {
                        bench_print_csv_perq(q, (bench_backend_id_t)b,
                                             per_backend[b], n_pos, tq);
                    }
                }
            } else if (opts.output == FMT_PGFPLOTS) {
                bench_print_pgfplots_perq(&fake_gate, &opts, q,
                                          per_backend, n_pos, tq);
            }

            for (int b = 0; b < BACKEND_COUNT; ++b)
                free(alloc_results[b]);
        }
    }

    if (!opts.per_qubit && opts.output == FMT_PGFPLOTS)
        bench_print_pgfplots_agg(&fake_gate, &opts, qubit_list, n_qubit_configs,
                                 (const bench_result_t (*)[BACKEND_COUNT])agg_all);

    return 0;
}
