/**
 * @file profile_channel_tiled.c
 * @brief Micro-profiling harness for channel_tiled_1q() vs h_tiled()
 *
 * Profiles the tiled 1-qubit channel (depolarizing, p=0.1) alongside the
 * Hadamard gate to identify the source of the performance gap.
 *
 * Usage:
 *   ./profile_channel_tiled [options]
 *
 * Options:
 *   --min-qubits N    Minimum qubit count (default: LOG_TILE_DIM)
 *   --max-qubits N    Maximum qubit count (default: 12)
 *   --warmup N        Warmup iterations (default: 200)
 *   --iters N         Measurement iterations (default: 2000)
 *   --csv FILE        Write CSV output to FILE
 *   --target N        Profile only target qubit N (default: all valid targets)
 *   -v, --verbose     Print extra info
 */

#include "profile.h"
#include <string.h>
#include <getopt.h>
#include <math.h>

/* Channel types and functions — declared here to avoid header conflict
 * (channel.h and gate.h both define insertBit0 with different idx types) */
typedef struct { qubit_t n_qubits; uint64_t n_terms; cplx_t* data; } kraus_t;
typedef struct { qubit_t n_qubits; cplx_t* data; } superop_t;
extern superop_t kraus_to_superop(const kraus_t *kraus);

/* Backend function declarations */
extern void channel_tiled_1q(state_t *state, const cplx_t * restrict sop, qubit_t target);
extern void channel_packed_1q(state_t *state, const cplx_t * restrict sop, qubit_t target);
extern void h_tiled(state_t *state, const qubit_t target);

/*
 * =====================================================================================================================
 * Depolarizing channel constructor (p=0.1)
 * =====================================================================================================================
 */

static kraus_t make_depolarizing(double p) {
    kraus_t k;
    k.n_qubits = 1;
    k.n_terms = 4;
    k.data = (cplx_t *)calloc(4 * 4, sizeof(cplx_t));

    double s0 = sqrt(1.0 - 3.0 * p / 4.0);
    double sp = sqrt(p / 4.0);

    /* K0 = s0 * I */
    k.data[0] = s0; k.data[3] = s0;
    /* K1 = sp * X */
    k.data[5] = sp; k.data[6] = sp;
    /* K2 = sp * Y */
    k.data[9] = -sp * I; k.data[10] = sp * I;
    /* K3 = sp * Z */
    k.data[12] = sp; k.data[15] = -sp;

    return k;
}

/*
 * =====================================================================================================================
 * Channel profiling runner
 * =====================================================================================================================
 */

typedef void (*channel_1q_fn)(state_t *, const cplx_t * restrict, qubit_t);

static prof_result_t prof_run_channel_1q(
    channel_1q_fn   fn,
    const cplx_t   *sop,
    const char     *func_name,
    state_type_t    stype,
    qubit_t         qubits,
    qubit_t         target,
    int             warmup,
    int             iters)
{
    prof_result_t r = {0};
    r.func_name = func_name;
    r.qubits = qubits;
    r.target = target;

    if (stype == MIXED_TILED) {
        r.path = (target < LOG_TILE_DIM) ? "within-tile" : "cross-tile";
        r.touch_bytes = prof_tiled_bytes(qubits);
    } else {
        r.path = "packed";
        r.touch_bytes = prof_packed_bytes(qubits);
    }

    state_t state = {0};
    state.type = stype;
    state_plus(&state, qubits);
    r.data_bytes = (size_t)state_len(&state) * sizeof(cplx_t);

    /* Warm up */
    for (int i = 0; i < warmup; i++) {
        fn(&state, sop, target);
    }

    /* Measure */
    double *samples = (double *)malloc((size_t)iters * sizeof(double));
    for (int i = 0; i < iters; i++) {
        uint64_t t0 = prof_time_ns();
        fn(&state, sop, target);
        uint64_t t1 = prof_time_ns();
        samples[i] = (double)(t1 - t0);
    }

    r.stats = prof_compute_stats(samples, iters);
    free(samples);
    state_free(&state);

    if (r.stats.median_ns > 0) {
        r.bw_gbs = (double)r.touch_bytes / r.stats.median_ns;
        uint64_t dim = (uint64_t)1 << qubits;
        uint64_t n_tiles = (dim + TILE_DIM - 1) / TILE_DIM;
        uint64_t total_elems = n_tiles * (n_tiles + 1) / 2 * TILE_SIZE;
        r.elem_per_ns = (double)total_elems / r.stats.median_ns;
    }

    return r;
}

/*
 * =====================================================================================================================
 * Command-line options
 * =====================================================================================================================
 */

typedef struct {
    qubit_t min_qubits;
    qubit_t max_qubits;
    int     warmup;
    int     iters;
    char    csv_path[256];
    int     specific_target;
    int     verbose;
} options_t;

static options_t parse_args(int argc, char *argv[]) {
    options_t o = {
        .min_qubits = (LOG_TILE_DIM < 2) ? 2 : LOG_TILE_DIM,
        .max_qubits = PROF_DEFAULT_MAX_QUBITS,
        .warmup     = PROF_DEFAULT_WARMUP,
        .iters      = PROF_DEFAULT_ITERS,
        .csv_path   = "",
        .specific_target = -1,
        .verbose    = 0,
    };

    static struct option long_opts[] = {
        {"min-qubits", required_argument, NULL, 'm'},
        {"max-qubits", required_argument, NULL, 'M'},
        {"warmup",     required_argument, NULL, 'w'},
        {"iters",      required_argument, NULL, 'i'},
        {"csv",        required_argument, NULL, 'c'},
        {"target",     required_argument, NULL, 't'},
        {"verbose",    no_argument,       NULL, 'v'},
        {"help",       no_argument,       NULL, 'h'},
        {NULL, 0, NULL, 0}
    };

    int ch;
    while ((ch = getopt_long(argc, argv, "m:M:w:i:c:t:vh", long_opts, NULL)) != -1) {
        switch (ch) {
            case 'm': o.min_qubits = (qubit_t)atoi(optarg); break;
            case 'M': o.max_qubits = (qubit_t)atoi(optarg); break;
            case 'w': o.warmup = atoi(optarg); break;
            case 'i': o.iters  = atoi(optarg); break;
            case 'c': strncpy(o.csv_path, optarg, sizeof(o.csv_path) - 1); break;
            case 't': o.specific_target = atoi(optarg); break;
            case 'v': o.verbose = 1; break;
            case 'h':
                printf("Usage: %s [--min-qubits N] [--max-qubits N] [--warmup N] "
                       "[--iters N] [--csv FILE] [--target N] [-v]\n", argv[0]);
                exit(0);
            default: break;
        }
    }

    if (o.min_qubits < 2) o.min_qubits = 2;
    return o;
}

/*
 * =====================================================================================================================
 * Main
 * =====================================================================================================================
 */

int main(int argc, char *argv[]) {
    options_t opts = parse_args(argc, argv);

    printf("=== channel_tiled_1q() vs h_tiled() profiling ===\n");
    printf("Tile config: TILE_DIM=%d, LOG_TILE_DIM=%d, TILE_SIZE=%d\n",
           TILE_DIM, LOG_TILE_DIM, TILE_SIZE);
    printf("Warmup: %d, Iterations: %d\n", opts.warmup, opts.iters);
    printf("Qubit range: %d to %d\n", opts.min_qubits, opts.max_qubits);
    printf("Channel: depolarizing (p=0.1)\n\n");

    FILE *csv = NULL;
    if (opts.csv_path[0]) {
        csv = fopen(opts.csv_path, "w");
        if (!csv) {
            fprintf(stderr, "Cannot open %s for writing\n", opts.csv_path);
            return 1;
        }
        prof_print_csv_header(csv);
    }

    /* Build depolarizing superoperator */
    kraus_t kraus = make_depolarizing(0.1);
    superop_t sop = kraus_to_superop(&kraus);

    /*
     * =====================================================================================================================
     * Per-target comparison: channel vs Hadamard
     * =====================================================================================================================
     */
    printf("=== Per-target comparison ===\n");
    prof_print_header();

    for (qubit_t n = opts.min_qubits; n <= opts.max_qubits; n++) {
        for (qubit_t t = 0; t < n; t++) {
            if (opts.specific_target >= 0 && t != (qubit_t)opts.specific_target)
                continue;

            prof_result_t rc = prof_run_channel_1q(
                channel_tiled_1q, sop.data, "ch_tiled", MIXED_TILED,
                n, t, opts.warmup, opts.iters);
            prof_print_result(&rc);
            if (csv) prof_print_csv_row(csv, &rc);

            prof_result_t rh = prof_run_1q(
                h_tiled, "h_tiled", MIXED_TILED,
                n, t, opts.warmup, opts.iters);
            prof_print_result(&rh);
            if (csv) prof_print_csv_row(csv, &rh);
        }
        printf("\n");
    }

    /*
     * =====================================================================================================================
     * Scaling analysis: fixed target, varying qubits
     * =====================================================================================================================
     */
    printf("\n=== Scaling: target=0 (within-tile) ===\n");
    prof_print_header();
    for (qubit_t n = opts.min_qubits; n <= opts.max_qubits; n++) {
        prof_result_t rc = prof_run_channel_1q(
            channel_tiled_1q, sop.data, "ch_tiled", MIXED_TILED,
            n, 0, opts.warmup, opts.iters);
        prof_print_result(&rc);

        prof_result_t rh = prof_run_1q(
            h_tiled, "h_tiled", MIXED_TILED,
            n, 0, opts.warmup, opts.iters);
        prof_print_result(&rh);
    }

    if (opts.max_qubits > LOG_TILE_DIM) {
        printf("\n=== Scaling: target=%d (cross-tile) ===\n", LOG_TILE_DIM);
        prof_print_header();
        for (qubit_t n = LOG_TILE_DIM + 1; n <= opts.max_qubits; n++) {
            prof_result_t rc = prof_run_channel_1q(
                channel_tiled_1q, sop.data, "ch_tiled", MIXED_TILED,
                n, LOG_TILE_DIM, opts.warmup, opts.iters);
            prof_print_result(&rc);

            prof_result_t rh = prof_run_1q(
                h_tiled, "h_tiled", MIXED_TILED,
                n, LOG_TILE_DIM, opts.warmup, opts.iters);
            prof_print_result(&rh);
        }
    }

    free(kraus.data);
    free(sop.data);

    if (csv) {
        fclose(csv);
        printf("\nCSV written to: %s\n", opts.csv_path);
    }

    printf("\n=== Profiling complete ===\n");
    return 0;
}
