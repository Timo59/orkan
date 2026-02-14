/**
 * @file profile_x_tiled.c
 * @brief Micro-profiling harness for x_tiled()
 *
 * Profiles the tiled Pauli-X gate across multiple qubit counts and target
 * qubits, separating within-tile (target < LOG_TILE_DIM) from cross-tile
 * (target >= LOG_TILE_DIM) code paths.
 *
 * Usage:
 *   ./profile_x_tiled [options]
 *
 * Options:
 *   --min-qubits N    Minimum qubit count (default: LOG_TILE_DIM)
 *   --max-qubits N    Maximum qubit count (default: 12)
 *   --warmup N        Warmup iterations (default: 200)
 *   --iters N         Measurement iterations (default: 2000)
 *   --csv FILE        Write CSV output to FILE
 *   --target N        Profile only target qubit N (default: all valid targets)
 *   --compare         Also profile x_packed() for comparison
 *   -v, --verbose     Print extra info (state sizes, tile config)
 */

#include "profile.h"
#include <string.h>
#include <getopt.h>

/* Forward declarations for backend functions (not in public gate.h) */
extern void x_tiled(state_t *state, const qubit_t target);
extern void x_packed(state_t *state, const qubit_t target);

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
    int     specific_target;  /* -1 = all targets */
    int     compare;
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
        .compare    = 0,
        .verbose    = 0,
    };

    static struct option long_opts[] = {
        {"min-qubits", required_argument, NULL, 'm'},
        {"max-qubits", required_argument, NULL, 'M'},
        {"warmup",     required_argument, NULL, 'w'},
        {"iters",      required_argument, NULL, 'i'},
        {"csv",        required_argument, NULL, 'c'},
        {"target",     required_argument, NULL, 't'},
        {"compare",    no_argument,       NULL, 'C'},
        {"verbose",    no_argument,       NULL, 'v'},
        {"help",       no_argument,       NULL, 'h'},
        {NULL, 0, NULL, 0}
    };

    int ch;
    while ((ch = getopt_long(argc, argv, "m:M:w:i:c:t:Cvh", long_opts, NULL)) != -1) {
        switch (ch) {
            case 'm': o.min_qubits = (qubit_t)atoi(optarg); break;
            case 'M': o.max_qubits = (qubit_t)atoi(optarg); break;
            case 'w': o.warmup = atoi(optarg); break;
            case 'i': o.iters  = atoi(optarg); break;
            case 'c': strncpy(o.csv_path, optarg, sizeof(o.csv_path) - 1); break;
            case 't': o.specific_target = atoi(optarg); break;
            case 'C': o.compare = 1; break;
            case 'v': o.verbose = 1; break;
            case 'h':
                printf("Usage: %s [--min-qubits N] [--max-qubits N] [--warmup N] "
                       "[--iters N] [--csv FILE] [--target N] [--compare] [-v]\n", argv[0]);
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

    printf("=== x_tiled() profiling ===\n");
    printf("Tile config: TILE_DIM=%d, LOG_TILE_DIM=%d, TILE_SIZE=%d\n",
           TILE_DIM, LOG_TILE_DIM, TILE_SIZE);
    printf("Warmup: %d, Iterations: %d\n", opts.warmup, opts.iters);
    printf("Qubit range: %d to %d\n\n", opts.min_qubits, opts.max_qubits);

    FILE *csv = NULL;
    if (opts.csv_path[0]) {
        csv = fopen(opts.csv_path, "w");
        if (!csv) {
            fprintf(stderr, "Cannot open %s for writing\n", opts.csv_path);
            return 1;
        }
        prof_print_csv_header(csv);
    }

    /* Print table header */
    prof_print_header();

    for (qubit_t n = opts.min_qubits; n <= opts.max_qubits; n++) {
        /* Tiled storage requires at least LOG_TILE_DIM qubits for a full tile.
         * For smaller n, dim < TILE_DIM and the within-tile path still works
         * but cross-tile paths don't exist. */
        qubit_t max_target = n - 1;

        if (opts.verbose) {
            uint64_t dim = (uint64_t)1 << n;
            uint64_t n_tiles = (dim + TILE_DIM - 1) / TILE_DIM;
            printf("\n--- %d qubits: dim=%llu, n_tiles=%llu, "
                   "data=%.1f KiB ---\n",
                   n, (unsigned long long)dim, (unsigned long long)n_tiles,
                   (double)prof_tiled_bytes(n) / 2048.0);
        }

        for (qubit_t t = 0; t <= max_target; t++) {
            if (opts.specific_target >= 0 && t != (qubit_t)opts.specific_target)
                continue;

            /* Profile x_tiled */
            prof_result_t r = prof_run_1q(
                x_tiled, "x_tiled", MIXED_TILED,
                n, t, opts.warmup, opts.iters);
            prof_print_result(&r);
            if (csv) prof_print_csv_row(csv, &r);

            /* Optionally compare with x_packed */
            if (opts.compare) {
                prof_result_t rp = prof_run_1q(
                    x_packed, "x_packed", MIXED_PACKED,
                    n, t, opts.warmup, opts.iters);
                prof_print_result(&rp);
                if (csv) prof_print_csv_row(csv, &rp);
            }
        }
    }

    if (csv) {
        fclose(csv);
        printf("\nCSV written to: %s\n", opts.csv_path);
    }

    /*
     * =====================================================================================================================
     * Scaling analysis: fix target, vary qubits
     * =====================================================================================================================
     */
    printf("\n\n=== Scaling analysis (fixed target, varying qubits) ===\n");
    prof_print_header();

    /* Within-tile target (target=0, always within-tile) */
    for (qubit_t n = opts.min_qubits; n <= opts.max_qubits; n++) {
        prof_result_t r = prof_run_1q(
            x_tiled, "x_tiled", MIXED_TILED,
            n, 0, opts.warmup, opts.iters);
        prof_print_result(&r);
        if (csv) prof_print_csv_row(csv, &r);
    }

    /* Cross-tile target (target=LOG_TILE_DIM, when n > LOG_TILE_DIM) */
    if (opts.max_qubits > LOG_TILE_DIM) {
        printf("\n");
        for (qubit_t n = LOG_TILE_DIM + 1; n <= opts.max_qubits; n++) {
            prof_result_t r = prof_run_1q(
                x_tiled, "x_tiled", MIXED_TILED,
                n, LOG_TILE_DIM, opts.warmup, opts.iters);
            prof_print_result(&r);
            if (csv) prof_print_csv_row(csv, &r);
        }
    }

    printf("\n=== Profiling complete ===\n");

    /*
     * =====================================================================================================================
     * Suggested next steps (printed as a mini-report)
     * =====================================================================================================================
     */
    printf("\n");
    printf("== Recommended profiling tools for deeper analysis ==\n\n");
    printf("1. Apple Instruments - Time Profiler\n");
    printf("   Record a trace to find hotspots within x_tiled():\n");
    printf("     xcrun xctrace record --template 'Time Profiler' \\\n");
    printf("       --launch -- ./profile_x_tiled --min-qubits 10 --max-qubits 10 --target 0 --iters 50000\n\n");
    printf("2. Apple Instruments - Counters (CPU Counters)\n");
    printf("   Measure cache misses, branch mispredictions, IPC on Apple Silicon:\n");
    printf("     xcrun xctrace record --template 'CPU Counters' \\\n");
    printf("       --launch -- ./profile_x_tiled --min-qubits 10 --max-qubits 10 --target 0 --iters 50000\n\n");
    printf("3. Compiler vectorization report\n");
    printf("   Rebuild with optimization remarks to see what the compiler vectorized:\n");
    printf("     cmake --build . -- CFLAGS=\"-Rpass=loop-vectorize -Rpass-missed=loop-vectorize "
           "-Rpass-analysis=loop-vectorize\"\n\n");
    printf("4. LLVM Machine Code Analyzer (llvm-mca)\n");
    printf("   Analyze the inner loop's instruction throughput and bottlenecks:\n");
    printf("     Extract assembly, then: llvm-mca -mcpu=apple-m1 inner_loop.s\n\n");
    printf("5. sample (quick stack sampling)\n");
    printf("   Quick check of where time is spent (no recompilation needed):\n");
    printf("     sample <pid> 5 -file sample_output.txt\n\n");

    return 0;
}
