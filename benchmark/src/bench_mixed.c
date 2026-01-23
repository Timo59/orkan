/**
 * @file bench_mixed.c
 * @brief Benchmark comparing qlib's mixed state gate implementations against alternatives
 *
 * Compares:
 * 1. qlib packed sparse - Production implementation with packed lower-triangular storage
 * 2. BLAS dense - Full matrix UρU† via zgemm
 * 3. Naive loop - Element-by-element computation
 */

#include "bench.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <math.h>

#if defined(__APPLE__)
    #include <vecLib/cblas_new.h>
#elif defined(__linux__)
    #include <cblas.h>
#endif

/*
 * =====================================================================================================================
 * Gate matrices (2x2, column-major)
 * =====================================================================================================================
 */

#define INVSQRT2 0.7071067811865475

static const cplx_t XMAT[4] = {0.0, 1.0, 1.0, 0.0};
static const cplx_t YMAT[4] = {0.0, -I, I, 0.0};
static const cplx_t ZMAT[4] = {1.0, 0.0, 0.0, -1.0};
static const cplx_t HMAT[4] = {INVSQRT2, INVSQRT2, INVSQRT2, -INVSQRT2};

/*
 * =====================================================================================================================
 * Helper functions
 * =====================================================================================================================
 */

/** @brief Initialize a random mixed state (density matrix) */
static void init_random_mixed_state(state_t *state, qubit_t qubits) {
    state->type = MIXED;
    state->qubits = qubits;

    dim_t dim = (dim_t)1 << qubits;
    dim_t len = dim * (dim + 1) / 2;

    state->data = malloc(len * sizeof(cplx_t));
    if (!state->data) {
        fprintf(stderr, "init_random_mixed_state: allocation failed\n");
        return;
    }

    /* Fill with random values (not a valid density matrix, but fine for benchmarking) */
    for (dim_t i = 0; i < len; ++i) {
        double re = (double)rand() / RAND_MAX - 0.5;
        double im = (double)rand() / RAND_MAX - 0.5;
        state->data[i] = re + im * I;
    }
}

/** @brief Kronecker product: C = A ⊗ B */
static cplx_t* kron(dim_t ma, dim_t na, const cplx_t *A,
                    dim_t mb, dim_t nb, const cplx_t *B) {
    dim_t mc = ma * mb;
    dim_t nc = na * nb;
    cplx_t *C = malloc(mc * nc * sizeof(cplx_t));
    if (!C) return NULL;

    for (dim_t ja = 0; ja < na; ++ja) {
        for (dim_t ia = 0; ia < ma; ++ia) {
            cplx_t aij = A[ia + ja * ma];
            for (dim_t jb = 0; jb < nb; ++jb) {
                for (dim_t ib = 0; ib < mb; ++ib) {
                    dim_t ic = ib + ia * mb;
                    dim_t jc = jb + ja * nb;
                    C[ic + jc * mc] = aij * B[ib + jb * mb];
                }
            }
        }
    }
    return C;
}

/** @brief Build full gate matrix I⊗...⊗U⊗...⊗I for single-qubit gate */
static cplx_t* build_full_gate_matrix(qubit_t qubits, const cplx_t gate[4], qubit_t target) {
    dim_t dim = (dim_t)1 << qubits;

    /* Start with identity for qubits above target */
    qubit_t pos = qubits - 1 - target;
    dim_t dim_left = (dim_t)1 << pos;

    /* Identity matrix for left qubits */
    cplx_t *id_left = calloc(dim_left * dim_left, sizeof(cplx_t));
    if (!id_left) return NULL;
    for (dim_t i = 0; i < dim_left; ++i) id_left[i + i * dim_left] = 1.0;

    /* Kronecker with gate */
    cplx_t *tmp = kron(dim_left, dim_left, id_left, 2, 2, gate);
    free(id_left);
    if (!tmp) return NULL;

    dim_left *= 2;

    /* Identity for right qubits */
    dim_t dim_right = (dim_t)1 << target;
    cplx_t *id_right = calloc(dim_right * dim_right, sizeof(cplx_t));
    if (!id_right) { free(tmp); return NULL; }
    for (dim_t i = 0; i < dim_right; ++i) id_right[i + i * dim_right] = 1.0;

    /* Final Kronecker product */
    cplx_t *U = kron(dim_left, dim_left, tmp, dim_right, dim_right, id_right);
    free(tmp);
    free(id_right);

    return U;
}

/** @brief Unpack lower-triangular to full Hermitian matrix */
static cplx_t* unpack_density(dim_t dim, const cplx_t *packed) {
    cplx_t *full = malloc(dim * dim * sizeof(cplx_t));
    if (!full) return NULL;

    dim_t k = 0;
    for (dim_t j = 0; j < dim; ++j) {
        for (dim_t i = j; i < dim; ++i) {
            full[i + j * dim] = packed[k];
            full[j + i * dim] = conj(packed[k]);
            ++k;
        }
    }
    return full;
}

/** @brief Pack full Hermitian matrix to lower-triangular */
static void pack_density(dim_t dim, const cplx_t *full, cplx_t *packed) {
    dim_t k = 0;
    for (dim_t j = 0; j < dim; ++j) {
        for (dim_t i = j; i < dim; ++i) {
            packed[k++] = full[i + j * dim];
        }
    }
}

/** @brief Compute UρU† using BLAS zgemm */
static void apply_uruh_blas(dim_t dim, const cplx_t *U, cplx_t *rho) {
    cplx_t *tmp = malloc(dim * dim * sizeof(cplx_t));
    if (!tmp) return;

    const cplx_t alpha = 1.0;
    const cplx_t beta = 0.0;

    /* tmp = U * rho */
    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                dim, dim, dim, &alpha, U, dim, rho, dim, &beta, tmp, dim);

    /* rho = tmp * U† */
    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasConjTrans,
                dim, dim, dim, &alpha, tmp, dim, U, dim, &beta, rho, dim);

    free(tmp);
}

/** @brief Compute UρU† using naive loops */
static void apply_uruh_naive(dim_t dim, const cplx_t *U, cplx_t *rho) {
    cplx_t *tmp = calloc(dim * dim, sizeof(cplx_t));
    cplx_t *out = calloc(dim * dim, sizeof(cplx_t));
    if (!tmp || !out) { free(tmp); free(out); return; }

    /* tmp = U * rho */
    for (dim_t j = 0; j < dim; ++j) {
        for (dim_t i = 0; i < dim; ++i) {
            cplx_t sum = 0.0;
            for (dim_t k = 0; k < dim; ++k) {
                sum += U[i + k * dim] * rho[k + j * dim];
            }
            tmp[i + j * dim] = sum;
        }
    }

    /* out = tmp * U† */
    for (dim_t j = 0; j < dim; ++j) {
        for (dim_t i = 0; i < dim; ++i) {
            cplx_t sum = 0.0;
            for (dim_t k = 0; k < dim; ++k) {
                sum += tmp[i + k * dim] * conj(U[j + k * dim]);
            }
            out[i + j * dim] = sum;
        }
    }

    memcpy(rho, out, dim * dim * sizeof(cplx_t));
    free(tmp);
    free(out);
}

/*
 * =====================================================================================================================
 * Benchmark implementations
 * =====================================================================================================================
 */

bench_result_t bench_qlib_packed(qubit_t qubits, const char *gate_name,
                                  qs_error_t (*gate_fn)(state_t*, qubit_t),
                                  int iterations, int warmup) {
    bench_result_t result = {0};
    result.qubits = qubits;
    result.dim = (dim_t)1 << qubits;
    result.gate_name = gate_name;
    result.method = "qlib_packed";
    result.iterations = iterations;

    /* Measure actual memory usage */
    size_t mem_before = bench_get_rss();
    state_t state = {0};
    init_random_mixed_state(&state, qubits);
    if (!state.data) return result;
    size_t mem_after = bench_get_rss();
    result.memory_bytes = (mem_after > mem_before) ? (mem_after - mem_before) : bench_packed_size(qubits);

    /* Warm-up */
    for (int i = 0; i < warmup; ++i) {
        gate_fn(&state, 0);
    }

    /* Timed iterations */
    uint64_t start = bench_time_ns();
    for (int i = 0; i < iterations; ++i) {
        for (qubit_t t = 0; t < qubits; ++t) {
            gate_fn(&state, t);
        }
    }
    uint64_t end = bench_time_ns();

    result.time_ms = bench_ns_to_ms(end - start);
    result.ops_per_sec = (double)(iterations * qubits) / (result.time_ms / 1000.0);

    state_free(&state);
    return result;
}

bench_result_t bench_blas_dense(qubit_t qubits, const char *gate_name,
                                 const cplx_t gate_mat[4],
                                 int iterations, int warmup) {
    bench_result_t result = {0};
    result.qubits = qubits;
    result.dim = (dim_t)1 << qubits;
    result.gate_name = gate_name;
    result.method = "blas_dense";
    result.iterations = iterations;

    dim_t dim = result.dim;

    /* Measure actual memory usage */
    size_t mem_before = bench_get_rss();

    /* Create full density matrix (random) */
    cplx_t *rho = malloc(dim * dim * sizeof(cplx_t));
    if (!rho) return result;
    for (dim_t i = 0; i < dim * dim; ++i) {
        rho[i] = (double)rand() / RAND_MAX + ((double)rand() / RAND_MAX) * I;
    }

    /* Build full gate matrices for each target */
    cplx_t **U = malloc(qubits * sizeof(cplx_t*));
    if (!U) { free(rho); return result; }
    for (qubit_t t = 0; t < qubits; ++t) {
        U[t] = build_full_gate_matrix(qubits, gate_mat, t);
    }

    size_t mem_after = bench_get_rss();
    result.memory_bytes = (mem_after > mem_before) ? (mem_after - mem_before) : bench_dense_size(qubits);

    /* Warm-up */
    for (int i = 0; i < warmup; ++i) {
        apply_uruh_blas(dim, U[0], rho);
    }

    /* Timed iterations */
    uint64_t start = bench_time_ns();
    for (int i = 0; i < iterations; ++i) {
        for (qubit_t t = 0; t < qubits; ++t) {
            apply_uruh_blas(dim, U[t], rho);
        }
    }
    uint64_t end = bench_time_ns();

    result.time_ms = bench_ns_to_ms(end - start);
    result.ops_per_sec = (double)(iterations * qubits) / (result.time_ms / 1000.0);

    /* Cleanup */
    free(rho);
    for (qubit_t t = 0; t < qubits; ++t) free(U[t]);
    free(U);

    return result;
}

bench_result_t bench_naive_loop(qubit_t qubits, const char *gate_name,
                                 const cplx_t gate_mat[4],
                                 int iterations, int warmup) {
    bench_result_t result = {0};
    result.qubits = qubits;
    result.dim = (dim_t)1 << qubits;
    result.gate_name = gate_name;
    result.method = "naive_loop";
    result.iterations = iterations;

    dim_t dim = result.dim;

    /* Measure actual memory usage */
    size_t mem_before = bench_get_rss();

    /* Create full density matrix (random) */
    cplx_t *rho = malloc(dim * dim * sizeof(cplx_t));
    if (!rho) return result;
    for (dim_t i = 0; i < dim * dim; ++i) {
        rho[i] = (double)rand() / RAND_MAX + ((double)rand() / RAND_MAX) * I;
    }

    /* Build full gate matrices for each target */
    cplx_t **U = malloc(qubits * sizeof(cplx_t*));
    if (!U) { free(rho); return result; }
    for (qubit_t t = 0; t < qubits; ++t) {
        U[t] = build_full_gate_matrix(qubits, gate_mat, t);
    }

    size_t mem_after = bench_get_rss();
    result.memory_bytes = (mem_after > mem_before) ? (mem_after - mem_before) : bench_dense_size(qubits);

    /* Warm-up */
    for (int i = 0; i < warmup; ++i) {
        apply_uruh_naive(dim, U[0], rho);
    }

    /* Timed iterations */
    uint64_t start = bench_time_ns();
    for (int i = 0; i < iterations; ++i) {
        for (qubit_t t = 0; t < qubits; ++t) {
            apply_uruh_naive(dim, U[t], rho);
        }
    }
    uint64_t end = bench_time_ns();

    result.time_ms = bench_ns_to_ms(end - start);
    result.ops_per_sec = (double)(iterations * qubits) / (result.time_ms / 1000.0);

    /* Cleanup */
    free(rho);
    for (qubit_t t = 0; t < qubits; ++t) free(U[t]);
    free(U);

    return result;
}

/*
 * =====================================================================================================================
 * Output functions
 * =====================================================================================================================
 */

void bench_print_result(const bench_result_t *result, int verbose) {
    printf("  %-14s %10.3f ms (%'12.0f ops/sec)",
           result->method, result->time_ms, result->ops_per_sec);
    if (verbose) {
        printf("  [%zu bytes]", result->memory_bytes);
    }
    printf("\n");
}

void bench_print_csv_header(void) {
    printf("qubits,dim,gate,method,time_ms,ops_per_sec,memory_bytes,iterations\n");
}

void bench_print_csv(const bench_result_t *result) {
    printf("%u,%d,%s,%s,%.6f,%.0f,%zu,%d\n",
           result->qubits, result->dim, result->gate_name, result->method,
           result->time_ms, result->ops_per_sec, result->memory_bytes, result->iterations);
}

/*
 * =====================================================================================================================
 * Command-line parsing
 * =====================================================================================================================
 */

static void print_usage(const char *prog) {
    printf("Usage: %s [options]\n\n", prog);
    printf("Options:\n");
    printf("  --min-qubits N   Minimum qubit count (default: %d)\n", BENCH_DEFAULT_MIN_QUBITS);
    printf("  --max-qubits N   Maximum qubit count (default: %d)\n", BENCH_DEFAULT_MAX_QUBITS);
    printf("  --iterations N   Number of iterations (default: %d)\n", BENCH_DEFAULT_ITERATIONS);
    printf("  --warmup N       Warm-up iterations (default: %d)\n", BENCH_DEFAULT_WARMUP);
    printf("  --csv            Output in CSV format\n");
    printf("  --pgfplots       Output in pgfplots-compatible .dat format\n");
    printf("  --verbose        Show additional details\n");
    printf("  --help           Show this help\n");
}

bench_options_t bench_parse_options(int argc, char *argv[]) {
    bench_options_t opts = {
        .min_qubits = BENCH_DEFAULT_MIN_QUBITS,
        .max_qubits = BENCH_DEFAULT_MAX_QUBITS,
        .iterations = BENCH_DEFAULT_ITERATIONS,
        .warmup = BENCH_DEFAULT_WARMUP,
        .csv_output = 0,
        .pgfplots_output = 0,
        .verbose = 0
    };

    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "--min-qubits") == 0 && i + 1 < argc) {
            opts.min_qubits = (qubit_t)atoi(argv[++i]);
        } else if (strcmp(argv[i], "--max-qubits") == 0 && i + 1 < argc) {
            opts.max_qubits = (qubit_t)atoi(argv[++i]);
        } else if (strcmp(argv[i], "--iterations") == 0 && i + 1 < argc) {
            opts.iterations = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--warmup") == 0 && i + 1 < argc) {
            opts.warmup = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--csv") == 0) {
            opts.csv_output = 1;
        } else if (strcmp(argv[i], "--pgfplots") == 0) {
            opts.pgfplots_output = 1;
        } else if (strcmp(argv[i], "--verbose") == 0) {
            opts.verbose = 1;
        } else if (strcmp(argv[i], "--help") == 0) {
            print_usage(argv[0]);
            exit(0);
        }
    }

    return opts;
}

/*
 * =====================================================================================================================
 * Gate definitions for benchmarking
 * =====================================================================================================================
 */

typedef struct gate_def {
    const char *name;
    qs_error_t (*fn)(state_t*, qubit_t);
    const cplx_t *mat;
} gate_def_t;

static const gate_def_t gates[] = {
    {"X", x, XMAT},
    {"H", h, HMAT},
    {"Z", z, ZMAT},
    {NULL, NULL, NULL}
};

#define NUM_GATES 3
#define MAX_QUBIT_CONFIGS 16
#define NUM_METHODS 5  /* qlib, blas, naive, quest, qpp */

/** @brief Storage for pgfplots output (collected during benchmark) */
typedef struct {
    qubit_t qubits[MAX_QUBIT_CONFIGS];
    int num_configs;
    double time_ms[NUM_GATES][MAX_QUBIT_CONFIGS][NUM_METHODS];
    size_t memory[NUM_GATES][MAX_QUBIT_CONFIGS][NUM_METHODS];
} pgfplots_data_t;

static pgfplots_data_t pgf_data = {0};

/** @brief Method indices for pgfplots data */
enum { M_QLIB = 0, M_BLAS = 1, M_NAIVE = 2, M_QUEST = 3, M_QPP = 4 };

/** @brief Print pgfplots-compatible .dat output */
static void print_pgfplots_output(void) {
    const char *method_names[] = {"qlib", "blas", "quest", "qpp"};
    int method_indices[] = {M_QLIB, M_BLAS, M_QUEST, M_QPP};
    int num_methods = 4;

    /* Print timing tables (one per gate) */
    for (int g = 0; g < NUM_GATES; ++g) {
        printf("# %s gate timing (ms)\n", gates[g].name);
        printf("%-8s", "qubits");
        for (int m = 0; m < num_methods; ++m) {
            printf("  %-12s", method_names[m]);
        }
        printf("\n");

        for (int q = 0; q < pgf_data.num_configs; ++q) {
            printf("%-8u", pgf_data.qubits[q]);
            for (int m = 0; m < num_methods; ++m) {
                double val = pgf_data.time_ms[g][q][method_indices[m]];
                if (val > 0) {
                    printf("  %-12.6f", val);
                } else {
                    printf("  %-12s", "nan");
                }
            }
            printf("\n");
        }
        printf("\n");
    }

    /* Print memory table (same for all gates, use first gate) */
    printf("# Memory usage (MB)\n");
    printf("%-8s", "qubits");
    for (int m = 0; m < num_methods; ++m) {
        printf("  %-12s", method_names[m]);
    }
    printf("\n");

    for (int q = 0; q < pgf_data.num_configs; ++q) {
        printf("%-8u", pgf_data.qubits[q]);
        for (int m = 0; m < num_methods; ++m) {
            size_t bytes = pgf_data.memory[0][q][method_indices[m]];
            if (bytes > 0) {
                printf("  %-12.3f", (double)bytes / (1024.0 * 1024.0));
            } else {
                printf("  %-12s", "nan");
            }
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

    /* Seed random number generator for reproducibility */
    srand(42);

    if (opts.csv_output) {
        bench_print_csv_header();
    } else if (!opts.pgfplots_output) {
        printf("\nMixed State Gate Benchmarks\n");
        printf("===========================\n");
        printf("Comparing: qlib packed sparse vs BLAS dense vs naive loop");
#ifdef WITH_QUEST
        printf(" vs QuEST");
#endif
#ifdef WITH_QPP
        printf(" vs Quantum++");
#endif
        printf("\n");
        printf("Iterations: %d, Warm-up: %d\n", opts.iterations, opts.warmup);
#ifdef WITH_QUEST
        printf("QuEST: enabled (https://quest.qtechtheory.org/)\n");
#endif
#ifdef WITH_QPP
        printf("Quantum++: enabled (https://github.com/softwareQinc/qpp)\n");
#endif
        printf("\n");
    }

    int qi = 0;  /* Qubit configuration index for pgfplots */

    for (qubit_t qubits = opts.min_qubits; qubits <= opts.max_qubits; qubits += 2) {
        dim_t dim = (dim_t)1 << qubits;

        /* Skip if memory would be too large (> 4GB for dense) */
        if (bench_dense_size(qubits) > 4ULL * 1024 * 1024 * 1024) {
            if (!opts.csv_output && !opts.pgfplots_output) {
                printf("Skipping %u qubits (dense matrix > 4GB)\n", qubits);
            }
            continue;
        }

        /* Store qubit count for pgfplots */
        if (opts.pgfplots_output && qi < MAX_QUBIT_CONFIGS) {
            pgf_data.qubits[qi] = qubits;
        }

        if (!opts.csv_output && !opts.pgfplots_output) {
            printf("qubits: %u, dim: %ld\n", qubits, dim);
        }

        for (int g = 0; gates[g].name != NULL; ++g) {
            if (!opts.csv_output && !opts.pgfplots_output) {
                printf("\nGate: %s\n", gates[g].name);
            }

            /* Benchmark qlib packed */
            bench_result_t r_packed = bench_qlib_packed(
                qubits, gates[g].name, gates[g].fn,
                opts.iterations, opts.warmup);

            /* Benchmark BLAS dense (skip for larger systems - O(N³) is too slow) */
            bench_result_t r_dense = {0};
            if (qubits <= 8) {
                r_dense = bench_blas_dense(
                    qubits, gates[g].name, gates[g].mat,
                    opts.iterations, opts.warmup);
            }

            /* Benchmark naive loop (skip for larger systems) */
            bench_result_t r_naive = {0};
            if (qubits <= 8) {
                r_naive = bench_naive_loop(
                    qubits, gates[g].name, gates[g].mat,
                    opts.iterations / 10, opts.warmup / 10);  /* Fewer iterations for slow naive */
            }

#ifdef WITH_QUEST
            /* Benchmark QuEST framework */
            bench_result_t r_quest = bench_quest(
                qubits, gates[g].name,
                opts.iterations, opts.warmup);
#endif

#ifdef WITH_QPP
            /* Benchmark Quantum++ framework */
            bench_result_t r_qpp = bench_qpp(
                qubits, gates[g].name,
                opts.iterations, opts.warmup);
#endif

            /* Store results for pgfplots output */
            if (opts.pgfplots_output && qi < MAX_QUBIT_CONFIGS) {
                pgf_data.time_ms[g][qi][M_QLIB] = r_packed.time_ms;
                pgf_data.time_ms[g][qi][M_BLAS] = r_dense.time_ms;
                pgf_data.time_ms[g][qi][M_NAIVE] = r_naive.time_ms;
                pgf_data.memory[g][qi][M_QLIB] = r_packed.memory_bytes;
                pgf_data.memory[g][qi][M_BLAS] = r_dense.memory_bytes;
                pgf_data.memory[g][qi][M_NAIVE] = r_naive.memory_bytes;
#ifdef WITH_QUEST
                pgf_data.time_ms[g][qi][M_QUEST] = r_quest.time_ms;
                pgf_data.memory[g][qi][M_QUEST] = r_quest.memory_bytes;
#endif
#ifdef WITH_QPP
                pgf_data.time_ms[g][qi][M_QPP] = r_qpp.time_ms;
                pgf_data.memory[g][qi][M_QPP] = r_qpp.memory_bytes;
#endif
            }

            if (opts.csv_output) {
                bench_print_csv(&r_packed);
                if (qubits <= 8) bench_print_csv(&r_dense);
                if (qubits <= 8) bench_print_csv(&r_naive);
#ifdef WITH_QUEST
                bench_print_csv(&r_quest);
#endif
#ifdef WITH_QPP
                bench_print_csv(&r_qpp);
#endif
            } else if (!opts.pgfplots_output) {
                bench_print_result(&r_packed, opts.verbose);
                if (qubits <= 8) bench_print_result(&r_dense, opts.verbose);
                if (qubits <= 8) bench_print_result(&r_naive, opts.verbose);
#ifdef WITH_QUEST
                bench_print_result(&r_quest, opts.verbose);
#endif
#ifdef WITH_QPP
                bench_print_result(&r_qpp, opts.verbose);
#endif

                /* Speedup */
                printf("  Speedup:");
                int has_speedup = 0;
                if (qubits <= 8 && r_dense.time_ms > 0) {
                    double speedup_dense = r_dense.time_ms / r_packed.time_ms;
                    printf(" %.1fx vs dense", speedup_dense);
                    has_speedup = 1;
                }
                if (qubits <= 8 && r_naive.time_ms > 0) {
                    double speedup_naive = r_naive.time_ms / r_packed.time_ms;
                    /* Scale for different iteration counts */
                    speedup_naive *= 10;
                    printf("%s~%.0fx vs naive", has_speedup ? ", " : " ", speedup_naive);
                    has_speedup = 1;
                }
#ifdef WITH_QUEST
                if (r_quest.time_ms > 0) {
                    double speedup_quest = r_quest.time_ms / r_packed.time_ms;
                    printf(", %.1fx vs QuEST", speedup_quest);
                }
#endif
#ifdef WITH_QPP
                if (r_qpp.time_ms > 0) {
                    double speedup_qpp = r_qpp.time_ms / r_packed.time_ms;
                    printf(", %.1fx vs Quantum++", speedup_qpp);
                }
#endif
                printf("\n");

                /* Memory comparison (actual RSS measured) */
                size_t max_mem = r_packed.memory_bytes;
                if (r_dense.memory_bytes > max_mem) max_mem = r_dense.memory_bytes;
                if (r_naive.memory_bytes > max_mem) max_mem = r_naive.memory_bytes;
#ifdef WITH_QUEST
                if (r_quest.memory_bytes > max_mem) max_mem = r_quest.memory_bytes;
#endif
#ifdef WITH_QPP
                if (r_qpp.memory_bytes > max_mem) max_mem = r_qpp.memory_bytes;
#endif

                /* Auto-scale units based on largest value */
                const char *unit;
                double scale;
                if (max_mem >= 1024 * 1024) {
                    unit = "MB";
                    scale = 1024.0 * 1024.0;
                } else {
                    unit = "KB";
                    scale = 1024.0;
                }

                printf("  Memory: qlib=%.1f %s", (double)r_packed.memory_bytes / scale, unit);
                if (qubits <= 8 && r_dense.memory_bytes > 0) {
                    printf(", dense=%.1f %s", (double)r_dense.memory_bytes / scale, unit);
                }
                if (qubits <= 8 && r_naive.memory_bytes > 0) {
                    printf(", naive=%.1f %s", (double)r_naive.memory_bytes / scale, unit);
                }
#ifdef WITH_QUEST
                if (r_quest.memory_bytes > 0) {
                    printf(", QuEST=%.1f %s", (double)r_quest.memory_bytes / scale, unit);
                }
#endif
#ifdef WITH_QPP
                if (r_qpp.memory_bytes > 0) {
                    printf(", Qpp=%.1f %s", (double)r_qpp.memory_bytes / scale, unit);
                }
#endif
                if (qubits <= 8 && r_packed.memory_bytes > 0 && r_dense.memory_bytes > 0) {
                    double savings = 100.0 * (1.0 - (double)r_packed.memory_bytes / (double)r_dense.memory_bytes);
                    printf(" (qlib saves %.0f%%)", savings);
                }
                printf("\n");
            }
        }

        if (!opts.csv_output && !opts.pgfplots_output) {
            printf("\n" "─────────────────────────────────────────────────\n");
        }

        /* Increment qubit config index for pgfplots */
        if (opts.pgfplots_output) {
            qi++;
        }
    }

    /* Output pgfplots data */
    if (opts.pgfplots_output) {
        pgf_data.num_configs = qi;
        print_pgfplots_output();
    }

#ifdef WITH_QUEST
    bench_quest_cleanup();
#endif

    return 0;
}
