/**
 * @file bench.h
 * @brief Shared types, constants, and function declarations for the gate benchmark.
 *
 * Defines the gate enum, backend vtable, CLI options, position descriptors,
 * timing utilities, and output formatter declarations.
 */

#ifndef BENCH_H
#define BENCH_H

#include "bench_result.h"
#include "q_types.h"
#include "state.h"
#include "gate.h"

#include <time.h>
#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/* =====================================================================
 * Gate identification
 * ===================================================================== */

typedef enum {
    BG_X, BG_Y, BG_Z, BG_H, BG_S, BG_T,
    BG_RX, BG_RY, BG_RZ,
    BG_CX, BG_CY, BG_CZ,
    BG_CCX,
    BG_COUNT
} bench_gate_id_t;

typedef struct {
    bench_gate_id_t id;
    const char     *name;         /**< Canonical display name (e.g. "X", "CNOT") */
    int             n_qubits;     /**< 1, 2, or 3 */
    int             parametrised; /**< 1 if --par is required */
} bench_gate_info_t;

/** Gate info table — defined in bench_main.c */
extern const bench_gate_info_t bench_gate_table[BG_COUNT];

/**
 * @brief Resolve a gate name string to a gate id.
 *
 * Case-insensitive.  Supports aliases: CX/CNOT, CCX/TOFFOLI.
 * Returns BG_COUNT if the name is not recognised.
 */
bench_gate_id_t bench_resolve_gate(const char *name);

/* =====================================================================
 * Backend identification
 * ===================================================================== */

typedef enum {
    BACKEND_QLIB_TILED,
    BACKEND_QLIB_PACKED,
    BACKEND_BLAS,
    BACKEND_QUEST,
    BACKEND_QULACS,
    BACKEND_AER,
    BACKEND_COUNT
} bench_backend_id_t;

/**
 * @brief Backend vtable.
 *
 * Each backend implements init/apply/cleanup/touch_bytes/available.
 * The timing loop in bench_run.c calls these through the vtable,
 * making it fully backend-agnostic.
 */
typedef struct bench_backend {
    const char *name;

    /** Allocate state, resolve gate function pointer into context. */
    void *(*init)(qubit_t qubits, bench_gate_id_t gate, double par);

    /** Apply one gate at the given position (gate already bound in context). */
    void  (*apply)(void *ctx, const qubit_t *pos);

    /** Free the state and context. */
    void  (*cleanup)(void *ctx);

    /** Bytes read+written per gate application (storage-specific). */
    size_t (*touch_bytes)(qubit_t qubits);

    /** Return 1 if this backend can run at the given qubit count. */
    int   (*available)(qubit_t qubits);
} bench_backend_t;

/** Backend vtable array — populated by init_backends() in bench_main.c */
extern bench_backend_t bench_backends[BACKEND_COUNT];

/** Individual backend vtables — defined in their respective source files */
extern const bench_backend_t qlib_tiled_backend;
extern const bench_backend_t qlib_packed_backend;
extern const bench_backend_t blas_backend;
#ifdef HAVE_QUEST
extern const bench_backend_t quest_backend;
#endif
#ifdef HAVE_QULACS
extern const bench_backend_t qulacs_backend;
#endif
#ifdef HAVE_AER
extern const bench_backend_t aer_backend;
#endif

/* =====================================================================
 * CLI options
 * ===================================================================== */

typedef enum {
    FMT_CONSOLE,
    FMT_CSV,
    FMT_PGFPLOTS
} bench_fmt_t;

typedef struct {
    bench_gate_id_t gate;
    double          par;
    int             par_provided;
    qubit_t         min_qubits;   /**< Default: 2 */
    qubit_t         max_qubits;   /**< Default: 8 */
    int             samples;      /**< Default: 10 */
    int             iterations;   /**< Default: 100; auto-tuned per qubit count if not set */
    int             iterations_explicit; /**< 1 if --iterations was given on CLI */
    int             warm_up;      /**< Default: 3 */
    int             per_qubit;    /**< 0 or 1 */
    bench_fmt_t     output;
} bench_options_t;

/** @brief Parse command-line arguments into options.  Exits on error. */
bench_options_t bench_parse_options(int argc, char *argv[]);

/* =====================================================================
 * Qubit positions
 * ===================================================================== */

/**
 * @brief Maximum number of qubit positions.
 *
 * 3Q gate on 16 qubits: 16*15*14 = 3360.  4096 covers up to ~16 qubits.
 */
#define BENCH_MAX_POSITIONS 4096

/**
 * @brief Generate all qubit positions for a gate on an n-qubit system.
 *
 * For 1Q gates: n positions (one per target).
 * For 2Q gates: n*(n-1) ordered pairs (control, target).
 * For 3Q gates: n*(n-1)*(n-2) ordered triples (c1, c2, target).
 *
 * Positions are stored as flat qubit_t arrays: pos[i*gate_qubits + j].
 *
 * @param n_qubits     System size.
 * @param gate_qubits  1, 2, or 3.
 * @param out          Output array of at least BENCH_MAX_POSITIONS * gate_qubits elements.
 * @return             Number of positions generated.
 */
int bench_gen_positions(qubit_t n_qubits, int gate_qubits, qubit_t *out);

/* =====================================================================
 * Per-position result (for per-qubit mode)
 * ===================================================================== */

typedef struct {
    bench_result_t result;
    qubit_t        pos[3];  /**< Qubit position (target / control,target / c1,c2,target) */
} bench_pos_result_t;

/* =====================================================================
 * Timing loop (bench_run.c)
 * ===================================================================== */

/**
 * @brief Run aggregate-mode benchmark for one backend at one qubit count.
 *
 * Times full sweeps across all positions, divides by positions to get
 * per-gate-application time.  Produces one bench_result_t.
 */
void bench_run_aggregate(const bench_backend_t *be, void *ctx,
                         qubit_t qubits,
                         const qubit_t *positions, int n_pos, int gate_qubits,
                         const bench_options_t *opts,
                         bench_result_t *out);

/**
 * @brief Run per-qubit-mode benchmark for one backend at one qubit count.
 *
 * Times each position independently.  Produces one bench_pos_result_t
 * per position.
 */
void bench_run_per_qubit(const bench_backend_t *be, void *ctx,
                         qubit_t qubits,
                         const qubit_t *positions, int n_pos, int gate_qubits,
                         const bench_options_t *opts,
                         bench_pos_result_t *out);

/* =====================================================================
 * Output formatters
 * ===================================================================== */

void bench_print_console_agg(qubit_t qubits, const bench_gate_info_t *gate,
                             const bench_options_t *opts,
                             const bench_result_t results[BACKEND_COUNT]);

void bench_print_console_perq(qubit_t qubits, const bench_gate_info_t *gate,
                              const bench_options_t *opts,
                              const bench_pos_result_t *per_backend[BACKEND_COUNT],
                              int n_pos, int gate_qubits);

void bench_print_csv_header(const bench_gate_info_t *gate, int per_qubit);

void bench_print_csv_agg(qubit_t qubits,
                         const bench_result_t results[BACKEND_COUNT]);

void bench_print_csv_perq(qubit_t qubits, bench_backend_id_t backend,
                          const bench_pos_result_t *results, int n_pos,
                          int gate_qubits);

void bench_print_pgfplots_agg(const bench_gate_info_t *gate,
                              const bench_options_t *opts,
                              const qubit_t *qubit_list, int n_qubits,
                              const bench_result_t (*results)[BACKEND_COUNT]);

void bench_print_pgfplots_perq(const bench_gate_info_t *gate,
                               const bench_options_t *opts,
                               qubit_t qubits,
                               const bench_pos_result_t *per_backend[BACKEND_COUNT],
                               int n_pos, int gate_qubits);

/* =====================================================================
 * Timing utilities
 * ===================================================================== */

/** @brief Get current time in nanoseconds (CLOCK_MONOTONIC). */
static inline uint64_t bench_time_ns(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (uint64_t)ts.tv_sec * 1000000000ULL + (uint64_t)ts.tv_nsec;
}

/** @brief CPU/compiler barrier to serialise timing measurements. */
static inline void bench_timing_barrier(void) {
#if defined(__aarch64__)
    __asm__ volatile("isb" ::: "memory");
#elif defined(__x86_64__)
    __asm__ volatile("lfence" ::: "memory");
#else
    __asm__ volatile("" ::: "memory");
#endif
}

/* =====================================================================
 * Random fill and Hermitian initialisation (bench_util.c)
 * ===================================================================== */

/**
 * @brief Fill n bytes with pseudo-random data seeded from OS entropy.
 *
 * Uses xoshiro256** expansion from a splitmix64-seeded state.
 */
void bench_fill_random(void *buf, size_t n);

/**
 * @brief Fill a flat dim x dim matrix as random Hermitian.
 *
 * Random-fills the full matrix, then enforces: upper triangle = conj(lower),
 * diagonal elements real.  Used by BLAS and competitor backends.
 */
void bench_init_hermitian_dense(cplx_t *rho, dim_t dim);

/**
 * @brief Fill packed lower-triangle storage as random Hermitian.
 *
 * Random-fills the packed array, then zeros the imaginary part of each
 * diagonal element.
 */
void bench_init_hermitian_packed(cplx_t *data, qubit_t qubits);

/**
 * @brief Fill tiled lower-triangle storage as random Hermitian.
 *
 * Random-fills all tile storage, then enforces Hermiticity within each
 * diagonal tile (upper = conj(lower), diagonal elements real).
 */
void bench_init_hermitian_tiled(cplx_t *data, qubit_t qubits);

/* =====================================================================
 * Huge-page memory utilities (bench_util.c)
 * ===================================================================== */

#define BENCH_HUGEPAGE_THRESHOLD ((size_t)1 << 30)

void *bench_alloc_huge(size_t n);
void  bench_free_huge(void *p, size_t n);

#ifdef __cplusplus
}
#endif

#endif /* BENCH_H */
