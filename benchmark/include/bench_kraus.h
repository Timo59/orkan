/**
 * @file bench_kraus.h
 * @brief Channel table and Kraus backend vtable for the Kraus benchmark.
 *
 * Reuses types and formatters from bench.h. Adds only:
 *   - Channel definitions (Kraus operator construction)
 *   - Kraus backend vtable (different init signature from gate backend)
 *   - Kraus-specific run functions
 */

#ifndef BENCH_KRAUS_H
#define BENCH_KRAUS_H

#include "bench.h"

#ifdef __cplusplus
extern "C" {
#endif

/* =====================================================================
 * Channel definitions
 * ===================================================================== */

typedef enum {
    BC_DEPOLARIZE,
    BC_COUNT
} bench_channel_id_t;

/**
 * @brief Build the Kraus operators for a channel.
 *
 * Writes n_ops matrices of size d x d (d = 2^tgt_qubits) into out,
 * stored contiguously in row-major order: out[k * d * d + i * d + j].
 */
typedef void (*bench_channel_build_fn)(double par, cplx_t *out);

typedef struct {
    bench_channel_id_t  id;
    const char         *name;
    int                 n_ops;
    int                 tgt_qubits;
    int                 parametrised;
    bench_channel_build_fn build;
} bench_channel_info_t;

extern const bench_channel_info_t bench_channel_table[BC_COUNT];

bench_channel_id_t bench_resolve_channel(const char *name);

/* =====================================================================
 * Kraus backend vtable
 * ===================================================================== */

typedef struct {
    const char *name;

    void *(*init)(qubit_t sys_qubits, int n_ops,
                  int tgt_qubits, const cplx_t *ops);

    void  (*apply)(void *ctx, const qubit_t *pos);
    void  (*cleanup)(void *ctx);
    size_t (*touch_bytes)(qubit_t sys_qubits);
    int   (*available)(qubit_t sys_qubits);
} bench_kraus_backend_t;

#define KRAUS_BACKEND_COUNT 6

extern const bench_kraus_backend_t kraus_qlib_tiled_backend;
extern const bench_kraus_backend_t kraus_qlib_packed_backend;
extern const bench_kraus_backend_t kraus_blas_backend;

#ifdef HAVE_QUEST
extern const bench_kraus_backend_t kraus_quest_backend;
#endif
#ifdef HAVE_QULACS
extern const bench_kraus_backend_t kraus_qulacs_backend;
#endif
#ifdef HAVE_AER
extern const bench_kraus_backend_t kraus_aer_backend;
#endif

/* =====================================================================
 * Position generation
 * ===================================================================== */

int kraus_gen_positions(qubit_t n_qubits, int tgt_qubits, qubit_t *out);

/* =====================================================================
 * Kraus-specific run functions
 * ===================================================================== */

void kraus_run_aggregate(const bench_kraus_backend_t *be, void *ctx,
                         qubit_t qubits,
                         const qubit_t *positions, int n_pos, int tgt_qubits,
                         const bench_options_t *opts,
                         bench_result_t *out);

void kraus_run_per_qubit(const bench_kraus_backend_t *be, void *ctx,
                         qubit_t qubits,
                         const qubit_t *positions, int n_pos, int tgt_qubits,
                         const bench_options_t *opts,
                         bench_pos_result_t *out);

#ifdef __cplusplus
}
#endif

#endif /* BENCH_KRAUS_H */
