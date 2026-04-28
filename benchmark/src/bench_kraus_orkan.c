/**
 * @file bench_kraus_orkan.c
 * @brief orkan_tiled and orkan_packed backend wrappers for the Kraus benchmark.
 *
 * Pre-builds the superoperator at init (not timed), then calls channel_1q()
 * in the hot path.  Mirrors the fairness model of QuEST and Aer which also
 * pre-build superoperators.
 */

#include "bench_kraus.h"

/* Forward-declare channel types and functions to avoid including channel.h,
 * which has a conflicting insertBit0 definition with gate.h (pulled in by
 * bench.h).  Only the public API needed here: kraus_t, superop_t,
 * kraus_to_superop(), and channel_1q(). */
typedef struct {
    qubit_t  n_qubits;
    uint64_t n_terms;
    cplx_t  *data;
} kraus_t;

typedef struct {
    qubit_t  n_qubits;
    cplx_t  *data;
} superop_t;

superop_t kraus_to_superop(const kraus_t *kraus);
void channel_1q(state_t *state, const superop_t *sop, qubit_t target);

/* =====================================================================
 * Context: state + pre-built superoperator
 * ===================================================================== */

typedef struct {
    state_t   *state;
    superop_t  sop;
} kraus_orkan_ctx_t;

/* =====================================================================
 * Apply — single call, zero branching
 * ===================================================================== */

static void kraus_orkan_apply(void *c, const qubit_t *pos) {
    kraus_orkan_ctx_t *ctx = (kraus_orkan_ctx_t *)c;
    channel_1q(ctx->state, &ctx->sop, pos[0]);
}

/* =====================================================================
 * Cleanup (shared)
 * ===================================================================== */

static void kraus_orkan_cleanup(void *c) {
    kraus_orkan_ctx_t *ctx = (kraus_orkan_ctx_t *)c;
    state_free(ctx->state);
    free(ctx->state);
    free(ctx->sop.data);
    free(ctx);
}

/* =====================================================================
 * Build superoperator from benchmark Kraus ops
 *
 * The benchmark channel builder stores operators in row-major order.
 * kraus_to_superop() expects the same layout in kraus_t.data.
 * ===================================================================== */

static superop_t build_superop(int n_ops, int tgt_qubits, const cplx_t *ops) {
    idx_t d = (idx_t)1 << tgt_qubits;
    kraus_t kraus = {
        .n_qubits = (qubit_t)tgt_qubits,
        .n_terms  = (uint64_t)n_ops,
        .data     = (cplx_t *)ops   /* kraus_to_superop reads but does not modify */
    };
    return kraus_to_superop(&kraus);
}

/* =====================================================================
 * orkan_tiled backend
 * ===================================================================== */

static void *kraus_orkan_tiled_init(qubit_t sys_qubits, int n_ops,
                                    int tgt_qubits, const cplx_t *ops) {
    if (tgt_qubits != 1) return NULL;  /* only 1-qubit channels for now */

    kraus_orkan_ctx_t *ctx = calloc(1, sizeof(kraus_orkan_ctx_t));
    if (!ctx) return NULL;

    ctx->state = malloc(sizeof(state_t));
    if (!ctx->state) { free(ctx); return NULL; }
    ctx->state->type = MIXED_TILED;
    state_init(ctx->state, sys_qubits, NULL);
    if (!ctx->state->data) { free(ctx->state); free(ctx); return NULL; }
    bench_init_hermitian_tiled(ctx->state->data, sys_qubits);

    ctx->sop = build_superop(n_ops, tgt_qubits, ops);
    return ctx;
}

static size_t kraus_orkan_tiled_touch_bytes(qubit_t qubits) {
    idx_t dim = (idx_t)1 << qubits;
    if (qubits < LOG_TILE_DIM) {
        return 2 * (size_t)dim * (size_t)dim * sizeof(cplx_t);
    }
    idx_t n_tiles = (dim + TILE_DIM - 1) / TILE_DIM;
    size_t n_tile_pairs = (size_t)n_tiles * (n_tiles + 1) / 2;
    return 2 * n_tile_pairs * TILE_SIZE * sizeof(cplx_t);
}

static int kraus_orkan_available(qubit_t qubits) {
    (void)qubits;
    return 1;
}

const bench_kraus_backend_t kraus_orkan_tiled_backend = {
    .name        = "orkan_tiled",
    .init        = kraus_orkan_tiled_init,
    .apply       = kraus_orkan_apply,
    .cleanup     = kraus_orkan_cleanup,
    .touch_bytes = kraus_orkan_tiled_touch_bytes,
    .available   = kraus_orkan_available
};

/* =====================================================================
 * orkan_packed backend
 * ===================================================================== */

static void *kraus_orkan_packed_init(qubit_t sys_qubits, int n_ops,
                                     int tgt_qubits, const cplx_t *ops) {
    if (tgt_qubits != 1) return NULL;

    kraus_orkan_ctx_t *ctx = calloc(1, sizeof(kraus_orkan_ctx_t));
    if (!ctx) return NULL;

    ctx->state = malloc(sizeof(state_t));
    if (!ctx->state) { free(ctx); return NULL; }
    ctx->state->type = MIXED_PACKED;
    state_init(ctx->state, sys_qubits, NULL);
    if (!ctx->state->data) { free(ctx->state); free(ctx); return NULL; }
    bench_init_hermitian_packed(ctx->state->data, sys_qubits);

    ctx->sop = build_superop(n_ops, tgt_qubits, ops);
    return ctx;
}

static size_t kraus_orkan_packed_touch_bytes(qubit_t qubits) {
    idx_t dim = (idx_t)1 << qubits;
    return 2 * ((size_t)dim * (dim + 1) / 2) * sizeof(cplx_t);
}

const bench_kraus_backend_t kraus_orkan_packed_backend = {
    .name        = "orkan_packed",
    .init        = kraus_orkan_packed_init,
    .apply       = kraus_orkan_apply,
    .cleanup     = kraus_orkan_cleanup,
    .touch_bytes = kraus_orkan_packed_touch_bytes,
    .available   = kraus_orkan_available
};
