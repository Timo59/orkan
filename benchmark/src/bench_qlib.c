/**
 * @file bench_qlib.c
 * @brief qlib_tiled and qlib_packed backend wrappers for the gate benchmark.
 */

#include "bench.h"

/* =====================================================================
 * Context: state + resolved gate function pointer
 * ===================================================================== */

typedef struct {
    state_t *state;
    union {
        void (*fn_1q)(state_t*, qubit_t);
        void (*fn_1q_par)(state_t*, qubit_t, double);
        void (*fn_2q)(state_t*, qubit_t, qubit_t);
        void (*fn_3q)(state_t*, qubit_t, qubit_t, qubit_t);
    };
    double par;
    int    gate_qubits;
    int    has_par;
} qlib_ctx_t;

/* =====================================================================
 * Gate resolution — called once at init
 * ===================================================================== */

static void qlib_resolve_gate(qlib_ctx_t *ctx, bench_gate_id_t gate, double par) {
    ctx->par = par;
    ctx->has_par = 0;
    ctx->gate_qubits = bench_gate_table[gate].n_qubits;

    switch (gate) {
    case BG_X:   ctx->fn_1q = x;  break;
    case BG_Y:   ctx->fn_1q = y;  break;
    case BG_Z:   ctx->fn_1q = z;  break;
    case BG_H:   ctx->fn_1q = h;  break;
    case BG_S:   ctx->fn_1q = s;  break;
    case BG_T:   ctx->fn_1q = t;  break;
    case BG_RX:  ctx->fn_1q_par = rx;  ctx->has_par = 1; break;
    case BG_RY:  ctx->fn_1q_par = ry;  ctx->has_par = 1; break;
    case BG_RZ:  ctx->fn_1q_par = rz;  ctx->has_par = 1; break;
    case BG_CX:  ctx->fn_2q = cx;  break;
    case BG_CY:  ctx->fn_2q = cy;  break;
    case BG_CZ:  ctx->fn_2q = cz;  break;
    case BG_CCX: ctx->fn_3q = ccx; break;
    default:
        fprintf(stderr, "bench_qlib: unknown gate id %d\n", (int)gate);
        abort();
    }
}

/* =====================================================================
 * Apply — branches on arity and parametrised flag (both constant
 * for the lifetime of a benchmark run, so branch prediction is perfect)
 * ===================================================================== */

static void qlib_apply(void *ctx, const qubit_t *pos) {
    qlib_ctx_t *c = (qlib_ctx_t *)ctx;
    switch (c->gate_qubits) {
    case 1:
        if (c->has_par) c->fn_1q_par(c->state, pos[0], c->par);
        else            c->fn_1q(c->state, pos[0]);
        break;
    case 2: c->fn_2q(c->state, pos[0], pos[1]); break;
    case 3: c->fn_3q(c->state, pos[0], pos[1], pos[2]); break;
    default: __builtin_unreachable();
    }
}

/* =====================================================================
 * Cleanup (shared)
 * ===================================================================== */

static void qlib_cleanup(void *ctx) {
    qlib_ctx_t *c = (qlib_ctx_t *)ctx;
    state_free(c->state);
    free(c->state);
    free(c);
}

/* =====================================================================
 * qlib_tiled backend
 * ===================================================================== */

static void *qlib_tiled_init(qubit_t qubits, bench_gate_id_t gate, double par) {
    qlib_ctx_t *c = calloc(1, sizeof(qlib_ctx_t));
    if (!c) return NULL;
    c->state = malloc(sizeof(state_t));
    if (!c->state) { free(c); return NULL; }
    c->state->type = MIXED_TILED;
    state_init(c->state, qubits, NULL);
    if (!c->state->data) { free(c->state); free(c); return NULL; }
    bench_init_hermitian_tiled(c->state->data, qubits);
    qlib_resolve_gate(c, gate, par);
    return c;
}

static size_t qlib_tiled_touch_bytes(qubit_t qubits) {
    dim_t dim = (dim_t)1 << qubits;
    if (qubits < LOG_TILE_DIM) {
        /* State fits in a single diagonal tile — gate touches all
         * dim x dim elements within the tile, not the full
         * TILE_DIM x TILE_DIM padded storage. */
        return 2 * (size_t)dim * (size_t)dim * sizeof(cplx_t);
    }
    dim_t n_tiles = (dim + TILE_DIM - 1) / TILE_DIM;
    size_t n_tile_pairs = (size_t)n_tiles * (n_tiles + 1) / 2;
    return 2 * n_tile_pairs * TILE_SIZE * sizeof(cplx_t);
}

static int qlib_available(qubit_t qubits) {
    (void)qubits;
    return 1;
}

const bench_backend_t qlib_tiled_backend = {
    .name        = "qlib_tiled",
    .init        = qlib_tiled_init,
    .apply       = qlib_apply,
    .cleanup     = qlib_cleanup,
    .touch_bytes = qlib_tiled_touch_bytes,
    .available   = qlib_available
};

/* =====================================================================
 * qlib_packed backend
 * ===================================================================== */

static void *qlib_packed_init(qubit_t qubits, bench_gate_id_t gate, double par) {
    qlib_ctx_t *c = calloc(1, sizeof(qlib_ctx_t));
    if (!c) return NULL;
    c->state = malloc(sizeof(state_t));
    if (!c->state) { free(c); return NULL; }
    c->state->type = MIXED_PACKED;
    state_init(c->state, qubits, NULL);
    if (!c->state->data) { free(c->state); free(c); return NULL; }
    bench_init_hermitian_packed(c->state->data, qubits);
    qlib_resolve_gate(c, gate, par);
    return c;
}

static size_t qlib_packed_touch_bytes(qubit_t qubits) {
    dim_t dim = (dim_t)1 << qubits;
    return 2 * ((size_t)dim * (dim + 1) / 2) * sizeof(cplx_t);
}

const bench_backend_t qlib_packed_backend = {
    .name        = "qlib_packed",
    .init        = qlib_packed_init,
    .apply       = qlib_apply,
    .cleanup     = qlib_cleanup,
    .touch_bytes = qlib_packed_touch_bytes,
    .available   = qlib_available
};
