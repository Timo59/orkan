/**
 * @file bench_qlib.c
 * @brief qlib_tiled and qlib_packed backend wrappers for the gate benchmark.
 */

#include "bench.h"

/* =====================================================================
 * Context: state + resolved apply function
 * ===================================================================== */

typedef void (*qlib_apply_fn)(void *ctx, const qubit_t *pos);

typedef struct {
    state_t       *state;
    qlib_apply_fn  apply;
    double         par;
} qlib_ctx_t;

/* =====================================================================
 * Per-gate apply functions — zero branching in hot path
 * ===================================================================== */

static void apply_x(void *c, const qubit_t *p)   { x(((qlib_ctx_t*)c)->state, p[0]); }
static void apply_y(void *c, const qubit_t *p)   { y(((qlib_ctx_t*)c)->state, p[0]); }
static void apply_z(void *c, const qubit_t *p)   { z(((qlib_ctx_t*)c)->state, p[0]); }
static void apply_h(void *c, const qubit_t *p)   { h(((qlib_ctx_t*)c)->state, p[0]); }
static void apply_s(void *c, const qubit_t *p)   { s(((qlib_ctx_t*)c)->state, p[0]); }
static void apply_t(void *c, const qubit_t *p)   { t(((qlib_ctx_t*)c)->state, p[0]); }
static void apply_rx(void *c, const qubit_t *p)  { qlib_ctx_t *q = (qlib_ctx_t*)c; rx(q->state, p[0], q->par); }
static void apply_ry(void *c, const qubit_t *p)  { qlib_ctx_t *q = (qlib_ctx_t*)c; ry(q->state, p[0], q->par); }
static void apply_rz(void *c, const qubit_t *p)  { qlib_ctx_t *q = (qlib_ctx_t*)c; rz(q->state, p[0], q->par); }
static void apply_cx(void *c, const qubit_t *p)  { cx(((qlib_ctx_t*)c)->state, p[0], p[1]); }
static void apply_cy(void *c, const qubit_t *p)  { cy(((qlib_ctx_t*)c)->state, p[0], p[1]); }
static void apply_cz(void *c, const qubit_t *p)  { cz(((qlib_ctx_t*)c)->state, p[0], p[1]); }
static void apply_ccx(void *c, const qubit_t *p) { ccx(((qlib_ctx_t*)c)->state, p[0], p[1], p[2]); }

/* =====================================================================
 * Gate resolution — called once at init
 * ===================================================================== */

static qlib_apply_fn qlib_resolve_gate(bench_gate_id_t gate) {
    switch (gate) {
    case BG_X:   return apply_x;
    case BG_Y:   return apply_y;
    case BG_Z:   return apply_z;
    case BG_H:   return apply_h;
    case BG_S:   return apply_s;
    case BG_T:   return apply_t;
    case BG_RX:  return apply_rx;
    case BG_RY:  return apply_ry;
    case BG_RZ:  return apply_rz;
    case BG_CX:  return apply_cx;
    case BG_CY:  return apply_cy;
    case BG_CZ:  return apply_cz;
    case BG_CCX: return apply_ccx;
    default:
        fprintf(stderr, "bench_qlib: unknown gate id %d\n", (int)gate);
        abort();
    }
}

/* =====================================================================
 * Apply — single indirect call, zero branching
 * ===================================================================== */

static void qlib_apply(void *ctx, const qubit_t *pos) {
    qlib_ctx_t *c = (qlib_ctx_t *)ctx;
    c->apply(ctx, pos);
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
    c->par = par;
    c->apply = qlib_resolve_gate(gate);
    return c;
}

static size_t qlib_tiled_touch_bytes(qubit_t qubits) {
    idx_t dim = (idx_t)1 << qubits;
    if (qubits < LOG_TILE_DIM) {
        /* State fits in a single diagonal tile — gate touches all
         * dim x dim elements within the tile, not the full
         * TILE_DIM x TILE_DIM padded storage. */
        return 2 * (size_t)dim * (size_t)dim * sizeof(cplx_t);
    }
    idx_t n_tiles = (dim + TILE_DIM - 1) / TILE_DIM;
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
    c->par = par;
    c->apply = qlib_resolve_gate(gate);
    return c;
}

static size_t qlib_packed_touch_bytes(qubit_t qubits) {
    idx_t dim = (idx_t)1 << qubits;
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
