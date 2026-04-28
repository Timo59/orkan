/**
 * @file bench_orkan.c
 * @brief orkan_tiled and orkan_packed backend wrappers for the gate benchmark.
 */

#include "bench.h"

/* =====================================================================
 * Context: state + resolved apply function
 * ===================================================================== */

typedef void (*orkan_apply_fn)(void *ctx, const qubit_t *pos);

typedef struct {
    state_t       *state;
    orkan_apply_fn  apply;
    double         par;
} orkan_ctx_t;

/* =====================================================================
 * Per-gate apply functions — zero branching in hot path
 * ===================================================================== */

static void apply_x(void *c, const qubit_t *p)   { x(((orkan_ctx_t*)c)->state, p[0]); }
static void apply_y(void *c, const qubit_t *p)   { y(((orkan_ctx_t*)c)->state, p[0]); }
static void apply_z(void *c, const qubit_t *p)   { z(((orkan_ctx_t*)c)->state, p[0]); }
static void apply_h(void *c, const qubit_t *p)   { h(((orkan_ctx_t*)c)->state, p[0]); }
static void apply_s(void *c, const qubit_t *p)   { s(((orkan_ctx_t*)c)->state, p[0]); }
static void apply_t(void *c, const qubit_t *p)   { t(((orkan_ctx_t*)c)->state, p[0]); }
static void apply_rx(void *c, const qubit_t *p)  { orkan_ctx_t *q = (orkan_ctx_t*)c; rx(q->state, p[0], q->par); }
static void apply_ry(void *c, const qubit_t *p)  { orkan_ctx_t *q = (orkan_ctx_t*)c; ry(q->state, p[0], q->par); }
static void apply_rz(void *c, const qubit_t *p)  { orkan_ctx_t *q = (orkan_ctx_t*)c; rz(q->state, p[0], q->par); }
static void apply_cx(void *c, const qubit_t *p)  { cx(((orkan_ctx_t*)c)->state, p[0], p[1]); }
static void apply_cy(void *c, const qubit_t *p)  { cy(((orkan_ctx_t*)c)->state, p[0], p[1]); }
static void apply_cz(void *c, const qubit_t *p)  { cz(((orkan_ctx_t*)c)->state, p[0], p[1]); }
static void apply_ccx(void *c, const qubit_t *p) { ccx(((orkan_ctx_t*)c)->state, p[0], p[1], p[2]); }

/* =====================================================================
 * Gate resolution — called once at init
 * ===================================================================== */

static orkan_apply_fn orkan_resolve_gate(bench_gate_id_t gate) {
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
        fprintf(stderr, "bench_orkan: unknown gate id %d\n", (int)gate);
        abort();
    }
}

/* =====================================================================
 * Apply — single indirect call, zero branching
 * ===================================================================== */

static void orkan_apply(void *ctx, const qubit_t *pos) {
    orkan_ctx_t *c = (orkan_ctx_t *)ctx;
    c->apply(ctx, pos);
}

/* =====================================================================
 * Cleanup (shared)
 * ===================================================================== */

static void orkan_cleanup(void *ctx) {
    orkan_ctx_t *c = (orkan_ctx_t *)ctx;
    state_free(c->state);
    free(c->state);
    free(c);
}

/* =====================================================================
 * orkan_tiled backend
 * ===================================================================== */

static void *orkan_tiled_init(qubit_t qubits, bench_gate_id_t gate, double par) {
    orkan_ctx_t *c = calloc(1, sizeof(orkan_ctx_t));
    if (!c) return NULL;
    c->state = malloc(sizeof(state_t));
    if (!c->state) { free(c); return NULL; }
    c->state->type = MIXED_TILED;
    state_init(c->state, qubits, NULL);
    if (!c->state->data) { free(c->state); free(c); return NULL; }
    bench_init_hermitian_tiled(c->state->data, qubits);
    c->par = par;
    c->apply = orkan_resolve_gate(gate);
    return c;
}

static size_t orkan_tiled_touch_bytes(qubit_t qubits) {
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

static int orkan_available(qubit_t qubits) {
    (void)qubits;
    return 1;
}

const bench_backend_t orkan_tiled_backend = {
    .name        = "orkan_tiled",
    .init        = orkan_tiled_init,
    .apply       = orkan_apply,
    .cleanup     = orkan_cleanup,
    .touch_bytes = orkan_tiled_touch_bytes,
    .available   = orkan_available
};

/* =====================================================================
 * orkan_packed backend
 * ===================================================================== */

static void *orkan_packed_init(qubit_t qubits, bench_gate_id_t gate, double par) {
    orkan_ctx_t *c = calloc(1, sizeof(orkan_ctx_t));
    if (!c) return NULL;
    c->state = malloc(sizeof(state_t));
    if (!c->state) { free(c); return NULL; }
    c->state->type = MIXED_PACKED;
    state_init(c->state, qubits, NULL);
    if (!c->state->data) { free(c->state); free(c); return NULL; }
    bench_init_hermitian_packed(c->state->data, qubits);
    c->par = par;
    c->apply = orkan_resolve_gate(gate);
    return c;
}

static size_t orkan_packed_touch_bytes(qubit_t qubits) {
    idx_t dim = (idx_t)1 << qubits;
    return 2 * ((size_t)dim * (dim + 1) / 2) * sizeof(cplx_t);
}

const bench_backend_t orkan_packed_backend = {
    .name        = "orkan_packed",
    .init        = orkan_packed_init,
    .apply       = orkan_apply,
    .cleanup     = orkan_cleanup,
    .touch_bytes = orkan_packed_touch_bytes,
    .available   = orkan_available
};
