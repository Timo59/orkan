/**
 * @file bench_quest.c
 * @brief QuEST density-matrix backend for the gate benchmark.
 */

#ifdef HAVE_QUEST

#include "quest/include/quest.h"
#include "bench.h"

_Static_assert(sizeof(qcomp) == sizeof(cplx_t),
    "QuEST qcomp size must match cplx_t (double _Complex)");

/* =====================================================================
 * QuEST environment — singleton init/cleanup
 * ===================================================================== */

static int quest_env_initialised = 0;

static void quest_ensure_env(void) {
    if (!quest_env_initialised) {
        initQuESTEnv();
        quest_env_initialised = 1;
    }
}

/* =====================================================================
 * Context: Qureg + resolved apply function
 * ===================================================================== */

typedef void (*quest_apply_fn)(void *ctx, const qubit_t *pos);

typedef struct {
    Qureg          qureg;
    quest_apply_fn apply;
    double         par;
} quest_ctx_t;

/* =====================================================================
 * Per-gate apply functions — zero branching in hot path
 * ===================================================================== */

static void apply_x(void *c, const qubit_t *p) {
    applyPauliX(((quest_ctx_t*)c)->qureg, (int)p[0]);
}
static void apply_y(void *c, const qubit_t *p) {
    applyPauliY(((quest_ctx_t*)c)->qureg, (int)p[0]);
}
static void apply_z(void *c, const qubit_t *p) {
    applyPauliZ(((quest_ctx_t*)c)->qureg, (int)p[0]);
}
static void apply_h(void *c, const qubit_t *p) {
    applyHadamard(((quest_ctx_t*)c)->qureg, (int)p[0]);
}
static void apply_s(void *c, const qubit_t *p) {
    applyS(((quest_ctx_t*)c)->qureg, (int)p[0]);
}
static void apply_t(void *c, const qubit_t *p) {
    applyT(((quest_ctx_t*)c)->qureg, (int)p[0]);
}
static void apply_rx(void *c, const qubit_t *p) {
    quest_ctx_t *q = (quest_ctx_t*)c;
    applyRotateX(q->qureg, (int)p[0], q->par);
}
static void apply_ry(void *c, const qubit_t *p) {
    quest_ctx_t *q = (quest_ctx_t*)c;
    applyRotateY(q->qureg, (int)p[0], q->par);
}
static void apply_rz(void *c, const qubit_t *p) {
    quest_ctx_t *q = (quest_ctx_t*)c;
    applyRotateZ(q->qureg, (int)p[0], q->par);
}
static void apply_cx(void *c, const qubit_t *p) {
    applyControlledPauliX(((quest_ctx_t*)c)->qureg, (int)p[0], (int)p[1]);
}
static void apply_cy(void *c, const qubit_t *p) {
    applyControlledPauliY(((quest_ctx_t*)c)->qureg, (int)p[0], (int)p[1]);
}
static void apply_cz(void *c, const qubit_t *p) {
    applyControlledPauliZ(((quest_ctx_t*)c)->qureg, (int)p[0], (int)p[1]);
}
static void apply_ccx(void *c, const qubit_t *p) {
    quest_ctx_t *q = (quest_ctx_t*)c;
    int ctrls[2] = {(int)p[0], (int)p[1]};
    applyMultiControlledPauliX(q->qureg, ctrls, 2, (int)p[2]);
}

/* =====================================================================
 * Gate resolution
 * ===================================================================== */

static quest_apply_fn quest_resolve_gate(bench_gate_id_t gate) {
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
    default:     return NULL;
    }
}

/* =====================================================================
 * Vtable implementation
 * ===================================================================== */

static void *quest_init(qubit_t qubits, bench_gate_id_t gate, double par) {
    quest_ensure_env();

    quest_apply_fn fn = quest_resolve_gate(gate);
    if (!fn) return NULL;

    quest_ctx_t *ctx = malloc(sizeof(quest_ctx_t));
    if (!ctx) return NULL;

    ctx->qureg = createDensityQureg((int)qubits);
    ctx->apply = fn;
    ctx->par   = par;

    /* Random Hermitian init */
    bench_init_hermitian_dense((cplx_t*)ctx->qureg.cpuAmps,
                               (idx_t)1 << qubits);

    return ctx;
}

static void quest_apply(void *c, const qubit_t *pos) {
    ((quest_ctx_t*)c)->apply(c, pos);
}

static void quest_cleanup(void *c) {
    quest_ctx_t *ctx = (quest_ctx_t*)c;
    destroyQureg(ctx->qureg);
    free(ctx);
}

static size_t quest_touch_bytes(qubit_t qubits) {
    size_t dim = (size_t)1 << qubits;
    return 2 * dim * dim * sizeof(cplx_t);
}

static int quest_available(qubit_t qubits) {
    (void)qubits;
    return 1;
}

const bench_backend_t quest_backend = {
    .name        = "QuEST",
    .init        = quest_init,
    .apply       = quest_apply,
    .cleanup     = quest_cleanup,
    .touch_bytes = quest_touch_bytes,
    .available   = quest_available
};

#endif /* HAVE_QUEST */
