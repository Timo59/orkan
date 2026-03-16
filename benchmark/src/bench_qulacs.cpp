/**
 * @file bench_qulacs.cpp
 * @brief Qulacs density-matrix backend for the gate benchmark.
 */

#ifdef HAVE_QULACS

#include <csim/update_ops_dm.hpp>
#include <csim/type.hpp>
#include <cstdlib>
#include <cstdio>

#include "bench.h"

static_assert(sizeof(CTYPE) == sizeof(cplx_t),
    "Qulacs CTYPE size must match cplx_t (double _Complex)");

/* =====================================================================
 * Context: flat DM array + resolved apply function
 * ===================================================================== */

typedef void (*qulacs_apply_fn)(void *ctx, const qubit_t *pos);

struct qulacs_ctx_t {
    CTYPE           *rho;
    ITYPE            dim;
    qulacs_apply_fn  apply;
    double           par;
    CTYPE            cy_mat[4];   /* Y gate matrix for CY via controlled gate */
    CTYPE            cx_mat[4];   /* X gate matrix for CCX via multi-controlled gate */
};

/* =====================================================================
 * Per-gate apply functions — zero branching in hot path
 * ===================================================================== */

static void apply_x(void *c, const qubit_t *p) {
    auto *q = static_cast<qulacs_ctx_t*>(c);
    dm_X_gate(static_cast<UINT>(p[0]), q->rho, q->dim);
}
static void apply_y(void *c, const qubit_t *p) {
    auto *q = static_cast<qulacs_ctx_t*>(c);
    dm_Y_gate(static_cast<UINT>(p[0]), q->rho, q->dim);
}
static void apply_z(void *c, const qubit_t *p) {
    auto *q = static_cast<qulacs_ctx_t*>(c);
    dm_Z_gate(static_cast<UINT>(p[0]), q->rho, q->dim);
}
static void apply_h(void *c, const qubit_t *p) {
    auto *q = static_cast<qulacs_ctx_t*>(c);
    dm_H_gate(static_cast<UINT>(p[0]), q->rho, q->dim);
}
static void apply_s(void *c, const qubit_t *p) {
    auto *q = static_cast<qulacs_ctx_t*>(c);
    dm_S_gate(static_cast<UINT>(p[0]), q->rho, q->dim);
}
static void apply_t(void *c, const qubit_t *p) {
    auto *q = static_cast<qulacs_ctx_t*>(c);
    dm_T_gate(static_cast<UINT>(p[0]), q->rho, q->dim);
}
static void apply_rx(void *c, const qubit_t *p) {
    auto *q = static_cast<qulacs_ctx_t*>(c);
    dm_RX_gate(static_cast<UINT>(p[0]), q->par, q->rho, q->dim);
}
static void apply_ry(void *c, const qubit_t *p) {
    auto *q = static_cast<qulacs_ctx_t*>(c);
    dm_RY_gate(static_cast<UINT>(p[0]), q->par, q->rho, q->dim);
}
static void apply_rz(void *c, const qubit_t *p) {
    auto *q = static_cast<qulacs_ctx_t*>(c);
    dm_RZ_gate(static_cast<UINT>(p[0]), q->par, q->rho, q->dim);
}
static void apply_cx(void *c, const qubit_t *p) {
    auto *q = static_cast<qulacs_ctx_t*>(c);
    dm_CNOT_gate(static_cast<UINT>(p[0]), static_cast<UINT>(p[1]),
                  q->rho, q->dim);
}
static void apply_cy(void *c, const qubit_t *p) {
    /* No dedicated dm_CY_gate — use controlled single-qubit dense matrix */
    auto *q = static_cast<qulacs_ctx_t*>(c);
    UINT ctrl = static_cast<UINT>(p[0]);
    UINT ctrl_val = 1;
    dm_multi_qubit_control_single_qubit_dense_matrix_gate(
        &ctrl, &ctrl_val, 1, static_cast<UINT>(p[1]),
        q->cy_mat, q->rho, q->dim);
}
static void apply_cz(void *c, const qubit_t *p) {
    auto *q = static_cast<qulacs_ctx_t*>(c);
    dm_CZ_gate(static_cast<UINT>(p[0]), static_cast<UINT>(p[1]),
                q->rho, q->dim);
}
static void apply_ccx(void *c, const qubit_t *p) {
    /* No dedicated dm_CCX_gate — use multi-controlled with X matrix */
    auto *q = static_cast<qulacs_ctx_t*>(c);
    UINT ctrls[2] = {static_cast<UINT>(p[0]), static_cast<UINT>(p[1])};
    UINT vals[2]  = {1, 1};
    dm_multi_qubit_control_single_qubit_dense_matrix_gate(
        ctrls, vals, 2, static_cast<UINT>(p[2]),
        q->cx_mat, q->rho, q->dim);
}

/* =====================================================================
 * Vtable implementation
 * ===================================================================== */

static void *qulacs_init(qubit_t qubits, bench_gate_id_t gate, double par) {
    auto *ctx = static_cast<qulacs_ctx_t*>(std::calloc(1, sizeof(qulacs_ctx_t)));
    if (!ctx) return nullptr;

    ctx->dim = static_cast<ITYPE>(1) << qubits;
    size_t n_elems = static_cast<size_t>(ctx->dim) * static_cast<size_t>(ctx->dim);
    ctx->rho = static_cast<CTYPE*>(std::malloc(n_elems * sizeof(CTYPE)));
    if (!ctx->rho) { std::free(ctx); return nullptr; }
    ctx->par = par;

    /* Random Hermitian init */
    bench_init_hermitian_dense(
        reinterpret_cast<cplx_t*>(ctx->rho), static_cast<dim_t>(ctx->dim));

    /* Pre-build gate matrices for CY and CCX */
    /* Y = [[0, -i], [i, 0]] in row-major (Qulacs convention): mat[r*2+c] */
    ctx->cy_mat[0] = CTYPE(0.0, 0.0);
    ctx->cy_mat[1] = CTYPE(0.0, -1.0);
    ctx->cy_mat[2] = CTYPE(0.0, 1.0);
    ctx->cy_mat[3] = CTYPE(0.0, 0.0);
    /* X = [[0, 1], [1, 0]] in row-major (Qulacs convention): mat[r*2+c] */
    ctx->cx_mat[0] = CTYPE(0.0, 0.0);
    ctx->cx_mat[1] = CTYPE(1.0, 0.0);
    ctx->cx_mat[2] = CTYPE(1.0, 0.0);
    ctx->cx_mat[3] = CTYPE(0.0, 0.0);

    /* Resolve apply function */
    switch (gate) {
    case BG_X:   ctx->apply = apply_x;   break;
    case BG_Y:   ctx->apply = apply_y;   break;
    case BG_Z:   ctx->apply = apply_z;   break;
    case BG_H:   ctx->apply = apply_h;   break;
    case BG_S:   ctx->apply = apply_s;   break;
    case BG_T:   ctx->apply = apply_t;   break;
    case BG_RX:  ctx->apply = apply_rx;  break;
    case BG_RY:  ctx->apply = apply_ry;  break;
    case BG_RZ:  ctx->apply = apply_rz;  break;
    case BG_CX:  ctx->apply = apply_cx;  break;
    case BG_CY:  ctx->apply = apply_cy;  break;
    case BG_CZ:  ctx->apply = apply_cz;  break;
    case BG_CCX: ctx->apply = apply_ccx; break;
    default:
        std::free(ctx->rho);
        std::free(ctx);
        return nullptr;
    }
    return ctx;
}

static void qulacs_apply(void *c, const qubit_t *pos) {
    static_cast<qulacs_ctx_t*>(c)->apply(c, pos);
}

static void qulacs_cleanup(void *c) {
    auto *ctx = static_cast<qulacs_ctx_t*>(c);
    std::free(ctx->rho);
    std::free(ctx);
}

static size_t qulacs_touch_bytes(qubit_t qubits) {
    size_t dim = static_cast<size_t>(1) << qubits;
    return 2 * dim * dim * sizeof(CTYPE);
}

static int qulacs_available(qubit_t qubits) {
    (void)qubits;
    return 1;
}

extern "C" {
const bench_backend_t qulacs_backend = {
    .name        = "Qulacs",
    .init        = qulacs_init,
    .apply       = qulacs_apply,
    .cleanup     = qulacs_cleanup,
    .touch_bytes = qulacs_touch_bytes,
    .available   = qulacs_available
};
}

#endif /* HAVE_QULACS */
