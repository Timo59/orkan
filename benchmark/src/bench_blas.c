/**
 * @file bench_blas.c
 * @brief BLAS baseline: ρ' = UρU† via zgemm, pre-computed Kronecker matrices.
 *
 * Capped at 10 qubits.  Above 10, returns unavailable.
 */

#include "bench.h"
#include <math.h>
#include <complex.h>

#if defined(__APPLE__)
    #include <vecLib/cblas_new.h>
#elif defined(__linux__)
    #include <cblas.h>
#endif

#define BLAS_MAX_QUBITS 10
#define S2 0.7071067811865475

/* =====================================================================
 * Gate matrices (column-major)
 * ===================================================================== */

static const cplx_t MAT_X[4] = {0, 1, 1, 0};
static const cplx_t MAT_Y[4] = {0, I, -I, 0};
static const cplx_t MAT_Z[4] = {1, 0, 0, -1};
static const cplx_t MAT_H[4] = {S2, S2, S2, -S2};
static const cplx_t MAT_S[4] = {1, 0, 0, I};
static const cplx_t MAT_T[4] = {1, 0, 0, S2 + S2*I};
static const cplx_t MAT_CX[16] = {1,0,0,0, 0,0,0,1, 0,0,1,0, 0,1,0,0};
static const cplx_t MAT_CY[16] = {1,0,0,0, 0,0,0,I, 0,0,1,0, 0,-I,0,0};
static const cplx_t MAT_CZ[16] = {1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,-1};

static void build_rx(double t, cplx_t o[4]) {
    double c=cos(t/2), s=sin(t/2);
    o[0]=c; o[1]=-s*I; o[2]=-s*I; o[3]=c;
}
static void build_ry(double t, cplx_t o[4]) {
    double c=cos(t/2), s=sin(t/2);
    o[0]=c; o[1]=s; o[2]=-s; o[3]=c;
}
static void build_rz(double t, cplx_t o[4]) {
    o[0]=cexp(-I*t/2); o[1]=0; o[2]=0; o[3]=cexp(I*t/2);
}

/* =====================================================================
 * Kronecker product: C = A ⊗ B (column-major)
 * ===================================================================== */

static cplx_t *kron(idx_t ma, idx_t na, const cplx_t *A,
                     idx_t mb, idx_t nb, const cplx_t *B) {
    idx_t mc = ma*mb, nc = na*nb;
    cplx_t *C = malloc((size_t)(mc*nc) * sizeof(cplx_t));
    if (!C) return NULL;
    for (idx_t ja = 0; ja < na; ++ja)
        for (idx_t ia = 0; ia < ma; ++ia) {
            cplx_t a = A[ia + ja*ma];
            for (idx_t jb = 0; jb < nb; ++jb)
                for (idx_t ib = 0; ib < mb; ++ib)
                    C[(ib+ia*mb) + (jb+ja*nb)*mc] = a * B[ib+jb*mb];
        }
    return C;
}

/* =====================================================================
 * Full gate matrix builders
 * ===================================================================== */

static cplx_t *build_full_1q(qubit_t qubits, const cplx_t g[4], qubit_t tgt) {
    qubit_t pos = qubits - 1 - tgt;
    idx_t dl = (idx_t)1 << pos;
    cplx_t *idl = calloc((size_t)(dl*dl), sizeof(cplx_t));
    if (!idl) return NULL;
    for (idx_t i = 0; i < dl; ++i) idl[i+i*dl] = 1.0;
    cplx_t *tmp = kron(dl, dl, idl, 2, 2, g);
    free(idl);
    if (!tmp) return NULL;
    dl *= 2;
    idx_t dr = (idx_t)1 << tgt;
    cplx_t *idr = calloc((size_t)(dr*dr), sizeof(cplx_t));
    if (!idr) { free(tmp); return NULL; }
    for (idx_t i = 0; i < dr; ++i) idr[i+i*dr] = 1.0;
    cplx_t *U = kron(dl, dl, tmp, dr, dr, idr);
    free(tmp); free(idr);
    return U;
}

static cplx_t *build_full_2q(qubit_t qubits, const cplx_t g[16],
                               qubit_t q1, qubit_t q2) {
    idx_t dim = (idx_t)1 << qubits;
    cplx_t *U = calloc((size_t)(dim*dim), sizeof(cplx_t));
    if (!U) return NULL;
    for (idx_t i = 0; i < dim; ++i) {
        int b1 = (i>>q1)&1, b2 = (i>>q2)&1;
        int gc = b2*2+b1;
        for (int gr = 0; gr < 4; ++gr) {
            cplx_t v = g[gr+gc*4];
            if (cabs(v) < 1e-15) continue;
            idx_t j = i;
            j = (j & ~((idx_t)1<<q1)) | ((idx_t)(gr&1)<<q1);
            j = (j & ~((idx_t)1<<q2)) | ((idx_t)((gr>>1)&1)<<q2);
            U[j+i*dim] = v;
        }
    }
    return U;
}

static cplx_t *build_full_ccx(qubit_t qubits, qubit_t c1, qubit_t c2, qubit_t tgt) {
    idx_t dim = (idx_t)1 << qubits;
    cplx_t *U = calloc((size_t)(dim*dim), sizeof(cplx_t));
    if (!U) return NULL;
    for (idx_t i = 0; i < dim; ++i) {
        idx_t j = ((i>>c1)&1) && ((i>>c2)&1) ? i ^ ((idx_t)1<<tgt) : i;
        U[j+i*dim] = 1.0;
    }
    return U;
}

/* =====================================================================
 * ρ' = UρU† via zgemm
 * ===================================================================== */

static void apply_uruh(idx_t dim, const cplx_t *U, cplx_t *rho, cplx_t *tmp) {
    const cplx_t one = 1.0, zero = 0.0;
    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                dim, dim, dim, &one, U, dim, rho, dim, &zero, tmp, dim);
    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasConjTrans,
                dim, dim, dim, &one, tmp, dim, U, dim, &zero, rho, dim);
}

/* =====================================================================
 * Context + per-arity apply functions
 * ===================================================================== */

typedef void (*blas_apply_fn)(void *ctx, const qubit_t *pos);

typedef struct {
    cplx_t       *rho;
    cplx_t       *tmp;
    cplx_t      **U_lookup;
    int           lookup_size;
    idx_t         dim;
    qubit_t       n_qubits;
    blas_apply_fn apply;
} blas_ctx_t;

static void blas_apply_1q(void *c, const qubit_t *pos) {
    blas_ctx_t *b = (blas_ctx_t *)c;
    apply_uruh(b->dim, b->U_lookup[pos[0]], b->rho, b->tmp);
}
static void blas_apply_2q(void *c, const qubit_t *pos) {
    blas_ctx_t *b = (blas_ctx_t *)c;
    apply_uruh(b->dim, b->U_lookup[(int)pos[0] * b->n_qubits + (int)pos[1]],
               b->rho, b->tmp);
}
static void blas_apply_3q(void *c, const qubit_t *pos) {
    blas_ctx_t *b = (blas_ctx_t *)c;
    int n = b->n_qubits;
    apply_uruh(b->dim, b->U_lookup[(int)pos[0]*n*n + (int)pos[1]*n + (int)pos[2]],
               b->rho, b->tmp);
}

/* =====================================================================
 * Vtable
 * ===================================================================== */

static void *blas_init(qubit_t qubits, bench_gate_id_t gate, double par) {
    if (qubits > BLAS_MAX_QUBITS) return NULL;
    const bench_gate_info_t *gi = &bench_gate_table[gate];
    idx_t dim = (idx_t)1 << qubits;
    int n = (int)qubits;

    blas_ctx_t *ctx = calloc(1, sizeof(blas_ctx_t));
    if (!ctx) return NULL;
    ctx->dim = dim;
    ctx->n_qubits = qubits;

    ctx->rho = malloc((size_t)(dim*dim) * sizeof(cplx_t));
    ctx->tmp = malloc((size_t)(dim*dim) * sizeof(cplx_t));
    if (!ctx->rho || !ctx->tmp) goto fail;
    bench_init_hermitian_dense(ctx->rho, dim);

    /* Resolve gate matrix */
    cplx_t par_mat[4];
    const cplx_t *m1q = NULL, *m2q = NULL;
    switch (gate) {
    case BG_X:  m1q=MAT_X; break;  case BG_Y:  m1q=MAT_Y; break;
    case BG_Z:  m1q=MAT_Z; break;  case BG_H:  m1q=MAT_H; break;
    case BG_S:  m1q=MAT_S; break;  case BG_T:  m1q=MAT_T; break;
    case BG_RX: build_rx(par,par_mat); m1q=par_mat; break;
    case BG_RY: build_ry(par,par_mat); m1q=par_mat; break;
    case BG_RZ: build_rz(par,par_mat); m1q=par_mat; break;
    case BG_CX: m2q=MAT_CX; break;
    case BG_CY: m2q=MAT_CY; break;
    case BG_CZ: m2q=MAT_CZ; break;
    case BG_CCX: break;
    default: goto fail;
    }

    /* Build lookup table of pre-computed full gate matrices */
    if (gi->n_qubits == 1) {
        ctx->apply = blas_apply_1q;
        ctx->lookup_size = n;
        ctx->U_lookup = calloc((size_t)n, sizeof(cplx_t*));
        if (!ctx->U_lookup) goto fail;
        for (int t = 0; t < n; ++t) {
            ctx->U_lookup[t] = build_full_1q(qubits, m1q, (qubit_t)t);
            if (!ctx->U_lookup[t]) goto fail;
        }
    } else if (gi->n_qubits == 2) {
        ctx->apply = blas_apply_2q;
        ctx->lookup_size = n * n;
        ctx->U_lookup = calloc((size_t)(n*n), sizeof(cplx_t*));
        if (!ctx->U_lookup) goto fail;
        for (int c = 0; c < n; ++c)
            for (int t = 0; t < n; ++t) {
                if (c == t) continue;
                ctx->U_lookup[c*n+t] = build_full_2q(qubits, m2q, (qubit_t)c, (qubit_t)t);
                if (!ctx->U_lookup[c*n+t]) goto fail;
            }
    } else {
        ctx->apply = blas_apply_3q;
        ctx->lookup_size = n * n * n;
        ctx->U_lookup = calloc((size_t)(n*n*n), sizeof(cplx_t*));
        if (!ctx->U_lookup) goto fail;
        for (int a = 0; a < n; ++a)
            for (int b = 0; b < n; ++b)
                for (int t = 0; t < n; ++t) {
                    if (a==b || a==t || b==t) continue;
                    ctx->U_lookup[a*n*n+b*n+t] = build_full_ccx(qubits,
                        (qubit_t)a, (qubit_t)b, (qubit_t)t);
                    if (!ctx->U_lookup[a*n*n+b*n+t]) goto fail;
                }
    }
    return ctx;

fail:
    if (ctx->U_lookup) {
        for (int i = 0; i < ctx->lookup_size; ++i) free(ctx->U_lookup[i]);
        free(ctx->U_lookup);
    }
    free(ctx->tmp); free(ctx->rho); free(ctx);
    return NULL;
}

static void blas_apply(void *c, const qubit_t *pos) {
    ((blas_ctx_t *)c)->apply(c, pos);
}

static void blas_cleanup(void *c) {
    blas_ctx_t *ctx = (blas_ctx_t *)c;
    if (ctx->U_lookup) {
        for (int i = 0; i < ctx->lookup_size; ++i) free(ctx->U_lookup[i]);
        free(ctx->U_lookup);
    }
    free(ctx->tmp); free(ctx->rho); free(ctx);
}

static size_t blas_touch_bytes(qubit_t qubits) {
    size_t dim = (size_t)1 << qubits;
    return 2 * dim * dim * sizeof(cplx_t);
}

static int blas_available(qubit_t qubits) {
    return qubits <= BLAS_MAX_QUBITS;
}

const bench_backend_t blas_backend = {
    .name        = "BLAS",
    .init        = blas_init,
    .apply       = blas_apply,
    .cleanup     = blas_cleanup,
    .touch_bytes = blas_touch_bytes,
    .available   = blas_available
};
